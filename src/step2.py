"""
Within this step we run the following :

    Take the resulting tables from oligodb, (usearch -search_oligodb)
        op_lib_dir/01-us_ops/{lib}_insbc.txt 
    which have the following columns:
        query, target, id, alnlen, mism, opens, qlo, qhi, tlo, thi, evalue
    Column names:
        query: Read name
        target: Oligo name (flanking region match)
        id: Percent identity
        alnlen: Alignment length
        mism: # of mismatch columns
        opens: # of gap open columns
        qlo: Query-start (1-based start) ## despite usearch manual online saying 0-index...
        qhi: Query-end (1-based start)
        tlo: Target-start (1-based start)
        thi: Target-end (1-based start)
        evalue: Karlin-Altschul E Value
    
    We filter hit tables using "id_cutoff" in cfg file and exclude reads with multiple 
    hits to the same flanking oligo (chimeras). We generate .tab files for barcode and 
    insert positions, separately, and use this to extract sequences into fasta and fq files.

    Extractions from this step are computed to two dataframes (tables)
        fns: lib + ['bcs_tab'/'ins_tab']
    Sequence extractions from this step are computed to fasta and fq files
        fns: lib + _bc.fasta or _ins.fq]
    in the directory:
        op_lib_dir/02-pos_ins_bc/

Key functions:
    run_step_2_singlelib
        load_and_filter_tables_from_oligosearch
        extract_insert_BC_pos
            fiveThreeRemoveAndMerge
            writePosBCandPosInsertToFile
        extract_by_pos

"""

import re
import sys
import os
import logging
import json
import shutil
import pandas as pd
import numpy as np
from util_mod import force_create_dir
from typing import List, Set, Tuple, Dict
from collections import Counter


def run_step_2_singlelib(op_lib_dir: str, lib_name: str, cfg_d: Dict, interactive=False) -> None:
    """
    Args:
        cfg_d (d): is the entire config_dict
        op_lib_dir (path to lib dir)
        lib_name (str) name of library

    Description:
        General extraction from files, write two dataframes (tables) out to files at
            op_lib_dir/02-pos_ins_bc
        Uses positions to extract sequences of barcode and inserts to fasta and fastq, respectively.
    """

    print("\nRunning step 2 for lib " + lib_name)
    logs_dir = os.path.join(os.path.dirname(op_lib_dir), "Logs")
    oligos_op_dir = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["1"])

    dir_fs: List[str] = os.listdir(oligos_op_dir)

    # A list of oligo names like ['5bc', '3bc', '5insert', '3insert']
    oligo_names: List[str] = get_oligo_names_from_oligos_db(cfg_d)

    s2_cfg = cfg_d["step_2"]

    # Part 1: filter search_oligodb output and obtain insert and barcode positions
    data_filtered, summary_d = load_and_filter_tables_from_oligosearch(
        oligos_op_dir, lib_name, oligo_names, cfg_d, id_cutoff=s2_cfg["oligo_id_cutoff"]
    )

    # Part 2: extract insert and barcode positions into tables
    pos_op_dir, nBadBarcode, bbCtr = extract_insert_BC_pos(
        data_filtered, lib_name, cfg_d, op_lib_dir
    )

    summary_d["num_barcodes_excluded"] = nBadBarcode
    summary_d["badBarcodeCounter"] = bbCtr

    # Write summary_d to log dir
    write_summary_to_log(logs_dir, summary_d, lib_name)

    # Part 3: extract barcode and insert sequences
    # Working on barcodes
    bc_fp = os.path.join(pos_op_dir, lib_name + cfg_d["d"]["fns"]["2"]["bcs_tab"])
    if not os.path.exists(bc_fp):
        raise Exception(f"Could not find expected barcode file at {bc_fp}")
    # Dict of read: tuple(start, end) which are both int
    bc_extract_d = import_id_to_extract(bc_fp, op_lib_dir=op_lib_dir)

    # Working on inserts
    ins_fp = os.path.join(pos_op_dir, lib_name + cfg_d["d"]["fns"]["2"]["ins_tab"])
    if not os.path.exists(ins_fp):
        raise Exception(f"Could not find expected insert file at {ins_fp}")
    ins_extract_d = import_id_to_extract(ins_fp, op_lib_dir=op_lib_dir)

    fq_fp: str = os.path.join(
        oligos_op_dir, lib_name + cfg_d["d"]["fns"]["1"]["trimmed"]
    )
    bc_log_list, ins_log_list = extract_by_pos(
        lib_name, fq_fp, bc_extract_d, ins_extract_d, pos_op_dir, cfg_d
    )
    with open(os.path.join(logs_dir, lib_name + "_step2_extraction.txt"), "w") as g:
        g.write("BARCODE INFO:\n")
        g.write("\n".join(bc_log_list))
        g.write("INSERT INFO:\n")
        g.write("\n".join(ins_log_list))

    print("Finished Step 2.")

    return None


def load_and_filter_tables_from_oligosearch(
    tables_dir, lib, oligo_names, cfg_d, id_cutoff=90.0
) -> Tuple[pd.DataFrame, Dict]:
    """
    Args:
        tables_dir (str): Path to directory with tables
        libs list<str>: Name of libraries
        oligo_names list<str>: Names of all oligos used in search_oligodb.
            4 oligos flanking the two regions we are looking to extract.
            e.g. ['5insert', '3insert', '3bc', '5bc']
        id_cutoff (float): Min percent ID match to oligos from usearch to keep reads
    Description:
        Imports output from usearch -search_oligodb
        filter_single_lib is called to filter tables for each library.
    """

    lib_fp = os.path.join(tables_dir, lib + cfg_d["d"]["fns"]["1"]["insbc"])

    if not os.path.isfile(lib_fp):
        raise Exception(f"Search oligodb table result at {lib_fp} not found.")

    lib_df, summary_d = filter_single_lib(lib_fp, oligo_names, lib, id_cutoff=id_cutoff)
    return lib_df, summary_d


def filter_single_lib(
    lib_fp: str, oligo_names: List[str], lib: str, id_cutoff=90
) -> Tuple[pd.DataFrame, Dict]:
    """
    Args:
        lib_fp (str): Path to output file from step 1's usearch oligodb (lib + cfg_d['d']['fns']['1']['insbc'])
        oligo_names (list<str>)  A list like ['5bc', '3bc', '5insert', '3insert']
        id_cutoff (float): Min percent ID match to oligos from usearch to keep reads
    Description:
        Filter by percent ID >= id_cutoff and ALL 4 of the oligos were found once per read.
        We remove queries with multiple hits to a single oligo - likely due to ligation during library prep.
    Returns:
        lib_df (pd.DataFrame):
                    query: read name
                    target: target oligo read name
                    id: Percent identity
                    alnlen: Alignment length
                    mism: # of mismatch columns
                    opens: # of gap open columns
                    qlo: Query-start (1-indexed)
                    qhi: Query-end (1-indexed)
                    tlo: Target-start (1-indexed)
                    thi: Target-end (1-indexed)
                    evalue: Karlin-Altschul E Value
        summary_d ():
                    "tot_reads >= ID_cutoff": # of reads after filtering for ID_cutoff
                    "concatenated_reads_excluded": # of concatenated, removed
                    "tot_reads_kept": final # of reads kept after all filters
                    "perc_reads_kept":
                    "oligo_info_d": info on each oligo search, stats on
                        <total_oligo_hits, total_reads_w_hits, reads_w_multiple_hits, perc_reads_w_multiple_hits>
    """
    dtypes = {
        "query": str,
        "target": str,
        "id": np.float64,
        "alnlen": np.int32,
        "mism": np.int32,
        "opens": np.int32,
        "qlo": np.int32,
        "qhi": np.int32,
        "tlo": np.int32,
        "thi": np.int32,
        "evalue": np.float64,
    }
    col_nms = [
        "query",
        "target",
        "id",
        "alnlen",
        "mism",
        "opens",
        "qlo",
        "qhi",
        "tlo",
        "thi",
        "evalue",
    ]

    # Import oligodb output file using dtypes above
    lib_df = pd.read_csv(lib_fp, sep="\t", names=col_nms, dtype=dtypes)

    # Stores read names as a list of sets to extract intersection
    unique_reads_list: List[Set[str]] = []

    # For each library we have 4 oligo dataframes stored in this list
    lib_dfs_list = []

    # This maintains all the read names which occur more than once in each oligo search output.
    concat_read_names = set()

    # This mintains data for each library and all oligos for output in Logs/
    oligo_info_d = {}

    # Below we need to find reads (queries) that hit all 4 oligos once.
    # A read needs to have all 4 matches in order to extract both barcode and insert sequences.
    for oligo in oligo_names:
        # e.g. oligo = '5bc'

        oligo_df = lib_df[lib_df["target"] == oligo]

        # filters by % ID cutoff
        oligo_df = oligo_df[oligo_df["id"] >= id_cutoff][
            ["query", "target", "qlo", "qhi"]
        ]
        # stores reads as a list of sets to extract intersection
        unique_reads_list.append(oligo_df["query"].unique())

        # compute stats for Logs/
        tot_oligo_hits = oligo_df.shape[0]
        tot_reads = oligo_df["query"].nunique()

        # reads with multiple hits to one oligo, concatenated amplicons
        concat_reads_df = oligo_df[oligo_df["query"].duplicated(keep="first")]
        uniq_concat_reads = concat_reads_df["query"].unique()
        nUniq_concat_reads = len(uniq_concat_reads)
        concat_read_names.update(set(uniq_concat_reads.tolist()))

        if tot_oligo_hits > 0:
            perc_concat_tot = round((nUniq_concat_reads / tot_reads) * 100, 2)
        else:
            perc_concat_tot = "NA"

        oligo_info_d[oligo] = {
            "total_oligo_hits": tot_oligo_hits,
            "total_reads_w_hits": tot_reads,
            "reads_w_multiple_hits": nUniq_concat_reads,
            "perc_reads_w_multiple_hits": perc_concat_tot,
        }
        lib_dfs_list.append(oligo_df)

    # Now we combine all dataframes for each oligo and exclude concatenated reads.
    lib_df = pd.concat(lib_dfs_list)
    del lib_dfs_list
    tot_hits = lib_df[
        "query"
    ].nunique()  # post-ID filtering, contains concatenated reads

    ls_remove = list(set(concat_read_names))
    print(
        f"Removing {len(ls_remove)} concatenated reads out of {tot_reads} reads (post-ID_cutoff filter) from lib {lib}."
    )
    lib_df = lib_df[~lib_df["query"].isin(ls_remove)]

    # extract common reads from all oligo tables: reads with single hits to all 4 oligos
    # Note: i goes from 1 -> 3 (inclusive), so it intersects the reads
    # from all oligo subsets.
    reads_common = set(unique_reads_list[0])
    for i in range(1, len(unique_reads_list)):  # 1,2,3
        reads_common = reads_common.intersection(set(unique_reads_list[i]))
    lib_df = lib_df[lib_df["query"].isin(reads_common)]

    if tot_hits > 0:
        perc_reads_kept = round((lib_df["query"].nunique() / tot_hits) * 100, 2)
    else:
        perc_reads_kept = "NA"

    summary_d = {
        "tot_reads >= ID_cutoff": tot_hits,
        "concatenated_reads_excluded": len(ls_remove),
        "tot_reads_kept": lib_df["query"].nunique(),
        "perc_reads_kept": perc_reads_kept,
        "oligo_info_d": oligo_info_d,
    }

    return lib_df, summary_d


def extract_insert_BC_pos(data_filtered, lib, cfg_d, op_lib_dir) -> Tuple[str, int, Dict[int, int]]:
    """
    Args:
        cfg_d (dict): contains the following important info:
            "primer_info" ->
                "flanking_names" -> flank_d
                    flank_d: Contains these keys:
                       'BC_flanking' -> {'5': '{5name}', '3': '{3name}'}
                       'insert_flanking' -> {'5': '{5name}', '3': '{3name}'}

    Description:
        Insert region: end of 5' insert oligo to  start of 3' insert oligo
        Barcode region: end of 5' BC oligo to start of 3' BC oligo
        We need to get related groupings of 5' and 3'

        We call writePosBCandPosInsertToFile to write files out to op_lib_dir/02-pos_ins_bc
        at file names 'lib + cfg_d['d']['fns']['2']['bcs_tab']' and 'lib + cfg_d['d']['fns']['2']['ins_tab']'

    Returns:
        nBadBarcode
        bbCtr
    """

    data_grouped = data_filtered.groupby(["target"])  # target = oligo name

    # The length of this list should be 2, and each sublist should also have length of 2
    # since there are 4 oligos
    flank_d: Dict[str, Dict[str, str]] = cfg_d["primer_info"]["flanking_names"]
    pos_bc = pd.DataFrame()
    pos_ins = pd.DataFrame()

    # Below, we update the two DataFrames 'pos_bc' and 'pos_ins'
    res: Tuple = fiveThreeRemoveAndMerge(lib, data_grouped, flank_d, pos_bc, pos_ins)
    pos_bc: pd.DataFrame = res[0]
    pos_ins: pd.DataFrame = res[1]
    nBadBarcode: int = res[2]
    # bb = badBarcode: bb counter and bb list, which is bbCtr in reverse for plotting
    bbCtr: Dict[int, int] = res[3]
    bbL: List[Tuple[int, int]] = res[4]

    pos_op_dir: str = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["2"])
    force_create_dir(pos_op_dir)
    writePosBCandPosInsertToFile(lib, pos_bc, pos_ins, pos_op_dir, cfg_d, bbL)

    # in case these aren't deleted
    del pos_bc, pos_ins

    return pos_op_dir, nBadBarcode, bbCtr


def fiveThreeRemoveAndMerge(lib, data_grouped, flank_d, pos_bc, pos_ins) -> Tuple[pd.DataFrame, pd.DataFrame, int, Dict[int,int], List[Tuple[int, int]] ]:
    """
    Args:
        flank_d: Contains these keys:
           'BC_flanking' -> {'5': '{5name}', '3': '{3name}'}
           'insert_flanking' -> {'5': '{5name}', '3': '{3name}'}
    Description:
        This appends values to pos_bc and pos_ins dataframes.
        No need to return them, since they have only a single copy.
    Returns:
        pos_bc: A merged dataframe on query to barcode positions (1-based start)
        pos_ins: A merged dataframe on query to insert positions (1-based start)
        nBadBarcode: Barcodes that are not 20 nt long
        badBarcodeCounter: A dictionary of bad barcode length to occurrence
        badBarcodeList: A reverse sorted list of tuples with (occurrence, bad barcode length), used for plotting
    """
    ## merge the two flanking positions of BC
    ## each table has: ["query", "target", "qlo", "qhi"], grouped on "target"
    data5 = data_grouped.get_group(flank_d["BC_flanking"]["5"])
    data3 = data_grouped.get_group(flank_d["BC_flanking"]["3"])
    df_merged_bc = data5.merge(data3, on="query", how="inner")

    # checks if barcode region is 20 nt
    ls_bc_remove = (
        df_merged_bc[df_merged_bc["qlo_y"] - df_merged_bc["qhi_x"] - 1 != 20]["query"]
        .unique()
        .tolist()
    )
    df_bc_keep = df_merged_bc[df_merged_bc["qlo_y"] - df_merged_bc["qhi_x"] - 1 == 20]

    # good barcodes into df
    pos_bc = pos_bc.append(df_bc_keep[["query", "qhi_x", "qlo_y"]])

    ## merge the two flanking regions of insert
    data5 = data_grouped.get_group(flank_d["insert_flanking"]["5"])
    data3 = data_grouped.get_group(flank_d["insert_flanking"]["3"])
    df_merged_ins = data5.merge(data3, on="query", how="inner")

    # remove reads from insert table where barcode is not 20 nt long
    df_ins_keep = df_merged_ins[~df_merged_ins["query"].isin(ls_bc_remove)]
    pos_ins = pos_ins.append(df_ins_keep[["query", "qhi_x", "qlo_y"]])
    out_cols = ["read", "start", "end"]
    pos_bc = pos_bc.set_axis(out_cols, axis=1)
    pos_ins = pos_ins.set_axis(out_cols, axis=1)

    # adjust to obtain positions of bc and ins regions (1-based starts)
    # recall that qlo_y = start of of 3' oligo hit; qhi_x = end of 5' oligo hit; 1-based starts
    pos_bc.loc[:, "start"] = pos_bc["start"] + 1
    pos_ins.loc[:, "start"] = pos_ins["start"] + 1
    pos_bc.loc[:, "end"] = pos_bc["end"] - 1
    pos_ins.loc[:, "end"] = pos_ins["end"] - 1

    # Computing lengths of bad barcodes
    print("Computing lengths of bad barcodes.")
    bc_rm_df = df_merged_bc[df_merged_bc["qlo_y"] - df_merged_bc["qhi_x"] - 1 != 20]
    a = bc_rm_df["qlo_y"].to_list()  # start of 3' oligo hit, 0-based start
    b = bc_rm_df["qhi_x"].to_list()  # end of 5' oligo hit
    # a-b +1 (from start)-1 (from end) +1 (for len) (positions are flanking barcode, not start and end of barcode)
    badBarcode_len: List[int] = [a[i] - b[i] - 1 for i in range(len(a))]

    badBarcodeCounter = Counter(badBarcode_len)  # len: freq
    badBarcodeList: List[Tuple[int, int]] = [
        (v, k) for k, v in badBarcodeCounter.items()
    ]
    badBarcodeList = sorted(
        badBarcodeList, reverse=True, key=lambda x: x[0]
    )  # freq: len
    nBadBarcode = len(ls_bc_remove)
    print(
        f"{len(ls_bc_remove)} reads were excluded from "
        + f"{lib} due to barcode length not = 20 nt."
    )
    if pos_bc.shape[0] != pos_ins.shape[0]:
        raise Exception(
            "WARNING: Expecting number of read with bc = number of reads with insert"
        )
    return pos_bc, pos_ins, nBadBarcode, badBarcodeCounter, badBarcodeList


def writePosBCandPosInsertToFile(lib, pos_bc, pos_ins, pos_op_dir, cfg_d, bbL) -> None:

    bc_path = os.path.join(pos_op_dir, lib + cfg_d["d"]["fns"]["2"]["bcs_tab"])
    pos_bc.to_csv(bc_path, sep="\t", header=False, index=False)
    ins_path = os.path.join(pos_op_dir, lib + cfg_d["d"]["fns"]["2"]["ins_tab"])
    pos_ins.to_csv(ins_path, sep="\t", header=False, index=False)
    print(
        "Files written to: "
        + lib
        + cfg_d["d"]["fns"]["2"]["bcs_tab"]
        + " and "
        + lib
        + cfg_d["d"]["fns"]["2"]["ins_tab"]
    )

    bad_bc_fp = os.path.join(pos_op_dir, lib + cfg_d["d"]["fns"]["2"]["bad_BC_lens"])
    with open(bad_bc_fp, "w") as g:
        g.write("nTimesSeen\tBarcodeLength\n")
        for x in bbL:
            g.write(f"{x[0]}\t{x[1]}\n")
    print("Wrote bad barcode lengths to " + lib + cfg_d["d"]["fns"]["2"]["bad_BC_lens"])


def write_summary_to_log(logs_dir, summary_d, lib) -> None:
    fp = os.path.join(logs_dir, lib + "_step2_oligo_hits.json")
    with open(fp, "w") as g:
        g.write(json.dumps(summary_d, indent=2))
    print("Wrote summary_d to file at " + fp)


def get_oligo_names_from_oligos_db(cfg_d) -> List[str]:
    oligo_names = []
    fl_nm = cfg_d["primer_info"]["flanking_names"]
    for k in fl_nm.keys():
        for n in ["5", "3"]:
            oligo_names.append(fl_nm[k][n])

    return oligo_names


def import_id_to_extract(in_file, op_lib_dir=".") -> Dict[str, Tuple[int, int]]:
    """
    Imports file and generates a dictionary mapping read to location info

    This is run on either '_bc.tab' or '_ins.tab':
        read_name start end
    """
    extract_d = {}

    with open(in_file) as id_handle:
        for line in id_handle:
            line_split = line.split("\t")
            read = line_split[0]
            start = int(line_split[1])
            end = int(line_split[2])
            if read in extract_d:
                print(f"WARNING: duplicate read found in {in_file}.")
            extract_d[read] = (start, end)
    print(f"Found {len(list(extract_d.keys()))} read identifiers in {in_file}")

    return extract_d


def extract_by_pos(lib, input_file, bc_extract_d, ins_extract_d, op_lib_dir, cfg_d) -> Tuple[List[str], List[str]]:
    """
    Args:
        input_file: fastq of trimmed reads at lib + cfg_d["d"]["fns"]["1"]["trimmed"]
        bc_extract_d (Dictionary mapping read name -> start and end of BC seq)
        ins_extract_d (Dictionary mapping read name -> start and end of insert seq)
        Note. positions are 1-based starts, need to convert to 0-index used by Python
    Description:
        We take .tab file of positions and extract regions from each read
        We write out to files: lib + _[bc.fasta/ins.fq]
        ["bcs_fasta"] is a FASTA file with read id and barcode sequence
        ["ins_fq"] is a FASTQ file of extracted insert sequences
    """

    # output file paths
    ins_op_fp = os.path.join(op_lib_dir, lib + cfg_d["d"]["fns"]["2"]["ins_fq"])
    bc_op_fp = os.path.join(op_lib_dir, lib + cfg_d["d"]["fns"]["2"]["bcs_fasta"])
    bc_count, ins_count = 0, 0
    line_num = 1
    bc_reads_dict = {"kept": set(), "excluded": set()}
    ins_reads_dict = {"kept": set(), "excluded": set()}
    in_handle = open(input_file, "r")
    bc_out_handle = open(bc_op_fp, "w")
    ins_out_handle = open(ins_op_fp, "w")

    line = in_handle.readline()  # fastq data
    while line != "":
        read_line = line.strip()
        seq = in_handle.readline().strip()
        separator = in_handle.readline().strip()
        qual = in_handle.readline().strip()
        if separator != "+":
            raise Exception(
                "Expecting '+' as separator, instead got:\n"
                + separator
                + "\nAt line number "
                + str(line_num + 2)
            )
        # Keeping the first term of read ID for output files
        # exclude first char from pacbio read IDs: @m64044_211202_122350/44/ccs
        read_id = read_line[1:]

        # extract sequences to fasta and fastq files
        bc_count = extract_bc_seq(
            read_id,
            bc_extract_d,
            bc_reads_dict,
            read_line,
            seq,
            qual,
            bc_out_handle,
            bc_count,
        )
        ins_count = extract_ins_seq(
            read_id,
            ins_extract_d,
            ins_reads_dict,
            read_line,
            seq,
            qual,
            ins_out_handle,
            ins_count,
        )
        line = in_handle.readline()
        line_num += 4

    in_handle.close()
    ins_out_handle.close()
    bc_out_handle.close()

    if ins_count < len(ins_extract_d):
        print(
            "WARNING: %i IDs from _ins.tab not found in %s"
            % (len(ins_extract_d) - ins_count, input_file)
        )
    if bc_count < len(bc_extract_d):
        print(
            "WARNING: %i IDs from _bc.tab not found in %s"
            % (len(bc_extract_d) - bc_count, input_file)
        )
    # Write out list of discarded read IDs (first term where bc or ins were not extracted
    bc_log_list: List[str] = report_and_write_discarded(
        lib, op_lib_dir, bc_reads_dict, bc_count, "bc", cfg_d
    )
    ins_log_list: List[str] = report_and_write_discarded(
        lib, op_lib_dir, ins_reads_dict, ins_count, "ins", cfg_d
    )

    return (bc_log_list, ins_log_list)


def extract_bc_seq(
    read_id,
    bc_extract_d,
    bc_reads_dict,
    read_line,
    seq,
    qual,
    bc_out_handle,
    bc_kept_count: int,
) -> int:
    """
    Extracts barcode sequences using positions and write to file.
    """
    if read_id in bc_extract_d:
        if read_id not in bc_reads_dict["kept"]:
            bc_reads_dict["kept"].add(read_id)

            # convert 1-based starts to 0-indexed with end being exclusive
            start = bc_extract_d[read_id][0] - 1
            end = bc_extract_d[read_id][1]
            seq_wanted = seq[start:end]

            bc_kept_count += 1
            # Write to fasta file
            bc_out_handle.write(f">{read_id}\n{seq_wanted}\n")
    else:
        bc_reads_dict["excluded"].add(read_id)

    return bc_kept_count


def extract_ins_seq(
    read_id,
    ins_extract_d,
    ins_reads_dict,
    read_line,
    seq,
    qual,
    ins_out_handle,
    ins_kept_count,
) -> int:
    """
    Extracts insert sequences using positions and write to file.
    """
    if read_id in ins_extract_d:
        if read_id not in ins_reads_dict["kept"]:
            ins_reads_dict["kept"].add(read_id)

            # convert 1-based starts to 0-indexed with end being exclusive
            start = ins_extract_d[read_id][0] - 1
            end = ins_extract_d[read_id][1]
            seq_wanted = seq[start:end]
            qual_wanted = qual[start:end]
            ins_kept_count += 1

            # Write to fastq file
            ins_out_handle.write(
                "@%s\n%s\n+\n%s\n" % (read_line[1:], seq_wanted, qual_wanted)
            )
    else:
        ins_reads_dict["excluded"].add(read_id)

    return ins_kept_count


def report_and_write_discarded(lib, op_lib_dir, reads_dict, count, seq_type, cfg_d) -> List[str]:
    """
    Note seq_type (string) is one of  'ins' or 'bc'
    """

    log_list: List[str] = []
    uniq_excluded = set(reads_dict["excluded"])

    log_list += [
        f"Reads kept: {len(reads_dict['kept'])};",
        f"Reads excluded: {len(reads_dict['excluded'])};",
        f"Duplicate reads excluded: {len(reads_dict['excluded']) -len(uniq_excluded)}",
        f"Extracted from {count} reads at {lib+cfg_d['d']['fns']['1']['trimmed']}.",
    ]
    # write read names not kept to _no_ins.txt or _no_bc.txt
    if len(reads_dict["excluded"]) != 0:
        crt_out_fp = os.path.join(op_lib_dir, lib + "_no_" + seq_type + ".txt")
        log_list += [
            f"{len(reads_dict['excluded'])} reads were not kept and wrote to {lib}_no_{seq_type}.txt\n"
        ]
        with open(crt_out_fp, "w") as out_f:
            for item in reads_dict["excluded"]:
                out_f.write("%s\n" % item)
        print("Wrote reads that weren't kept to file " + crt_out_fp)

    return log_list


def main():

    help_str = "python3 src/step2.py cfg_json inp_dir op_dir(tmp) 1"
    help_str = "OR\n"
    help_str = "python3 src/step2.py inp_dir oligos_dir 2"
    args = sys.argv
    if args[-1] not in ["1", "2"]:
        print(help_str)
        sys.exit(1)
    elif args[-1] == "1":
        test(args)
    else:
        intermediate_tests(args, tp=3)

    return None


if __name__ == "__main__":
    main()
