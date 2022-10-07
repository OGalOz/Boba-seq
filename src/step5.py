"""
In step 5, we combine barcodes to .paf output from minimap2 
    (insert hits against a reference genome).
    A series of filters are then applied to identify the best 
    position(s) for each barcode.
    Redundant barcodes in empty vectors can result in multiple inserts
    per barcodes, so barcodes that map to multiple locations are kept.
    This is then used to map insert fragments to genes based on
    coverage of complete ORFs.
    
Many key functions are in:
    collapse_bcs.py
"""

import os
import sys
import pandas as pd
import datetime
import json
from collections import defaultdict, Counter
from typing import Dict, List, TypeVar, Tuple
from parse_paf import parse_paf_file, get_read_name_to_info_dict
from collapse_bcs import bc_df_collapse_steps, export_collapsed_bc_files
from import_gff import DubSeq_import_gff
from validate import load_entire_cfg, verify_cfg_d, validate_collapse_params
import contig_collider


# Variable type - just for documentation purposes
T = TypeVar("T")


def run_step_5_singlelib(
    op_lib_dir, lib_name, cfg_d, pre_created_mid=False, export_bool=False
) -> str:
    """
    Args:
        op_lib_dir (str): tmp_dir/lib_name
        lib_name (str): name of the library
        cfg_d (d): Config dictionary (active, not filepath)
        pre_created_mid (bool) If debugging this file alone, you can run
                            from halfway if the halfway results have already
                            been created_
        export_bool (bool): If debugging this file alone, you can export
                            the results that pre_created_mid would take

    Key script files used: parse_paf.py, collapse_bcs.py, contig_collider.py
    """
    print("\nRunning step 5 for lib " + lib_name)
    log_list: List[str] = []

    minimap_dir = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["4"])
    if not os.path.exists(minimap_dir):
        raise RuntimeError("minimap2 dir not found in expected dir: " + op_lib_dir)

    # paf_fp: output from minimap2; barcode_fasta_fp: list of barcodes in fasta from step 2
    paf_fp = os.path.join(minimap_dir, lib_name + cfg_d["d"]["fns"]["4"]["minimap_op"])
    # converts .paf to df, filters using cfg_d["minimap_qual_min"]
    # keep first hit for reads with multiple hits (longer aln_len on top)
    # ignores info output, already wrote to Logs in step4
    paf_df, paf_multi_hits, info = parse_paf_file(paf_fp, cfg_d)
    del paf_multi_hits, info  # used in step4

    print("Getting read name to info dict")
    read2info: Dict = get_read_name_to_info_dict(paf_df)
    barcode_fasta_fp = os.path.join(
        op_lib_dir,
        cfg_d["d"]["steps2dirs"]["2"],
        lib_name + cfg_d["d"]["fns"]["2"]["bcs_fasta"],
    )
    read2bc: Dict[str, str] = create_read_to_barcode_dict(barcode_fasta_fp)
    log_list += ["Total barcode num from step 2: " + str(len(read2bc))]
    cols = ["bc"] + list(paf_df.columns)

    # Combines read2bc dictionary with paf output from minimap2, saves BCs without hits to genomes separately
    comb_df, missing_barcodes = combine_barcodes_with_other_data(
        read2info, read2bc, paf_df, cols
    )
    log_list += [
        "Barcodes excluded where insert has e.e. of >10 or has no hit from minimap2: "
        + str(len(missing_barcodes))
    ]
    log_list += [
        f"Fraction of barcodes with hits: {round(1-(len(missing_barcodes)/len(read2bc)), 3)}"
    ]
    del paf_df

    # ------------------------------------------------
    dfs_dir = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["5"])
    if os.path.exists(dfs_dir):
        print("Dataframe output directory found.")
    else:
        os.mkdir(dfs_dir)

    cp = validate_collapse_params(cfg_d)

    # input comb_df: filtered for perc_match, perc_cov only
    # mapping algorithm, generates various groups of barcodes
    # best_mappings_d: filtered for above, split by regions, filtered for BC_ratio, diff_frag_len to select best positions
    # best_mappings_df: cols ["bc", "tstart", "tend", "contig", "multi_frag", "strand"]
    (
        best_mappings_df,
        best_mappings_d,
        ls_bc_multi_mappings,
        multi_contig_bcs,
        failed_bcs,
        new_log,
    ) = bc_df_collapse_steps(comb_df, cp, dfs_dir, cfg_d, lib_name)
    log_list += new_log

    # Write to files: failed barcodes, barcodes mapped to multiple contigs, best_mappings_df
    # Plots histogram of mapped location to frequency (after all filters)
    new_log = export_collapsed_bc_files(
        best_mappings_df,
        ls_bc_multi_mappings,
        multi_contig_bcs,
        failed_bcs,
        dfs_dir,
        op_lib_dir,
        cfg_d,
        lib_name,
    )
    log_list += new_log

    # Imports genome annotation files
    gff_fp: str = get_gff_fp(cfg_d, lib_name)
    # Note that gff_df is already sorted
    gff_df: pd.DataFrame = DubSeq_import_gff(gff_fp)

    # Here we check that contig names match within best_mappings_df and gff_df
    gff_df, bc_df = contig_collider.match_contig_names(
        gff_df, best_mappings_df, debug=True
    )

    # MOST INVOLVED PART OF THE PROGRAM: <<<---
    gene_fragment_overlap(gff_df, bc_df)
    # --->>>

    # splits out locus_tags mapped to each BC into separate rows
    bc_df = reshape_dataframe(bc_df)

    # WRITING OUT BC_LOC_DF
    bc_df.to_csv(
        os.path.join(dfs_dir, cfg_d["d"]["fns"]["5"]["bc_loc_df"]),
        sep="\t",
        index=False,
    )
    # WRITING OUT GENES_COUNT_DF
    gff_df.to_csv(
        os.path.join(dfs_dir, cfg_d["d"]["fns"]["5"]["genes_count_df"]),
        sep="\t",
        index=False,
    )

    with open(
        os.path.join(os.path.dirname(op_lib_dir), "Logs", lib_name + "_step5_log.txt"),
        "w",
    ) as g:
        g.write("\n".join(log_list))

    print("Finished step 5.")
    return dfs_dir


def midway_run1(
    op_lib_dir, lib_name, comb_df_fp: str, cfg_d, bc_to_loc_dicts=None
) -> str:

    print("Starting step 5 midway")
    log_list: List[str] = []
    cp = validate_collapse_params(cfg_d)

    comb_df = pd.read_table(comb_df_fp)

    dfs_dir = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["5"])
    if os.path.exists(dfs_dir):
        print("Step 5 output directory found at ", dfs_dir)
    else:
        os.mkdir(dfs_dir)

    # input comb_df: filtered for perc_match, perc_cov only
    # mapping algorithm, generates various groups of barcodes
    # best_mappings_d: filtered for above, split by regions, filtered for BC_ratio, diff_frag_len to select best positions
    # best_mappings_df: cols ["bc", "tstart", "tend", "contig", "multi_frag", "strand"]
    (
        best_mappings_df,
        best_mappings_d,
        ls_bc_multi_mappings,
        multi_contig_bcs,
        failed_bcs,
        new_log,
    ) = bc_df_collapse_steps(
        comb_df,
        cp,
        dfs_dir,
        cfg_d,
        lib_name,
        print_dbg=True,
        bc_to_loc_dicts=bc_to_loc_dicts,
    )
    log_list += new_log

    # Write to file: failed barcodes, barcodes mapped to multiple contigs, best_mappings_df
    # Plots histogram of mapped location to frequency (after all filters)
    new_log = export_collapsed_bc_files(
        best_mappings_df,
        ls_bc_multi_mappings,
        multi_contig_bcs,
        failed_bcs,
        dfs_dir,
        op_lib_dir,
        cfg_d,
        lib_name,
    )
    log_list += new_log

    # Imports genome annotation files
    gff_fp: str = get_gff_fp(cfg_d, lib_name)
    # Note that gff_df is already sorted
    gff_df: pd.DataFrame = DubSeq_import_gff(gff_fp)

    # Here we check that contig names match within best_mappings_df and gff_df
    gff_df, bc_df = contig_collider.match_contig_names(
        gff_df, best_mappings_df, debug=True
    )

    # MOST INVOLVED PART OF THE PROGRAM: <<<---
    gene_fragment_overlap(gff_df, bc_df)
    # --->>>

    # splits out locus_tags mapped to each BC into separate rows
    bc_df = reshape_dataframe(bc_df)

    # WRITING OUT BC_LOC_DF
    bc_df.to_csv(
        os.path.join(dfs_dir, cfg_d["d"]["fns"]["5"]["bc_loc_df"]),
        sep="\t",
        index=False,
    )
    # WRITING OUT GENES_COUNT_DF
    gff_df.to_csv(
        os.path.join(dfs_dir, cfg_d["d"]["fns"]["5"]["genes_count_df"]),
        sep="\t",
        index=False,
    )

    with open(
        os.path.join(os.path.dirname(op_lib_dir), "Logs", lib_name + "_step5_log.txt"),
        "w",
    ) as g:
        g.write("\n".join(log_list))

    print("Finished step 5.")
    return dfs_dir


def reshape_dataframe(bc_df: pd.DataFrame) -> pd.DataFrame:
    """
    Description:
        Initially, within the bc dataframe, we have one row
        per fragment, and that row contains, within the column "locus_tag",
        both encoded genes and whether the direction is the same.
        We want to split each row into many rows, where each new row is
        a fragment with a different locus tag, and there being a
        new column called 'directionality'. If 'directionality' is 'same',
        then the fragment and gene are on the same strand, otherwise
        'directionality' will be set to 'opp', for (opposite).
    """
    print("Starting to reshape dataframe.")

    # We prepare the new dataframe by keeping a list of series
    new_dataframe: List[pd.Series] = []

    # We initialize the 'directionality' column
    bc_df["directionality"] = [""] * bc_df.shape[0]

    for ix, row in bc_df.iterrows():
        orig_locus_tags: str = row["locus_tag"]
        if not isinstance(orig_locus_tags, str):
            continue
        if "|" in orig_locus_tags:
            locus_tags = orig_locus_tags.split("|")
            for lct in locus_tags:
                if ",*" in lct:
                    lc, d = lct.split(",*")
                    new_row = row.copy(deep=True)
                    new_row.at["directionality"] = d
                    new_row.at["locus_tag"] = lc
                    new_dataframe.append(new_row)
        else:
            new_dataframe.append(row)

    print("Finished reshaping dataframe.")

    genes_on_single_rows = pd.DataFrame(new_dataframe)
    genes_on_single_rows.drop_duplicates(inplace=True, ignore_index=True)

    return genes_on_single_rows


def gene_fragment_overlap(
    genes: pd.DataFrame, fragments: pd.DataFrame, dbg=False
) -> None:
    """
    Desc: First we initialize pandas Series that are the same length as genes
    and fragments dataframes as part of these dataframes.
    Note this function uses DataFrame.at() to update values to df, does not explicitly return output
    genes has to have the following columns:
        pos_from
        pos_to
        contig
        strand
        locus_tag
    fragments has to have the following columns:
        tstart
        tend
        contig
        multi_frag
        strand
    """

    genes.sort_values(by="pos_from", inplace=True, ignore_index=True)
    fragments.sort_values(by="tstart", inplace=True, ignore_index=True)
    # Adding columns:
    genes["fragment_count"] = pd.Series(
        [0] * genes.shape[0], index=genes.index, dtype="int32"
    )
    genes["nSame_direc_frag"] = pd.Series(
        [0] * genes.shape[0], index=genes.index, dtype="int32"
    )
    fragments["gene_count"] = pd.Series(
        [0] * fragments.shape[0], index=fragments.index, dtype="int32"
    )
    fragments["locus_tag"] = pd.Series(
        [""] * fragments.shape[0], index=fragments.index, dtype="str"
    )
    tmp_frag = fragments[
        ["tstart", "tend", "gene_count", "locus_tag", "contig", "strand"]
    ]

    genes_pos_from = genes["pos_from"]
    nFragRows = fragments.shape[0]
    nGeneRows = genes.shape[0]
    print("Beginning fragment to gene mapping...")
    print("Number of fragments to iterate through:", nFragRows)
    print("Number of genes in this genome:", nGeneRows)

    # The following variables 'frac', 'nLoops', 'start_time'
    # are all used just to measure the time of the program,
    # they aren't essential for its functioning.
    frac = 1000
    nLoops = nFragRows // frac
    start_time = datetime.datetime.now()
    g_ix = 0
    gc_ix = tmp_frag.columns.get_loc("gene_count")
    lt_ix = tmp_frag.columns.get_loc("locus_tag")
    fc_ix = genes.columns.get_loc("fragment_count")
    nSame_ix = genes.columns.get_loc("nSame_direc_frag")

    for i in range(nLoops):
        for fi in range(i * frac, (i + 1) * frac):
            g_ix, break_bool = mid_gene_fragment_overlap_count(
                fi,
                tmp_frag,
                gc_ix,
                lt_ix,
                nSame_ix,
                fc_ix,
                g_ix,
                genes_pos_from,
                nGeneRows,
                nFragRows,
                nLoops,
                genes,
                start_time,
                dbg=dbg,
            )
            if break_bool:
                break
        if g_ix == nGeneRows:
            print(f"Reached final gene at fi {fi}, breaking.")
            break

    # This part of the program is the remaining lines that are outside
    # of the loops that were used to measure the timing of the program
    # It's like if you take 32 // 5 = 6, but the remainder is 2, so this
    # part runs the last 2. Except for in the case of this program,
    # it would be more like 35,853 // 1000 = 35, and remainder would be 853.
    # nLoops = 35, frac = 1000, nRows = 35,853
    if g_ix < nGeneRows:
        for fi in range(nLoops * frac, nFragRows):
            g_ix, break_bool = mid_gene_fragment_overlap_count(
                fi,
                tmp_frag,
                gc_ix,
                lt_ix,
                nSame_ix,
                fc_ix,
                g_ix,
                genes_pos_from,
                nGeneRows,
                nFragRows,
                nLoops,
                genes,
                start_time,
                dbg=dbg,
            )
            if break_bool:
                break
    else:
        print("Skipped computing last segment since fragments are past last gene")

    fragments["gene_count"] = tmp_frag["gene_count"]
    fragments["locus_tag"] = tmp_frag["locus_tag"]

    print("Completed fragment-to-gene counts")


def mid_gene_fragment_overlap_count(
    fi,
    tmp_frag,
    gc_ix,
    lt_ix,
    nSame_ix,
    fc_ix,
    g_ix,
    genes_pos_from,
    nGeneRows,
    nFragRows,
    nLoops,
    genes,
    start_time,
    dbg=False,
) -> Tuple[int, bool]:
    """
    Args:
        fi : frag index (row number)
        gc_ix: frag df column index of "gene_count"
        lt_ix = frag df column index of "locus_tag"
        nSame_ix = genes df column index of  "nSame_direc_frag"
        fc_ix = genes df column index of "fragment_count"
    """

    f = tmp_frag.iloc[fi]
    # This is the starting gene location within the sorted gene dataframe
    g_ix = update_gix(g_ix, genes_pos_from, f["tstart"], nGeneRows)
    for gi in range(g_ix, nGeneRows):
        # Now we go down the gene list until the gene starting point
        # is after the fragment's ending point, and then we move
        # to the next fragment, but keep the starting gene
        g: pd.Series = genes.iloc[gi]
        if (
            f["tstart"] <= g.pos_from
            and f["tend"] >= g.pos_to
            and f["contig"] == g.contig
        ):
            direc: str = "opp"
            if f["strand"] == g.strand:
                direc = "same"
            direc_bool = 0 if direc == "opp" else 1
            # gene_count +1 regardless of directionality
            tmp_frag.iat[fi, gc_ix] += 1
            if "|" in g.locus_tag:
                raise Exception(
                    "The character '|' cannot be in locus tags."
                    + f" Note locus tag: {g.locus_tag}, genes_df row {gi}."
                )
            tmp_frag.iat[fi, lt_ix] += g.locus_tag + ",*" + direc + "|"
            genes.iat[gi, nSame_ix] += direc_bool
            genes.iat[gi, fc_ix] += 1
        elif g.pos_from > f["tend"]:
            break
    break_bool = False
    if g_ix == nGeneRows:
        print(f"reached final gene at fi {fi}, breaking.")
        break_bool = True
    return g_ix, break_bool


def update_gix(g_ix, genes_pos_from, f_start, N) -> int:
    # N is the total number of values in genes_pos_from (nRows of g_df)
    # Gene index moves down whenever starting position is less than
    # the start of the fragment. As soon as the starting position is after
    # the start of the fragment the loop stops and the gene index no longer
    # increases.
    while g_ix < N and genes_pos_from[g_ix] < f_start:
        g_ix += 1

    return g_ix


def get_gff_fp(cfg_d, lib_name) -> str:
    if "gff_fp" in cfg_d:
        return cfg_d["gff_fp"]
    else:
        lib_ix = cfg_d["lib_names"].index(lib_name)
        return os.path.join(cfg_d["lib_genome_dir"], cfg_d["lib_genome_gffs"][lib_ix])


def filter_gff(gff_df) -> pd.DataFrame:
    # This only returns the rows that represent important genes (CDS)
    return gff_df[gff_df["type"] == "CDS"]


def combine_barcodes_with_other_data(
    read2info: Dict, read2bc: Dict, paf_df: pd.DataFrame, cols: List[str]
) -> Tuple[pd.DataFrame, List[Tuple[str, str]]]:
    # columns: bc, read_name, qlen, qstart, qend, strand, target,
    # tlen, tstart, tend, matches, aln_len, qual, perc_match_cov, perc_match

    print("Beginning to combine barcodes with .paf data.")
    op_l: List[List[T]] = []
    # Each tuple is (read_name, barcode)
    missing_barcode_reads: List[Tuple[str, str]] = []
    read_names = list(read2bc.keys())
    nBC_reads = len(read_names)
    L = 10**4
    for i in range(len(read_names)):
        read_name = read_names[i]
        if read_name in read2info:
            new_l = [read2bc[read_name]] + paf_df.iloc[read2info[read_name]].to_list()
            op_l.append(new_l)
        else:
            missing_barcode_reads.append((read_name, read2bc[read_name]))

    op_df = pd.DataFrame(op_l, columns=cols)
    return (op_df, missing_barcode_reads)


def create_read_to_barcode_dict(barcode_fasta_fp: str) -> Dict[str, str]:
    """
    Fasta file looks like
    >Read name
    BC
    """

    print("Starting to parse barcode fasta file")
    op_d: Dict[str, str] = {}
    FH = open(barcode_fasta_fp)
    c_line = FH.readline()
    nL = 1
    while c_line != "":
        if c_line[0] == ">":
            read_name = c_line.rstrip().split(" ")[0][1:]
        else:
            raise Exception("Expecting read line at line # " + str(nL))
        bc = FH.readline().rstrip()
        nL += 1
        if read_name in op_d:
            raise RuntimeError(
                "Was not expecting repeat read name for "
                + "barcode: read name "
                + read_name
                + ", line: "
                + str(nL)
            )
        op_d[read_name] = bc
        c_line = FH.readline()
        nL += 1

    FH.close()
    return op_d


def main():
    args = sys.argv
    if args[-1] == "1":
        _, op_lib_dir, comb_df_fp, cfg_fp, bc_to_loc_fp = args[:-1]
        lib_name = op_lib_dir.split("/")[-1]
        cfg_d = load_entire_cfg(cfg_fp)
        with open(bc_to_loc_fp, "r") as f:
            bc_to_loc_dicts = json.loads(f.read())
            bc_to_loc_dicts = {
                x: {
                    tuple(
                        json.loads(
                            k.replace("(", "[").replace(")", "]").replace("'", '"')
                        )
                    ): v
                    for k, v in y.items()
                }
                for x, y in bc_to_loc_dicts.items()
            }
        midway_run1(op_lib_dir, lib_name, comb_df_fp, cfg_d, bc_to_loc_dicts)
    else:
        print("python3 src/step5.py op_lib_dir, comb_df_fp, cfg_fp, 1")
        raise Exception("Could not recognize run style flag")


if __name__ == "__main__":
    main()
