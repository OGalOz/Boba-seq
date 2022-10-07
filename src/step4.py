"""
Within this step, we run Minimap2 to map inserts to reference genome
    Executes the following:
    minimap2 -x map-hifi {inp_gnm} {insert_fp} --paf-no-hit > {op_fp} &

    where 
        inp_gnm = FROM CONFIG lib_genome_dir/lib_genome.fasta
        insert_fp = op_lib_dir/03-vs_ffee/{lib_name}_insert.fasta
        op_fp = op_lib_dir/04-minimap2/{lib_name}.paf

Then it computes a summary report on .paf output for Logs/
    liness_in_paf: # of lines in .paf
    reads_in_paf: # of reads in .paf
    reads_not_mapped: # of reads with no hit
    reads_mapped: # of reads with hits
    hits_high_qual: # of lines after quality filtering
    reads_high_qual: # of reads after quality filtering
    reads_high_qual/reads_mapped, perc: % reads with hits that passed filter
    perc_reads_w_one_hit: similar to above, but with single hit
    perc_reads_w_more_hits: similar to above, but with multiple hits, wrote to file
"""

import os
import sys
import subprocess
import shutil
import logging
import json
import pandas as pd
from collections import Counter
from typing import Dict, List, Tuple
from get_bc_to_ins_file import get_read_to_seq_dict_from_fa
from parse_paf import parse_paf_file


def run_step_4_singlelib(op_lib_dir, lib_name, cfg_d) -> None:
    """
    Args:
        op_lib_dir: Central library dir
            Should contain subdirectories '02-pos_ins_bc', '01-us_ogdb', etc.

    cfg_d must have mapping of lib names to related reference genomes
    """
    print("\nRunning step 4 for lib " + lib_name)
    step4_log_d, minimap_dir = run_minimap_executable(op_lib_dir, lib_name, cfg_d)

    new_stats_d: Dict = compute_other_paf_stats(op_lib_dir, lib_name, cfg_d)

    step4_log_d.update(new_stats_d)
    log_fp = os.path.join(
        os.path.dirname(op_lib_dir), "Logs", lib_name + "_step4_log.json"
    )
    with open(log_fp, "w") as g:
        g.write(json.dumps(step4_log_d, indent=2))
    print("Wrote log to " + log_fp)

    print("Finished step 4.")


def run_minimap_executable(op_lib_dir, lib_name, cfg_d) -> Tuple[Dict, str]:
    insert_fp = get_step3_output(op_lib_dir, lib_name, cfg_d)
    minimap_exec = cfg_d["minimap2_exec_path"]
    genome_fp = get_genome_fp(lib_name, cfg_d)
    # minimap_dir: op_lib_dir/minimap
    minimap_dir: str = create_minimap_dir(op_lib_dir, cfg_d)
    op_fp = os.path.join(minimap_dir, lib_name + cfg_d["d"]["fns"]["4"]["minimap_op"])
    # log_d: stderr, file length of .paf output, percent hit
    log_d: Dict = run_minimap2(minimap_exec, genome_fp, insert_fp, op_fp)

    return log_d, minimap_dir


def run_minimap2(minimap_exec, inp_gnm, insert_fp, op_fp) -> Dict:
    """
    minimap2 -x map-hifi {inp_gnm} {insert_fp} --paf-no-hit > {op_fp} &
    """
    # Two options for the command: asm20/map-hifi
    command_args = [minimap_exec]
    command_args += ["-x", "map-hifi", inp_gnm, insert_fp, "-paf-no-hit"]

    print("Running minimap2 with the following command: " + ", ".join(command_args))
    with open(op_fp, "w") as outfile:
        res = subprocess.run(command_args, stdout=outfile, stderr=subprocess.PIPE)
    print(
        "Finished running minimap2 with the above command." + " Wrote to file " + op_fp
    )
    res_str = res.stderr.decode("utf-8")
    paf_length = get_file_length(op_fp)  # # of rows in .paf output
    log_d = {"stderr": res_str, "lines_in_paf": paf_length}

    return log_d


def get_file_length(fp: str) -> int:
    cmds = ["wc", "-l", fp]
    result = subprocess.run(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    x = result.stdout.decode("utf-8")
    fl = int(x.split(" ")[0])
    return fl


def compute_other_paf_stats(op_lib_dir, lib_name, cfg_d) -> Dict:
    """
    Computs stats for log.
    Uses parse_paf_file to convert .paf to df, remove reads w/ no hits,
    filters by quality, remove reads w/ multiple hits into separate output.
    """
    # read output fasta from step3 and output .paf from Minimap2
    insert_fp = get_step3_output(op_lib_dir, lib_name, cfg_d)
    minimap_dir = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["4"])
    op_fp = os.path.join(minimap_dir, lib_name + cfg_d["d"]["fns"]["4"]["minimap_op"])

    # convert inserts fasta to dictionary of read name:sequence, extract a list of read names
    seq_name_to_seq_d: Dict = get_read_to_seq_dict_from_fa(insert_fp)
    fasta_read_names: List[str] = list(seq_name_to_seq_d.keys())
    fasta_read_names = [x.split(" ")[0] for x in fasta_read_names]
    del seq_name_to_seq_d

    # Filter paf file to generate summary stats for Logs/
    # Filtered paf table is kept and used in step 5
    paf_df, paf_multi_hits, paf_info = parse_paf_file(op_fp, cfg_d)

    # Write reads with multiple hits to file in subdir 04
    op_fp = os.path.join(minimap_dir, lib_name + "_paf_multi_hits.tsv")
    print(f"Wrote reads with multiple hits to {op_fp}")
    paf_multi_hits.to_csv(op_fp, index=False, sep="\t")

    return paf_info


def get_step3_output(op_lib_dir, lib_name, cfg_d) -> str:
    # insert_fp is something like op_lib_dir/03-vs_ffee/{lib_name}_insert.fasta
    vsearch_op_dir = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["3"])
    if not os.path.exists(vsearch_op_dir):
        raise Exception(
            "Could not find previous vsearch output directory" + " at " + vsearch_op_dir
        )
    insert_fp: str = os.path.join(
        vsearch_op_dir, lib_name + cfg_d["d"]["fns"]["3"]["inserts"]
    )
    if not os.path.exists(insert_fp):
        raise Exception(
            "Could not find previous vsearch output 'insert' fasta" + " at " + insert_fp
        )
    return insert_fp


def create_minimap_dir(op_lib_dir, cfg_d, force_create=True) -> str:
    minimap_dir = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["4"])
    if os.path.exists(minimap_dir):
        if force_create:
            print("Removing minimap directory at " + minimap_dir)
            shutil.rmtree(minimap_dir)
        else:
            print(
                f"minimap directory already exists at {minimap_dir}"
                + ". Not doing anything."
            )
            return minimap_dir
    os.mkdir(minimap_dir)

    return minimap_dir


def get_genome_fp(lib_name, cfg_d) -> str:
    """
    Description:
        Get path to genome fp using config dict

    Examples of values for the following keys:
        lib_names, lib_genome_dir, lib_genome_filenames

        "lib_names": ["BC1127", "BC1124", "BC1130"],
        "lib_genome_dir": "/ref_genomes",
        "lib_genome_filenames": ["Bvulgatus_ATCC8482.fasta",
                                "Bovatus_ATCC8483.fasta", "Bstercoris_CC31F.fasta"],
    """
    lib_names = cfg_d["lib_names"]
    lib_genome_dir = cfg_d["lib_genome_dir"]
    lib_genome_filenames = cfg_d["lib_genome_filenames"]
    lib_ix = lib_names.index(lib_name)
    gnm_fn = lib_genome_filenames[lib_ix]
    gnm_fp = os.path.join(lib_genome_dir, gnm_fn)

    return gnm_fp


def main():
    args = sys.argv
    logging.basicConfig(level=logging.DEBUG)
    help_str = "python3 src/step4.py lib_dir lib_name cfg_fp 1\n"
    help_str += "OR\n"
    help_str += "python3 src/step4.py lib_dir lib_name cfg_fp 2\n"
    if args[-1] == "1":
        py_fn, lib_dir, lib_name, cfg_fp, indicator = args
        with open(cfg_fp, "r") as f:
            cfg_d = json.loads(f.read())
        run_step_4_singlelib(lib_dir, lib_name, cfg_d)
    elif args[-1] == "2":
        # checking statistics
        py_fn, lib_dir, lib_name, cfg_fp, indicator = args
        with open(cfg_fp, "r") as f:
            cfg_d = json.loads(f.read())
    else:
        print(help_str)

    sys.exit(0)


if __name__ == "__main__":
    main()
