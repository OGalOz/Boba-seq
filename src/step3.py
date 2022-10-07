"""
Within this step we apply the vsearch command -fastq filter
    on _ins.fq files from step 2 to filter reads by quality.
    Specifically, we filter by a max expected error of 10.

e.g.:
    for f in lib1, lib2; do
    vsearch -fastq_filter $f\_ins.fq -fastq_qmin 0 -fastq_qmax 93 
            -fastq_maxee 10 -fastaout $f\_insert.fasta
    done &
    fastq_qmax is set to 93 for long reads
    
We generate two output files to subdir 03:
    lib_insert.fasta -> Filtered insert sequences from insert fastq file input.
    lib_ins_uniq.fna -> A dereplicated version of the above fasta file.
                        The read names are changed to contain only the first term.
                        This information is used to generate histogram / table to
                        visualize sequencing depth of each insert sequence.
                        This .fna file is removed after plot generation.
We generate two output files to subdir Plots:
    .pdf of histogram to visualize sequencing depth
    .tsv of table of counts used to plot histogram
    
"""

import sys
import os
import subprocess
import shutil
import logging
import json
from prepare_plots import create_plots_dir, create_read_to_ins_histogram
from typing import List, Dict, Tuple, Union


def run_step_3_singlelib(op_lib_dir, lib_name, cfg_d) -> None:
    """
    Args:
        op_lib_dir: library output dir
        cfg_d must have key
        'vsearch' executable
    """

    log_list = []
    print("\nRunning step 3 for lib " + lib_name)

    # input config parameters and file path
    vs_exc = cfg_d["vsearch_exec_path"]
    max_ee = cfg_d["step_3"]["max_expected_error"]
    inp_fp = get_inp_fp(op_lib_dir, lib_name, cfg_d)

    # vs_dir is op_lib_dir/03-vs_ffee
    vs_dir = make_op_dir(op_lib_dir, lib_name, cfg_d)
    op_fp_ins = os.path.join(vs_dir, lib_name + cfg_d["d"]["fns"]["3"]["inserts"])
    op_fp_derep = os.path.join(vs_dir, lib_name + "_uniq_inserts.fna")

    # running vsearch commands
    log_list += ["vsearch -fastq_filter\n"]
    log_list += run_vsearch_filter_command(vs_exc, inp_fp, op_fp_ins, max_ee)
    log_list += ["\n-----------------------------------"]
    log_list += ["\nvsearch --derep_fulllength\n"]
    log_list += run_vsearch_dereplicate_command(vs_exc, op_fp_ins, op_fp_derep)

    # We generate histogram and table to Plots/
    plots_dir = create_plots_dir(op_lib_dir)
    plot_fp = os.path.join(
        plots_dir, lib_name + cfg_d["d"]["fns"]["P"]["uniq_ins_hist"]
    )
    counts_table_fp = os.path.join(
        plots_dir, lib_name + cfg_d["d"]["fns"]["P"]["uniq_ins_tbl"]
    )
    create_read_to_ins_histogram(op_fp_derep, plot_fp, counts_table_fp, lib_name)
    # remove unique inserts fasta
    os.unlink(op_fp_derep)

    # Write to Logs/
    log_dir = get_log_dir(op_lib_dir)
    with open(
        os.path.join(os.path.dirname(op_lib_dir), "Logs", lib_name + "_step3_log.txt"),
        "w",
    ) as g:
        g.write("\n".join(log_list))

    print("Finished step 3.")

    return None


def make_op_dir(op_lib_dir, lib_name, cfg_d) -> str:
    vs_op_dir = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["3"])
    if os.path.exists(vs_op_dir):
        shutil.rmtree(vs_op_dir, ignore_errors=True)
    os.mkdir(vs_op_dir)
    print("Made vsearch output dir at " + vs_op_dir)
    return vs_op_dir


def get_inp_fp(op_lib_dir, lib_name, cfg_d) -> str:
    ext_dir = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["2"])
    expected_fp = os.path.join(ext_dir, lib_name + cfg_d["d"]["fns"]["2"]["ins_fq"])
    if not os.path.exists(expected_fp):
        raise Exception(f"Could not find expected file at {expected_fp}.")
    else:
        return expected_fp


def run_vsearch_filter_command(vs_exc, inp_fp, op_fp, max_ee) -> List[str]:
    """
    vsearch -fastq_filter $f\_ins.fq -fastq_qmin 0 -fastq_qmax 93 -fastq_maxee 10
        -fastaout $f\_insert.fasta
    """

    command_args = [
        vs_exc,
        "-fastq_filter",
        inp_fp,
        "-fastq_qmin",
        "0",
        "-fastq_qmax",
        "93",
        "-fastq_maxee",
        str(max_ee),
        "-fastaout",
        op_fp,
    ]
    print(
        "Running vsearch fastq_filter with the following commands: "
        + ", ".join(command_args)
    )
    res = subprocess.run(command_args, capture_output=True)
    print(res)
    print("Finished running vsearch with the above results.")
    full_stderr: str = res.stderr.decode("utf-8")
    n_remaining, n_discarded = get_important_nums(full_stderr, "filter")
    tot_reads = n_remaining + n_discarded
    perc_kept = round(n_remaining / tot_reads * 100, 2)
    log_list = [
        "stderr: " + full_stderr,
        "n_remaining: " + str(n_remaining),
        "n_discarded: " + str(n_discarded),
        "perc_kept: " + str(perc_kept),
    ]
    return log_list


def run_vsearch_dereplicate_command(vs_exc, inp_fp, op_fp) -> List[str]:
    """
    vsearch --derep_fulllength inp_fp --output op_fp --sizeout
    """

    command_args = [
        vs_exc,
        "--derep_fulllength",
        inp_fp,
        "--output",
        op_fp,
        "--sizeout",
    ]
    print(
        "Running vsearch dereplicate with the following commands: "
        + ", ".join(command_args)
    )
    res = subprocess.run(command_args, capture_output=True)
    res_str = res.stderr.decode("utf-8")
    print(res_str)
    print("Finished running vsearch with the above results.")

    # We store the number n_unique in case we want it in the future - 
    # it isn't currently in use.
    n_unique = get_important_nums(res_str, "derep")

    log_list = ["stderr: " + res_str]

    return log_list


def get_important_nums(inp: str, typ: str) -> Union[int, Tuple[int,int]]:
    # in this function we parse the standard error outputs of vsearch filter/derep
    # so typ is one of 'filter' or 'derep'
    # inp is standard error string
    try:
        res = get_imp1(inp, typ)
    except RuntimeError:
        print("Failed to parse vsearch output using method 1.", "Trying second method.")
        res = get_imp2(inp, typ)
    return res


def get_imp2(inp: str, typ: str) -> Union[int, Tuple[int,int]]:
    # in this function we parse the standard error outputs of vsearch filter/derep
    # so typ is one of 'filter' or 'derep'
    # inp is standard error string
    if typ == "filter":
        lines = inp.split("\n")
        ix = 0
        while ix < len(lines):
            test_line = lines[ix]
            if test_line.startswith("Reading"):
                line = lines[ix + 1]
                words = line.split(" ")
                n1Str, n2Str = words[0], words[7]
                return int(n1Str), int(n2Str)
            ix += 1
        raise RuntimeError("Could not find line 'Reading' in vsearch output.")
    elif typ == "derep":
        lines = inp.split("\n")
        ix = 0
        while ix < len(lines):
            test_line = lines[ix]
            if test_line.startswith("Sorting"):
                line = lines[ix + 1]
                words = line.split(" ")
                return int(words[0])
            ix += 1
        raise RuntimeError("Could not find line 'Sorting' in vsearch output.")


def get_imp1(inp: str, typ: str) -> Union[int, Tuple[int,int]]:
    # in this function we parse the standard error outputs of vsearch filter/derep
    # so typ is one of 'filter' or 'derep'
    # inp is standard error string
    if typ == "filter":
        lines = inp.split("\n")
        test_line = lines[3]
        if not test_line.startswith("Reading"):
            raise RuntimeError("Didn't expect this output: " + ",\n".join(lines))
        line = lines[4]
        words = line.split(" ")
        n1Str, n2Str = words[0], words[7]
        return int(n1Str), int(n2Str)
    elif typ == "derep":
        lines = inp.split("\n")
        test_line = lines[5]
        if not test_line.startswith("Sorting"):
            raise RuntimeError("Didn't expect this output: " + ",\n".join(lines))
        line = lines[6]
        words = line.split(" ")
        return int(words[0])
    else:
        raise RuntimeError(
            f"Could not recognize typ: `{typ}`. Should be one of 'filter' or 'derep'."
        )


def get_log_dir(op_lib_dir) -> str:
    return os.path.join(os.path.dirname(op_lib_dir), "Logs")
