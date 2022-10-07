"""
    Before this step you can assume the files are split by library into small enough sizes (each file < 2.5 GB) for the free version of usearch.

    Within this step we run the following two usearch commands: 'search_pcr2' and 
                                                                'search_oligodb'
    For fastq file $f

    usearch -search_pcr2 $f\.fq -fwdprimer $fwd -revprimer $rev \
            -maxdiffs 2 -minamp 300 -maxamp 15000 -strand both \
            -fastqout $f\_trimmed1.fq -notmatchedfq $f\_noprimer.fq

    With the above usearch command we 'trim' the forward and reverse primers.

    Then we run:
    
        usearch -threads 20 -search_oligodb $f\_trimmed1.fq -db oligos.fa -strand plus \
                -maxdiffs 1 -userout $f\_ins_bc.txt \
                -notmatchedfq discarded_reads/$f\_no_insbc.fq \
                -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue

    With the above usearch command (search_oligodb) we search for short sequences flanking barcode and insert regions.
    Note. usearch positions are 0-indexed.

    Then we combine these output files into a single file per library and place them in subdir 01:
        01-us_ops/{lib}_trimmed.fq - Reads where primers are found and trimmed from search_pcr2 (Used later on)
        01-us_ops/{lib}_noprimer.fq - Reads with no primer (for debug/log)
        01-us_ops/{lib}_insbc.txt - Main output of search_oligodb (Used later on)
        01-us_ops/{lib}_no_insbc.fq - Reads where flanking regions were not found (for debug/log)
"""

import os
import sys
import logging
import json
import subprocess
import shutil
import time
from validate import verify_cfg_d
from typing import List, Dict, Tuple

def run_step_1_singlelib(op_lib_dir: str, lib_name: str, cfg_d: Dict, fq_fps: List[str]) -> None:
    """
    Args:
        op_lib_dir is a directory main_output/lib_name
        lib_name (str) Name of the library
        cfg_d keys must include:
            "step_1", whose inner dict must include
               "usearch_exec_path": (str) Path to executable file
               "search_pcr2": dict must include
                    "fwd": str (sequence)
                    "rev": str (sequence)
                    "maxdiffs": positive_int
                    "minamp": positive_int
                    "maxamp": positive_int
                "oligo_db_fp" (str): file path of oligo database file
        fq_fps (list<str>): FQ fps associated w this lib_name
    """
    print("\nRunning step 1 for lib " + lib_name)
    print("Will run usearch on the following files:\n" + ", ".join(fq_fps))

    # Set up log dir
    log_list = []
    logs_dir = os.path.join(os.path.dirname(op_lib_dir), "Logs")
    if not os.path.exists(logs_dir):
        os.mkdir(logs_dir)
        print(f"Made Logs directory at {logs_dir}")

    # write output dir for step 1
    usearch_op_dir: str = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["1"])
    if os.path.exists(usearch_op_dir):
        shutil.rmtree(usearch_op_dir, ignore_errors=True)
    os.mkdir(usearch_op_dir)
    print(f"Made output dir for step 1 at {usearch_op_dir}")

    # f_basenames is either lib or lib_X for files already split
    f_basenames = [extract_file_name(f) for f in fq_fps]

    usearch_pcr2_results, new_log = run_many_usearch_pcr2(
        fq_fps, f_basenames, cfg_d, usearch_op_dir
    )
    log_list += ["usearch -search_pcr2 step:"]
    log_list += new_log

    oligo_db_fp = cfg_d["primer_info"]["oligo_db_fp"]
    log_oligodb = run_many_usearch_search_oligodbs(
        usearch_pcr2_results, cfg_d, usearch_op_dir, oligo_db_fp
    )
    log_list += ["\nusearch -search_oligodb step:"]
    log_list += log_oligodb

    ## Concatenate output files if there were multiple inputs
    # skip this if there was a single input file <-> single output file
    if len(f_basenames) > 1:
        rmv_bool = cfg_d["step_1"]["remove_non_concatenated_oligo_ops"]
        concat_oligo_outputs(
            usearch_op_dir, f_basenames, lib_name, cfg_d, remove_old=rmv_bool
        )
        concatenate_usearch_pcr2_results(usearch_op_dir, lib_name, cfg_d)
    else:
        print("Single input fastq file, skipping concatenation step.")

    # Get number of reads in output after concatenation
    ls_fp = os.listdir(usearch_op_dir)
    log_list += ["\n"]
    for file in ls_fp:
        out_fp = os.path.join(usearch_op_dir, file)
        if any(x in file for x in ["fq", "fastq"]):
            num_reads = get_file_length(out_fp) / 4
        else:
            num_reads = get_file_length(out_fp)
        log_list += [f"Reads in {out_fp}: " + str(num_reads)]

    with open(
        os.path.join(logs_dir, lib_name + "_step1_log.txt"),
        "w",
    ) as g:
        g.write("\n".join(log_list))

    print("Finished step 1.")

    return None


def run_many_usearch_pcr2(files_to_process, f_basenames, cfg_d, usearch_op_dir) -> Tuple[Dict, List[str]]:
    """
    Args:
        files_to_process list(fp)
        f_basenames list(str)
    usearch -search_pcr2 $f\.fq -fwdprimer $fwd -revprimer $rev \
    -maxdiffs 2 -minamp 300 -maxamp 15000 -strand both \
    -fastqout $f\_trimmed1.fq -notmatchedfq discarded_reads/$f\_noprimer.fq

    Description:
        usearch -search_pcr2, https://www.drive5.com/usearch/manual/cmd_search_pcr2.html.
        'f_basenames' is a list which has a 1-1 correspondence with
        'files_to_process', in which the value in f_basenames is filename w/o extension 
        '.fq' or '.fastq' from files_to_process.

    Returns:
        usearch_pcr2_results (dict): {fq (str) -> {X}} Where X is
        ret_d = {
                "fq_out": fq_out, (STRING)
                "discarded_out": not_matched_out, (STRING)
                "stderr": res.stderr.decode('utf-8'), (STRING)
                "stdout": res.stdout.decode('utf-8') (STRING)
        }
        usearch_pcr2_op_dir (str): Path to usearch results directory
    """

    s1 = cfg_d["step_1"]
    usearch_exec_path = s1["usearch_exec_path"]
    arg_d = s1["search_pcr2"]
    usearch_pcr2_results = {}

    for i in range(len(files_to_process)):
        fq_fp = files_to_process[i]
        print(
            f"Running usearch search_pcr2 on file {fq_fp}, which is #{i + 1}/{len(files_to_process)}.\n"
        )
        f_basename = f_basenames[i]
        result, log_list = run_usearch_search_pcr2_command(
            fq_fp, f_basename, usearch_op_dir, usearch_exec_path, arg_d, cfg_d
        )
        usearch_pcr2_results[fq_fp] = result

    print("Finished running usearch search_pcr2, now on to running usearch oligodb")

    return (usearch_pcr2_results, log_list)


def run_usearch_search_pcr2_command(
    fq_fp: str, f_basename: str, op_lib_dir: str, usearch_exec_path, arg_d, cfg_d
) -> Tuple[Dict, List[str]]:
    """
    Args:
        op_lib_dir (str): Path to directory to write
                      'tabbed_out' and 'fastq_out'
        usearch_exec_path (str): Path to 'usearch' executable
        arg_d: (Comes from config file 'step1' - search_pcr2)
            "fwd": str (sequence)
            "rev": str (sequence)
            "maxdiffs": positive_int
            "minamp": positive_int
            "maxamp": positive_int
            "strand": "both", or "plus"

        Example:
            "fwd": "GTTCTTATCTTTGCAGTCTC",
            "rev": "GAGATTTACGCTTTGGTAAAAGTTGG",
            "maxdiffs": 2,
            "minamp": 300,
            "maxamp": 15000

    Description:
        usearch search_pcr2 looks at every read and takes the reads that contain
        the forward and reverse primer and at least a length of 'minamp' and
        at most a length of 'maxamp' and outputs it to the file listed at -fastqout

    Returns:
        ret_d = {
                "fq_out": fq_out, (STRING)
                "discarded_out": not_matched_out, (STRING)
                "stderr": res.stderr.decode('utf-8'), (STRING)
                "stdout": res.stdout.decode('utf-8') (STRING)
        }
    """

    fq_out: str = os.path.join(
        op_lib_dir, f_basename + cfg_d["d"]["fns"]["1"]["trimmed"]
    )
    not_matched_out: str = os.path.join(
        op_lib_dir, f_basename + cfg_d["d"]["fns"]["1"]["no_primer"]
    )

    command_args = [
        usearch_exec_path,
        "-search_pcr2",
        fq_fp,
        "-fwdprimer",
        arg_d["fwd"],
        "-revprimer",
        arg_d["rev"],
        "-maxdiffs",
        str(arg_d["maxdiffs"]),
        "-minamp",
        str(arg_d["minamp"]),
        "-maxamp",
        str(arg_d["maxamp"]),
        "-strand",
        "both",
        "-fastqout",
        fq_out,
        "-notmatchedfq",
        not_matched_out,
    ]
    # The output from the file is stored in res.stdout
    res = subprocess.run(command_args, capture_output=True)
    log_list: List[str] = ["stdout: " + res.stdout.decode("utf-8"),
            "stderr: " + res.stderr.decode("utf-8")
            ]
    ret_d = {"fq_out": fq_out}
    return ret_d, log_list


def get_file_length(fp: str) -> int:
    cmds = ["wc", "-l", fp]
    result = subprocess.run(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    x = result.stdout.decode("utf-8")
    fl = int(x.split(" ")[0])
    return fl


def run_many_usearch_search_oligodbs(
    usearch_pcr2_results, cfg_d, usearch_op_dir, oligo_db_fp
) -> List[str]:
    """
    Description:
        On every single file outputted from usearch search_pcr2 (which may be
        subsets of the input fastq files), we run usearch oligodb
    Returns:
       Dictionary containing the usearch oligodb running output
       for every file
    """

    print("Preparing to run usearch_oligodb on all output files from before.")
    us_oligodb_f2p = []
    us_oligodb_basenames = []
    for k, v in usearch_pcr2_results.items():
        us_oligodb_f2p.append(v["fq_out"])
        ## First, gets lib_X_trimmed as file basename, then exclude "_trimmed" for oligodb outputs
        trim = cfg_d["d"]["fns"]["1"]["trimmed"].split(".")[
            0
        ]  # something like '_trimmed'
        file_basename = extract_file_name(v["fq_out"]).split(trim)[0]
        us_oligodb_basenames.append(file_basename)

    log_list: List[str] = []
    for i in range(len(us_oligodb_f2p)):
        fq_fp = us_oligodb_f2p[i]
        f_basename = us_oligodb_basenames[i]
        log_list += run_usearch_search_oligodb_command(
            fq_fp,
            oligo_db_fp,
            cfg_d["step_1"]["usearch_exec_path"],
            cfg_d["step_1"]["search_oligodb"],
            f_basename,
            usearch_op_dir,
            cfg_d,
        )

    return log_list


def run_usearch_search_oligodb_command(
    fq_fp,
    oligo_db_fp,
    usearch_exec_path,
    arg_d,
    f_basename,
    usearch_op_dir,
    cfg_d,
    nThreads=None,
) -> List[str]:
    """
    Description: 
    Note that the input 'fq_fp' is the output of usearch -search_pcr2, and not the 
    original input fastq file. Reads from -search_pcr2 output are oriented based on fwd and rev primers.
    usearch -threads 20 -search_oligodb $f\_trimmed1.fq -db ~/pacbio/oligos_blunt-end.fa -strand plus \
            -maxdiffs 1 -userout $f\_ins_bc.txt \
            -notmatchedfq discarded/$f\_no_insbc.fq \
            -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue

    Writes all output to 'usearch_op_dir'.

    Returns:
        dictionary containing usearch results 
    """

    command_args: List[str] = [usearch_exec_path]
    if nThreads is not None:
        command_args += ["-threads", str(nThreads)]
    fp_out = os.path.join(usearch_op_dir, f_basename + cfg_d["d"]["fns"]["1"]["insbc"])
    no_match_out = os.path.join(
        usearch_op_dir, f_basename + cfg_d["d"]["fns"]["1"]["no_insbc"]
    )

    command_args += [
        "-search_oligodb",
        fq_fp,
        "-db",
        oligo_db_fp,
        "-strand",
        "plus",
        "-maxdiffs",
        str(arg_d["maxdiffs"]),
        "-userout",
        fp_out,
        "-notmatchedfq",
        no_match_out,
        "-userfields",
        "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue",
    ]

    print(
        "Running usearch search_oligodb with the following commands: "
        + ", ".join(command_args)
    )
    res = subprocess.run(command_args, capture_output=True)
    print(res)
    print("Finished running usearch with the above commands.")

    tmp_log_list: List[str] = ["stdout: " + res.stdout.decode("utf-8"),
                "stderr: " + res.stderr.decode("utf-8")
            ]

    return tmp_log_list


def concatenate_usearch_pcr2_results(usearch_op_dir: str, lib, cfg_d, remove_old=True) -> str:
    """
    Description:
        We combine multiple fastq files into a single fastq file
    Args:
        usearch_op_dir (str) : Path to directory with
                               outputs from usearch pcr2
                               command.
        lib_names list<str> names of libs.
    """

    op_d: str = usearch_op_dir
    all_lib_fs = os.listdir(op_d)

    print("Concatenating usearch pcr2 '.fq' files for library " + lib)

    lib_ins = [
        os.path.join(op_d, x)
        for x in all_lib_fs
        if x.endswith(cfg_d["d"]["fns"]["1"]["trimmed"])
    ]

    print(cfg_d["d"]["fns"]["1"]["trimmed"] + " files:\n" + ", ".join(lib_ins))

    concat_files_to_output(
        lib_ins, os.path.join(op_d, lib + cfg_d["d"]["fns"]["1"]["trimmed"])
    )

    lib_no_primers = [
        os.path.join(op_d, x)
        for x in all_lib_fs
        if x.endswith(cfg_d["d"]["fns"]["1"]["no_primer"])
    ]
    print(cfg_d["d"]["fns"]["1"]["no_primer"] + " files:\n" + ", ".join(lib_no_primers))
    concat_files_to_output(
        lib_no_primers, os.path.join(op_d, lib + cfg_d["d"]["fns"]["1"]["no_primer"])
    )

    if remove_old:
        print(f"Removing segmented outputs for lib {lib}")
        for f in lib_no_primers + lib_ins:
            os.unlink(f)

    print("Wrote concatenated files to " + op_d)

    return op_d


def concat_oligo_outputs(op_dir, f_basenames, lib, cfg_d, remove_old=False) -> None:
    """
    Args:
        op_dir: The directory in which the outputs of the usearch
                command search_oligodb are placed.
        lib_names (list<str>): The names of the original FASTQ files
                               without the ".fq" extension. If files
                               were previously split, then it takes
                               the split_files original name

        remove_old (bool): A flag that decides whether we remove
                           the old outputs after we concatenate them
                           and place the concatenated in a directory.
    """
    all_lib_fs = os.listdir(op_dir)

    print("Concatenating oligo '.txt' files for library: " + lib)

    lib_ins = [
        os.path.join(op_dir, x)
        for x in all_lib_fs
        if x.endswith(cfg_d["d"]["fns"]["1"]["insbc"])
    ]
    print(cfg_d["d"]["fns"]["1"]["insbc"] + " files:\n" + ", ".join(lib_ins))

    concat_files_to_output(
        lib_ins, os.path.join(op_dir, lib + cfg_d["d"]["fns"]["1"]["insbc"])
    )

    lib_no_ins = [
        os.path.join(op_dir, x)
        for x in all_lib_fs
        if x.endswith(cfg_d["d"]["fns"]["1"]["no_insbc"])
    ]
    print(cfg_d["d"]["fns"]["1"]["no_insbc"] + " files:\n" + ", ".join(lib_no_ins))
    concat_files_to_output(
        lib_no_ins, os.path.join(op_dir, lib + cfg_d["d"]["fns"]["1"]["no_insbc"])
    )

    if remove_old:
        print(f"Removing segmented outputs for lib {lib}")
        for f in lib_no_ins + lib_ins:
            os.unlink(f)

    print("Wrote concatenated files to " + op_dir)

    return None


def concat_files_to_output(inp_fs: List[str], op_f: str) -> None:
    """
    Concatenates files into a single output
    inp_fs List[str]: input files to concatenate.
    Note: If length of inp_fs is 0, nothing is written
    """
    ## Becomes an infinite loop if inp_fs has a value = op_f!
    if op_f in inp_fs:
        line_out = (
            "Output filename "
            + op_f
            + " already exists in subdir 01. Check redundancy in input fq."
        )
        raise Exception(line_out)
    else:
        with open(op_f, "wb") as wfd:
            for f in inp_fs:
                with open(f, "rb") as fd:
                    shutil.copyfileobj(fd, wfd)

    return None


def extract_file_name(full_path: str, already_split=False) -> str:
    """
    Extracts file name before extension .fq or .fastq,
    e.g. BC112.fq -> BC112
    If the files are already split, then value becomes 'lib_X'
    where X is the split file number.
    If already_split = True, returns lib instead of lib_X.
    """

    file_basename = os.path.basename(full_path)
    split_file_name = file_basename.split(".")
    if len(split_file_name) > 2:
        raise Exception(
            "Each fastq file in the input directory should only "
            "have 1 period within it, e.g. 'ABC.fq', not 'A.BC.fq'."
            " File looks like: " + file_basename
        )
    elif split_file_name[-1] not in ["fq", "fastq"]:
        raise Exception(
            "Each fastq file in the input directory should end "
            "with '.fq' or '.fastq'. Instead"
            " file named: " + file_basename
        )

    if split_file_name[0] == "Logs":
        raise Exception("Library name cannot be logs. f: " + full_path)
    elif split_file_name[0][0].isdigit():
        raise Exception(
            "Library name cannot start with a " + "numerical digit: " + full_path
        )
    if already_split == True:
        split_file_name = split_file_name.rsplit("_", 1)[0]

    return split_file_name[0]


def main():

    help_str = "python3 src/step1.py cfg_json inp_dir op_dir(tmp) 1"
    help_str = "OR\n"
    help_str = "python3 src/step1.py inp_dir oligos_dir 2"
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
