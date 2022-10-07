import os
import logging
import sys
import json
import pandas as pd
import matplotlib.pyplot as plt
import importlib
from util_mod import initialize_output_directory, load_json

from step1 import run_step_1_singlelib
from step2 import run_step_2_singlelib
from step3 import run_step_3_singlelib
from step4 import run_step_4_singlelib
from step5 import run_step_5_singlelib
from step6 import run_step_6_singlelib

from validate import load_entire_cfg, verify_cfg_d
from split_file_mod import split_files_and_get_lib2fq_fps
import traceback
import datetime


def run_all_steps(cfg_fp, inp_dir, op_dir):

    cfg_d = load_entire_cfg(cfg_fp)

    to_split_files = not cfg_d["step_1"]["files_pre_split"]
    # The following part is not run in parallel for Memory saving purposes
    lib2fq_fps = split_files_and_get_lib2fq_fps(
        inp_dir, cfg_d, to_split_files=to_split_files, debug=False
    )
    lib_names = cfg_d["lib_names"]
    logs_dir = initialize_output_directory(op_dir)

    res_str = sequential_runs(lib_names, lib2fq_fps, cfg_d, op_dir, logs_dir)
    return res_str


def sequential_runs(lib_names, lib2fq_fps, cfg_d, op_dir, logs_dir):
    """
    For each library in the list lib_names, we perform a
    complete run through.
    Args:
        lib_names: list<str>
        lib2fq_fps: dict mapping lib_name (str) -> list[fq_fp]
        cfg_d: main config dict
        op_dir (str): Path to output directory containing all sub outputs
        logs_dir (str): Path to output directory containing logs; op_dir/Logs
    """
    res_str = "successfully"
    # Sequential runs:
    for lib in lib_names:
        clean_between_libs()
        print("\nStarting sequential run on library: " + lib)
        try:
            complete_run(lib, lib2fq_fps, cfg_d, op_dir, logs_dir)
        except Exception as inst:
            print(f"Failed to run on library {lib}")
            print("Error message: " + str(inst))
            traceback.print_tb(inst.__traceback__)
            print(str(inst))
            res_str = "un" + res_str
    return res_str


def complete_run(lib, lib2fq_fps, cfg_d, op_dir, logs_dir):
    """
    This function clears out old run within the output directory,
    then starts a new run all the way from step 1.
    Args:
        lib (str) Name of library
        lib2fq_fps (d): Maps lib name to list of filepaths which make up lib
        cfg_d (d):
        op_dir (str): Path to output (tmp) directory
        logs_dir (str): Path to logs dir within op_dir (op_dir/Logs)
    """
    lib_op_dir = os.path.join(op_dir, lib)
    if os.path.exists(lib_op_dir):
        if os.path.isdir(lib_op_dir):
            shutil.rmtree(lib_op_dir)
        else:
            os.unlink(lib_op_dir)
    os.mkdir(lib_op_dir)
    fq_fps = lib2fq_fps[lib]
    run_step_1_singlelib(lib_op_dir, lib, cfg_d, fq_fps)
    run_step_2_singlelib(lib_op_dir, lib, cfg_d, interactive=False)
    run_step_3_singlelib(lib_op_dir, lib, cfg_d)
    # We only run the following code if we are not running on a metagenome
    if "metagenome_bool" not in cfg_d or not cfg_d["metagenome_bool"]:
        run_step_4_singlelib(lib_op_dir, lib, cfg_d)
        run_step_5_singlelib(lib_op_dir, lib, cfg_d)
        run_step_6_singlelib(lib_op_dir, lib, cfg_d)
    return None


def run_mid(cfg_fp, op_dir, step_num) -> str:
    """
    Desc:
        Run from the middle of a section (using previous step directories)
    Args:
        cfg_fp (str): Path to config file
        op_dir (str): Path to output directory containing the lib output
                        directories within it
        step_num (int): Represents which step we are in
    """
    if step_num < 2:
        raise Exception("Step must be 2 or greater")
    if not os.path.exists(op_dir):
        raise RuntimeError(
            f"Output directory not found at {op_dir},\n" + "try running from the start."
        )
    cfg_d = load_json(cfg_fp)
    verify_cfg_d(cfg_d)
    print("Succesfully validated config and inputs.")
    lib_names = cfg_d["lib_names"]
    ret_str = ""
    for lib in lib_names:
        clean_between_libs()
        lib_op_dir = os.path.join(op_dir, lib)
        try:
            if step_num <= 2:
                run_step_2_singlelib(lib_op_dir, lib, cfg_d, interactive=False)
            if step_num <= 3:
                run_step_3_singlelib(lib_op_dir, lib, cfg_d)
            if step_num <= 4:
                run_step_4_singlelib(lib_op_dir, lib, cfg_d)
            if step_num <= 5:
                run_step_5_singlelib(lib_op_dir, lib, cfg_d)
            if step_num <= 6:
                run_step_6_singlelib(lib_op_dir, lib, cfg_d)
        except Exception as inst:
            print("Failed to run for library " + lib)
            print("Error message:")
            print(inst)
            traceback.print_tb(inst.__traceback__)
            print("time now: " + str(datetime.datetime.now()))
            ret_str += f" Lib {lib} ran with error.\n"
            return "with error"
        ret_str += f"Lib {lib} ran succesfully.\n"

    return ret_str


def clean_between_libs():
    """
    Certain errors occur when you run multiple libs in a row
    """
    importlib.reload(pd)
    importlib.reload(plt)

    pass


def get_help_str():

    help_str = "To run the entire program from the beginning, use:\n"
    help_str += "python3 src/run_steps.py cfg_json inp_dir op_dir(tmp) 1\n"
    help_str += "To start the program at a certain point, use:\n"
    help_str += "python3 src/run_steps.py cfg_json inp_dir op_dir(tmp) [int=2,3,...]\n"
    return help_str


def main():
    start_time_str: str = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print(("*" * 10 + "\n") * 5 + f"NEW RUN {start_time_str}" + ("*" * 10 + "\n") * 5)
    args = sys.argv
    help_str = get_help_str()
    if args[-1] not in [str(x) for x in range(1, 7)]:
        print(help_str)
        sys.exit(1)
    else:
        cfg_fp = args[1]
        inp_dir = args[2]
        op_dir = args[3]
        res_str = "default"
        if args[-1] == "1":
            res_str = run_all_steps(cfg_fp, inp_dir, op_dir)
        elif args[-1] == "2":
            res_str = run_mid(cfg_fp, op_dir, 2)
        elif args[-1] == "3":
            res_str = run_mid(cfg_fp, op_dir, 3)
        elif args[-1] == "4":
            res_str = run_mid(cfg_fp, op_dir, 4)
        elif args[-1] == "5":
            res_str = run_mid(cfg_fp, op_dir, 5)
        elif args[-1] == "6":
            res_str = run_mid(cfg_fp, op_dir, 6)
        else:
            raise Exception("Last arg not parsed properly?: " + args[-1])

    print(
        ("*" * 10 + "\n") * 5
        + f"Finished run {res_str}. \nStarted at time {start_time_str}"
        + f"\n Ended at time {str(datetime.datetime.now())}"
        + ("*" * 10 + "\n") * 5
    )

    return None


if __name__ == "__main__":
    main()
