import os
import sys
import shutil
import logging
import json
from typing import List


def initialize_output_directory(op_dir, interactive=False):
    """
    We clear the output directory if it already exists (remov
    Some directories we need to add to tmp dir:
        'Logs'
    """
    if os.path.exists(op_dir):
        if interactive:
            print(
                f"Output directory already exists at {op_dir}."
                + " The program will delete the directory and"
                + " recreate it unless indicated otherwise."
                + " The process will continue in 5 seconds."
                + " Press Ctrl+C to escape the process."
            )
            time.sleep(5)
        print(f"Removing {op_dir}.")
        shutil.rmtree(op_dir, ignore_errors=True)
        if os.path.exists(op_dir):
            raise Exception(
                "Failed to clear output directory at {op_dir}. "
                + "Please remove the directory and rerun program."
            )
    print(f"Creating output directory at {op_dir}")
    os.mkdir(op_dir)

    logs_dir = os.path.join(op_dir, "Logs")
    os.mkdir(logs_dir)
    return logs_dir


def force_create_dir(op_dir):
    if os.path.exists(op_dir):
        shutil.rmtree(op_dir, ignore_errors=True)
    os.mkdir(op_dir)


def prepare_for_steps(op_dir, step_num, lib, cfg_d=None):

    if step_num == 2:
        lib_dir = os.path.join(op_dir, lib)
        oligo_concat_dir = os.path.join(lib_dir, "us_ogdb", "concat")
        return [oligo_concat_dir, lib_dir]
    elif step_num == 3:
        lib_dir = os.path.join(op_dir, lib)
        fq_dir = os.path.join(lib_dir, "us_spcr", "concat")
        fq_fp = get_trimmed_fq_fp_from_fq_dir(fq_dir)
        pos_op_dir = os.path.join(lib_dir, "pos_ins_bc")
        return [fq_fp, pos_op_dir]
    elif step_num == 4:
        return None

    return None


def get_trimmed_fq_fp_from_fq_dir(fq_dir):
    fs = os.listdir(fq_dir)
    fq_fp = None
    for f in fs:
        if f.endswith("trimmed1.fq"):
            fq_fp = os.path.join(fq_dir, f)
    if fq_fp is None:
        raise Exception("Could not find target file X_trimmed1.fq")
    else:
        return fq_fp


def load_json(fp):
    with open(fp, "r") as f:
        x = json.loads(f.read())
    return x


def replace_all_instances_of_x_with_y_dir(
    dir_path: str, to_replace: str, to_insert: str, ignore_files: List[str], debug=True
) -> None:

    fs = os.listdir(dir_path)
    for x in ignore_files:
        fs.remove(x)
    fs = [f for f in fs if not os.path.isdir(os.path.join(dir_path, f))]
    print("REMAINING FS: ", fs)
    for f in fs:
        replace_all_instances_of_x_with_y_file(
            os.path.join(dir_path, f), to_replace, to_insert, debug=debug
        )

    return None


def replace_all_instances_of_x_with_y_file(
    fp: str, to_replace: str, to_insert: str, debug=True
):
    swap_file = fp + ".swp"
    OPFH = open(swap_file, "w")
    with open(fp, "r") as f:
        for line in f:
            new_line = line.replace(to_replace, to_insert)
            OPFH.write(new_line)

    OPFH.close()

    shutil.move(fp, fp + ".old")
    shutil.move(swap_file, fp)
    if debug:
        print(f"Replaced instances of `{to_replace}` with `{to_insert}` in file {fp}.")


def remove_old(p: str):
    for x in os.listdir(p):
        if x.endswith(".old"):
            shutil.move(os.path.join(p, x), os.path.join(p, x[:-4]))


if __name__ == "__main__":
    args = sys.argv
    if args[-1] != "33":
        print("Wrong way to run the file.")
        sys.exit(0)
    else:
        _, df_cfg_fp, src_dir, tkn = args
        with open(df_cfg_fp, "r") as f:
            dflt = json.loads(f.read())
        d1 = dflt["fns"]
        ignore = []
        current_dir = os.getcwd()
        src_dir = os.path.join(current_dir, src_dir)
        for x in ["1", "2", "3", "4", "5", "P"]:
            print("Running for " + x)
            new_d = d1[x]
            for k, v in new_d.items():
                to_replace = f'"{v}"'
                to_insert = f"cfg_d['d']['fns']['{x}']['{k}']"
                replace_all_instances_of_x_with_y_dir(
                    src_dir, to_replace, to_insert, ignore, debug=True
                )
