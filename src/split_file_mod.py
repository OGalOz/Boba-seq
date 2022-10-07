"""
In this file we get each library name,
    and split the input files into
    sizes of less than 2.5 GB (or whatever
    number is set).
    We return a dictionary with each
    library name pointing to the file(s)
    which make up the input.
Only the input dir is involved.
"""
import logging, os, sys
import shutil


def split_files_and_get_lib2fq_fps(
    inp_dir,
    cfg_d,
    to_split_files=True,
    debug=False,
    interactive=False,
    MAX_FILE_SIZE=2.5 * (2**30),
):
    """
    If files already split, the directory should be in good shape,
    and we should be pointed at directory with split files.
    "to_split_files" means we will be splitting files
    """

    lib_names = cfg_d["lib_names"]
    print("Given lib names: " + ", ".join(lib_names))

    input_files = os.listdir(inp_dir)
    print("Will work on the following files: " + ", ".join(input_files))

    if to_split_files:
        if "split_files" in input_files:
            print(("*" * 30 + "\n") * 4)
            if interactive:
                res = input(
                    "Looks like 'split_files' directory is already in input_dir, "
                    + "do you want to use split_files directory as input and avoid "
                    + "splitting files again? y/n"
                )
                if res.upper() not in ["Y", "N"]:
                    raise Exception("Result not one of 'y' or 'n', exiting program.")
                elif res.upper() == "N":
                    print("Please remove split_files directory and try again.")
                    sys.exit(0)

            sf_dir = os.path.join(inp_dir, "split_files")
            to_split_files = False

    if to_split_files:
        # We will split the files
        sf_dir = os.path.join(inp_dir, "split_files")
        print("Creating split file directory at " + sf_dir)
        os.mkdir(sf_dir)
        fs2prcs = []
        for f in input_files:
            print("Beginning work on splitting file " + f + " to small enough sizes.")
            fs2prcs += split_file(
                os.path.join(inp_dir, f), sf_dir, max_size=MAX_FILE_SIZE
            )
        f_basenames = [check_file_name(f) for f in fs2prcs]
    else:
        sf_dir = inp_dir

    # The files are already split
    fs2prcs = check_split_files(sf_dir, MAX_FILE_SIZE)
    lib2fs = get_lib2fs(fs2prcs, lib_names)
    return lib2fs


def check_split_files(sf_dir, MAX_FILE_SIZE):
    # sf_dir: split files dir

    fs = os.listdir(sf_dir)
    fs2prcs = [os.path.join(sf_dir, x) for x in fs]
    for f in fs2prcs:
        filesize = os.path.getsize(f)
        if filesize > MAX_FILE_SIZE:
            raise Exception(
                f"File {f} too large, "
                + f"are you sure the dir {sf_dir} "
                + "contains only previously split files?"
            )
    return fs2prcs


def get_lib2fs(fs2prcs, lib_names):
    lib2fs = {}

    # Sort lib names in reverse length order
    srtd_l = sorted(lib_names, reverse=True, key=lambda l: len(l))

    # We make a copy
    new_fs = fs2prcs[:]
    for lib in srtd_l:
        fs = [f for f in new_fs if validate_f_name(f, lib)]
        lib2fs[lib] = fs
        for f in fs:
            new_fs.remove(f)

    return lib2fs


def validate_f_name(f, lib):
    basename = os.path.basename(f)
    if basename.startswith(lib):
        if "." in basename:
            if basename.split(".")[-1] in ["fq", "fastq"]:
                return True
    return False


def split_file(inp_file, op_dir, max_size=2.5 * (2**30), remove_orig=False):
    """
    max_size (int): Max number of bytes in a file, 2**30 is a GB,
                    2.5*(2**30) would be 2.5 GB
    remove_orig (bool): Should we remove the input file after splitting it
    """
    # Checking file name:
    file_base = check_file_name(inp_file)

    # file_size should be bytes (int)
    file_size = os.path.getsize(inp_file)
    if file_size < max_size:
        print(
            f"Not splitting file {inp_file} since it's size is under"
            + f"the file size max limit ({file_size}<{max_size}). "
            + f"But the file will be copied to {op_dir}/{file_base}_1.fq"
        )
        file_op_loc = os.path.join(op_dir, file_base + "_1.fq")
        shutil.copy(inp_file, file_op_loc)
        return [file_op_loc]
    else:
        print(
            f"Splitting file {inp_file} since it's size is over the file size max limit."
        )

    FH = open(inp_file)

    files_to_process = []

    new_file_list = []
    crt_new_file_size = 0
    total_bytes_processed = 0
    read_name_line = FH.readline()
    file_num = 1
    line_num = 1
    while read_name_line != "":
        additional_lines_size = len(read_name_line)
        read_seq = FH.readline()
        additional_lines_size += len(read_seq)
        separator = FH.readline()
        if separator[0] != "+":
            raise Exception(
                f"Fastq file '{inp_file}' doesn't have "
                "separator '+' at expected loc at line " + str(line_num + 2)
            )
        else:
            additional_lines_size += len(separator)
        quality = FH.readline()
        additional_lines_size += len(quality)

        if crt_new_file_size + additional_lines_size > max_size:
            new_file_name = file_base + "_" + str(file_num) + ".fq"
            out_fp = os.path.join(op_dir, new_file_name)
            print("Writing out split file at " + out_fp)
            # If file size passes 1 GB, we print out sub_file
            with open(out_fp, "w") as g:
                g.write("".join(new_file_list))
            files_to_process.append(out_fp)
            new_file_list = []
            file_num += 1
            crt_new_file_size = 0

        crt_new_file_size += additional_lines_size
        total_bytes_processed += additional_lines_size
        new_file_list += [read_name_line, read_seq, separator, quality]
        line_num += 4
        read_name_line = FH.readline()
        if (line_num - 1) % 10**5 == 0:
            print(
                "Just processed line "
                + str(line_num)
                + "."
                + f"\nCompleted {round(((total_bytes_processed/file_size) * 100),3)}% "
                + "of file."
            )

    FH.close()

    if file_num > 1:
        if remove_orig:
            os.unlink(inp_file)
        new_file_name = file_base + "_" + str(file_num) + ".fq"
        out_fp = os.path.join(op_dir, new_file_name)
        print(
            "Completed 100% of file '"
            + file_base
            + "'. Writing out final split file at "
            + out_fp
        )
        # If file size passes max GB, we print out sub_file
        with open(out_fp, "w") as g:
            g.write("".join(new_file_list))
    else:
        files_to_process = [inp_file]

    return files_to_process


def check_file_name(full_path, already_split=False):
    """
    Normally gives the file name before the period,
    e.g. BC112.fq -> BC112
    If the files are already split, then we take the library
    names by getting rid of the last '_X' where X is the split
    file number.
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
            " file looks like: " + file_basename
        )

    if split_file_name[0] == "Logs":
        raise Exception("Library name cannot be logs. f: " + full_path)
    elif split_file_name[0][0].isdigit():
        raise Exception(
            "Library name cannot start with a " + "numerical digit: " + full_path
        )

    return split_file_name[0]
