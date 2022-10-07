"""
In this file we create PDFs of plots
"""
import os
import sys
import logging
import matplotlib.pyplot as plt

plt.switch_backend("agg")
import shutil


def create_plots_dir(lib_op_dir, force_create=True):

    plots_dir = os.path.join(lib_op_dir, "Plots")
    if os.path.exists(plots_dir):
        if force_create:
            print("Removing Plots directory at " + plots_dir)
            shutil.rmtree(plots_dir)
        else:
            print(
                f"Plots directory already exists at {plots_dir}"
                + ". Not doing anything."
            )
            return plots_dir
    os.mkdir(plots_dir)

    return plots_dir


def create_read_to_ins_histogram(
    lib_inserts_uniq_fp, plot_op_fp, table_op_fp, lib_name
):
    """
    Description:
        This plots the # of reads per insert seq against occurrence to visualize sequencing depth

        The input file is the output of vsearch from step3, ins_uniq.fna
        At the end of each read name (line starting with ">")
        there is ";size=X" where X is a positive integer that
        counts the number of times the following sequence appeared
        in the original file.
        This function tallies the number of sequence duplicates,
        so for example, if 7 sequences appeared 5 times, our output
        histogram would have the y-value 7 for the x-value 5.
    """

    # List of all the sizes (repeat values occur, e.g.
    # [3,3,3,5,5,9,9,9,9,9,2], except we ignore "1"
    sizes_arr = []
    one_count = 0
    # Below two are kept to create histogram
    min_num = float("inf")
    max_num = float("-inf")
    error_str = "Failed to create histogram, line #"
    print("Beginning to gather counts for histogram")
    FH = open(lib_inserts_uniq_fp, "r")
    c_line = FH.readline()
    line_num = 1
    while c_line != "":
        if c_line[0] != ">":
            raise Exception(
                "Expecting line to have '>' as first character."
                + error_str
                + str(line_num)
            )
        last_part = c_line.rstrip().split(";")[-1]
        if last_part[:4] != "size":
            raise Exception(
                "Expecting read name to end with size=X, "
                + "instead "
                + last_part
                + ". "
                + error_str
                + str(line_num)
            )
        pos_int = last_part.split("=")[-1]
        try:
            nReplicates = int(pos_int)
        except Exception as inst:
            logging.error("Failed to parse integer at line " + str(line_num))
            raise Exception(error_str + str(line_num))
        if nReplicates != 1:
            sizes_arr.append(nReplicates)
        else:
            one_count += 1

        # Below we skip all the non-read name lines
        c_line = FH.readline()
        line_num += 1
        while c_line != "" and c_line[0] != ">":
            c_line = FH.readline()
            line_num += 1
    FH.close()

    write_counts_table(sizes_arr, one_count, table_op_fp)

    # Creating plot
    tot_reads_more_than_1 = sum(sizes_arr)
    tot_inserts_more_than_1 = len(sizes_arr)
    min_count = min(sizes_arr)
    max_count = max(sizes_arr)
    bins = prepare_bins(min_count, max_count)
    try:
        plt.hist(sizes_arr, bins=bins)
        plt.title("Sequencing depth, " + lib_name)
        plt.xlabel("# of reads per insert sequence")
        plt.ylabel("Occurrence")
        plt.figtext(
            0.4, 0.7, f"Tot reads of inserts w/ >1 read = {tot_reads_more_than_1}"
        )
        plt.figtext(0.4, 0.8, f"Inserts w/ >1 read = {tot_inserts_more_than_1}")
        print("Writing file out to " + plot_op_fp)
        plt.savefig(plot_op_fp)
        plt.clf()
    except Exception as e:
        print("Failed to create histogram: ", str(e))
    return None


def prepare_bins(min_num, max_num):
    if min_num == float("inf") or max_num == float("-inf"):
        return None
    else:
        return list(range(min_num, max_num + 2))


def write_counts_table(counts_list, one_count, op_fp):
    count_d = create_count_dict(counts_list)
    count_d[1] = one_count
    keys = sorted(list(count_d.keys()))
    with open(op_fp, "w") as g:
        for k in keys:
            g.write(f"{k}\t{count_d[k]}\n")
    print("Wrote read counts per insert sequence and frequency to " + op_fp)


def create_count_dict(counts_list):
    count_d = {}
    for c in counts_list:
        if c in count_d:
            count_d[c] += 1
        else:
            count_d[c] = 1
    return count_d


def test_libuniq_hst(inp_fp, op_fp):
    logging.basicConfig(level=logging.DEBUG)

    create_lib_inserts_histogram(inp_fp, op_fp)


def main():
    args = sys.argv
    help_str = "python3 src/prepare_plots.py lib_inserts_fp op_plots_fp 1"
    if args[-1] == "1":
        py_fn, lib_uniq_inserts_fp, op_fp, indicator = args
        test_libuniq_hst(lib_uniq_inserts_fp, op_fp)
    else:
        print(help_str)

    sys.exit(0)


if __name__ == "__main__":
    main()
