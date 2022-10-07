"""
Within this step, we generate plots to visualize fragment lengths and genome coverage.

Note within histograms you can control the bin locations as such:
e.g. plt.hist(data, bins=[0, 4, 8, 12, 16, 20])

Columns for bc_mapped_df:

    bc, tstart, tend, orig_contig, multi_frag_x, contig, multi_frag_y, 
    strand, gene_count, locus_tag, directionality
"""
import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging

logging.getLogger("matplotlib").setLevel(logging.WARNING)
import sys
import statistics
from typing import List, Dict, Tuple
from collections import Counter
from get_bc_to_ins_file import get_read_to_seq_dict_from_fa
from contig_collider import special_match_contig_names


def run_step_6_singlelib(
    op_lib_dir, lib_name, cfg_d, plots_dirname="Plots", step_num=6
) -> None:
    """
    Within Step 6, we primarily create plots.

    We import bc_loc_df from step5 and exclude barcodes with multiple mappings.
    The fraction of barcodes mapped to multiple positions should be very low.
    This dataset is used for fragment size, fragment coverage, and fragment distribution plots.

    We import genes_count_df from step5, which includes barcodes with multiple mappings.
    This dataset is used for all gene coverage plots and gene coverage statistics.

    Specifically the following files are generated:
        loc_count_table.tsv
        _cumulative_gene_cov.pdf
        _fragment_cov.pdf
        _gene_cov.pdf
        _fragment_len.pdf
        _CONTIGNAME_map.pdf
        _gene_coverage.txt
    """
    print("\nRunning step 6 for lib " + lib_name)
    if plots_dirname == "Plots":
        plots_dir = os.path.join(op_lib_dir, "Plots")
    else:
        plots_dir = os.path.join(op_lib_dir, "Plots")
        plots_dir = os.path.join(plots_dir, plots_dirname)
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)
    print(f"All plots will be saved at {plots_dir}")

    # file paths and import datasets
    inp_dfs_dir = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["5"])
    bc_map_df_fp = os.path.join(inp_dfs_dir, cfg_d["d"]["fns"]["5"]["bc_loc_df"])
    redundant_bc_fp = os.path.join(inp_dfs_dir, cfg_d["d"]["fns"]["5"]["redundant_bc"])
    genes_df_fp = os.path.join(inp_dfs_dir, cfg_d["d"]["fns"]["5"]["genes_count_df"])

    bc_mapped_df = pd.read_table(bc_map_df_fp)

    # exclude barcodes mapped to multiple positions for plots and statistics
    with open(redundant_bc_fp, "r") as f:
        bc_multi_dict = json.load(f)
    ls_bc_multi = bc_multi_dict["bc_w_multi_mappings"]
    bc_mapped_noMulti_df = bc_mapped_df[~bc_mapped_df["bc"].isin(ls_bc_multi)]
    print(bc_mapped_noMulti_df.head(5))

    genes_df = pd.read_table(genes_df_fp)
    nGenes = genes_df.shape[0]
    log_list = ["Total number of genes in gff: " + str(nGenes)]

    fragment_size_distribution_plot(
        bc_mapped_noMulti_df, plots_dir, cfg_d, lib_name, log_list
    )
    gene_coverage_plot(genes_df, plots_dir, cfg_d, lib_name)
    fragment_coverage_plot(bc_mapped_noMulti_df, plots_dir, cfg_d, lib_name)
    cumulative_gene_coverage_plot(genes_df, plots_dir, cfg_d, lib_name)

    # These contigs come from genome.fna file
    genome_contig2len: Dict[str,int] = get_contig_to_len(cfg_d, lib_name)
    # These contigs come from 'gff'
    contigs: List[str] = list(bc_mapped_noMulti_df["contig"].unique())
    contig2len = reconcile_contig2len(genome_contig2len, contigs)
    create_all_contigs_fragment_maps(
        bc_mapped_noMulti_df, plots_dir, cfg_d, lib_name, contigs, contig2len
    )
    get_frag_cov_per_window(
        bc_mapped_noMulti_df, plots_dir, cfg_d, lib_name, contigs, contig2len
    )

    ## Gene stats calculations for log file
    n_genes_not_hit = get_genes_not_hit(genes_df)
    n_genes_covered = get_genes_covered(genes_df)
    n_genes_opp_only = get_genes_opp_dir_only(genes_df)
    n_genes_same_only = get_genes_same_dir_only(genes_df)
    n_genes_opp_tot = get_counts_opp_dir(genes_df)
    n_genes_same_tot = get_counts_same_dir(genes_df)
    perc_genes_not_hit = round(n_genes_not_hit / nGenes * 100, 2)
    perc_genes_opp_only = round(n_genes_opp_only / nGenes * 100, 2)
    perc_genes_same_only = round(n_genes_same_only / nGenes * 100, 2)
    perc_genes_opp_tot = round(
        n_genes_opp_tot / (n_genes_opp_tot + n_genes_same_tot) * 100, 2
    )
    perc_genes_same_tot = round(
        n_genes_same_tot / (n_genes_opp_tot + n_genes_same_tot) * 100, 2
    )
    ratio_opp_to_same = round(n_genes_opp_tot / n_genes_same_tot, 2)

    log_list += [
        str(nGenes - n_genes_not_hit)
        + " genes covered, "
        + str(100 - (perc_genes_not_hit))
        + "%."
    ]
    log_list += [
        str(n_genes_not_hit) + " genes not covered, " + str(perc_genes_not_hit) + "%."
    ]
    log_list += [
        str(n_genes_opp_only)
        + " genes always oriented opposite to promoter, "
        + str(perc_genes_opp_only)
        + "%."
    ]
    log_list += [
        str(n_genes_same_only)
        + " genes always same strand as promoter, "
        + str(perc_genes_same_only)
        + "%."
    ]
    log_list += [
        str(n_genes_opp_tot)
        + " occurrence of gene oriented opposite to promoter, "
        + str(perc_genes_opp_tot)
        + "%."
    ]
    log_list += [
        str(n_genes_same_tot)
        + " occurrence of gene on same strand as promoter, "
        + str(perc_genes_same_tot)
        + "%."
    ]
    log_list += [
        "Ratio of genes oriented opposite / same as promoter: " + str(ratio_opp_to_same)
    ]

    ls_frag_per_gene = genes_df["fragment_count"].to_list()
    log_list += ["Fragment count per gene: "]
    log_list += ["Median: " + str(statistics.median(ls_frag_per_gene))]
    log_list += [
        "Range: " + str(min(ls_frag_per_gene)) + " to " + str(max(ls_frag_per_gene))
    ]

    bc_gene_count_uniq = bc_mapped_noMulti_df.drop_duplicates(
        subset=["bc", "gene_count"]
    )
    ls_gene_per_frag = bc_gene_count_uniq["gene_count"].to_list()
    log_list += ["Gene count per fragment: "]
    log_list += ["Median: " + str(statistics.median(ls_gene_per_frag))]
    log_list += [
        "Range: " + str(min(ls_gene_per_frag)) + " to " + str(max(ls_gene_per_frag))
    ]

    g_cov_log_fp = os.path.join(
        os.path.dirname(op_lib_dir), "Logs", lib_name + "_lib_stats.txt"
    )
    with open(g_cov_log_fp, "w") as g:
        g.write("\n".join(log_list))
        print("Wrote genome coverage stats to " + g_cov_log_fp)

    print(f"Finished step 6.")


def get_frag_cov_per_window(
    bc_map_df, plots_dir, cfg_d, lib_name, contigs, contig2len, window_size=20000
) -> None:
    """
    Desc:
        Within this function, we create windows of 20KB with the median fragment coverage
        per base pair, this is done per contig.
    """
    print("Binning genome into windows of ", window_size, "bp")
    if window_size % 1000 != 0:
        raise Exception(
            "Bin size for genome must be divisible by 1000, instead: "
            + str(window_size)
        )
    cntg_2_windows = {}
    for ctg in contigs:
        print("Binning contig ", ctg)
        contig_len = contig2len[ctg]

        if contig_len == 10**6 or contig_len < 6*(10**4):
            continue
        nFullWindows = contig_len // window_size
        last_window = contig_len % window_size
        bp_d = {x: 0 for x in range(contig_len)}
        # all fragments in contig
        contig_df = bc_map_df[bc_map_df["contig"] == ctg]
        contig_df.sort_values(by="tstart", inplace=True, ignore_index=True)
        for i in range(contig_df.shape[0]):
            start, end = contig_df.at[i, "tstart"], contig_df.at[i, "tend"]
            for j in range(start, end):
                bp_d[j] += 1  # counter for fragment coverage at each position
        window2median: Dict[Tuple[int, int], float] = {}

        plt.title(f"Fragment coverage in {window_size / 1000}kb windows, {lib_name}")
        plt.xlabel("Position in contig " + ctg + ", kb")
        plt.ylabel("Median fragment coverage per base")
        plt.grid(visible=True)
        plt.xlim(0, contig_len / 1000.0)
        max_y = 0
        for i in range(nFullWindows):
            rng = (i * window_size, (i + 1) * window_size)  # start, end
            halfway_point = (rng[0] + ((rng[1] - rng[0]) / 2)) / 1000.0
            ls_cov_in_window = [bp_d[j] for j in range(rng[0], rng[1])]
            median = statistics.median(ls_cov_in_window)
            max_y = max(max_y, median)
            xs = [halfway_point, halfway_point]
            ys = [0, median]
            plt.plot(xs, ys, color="grey")
            # Storing data:
            window2median[rng] = median
        plt.ylim = (0, max_y + 1)
        # last window is <20kb in length, only plot if >5kb
        last_rng = (nFullWindows * window_size, contig_len)
        last_win_len = last_rng[1] - last_rng[0]
        if last_win_len > 5000:
            ls_cov_in_window = [bp_d[j] for j in range(last_rng[0], last_rng[1])]
            median = statistics.median(ls_cov_in_window)
            max_y = max(max_y, median)
            xs = [last_rng[0], last_rng[1]]
            ys = [0, median]
            plt.plot(xs, ys, color="grey")
            window2median[last_rng] = median

        op_fp = os.path.join(
            plots_dir, lib_name + "_" + ctg + cfg_d["d"]["fns"]["P"]["contig_frag_cov"]
        )
        plt.savefig(op_fp)
        plt.clf()
        print("Wrote fragment coverage per genomic windows plot to " + op_fp + ".")
        # Here we store the results
        cntg_2_windows[ctg] = window2median


def get_contig_to_len(cfg_d, lib_name) -> Dict[str, int]:
    # Returns contig name mapped to length
    ix = cfg_d["lib_names"].index(lib_name)
    genomes_dir = cfg_d["lib_genome_dir"]
    genome_fna_file = os.path.join(genomes_dir, cfg_d["lib_genome_filenames"][ix])
    contig2seq = get_read_to_seq_dict_from_fa(genome_fna_file)
    contig2len = {}
    for ctg in contig2seq:
        contig2len[ctg.split(" ")[0]] = len(contig2seq[ctg])
    return contig2len


def get_genes_not_hit(genes_df) -> int:
    x = genes_df[genes_df["fragment_count"] == 0]
    return x.shape[0]


def get_genes_covered(genes_df) -> int:
    x = genes_df[genes_df["fragment_count"] != 0]
    return x.shape[0]


def get_genes_opp_dir_only(genes_df) -> int:
    x = genes_df[genes_df["nSame_direc_frag"] == 0]
    return x.shape[0]


def get_genes_same_dir_only(genes_df) -> int:
    x = genes_df[genes_df["nSame_direc_frag"] == genes_df["fragment_count"]]
    return x.shape[0]


def get_counts_opp_dir(genes_df) -> int:
    x = genes_df["fragment_count"].sum() - genes_df["nSame_direc_frag"].sum()
    return x


def get_counts_same_dir(genes_df) -> int:
    x = genes_df["nSame_direc_frag"].sum()
    return x


def fragment_size_distribution_plot(
    bc_mapped_df: pd.DataFrame, plots_dir: str, cfg_d, lib_name, log_list, bin_width=250
) -> None:
    f_lens = bc_mapped_df["frag_len"].values
    mx = max(f_lens)
    bins = [i * bin_width for i in range(mx // bin_width)] + [mx]
    plt.hist(f_lens, bins=bins, color="teal")
    plt.title("Distribution of fragment size, " + lib_name)
    plt.xlabel("Fragment size (bp)")
    plt.grid(visible=False)
    op_fp = os.path.join(plots_dir, lib_name + cfg_d["d"]["fns"]["P"]["frag_len"])
    plt.savefig(op_fp)
    plt.clf()
    print(
        "Wrote fragment size histogram to "
        + lib_name
        + cfg_d["d"]["fns"]["P"]["frag_len"]
    )

    # output stats
    log_list += ["Median fragment length: " + str(statistics.median(f_lens))]
    log_list += ["Min fragment length: " + str(min(f_lens))]
    log_list += ["Max fragment length: " + str(max(f_lens))]

    return log_list


def gene_coverage_plot(
    genes_df: pd.DataFrame, plots_dir: str, cfg_d, lib_name, bin_width=10
) -> None:
    fragment_counts = genes_df["fragment_count"]
    mx = max(fragment_counts)
    bins = [i * bin_width for i in range(mx // bin_width)] + [mx]
    plt.hist(fragment_counts, bins=bins, color="teal")
    plt.title("Gene coverage, " + lib_name)
    plt.xlabel("Number of fragments")
    plt.ylabel("Number of genes")
    plt.grid(visible=False)
    op_fp = os.path.join(plots_dir, lib_name + cfg_d["d"]["fns"]["P"]["gene_cov"])
    plt.savefig(op_fp)
    plt.clf()
    print(
        "Wrote gene coverage histogram to "
        + lib_name
        + cfg_d["d"]["fns"]["P"]["gene_cov"]
    )
    return None


def fragment_coverage_plot(
    bc_map_df: pd.DataFrame, plots_dir: str, cfg_d, lib_name, bin_width=1
) -> None:
    gcs = bc_map_df["gene_count"]
    mx = max(gcs)
    bins = [i * bin_width for i in range(mx // bin_width)] + [mx]
    plt.hist(gcs, bins=bins, color="teal")
    plt.title("Fragment coverage, " + lib_name)
    plt.xlabel("Number of genes")
    plt.ylabel("Number of fragments")
    plt.grid(visible=False)
    op_fp = os.path.join(plots_dir, lib_name + cfg_d["d"]["fns"]["P"]["frag_cov"])
    plt.savefig(op_fp)
    plt.clf()
    print(
        "Wrote fragment coverage histogram to "
        + lib_name
        + cfg_d["d"]["fns"]["P"]["frag_cov"]
    )
    return None


def cumulative_gene_coverage_plot(
    genes_df: pd.DataFrame, plots_dir: str, cfg_d, lib_name
) -> None:
    x = np.arange(0, genes_df.fragment_count.max() + 1, 1)
    y = np.array([0] * len(x))
    for gi, gene in genes_df.iterrows():
        y[gene.fragment_count] += 1
    for yi in range(len(y) - 2, -1, -1):
        y[yi] += y[yi + 1]
    y = y / genes_df.shape[0]
    plt.plot(x, y, color="teal")
    plt.grid(True)
    plt.title("Cumulative gene coverage, " + lib_name)
    plt.ylabel("Fraction of genes covered by X fragments")
    plt.xlabel("Number of fragments")
    op_fp = os.path.join(plots_dir, lib_name + cfg_d["d"]["fns"]["P"]["cu_gene_cov"])
    plt.savefig(op_fp)
    plt.clf()
    print(
        "Wrote cumulative gene coverage line plot to "
        + lib_name
        + cfg_d["d"]["fns"]["P"]["cu_gene_cov"]
    )
    return None


def create_all_contigs_fragment_maps(
    bc_map_df: pd.DataFrame, plots_dir: str, cfg_d, lib_name, contigs, contig2len: Dict
) -> None:
    """
    For each contig in the bc_map_df (column "contig"),
    we visualize distribution of insert fragments.
    We number each fragment for a sub dataframe (reset index), then the y
    axis will have 1 value per fragment.
    We sort them by tstart (starting location on the contig),
    then the amount of x axis they receive is relative to the length
    of the contig (which is 'tlen')
    """

    plt.style.use("seaborn-whitegrid")
    print("Total number of fragments:", bc_map_df.shape[0])
    print("Creating contig fragment map")
    print("All contigs:", contigs)
    size_sum = 0
    for contig in contigs:
        print("running on contig " + contig)
        contig_df = bc_map_df[bc_map_df["contig"] == contig]
        contig_df.sort_values(by="tstart", inplace=True, ignore_index=True)
        op_fp = os.path.join(
            plots_dir,
            lib_name
            + "_"
            + contig_df.at[0, "contig"]
            + cfg_d["d"]["fns"]["P"]["contig_frag_map"],
        )
        create_limited_contig_fragment_plot(
            contig_df, op_fp, lib_name, contig2len[contig]
        )
        print("Created contig fragment plot for contig " + contig)
        size_sum += contig_df.shape[0]


def create_limited_contig_fragment_plot(
    contig_df: pd.DataFrame,
    output_fp: str,
    lib_name,
    contig_length: int,
    fragment_start: int = 0,
    fragment_end: int = float("inf"),
    contig_start: int = 0,
    contig_end: int = float("inf"),
    contig_limit: bool = False,
    per_run: int = 10000,
) -> None:
    """
    Cols for contig_df:
    bc, tstart, tend, orig_contig, multi_frag_x, contig, multi_frag_y,
    strand, gene_count, locus_tag, directionality
    Desc:
        We create a single plot for each contig.
        There are optionally limitations on which fragment # you start on,
        and which locations within the contig the fragment starts between.
        TD: Error is with the ending location being frag_len of first one??
    Args:
        contig_df is a DataFrame which only holds the fragments that map
            to this contig. Column names used:
            "frag_len", "tstart", "tend", "contig",

        output_fp (str) Path to file where the PDF is saved
        fragment_start: Within this dataframe, choose the first fragment
                        to map
        fragment_end: Within this dataframe, choose the last fragment to map
        contig_start: Withint the contig, choose the lowest possible starting
                      base pair number.
        contig_end: Withint the contig, choose the highest possible starting
                      base pair number.
        contig_limit (bool): Whether we apply the contig_start and contig_end
                            numbers. If this is set to False we just ignore
                            those numbers entirely.
        per_run (int): Unimportant integer which defines how often we report
                        the number of fragments plotted thus far.
    """
    print(contig_df.columns)
    if fragment_end == float("inf"):
        fragment_end = contig_df.shape[0]
    if contig_end == float("inf"):
        contig_end = contig_length
    new_df = contig_df[fragment_start:fragment_end]
    # We make sure the fragment start is between the two given locations.
    if contig_limit:
        new_df = new_df[new_df["tstart"] >= contig_start][
            new_df["tstart"] <= contig_end
        ]

    nFrag: int = int(new_df.shape[0])
    contig_len: int = contig_end - contig_start

    plt.ylabel("Fragment Number")
    contig_name = new_df.at[0, "contig"]
    plt.xlabel(contig_name + ", " + lib_name)
    plt.xlim(contig_start, contig_end)
    # The following line can cause an error - restart kernel (?)
    plt.ylim(0, nFrag)

    nRuns = nFrag // per_run
    print("nRuns:", nRuns)
    max_start = max(new_df["tstart"].to_list())
    max_end = max(new_df["tend"].to_list())
    print("max start", max_start, "max end", max_end)

    for i in range(nRuns):
        for j in range(i * per_run, (i + 1) * per_run):
            # k = fragment_start + j
            k = j
            xs = [int(new_df.at[k, "tstart"]), int(new_df.at[k, "tend"])]
            ys = [k, k]
            plt.plot(xs, ys, color=str((j % per_run) / per_run))
        print(
            "plotted",
            (i + 1) * per_run,
            "/",
            nFrag,
            "for contig",
            contig_name,
            "so far.",
        )
        print(
            "tstart, tend here:",
            new_df.at[(i + 1) * per_run, "tstart"],
            new_df.at[(i + 1) * per_run, "tend"],
        )

    print("Wrapping up - computing remainder of fragments")
    for j in range(nRuns * per_run, nFrag):
        # k = fragment_start + j
        k = j
        xs = [new_df.at[k, "tstart"], new_df.at[k, "tend"]]
        ys = [k, k]
        plt.plot(xs, ys, color=str((j % per_run) / per_run))
    print("Done plotting for contig", contig_name)
    plt.savefig(output_fp)
    print(f"Wrote figure for contig {contig_name} to " + output_fp)
    plt.clf()


def reconcile_contig2len(contig2len: Dict[str, int], contigs: List[str]) -> Dict[str, int]:
    # This function exists because there may be separate contig names
    # like X.1 vs X
    origs = sorted(list(contig2len.keys()))
    contigs = sorted(contigs)

    if len(contigs) > len(origs):
        raise Exception(
            "More contig names in gff file than genome fna file, cannot continue."
        )

    # This dict matches original contig name to possible new contig name with suffix
    match_d: Dict[str,str] = special_match_contig_names(origs, contigs)

    new_ctg2_orig = {}
    for i in range(len(origs)):
        if match_d[origs[i]] in contigs:
            new_ctg2_orig[match_d[origs[i]]] = contig2len[origs[i]]

    return new_ctg2_orig


if __name__ == "__main__":
    _, op_lib_dir = sys.argv
    lib_name = op_lib_dir.split("/")[-1]
    run_step_6_singlelib(op_lib_dir, lib_name, {})
