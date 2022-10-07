"""
In this file we attempt to match each barcode to one or two regions;
Functions in which barcodes are decided are:
    get_best_mapping_for_barcodes
    
Two functions called under main:
    bc_df_collapse_steps
    export_collapsed_bc_files
    
Currently imported by step 5.
"""

import pandas as pd
import os
from typing import List, Dict, Tuple, Union
from collections import defaultdict, Counter
import sys
import json
import statistics
import matplotlib.pyplot as plt


def bc_df_collapse_steps(
    bc_df, cp, dfs_dir, cfg_d, lib_name, print_dbg=False, bc_to_loc_dicts=None
):
    """
    Filter df by perc_match and perc_cov
    Additional filters using get_best_mapping_for_barcodes
    Barcodes mapped to multiple pos and contigs are exported to ["redundant_bc"] --> two dictionaries, one with all multi mappings and one for BCs across diff contigs
    Barcodes with best mapped positions to ["best_bc_pos"]
    """

    log_list: List[str] = []
    init_BC_num = str(bc_df.shape[0])
    init_uniq_BC = str(bc_df["bc"].unique().shape[0])

    bc_df = filter_out_low_percent_match(bc_df, min_perc_match=cp["min_perc_match"])
    BC_num_qual_filter = str(bc_df.shape[0])
    post_qual_filter_uniq_BC = str(bc_df["bc"].unique().shape[0])

    bc_df = filter_out_low_perc_cov(bc_df, min_cov=cp["min_perc_cov"])
    bc_df.reset_index(inplace=True, drop=True)
    BC_num_perc_cov_filter = str(bc_df.shape[0])
    post_perc_cov_uniq_BC = str(bc_df["bc"].unique().shape[0])

    log_list.append("Total BCs: init, after perc_match filter, after perc_cov filter")
    log_list.append(
        ", ".join([init_BC_num, BC_num_qual_filter, BC_num_perc_cov_filter])
    )
    log_list.append("Unique BCs: init, after perc_match filter, after perc_cov filter")
    log_list.append(
        ", ".join([init_uniq_BC, post_qual_filter_uniq_BC, post_perc_cov_uniq_BC])
    )

    bc2loc_fp = os.path.join(dfs_dir, cfg_d["d"]["fns"]["5"]["bc2locs"])
    # convert df to dict of bc: {(locations): {count}}
    if not bc_to_loc_dicts:
        bc_to_loc_dicts, new_log = get_bc_to_locs(bc_df, bc2loc_fp, write_to_file=True)
        # bc_to_loc_dicts: Dict[str, Dict[Tuple[int, int, str, str], int]]
    log_list += new_log

    best_mappings_d, failed_bcs, multi_contig_bcs = get_best_mapping_for_barcodes(
        bc_to_loc_dicts, cp, print_dbg=print_dbg
    )

    # best_mapping_d: dict of bc: [Tuple() for each best location]
    # converts this dict into a form that can be easily converted to df
    pre_best_mappings_df = []
    for k, v in best_mappings_d.items():
        for z in v:
            z = list(z)
            pre_best_mappings_df.append([k] + z[:3] + ["T" if z[3] else "F"] + z[4:])
    df_cols = ["bc", "tstart", "tend", "contig", "multi_frag", "strand"]
    best_mappings_df = pd.DataFrame(pre_best_mappings_df, columns=df_cols)
    best_mappings_df["frag_len"] = (
        best_mappings_df["tend"] - best_mappings_df["tstart"] + 1
    )
    del pre_best_mappings_df

    # We convert bc_to_loc_dicts to df and merge onto best_mappings_df to generate count col
    ls_bc_to_loc = []
    for bc, bc_dict in bc_to_loc_dicts.items():
        for k, count in bc_dict.items():
            k = list(k)  # Tuple[int, int, str, str]
            ls_bc_to_loc.append([bc, k[0], k[1], k[2], k[3], count])
    df_cols = ["bc", "tstart", "tend", "contig", "strand", "count"]
    bc_to_loc_df = pd.DataFrame(ls_bc_to_loc, columns=df_cols)

    # merge best_mappings_df and bc_to_loc_df and get count column
    nRows_before = best_mappings_df.shape[0]
    best_mappings_df = best_mappings_df.merge(
        bc_to_loc_df, on=["bc", "tstart", "tend", "contig", "strand"], how="inner"
    )
    nRows_after = best_mappings_df.shape[0]
    del ls_bc_to_loc, bc_to_loc_df

    # barcodes with multiple mappings
    bc_multi_mappings = best_mappings_df[
        best_mappings_df.duplicated(subset="bc", keep="first")
    ]
    ls_bc_multi_mappings = bc_multi_mappings["bc"].unique().tolist()

    return [
        best_mappings_df,
        best_mappings_d,
        ls_bc_multi_mappings,
        multi_contig_bcs,
        failed_bcs,
        log_list,
    ]


def export_collapsed_bc_files(
    best_mappings_df,
    ls_bc_multi_mappings,
    multi_contig_bcs,
    failed_bcs,
    dfs_dir,
    op_lib_dir,
    cfg_d,
    lib_name,
):
    """
    Write to file: failed barcodes, multi_mappings, good barcodes
    Plots histogram via create_bc_to_loc_count_histogram
    """
    log_list: List[str] = []
    nUniq_good_bc = best_mappings_df["bc"].unique().shape[0]
    nUniq_multi_mapped = len(ls_bc_multi_mappings)
    nUniq_multi_contig = len(multi_contig_bcs)
    nUniq_failed = len(failed_bcs)
    perc_multi = round(nUniq_multi_mapped / nUniq_good_bc * 100, 3)
    perc_bad = round(nUniq_failed / (nUniq_failed + nUniq_good_bc) * 100, 3)

    # stats of mapped, failed, multi
    log_list.append(
        "Found "
        + str(best_mappings_df.shape[0])
        + " best mappings for "
        + str(nUniq_good_bc)
        + " barcodes."
    )
    log_list.append("nUniq barcodes that were mapped: " + str(nUniq_good_bc))
    log_list.append(
        "nUniq barcodes mapped to multiple contigs: " + str(nUniq_multi_contig)
    )
    log_list.append(
        "nUniq barcodes mapped to multiple positions: " + str(nUniq_multi_mapped)
    )
    log_list.append("Perc of barcodes mapped to multiple positions: " + str(perc_multi))
    log_list.append("nUniq barcodes that failed mapping: " + str(nUniq_failed))
    log_list.append("Perc of barcodes that failed mapping: " + str(perc_bad))

    # Write barcodes w/ multiple mappings
    multi_mapping_fp = os.path.join(dfs_dir, cfg_d["d"]["fns"]["5"]["redundant_bc"])
    out_multi = {
        "bc_w_multi_mappings": ls_bc_multi_mappings,
        "multi_contig_bcs": multi_contig_bcs,
    }
    with open(multi_mapping_fp, "w") as g:
        g.write(json.dumps(out_multi, indent=4))

    print(
        "Wrote barcodes w/ multiple mappings to "
        + cfg_d["d"]["fns"]["5"]["redundant_bc"]
    )

    # Write list of failed BCs
    failed_mappings_fp = os.path.join(dfs_dir, cfg_d["d"]["fns"]["5"]["failed_bcs"])
    with open(failed_mappings_fp, "w") as g:
        failed_d = {x[0]: x[1] for x in failed_bcs}
        printed_fail_d = {
            k: {str(m): n for m, n in v.items()} for k, v in failed_d.items()
        }
        g.write(json.dumps(printed_fail_d, indent=4))
        print("Wrote failed mappings to " + cfg_d["d"]["fns"]["5"]["failed_bcs"])

    # Write best_mappings_df to ["best_bc_pos"]
    best_mappings_fp = os.path.join(dfs_dir, cfg_d["d"]["fns"]["5"]["best_bc_pos"])
    best_mappings_df.to_csv(best_mappings_fp, sep="\t", index=False)
    print(
        "Wrote best locations for each barcode to "
        + cfg_d["d"]["fns"]["5"]["best_bc_pos"]
    )

    # Table of # of positions per BC and freq to log
    bc_pos_counts = best_mappings_df.groupby(["bc"]).size().to_frame("num_positions")
    bc_pos_counts = bc_pos_counts.reset_index()
    # list BCs with high # of mappings:
    n_cutoff = 2
    ls_bc_high_mappings = bc_pos_counts[bc_pos_counts["num_positions"] > n_cutoff][
        "bc"
    ].tolist()
    # convert to # pos per BC to freq:
    num_pos_per_BC_freq = (
        bc_pos_counts.groupby(["num_positions"]).size().to_frame("freq")
    )
    num_pos_per_BC_freq = num_pos_per_BC_freq.reset_index().sort_values(
        by="num_positions"
    )
    # write to log:
    log_list.append("Number of positions mapped per barcode to frequency:")
    log_list.append(str(num_pos_per_BC_freq.values.tolist()))
    log_list.append(f"Barcodes mapped to >{n_cutoff} locations:")
    log_list.append(", ".join(ls_bc_high_mappings))

    create_bc_to_loc_count_histogram(
        op_lib_dir, best_mappings_df, dfs_dir, cfg_d, lib_name
    )

    return log_list


def create_bc_to_loc_count_histogram(op_lib_dir, final_bc_df, dfs_dir, cfg_d, lib_name):
    """
    df input: after all filters, good BCs, including ones with multi locations.
    Writes to "nHits2LocFreq" as json: # of times a position is observed --> freq of this
    Generates a histogram of this to in op_lib_dir/Plots
    This is used to visualize redundancy in insert fragments from PCR amplification steps during cloning.
    """
    print("Creating location by barcode histogram")
    plots_dir = os.path.join(op_lib_dir, "Plots")
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)

    # First we count number of times a precise location is seen and write to dict
    loc2count_d = defaultdict(int)
    for i in range(final_bc_df.shape[0]):
        loc2count_d[(final_bc_df.at[i, "tstart"], final_bc_df.at[i, "tend"])] += 1

    # Next we take this dictionary with location -> number of times seen
    # and computes the frequency of each location seen
    # Ex. (number of times a location is seen) -> frequency of this
    reverse_d = defaultdict(int)
    for loc in loc2count_d:
        reverse_d[loc2count_d[loc]] += 1

    with open(os.path.join(dfs_dir, cfg_d["d"]["fns"]["P"]["nHits2LocFreq"]), "w") as g:
        g.write(json.dumps(reverse_d, indent=2, sort_keys=True))

    # We use matplotlib to make a histogram using the values from the above dict
    hist_points = list(loc2count_d.values())
    bins = list(range(min(hist_points), max(hist_points) + 2))
    #     plt.hist(hist_points, bins=max(list(reverse_d.keys())))
    plt.hist(hist_points, bins=bins)
    plt.title("BC coverage per location, " + lib_name)
    plt.xlabel("Number of BCs mapped to a precise location")
    plt.ylabel("Occurrence")
    op_fp = os.path.join(plots_dir, lib_name + cfg_d["d"]["fns"]["P"]["good_bc_hist"])
    plt.savefig(op_fp)
    plt.clf()

    print("Wrote number of BC at each location histogram to " + op_fp)
    return None


def get_best_mapping_for_barcodes(bc_to_loc_dicts, cp, print_dbg=False):
    """
    bc_to_loc_dicts is a dictionary mapping barcode to a dictionary which
    # holds tuple of location (start, end, contig) -> number of times that location
    # appeared
    e.g.
        {
            "ACC..": {
                (12,20, contigA, +) : 5
                (12,19, contigA, +): 1
            },
            "TACC": {
                (1224, 3000, contigC, -): 24,
                (1224, 2778, contigD, +): 8
            }
        }
    cp is a config dictionary that holds vars to help find best mapping
        out of a group of mappings
    e.g.
        "collapse_params": {
        "min_perc_cov": 0.9,
        "min_perc_match": 0.9,
        "min_BC_support_ratio": 3,
        "max_frag_diff_to_len": 0.25}

    Returns:
        best_mapping_d (dict)
            Maps barcode to tuple (start,end) of best location. Multi_frag_bool
            is True when there are multiple positions assigned to the same Barcode.
            "ACCT" : [(224, 351, contig, T (multi_frag bool) ),(3355, 4466, contig, T), ..]
            "ATTC": [(3355, 4466, contig, F)]
            .
        failed_bcs list<dict>: List of tuples,(bc, loc_ds) that failed
        Note that multi_frag bool will be True if there are many locations for a single
        barcode, and False when there is only one.

    """
    print("Beginning process of finding best locations for each barcode.")
    best_mapping_d = {}
    failed_bcs = []
    multi_contig_bcs = []
    for bc in bc_to_loc_dicts:
        # barcodes with a single position
        if len(bc_to_loc_dicts[bc]) == 1:
            # Note that in this case k is a tuple (start, end, contig, strand)
            # count is combined into tuple
            k = list(bc_to_loc_dicts[bc].keys())[0]
            #             count = bc_to_loc_dicts[bc][k]
            #             best_mapping_d[bc] = [(k[0], k[1], k[2], False, k[3], count)]
            best_mapping_d[bc] = [(k[0], k[1], k[2], False, k[3])]
        else:
            # for each bc, and a dict of bc2loc
            best_mappings: Union[List, bool] = compute_best_mappings(
                bc, bc_to_loc_dicts[bc], cp
            )
            if not best_mappings:
                # best_mapping returns False, failed to pass both filters
                failed_bcs.append((bc, bc_to_loc_dicts[bc]))
            else:
                best_mapping_d[bc] = best_mappings

                # If across multiple contigs, add to multi_contig_bcs
                if best_mappings[0][3]:  # if True
                    # Contig of first option
                    crt_contig = best_mappings[0][2]
                    for x in best_mappings:
                        if x[2] != crt_contig:
                            multi_contig_bcs.append(bc)
    multi_contig_bcs = list(
        set(multi_contig_bcs)
    )  # remove duplicates in case of >2 contigs

    # writes best_mapping_d for debugging
    if print_dbg:
        out_fp = os.path.join(dfs_dir, "best_mappings_debug.json")
        print("Printing debug best_mappings to", out_fp)
        with open(out_fp, "w") as g:
            g.write(json.dumps(best_mapping_d, indent=2))

    return (best_mapping_d, failed_bcs, multi_contig_bcs)


def compute_best_mappings(
    bc: str, loc2num_d: Dict, cp: Dict
) -> Union[List[Tuple[int, int, str, bool]], bool]:
    """
    loc2num_d looks like:
    {
        (12,20, contigA, +) : 5,
        (12,19, contigA, +): 1,
        (25,42, contigA. +): 3,
        (120,150, contigB, -): 13,
        .
    }
    cp is a config dict (as listed above)
    Uses min_BC_support_ratio and max_frag_diff_to_len to identify best positions.
    Returns:
        List of all good locations for this barcode
        e.g.
            [(12,20, contigA, True, +), (25,42, contigA, True, +), (120,150, contigB, True, -)]
    """
    # We split locations into groups by non-overlapping regions
    split_by_ranges: List[Dict] = split_locs_by_ranges(loc2num_d)

    # compute best position for each list of overlapping fragments
    final_loc_list = []
    for split_range_d in split_by_ranges:
        # split_range_d should always have at least a single key
        final_loc_list += compute_best_mappings_by_split_range(bc, split_range_d, cp)

    if len(final_loc_list) > 0:
        if len(final_loc_list) > 1:
            # multi frag for one barcode
            final_loc_list = [(x[0], x[1], x[2], True, x[3]) for x in final_loc_list]
            return final_loc_list
        else:
            # Not multi frag - single frag
            x = final_loc_list[0]
            final_loc_list = [(x[0], x[1], x[2], False, x[3])]
            return final_loc_list
    else:
        return False


def split_locs_by_ranges(loc2num_d: Dict) -> List[Dict]:
    """
    loc2num_d looks like:
    {
        (12,20, contigA, '+') : 5,
        (12,19, contigA, '+'): 1,
        (25,42, contigA. '+'): 3,
        (120,150, contigB, '-'): 13,
        .
    }
    Description:
        We want to split it into a list of dictionaries in which each
        subdictionary is a different 'range' of non-overlapping locations, including
        on separate contigs.
    Algorithm:
        First you split up by contigs into different bins (a list of dictionaries).
        Then you take each contig bin and split by range (non-overlapping regions)
        into a list of dictionaries where each list contains overlapping fragments.
    Returns:
        A list of dictionaries where each list contains overlapping fragments on the same contig.
    """
    bins = []
    contig_bins: List[Dict] = split_by_contigs(loc2num_d)
    for cbin in contig_bins:
        bins += split_by_range_bins(cbin)

    return bins


def split_by_contigs(loc2num_d) -> List[Dict]:
    """
    {
        (12,20, contigA, '+') : 5,
        (12,19, contigA, '+'): 1,
        (25,42, contigA. '+'): 3,
        (120,150, contigB, '-'): 13,
        .
        .
    }
    Description:
        We want to split it into a list of dictionaries in which each
        subdictionary is a different 'range' of non-overlapping locations, including
        on separate contigs.
    Algorithm:
        Sort by contig name. Maintain last contig name. Once it shifts, save the collection
        and start a new empty collection.
    Returns:
        List of dictionaries that share the same contig within each list
            {
        (12,20, contigA, '+'): 5,
        (12,19, contigA, '+'): 1,
        (25,42, contigA. '+'): 3,

        (120,150, contigB, '-'): 13,
        .
        .
    }
    """
    bins = []
    # sort ascending on third key value (contig)
    sorted_keys = sorted(list(loc2num_d.keys()), key=lambda x: x[2])
    crt_contig = sorted_keys[0][2]
    new_bin = []
    i = 0
    while i < len(sorted_keys):
        crt_k = sorted_keys[i]
        if crt_k[2] == crt_contig:
            new_bin.append(crt_k)
        else:
            # converts list keys to a subdict of loc2num_d using the list of keys previously stored; then stores this dict as a list in bin output
            new_bin = {k: loc2num_d[k] for k in new_bin}
            bins.append(new_bin)
            # starts a new collection of keys
            new_bin = [crt_k]
            crt_contig = crt_k[2]
        i += 1

    new_bin = {k: loc2num_d[k] for k in new_bin}
    bins.append(new_bin)

    return bins


def split_by_range_bins(loc2num_d) -> List[Dict]:
    """
    {
        (12,20, contigA, '+') : 5,
        (12,19, contigA, '+'): 1,
        (25,42, contigA. '+'): 3,
        (120,150, contigB, '-'): 13,
        .
    }
    Description:
        We want to split it into a list of dictionaries in which each
        subdictionary is a different 'range' of non-overlapping locations.
    Algorithm:
        Sort by starting value. Keep a fixed start and end range for current bin.
        First you have the first key, e.g. start = 12, end =20, and you initialize a new bin.
        Then you check the next value. If there is an overlap, you might expand the end, so
        now start = 12, end = 22. As long as there is overlap, you add to the current bin.
        Once a new value occurs that has no overlap, you add the bin to bins, and initialize
        a new empty bin.
    Return:
        List of dictionaries that share overlapping fragments within one list.
    """
    bins = []
    # sort ascending on first key value (start position)
    sorted_keys = sorted(list(loc2num_d.keys()), key=lambda x: x[0])
    end = sorted_keys[0][1]
    new_bin = []
    i = 0
    while i < len(sorted_keys):
        crt_k = sorted_keys[i]
        if crt_k[0] <= end:
            new_bin.append(crt_k)
            end = max(end, crt_k[1])
        else:
            # converts list keys to a subdict of loc2num_d using the list of keys previously stored; then stores this dict as a list in bin output
            new_bin = {k: loc2num_d[k] for k in new_bin}
            bins.append(new_bin)
            new_bin = [crt_k]
            end = crt_k[1]
        i += 1

    new_bin = {k: loc2num_d[k] for k in new_bin}
    bins.append(new_bin)

    return bins


def compute_best_mappings_by_split_range(
    bc: str, loc2num_d: Dict, cp: Dict
) -> List[Tuple[int, int, str]]:
    """
    At this point, the barcodes are already split up into a list of non-overlapping ranges.
    All the locations within one list should be overlapping.
    bc is a single barcode
    loc2num_d or split_range_d looks like:
    {
        (12,20, contigA, '+') : 5,
        (12,19, contigA, '+'): 1,
    }
    cp is a config dict (as listed above)
    Applies min_BC_support_ratio and max_frag_diff_to_len filters to identify a best position.
    Even if one filter passes, location with top count is kept.
    If neither filter passes, returns an empty list.
    Returns:
        List of all good locations for this barcode (keys in loc2num_d)
        Could be an empty list []
    """
    # list of barcodes that pass each filter:
    good_locs_by_ratio = []
    good_locs_by_len = []

    locs: List[Tuple[int, int, str]] = list(loc2num_d.keys())
    # sort by counts per position: highest to lowest
    mapping_vals = sorted(list(loc2num_d.values()), reverse=True)

    # single position for this region, return
    if len(mapping_vals) == 1:
        #         count = list(loc2num_d.values())[0]
        #         return [list(loc2num_d.keys())[0] + (count,)]
        return [list(loc2num_d.keys())[0]]

    ## mapping_vals is a list with length > 1
    # Filtering for min_BC_support_ratio parameter:
    if mapping_vals[0] / mapping_vals[1] >= cp["min_BC_support_ratio"]:
        for k in loc2num_d.keys():
            if loc2num_d[k] == mapping_vals[0]:
                good_locs_by_ratio.append(k)
                break

        if len(good_locs_by_ratio) == 0:
            raise RuntimeError(
                "Key for best location not found somehow at min_BC_support_ratio filter step?\n"
                + f"Key: {k}, Value: {mapping_vals[0]}, bc: {bc}"
            )

    # Filtering for max_frag_diff_to_len parameter:
    min_start, max_start, min_end, max_end = compute_min_start_and_end(locs)
    start_diff = max_start - min_start
    end_diff = max_end - min_end
    total_diff = start_diff + end_diff
    average_fragment_len = statistics.mean([x[1] - x[0] + 1 for x in locs])

    if total_diff / average_fragment_len < cp["max_frag_diff_to_len"]:
        for k in locs:
            if loc2num_d[k] == mapping_vals[0]:  # take highest count
                #                 count = loc2num_d[k]
                #                 k = k + (count,)
                good_locs_by_len.append(k)
        if len(good_locs_by_len) == 0:
            raise RuntimeError(
                "Key for best location not found somehow at max_frag_diff_to_len filter step?\n"
                + f"Key: {k}, Value: {mapping_vals[0]}, bc: {bc}"
            )

    # if loc passed both filters, both lists should be the same
    # if only one parameter passed, also okay
    good_locs_list = good_locs_by_len + good_locs_by_ratio
    good_locs_list = list(set(good_locs_list))  # remove duplicates

    return good_locs_list


def compute_min_start_and_end(locs: List[Tuple[int, int, str]]):
    # In this situation we are simply computing the min_starts and min_ends
    # Below is initialization of first tuple:
    min_start = locs[0][0]
    max_start = locs[0][0]
    min_end = locs[0][1]
    max_end = locs[0][1]
    for l in locs[1:]:
        if l[0] < min_start:
            min_start = l[0]
        elif l[0] > max_start:
            max_start = l[0]
        if l[1] < min_end:
            min_end = l[1]
        elif l[1] > max_end:
            max_end = l[1]
    return min_start, max_start, min_end, max_end


def get_bc_to_locs(bc_df, op_fp, write_to_file):
    """
    Returns a dictionary:
        {bc -> {
                (start, end, cpntig): Ntimes seen,
                (start, end, contig): Ntimes seen,
                .
                }
        bc ->
        .
        }
    """

    print("Extracting BC to location info...")

    bc_to_locs: Dict[str, List[str]] = defaultdict(list)
    x = 10**4
    y = bc_df.shape[0]
    for i in range(bc_df.shape[0]):
        bc_to_locs[bc_df.at[i, "bc"]].append(
            (
                int(bc_df.at[i, "tstart"]),
                int(bc_df.at[i, "tend"]),
                bc_df.at[i, "target"],
                bc_df.at[i, "strand"],
            )
        )
    #         if i % x == 0:
    #             print(f"Finished {i}/{y} so far")

    bc_to_loc_dicts = {}
    n_single_loc = 0
    for bc in bc_to_locs:
        bc_to_loc_dicts[bc] = Counter(bc_to_locs[bc])
        if len(bc_to_loc_dicts[bc]) == 1:
            n_single_loc += 1

    log_list: List[str] = []
    log_list.append(
        "BCs mapped to a single location (no need to filter): " + str(n_single_loc)
    )

    printable_bc2locdicts = {}
    for bc in bc_to_loc_dicts:
        printable_bc2locdicts[bc] = {str(k): v for k, v in bc_to_loc_dicts[bc].items()}

    if write_to_file == True:
        with open(op_fp, "w") as g:
            g.write(json.dumps(printable_bc2locdicts, indent=2))

    return bc_to_loc_dicts, log_list


def filter_out_low_percent_match(bc_df: pd.DataFrame, min_perc_match=1):
    # Note default 'min_perc_match' normally not used.
    new_df = bc_df[bc_df["perc_match"] >= min_perc_match]
    return new_df


def filter_out_low_perc_cov(bc_df, min_cov=0.92):
    # Note default 'min_cov' normally not used.
    new_df = bc_df[bc_df["perc_match_cov"] >= min_cov]
    return new_df


def import_barcodes_file(bc_fp: str) -> pd.DataFrame:
    """
    Should have the following columns:
        bc	read_name	qlen	qstart	qend	strand	target
        tlen	tstart	tend	matches	aln_len	qual	perc_match_cov
        perc_match	gene_count	locus_tags
    """
    df = pd.read_table(bc_fp)
    return df


def main_collapse_barcodes(op_lib_dir, lib_name, cfg_d):
    """
    cfg_d has to have the following params:
        "collapse_params":
            "min_perc_cov": float (Fraction of query overlap
                                    with alignment)
            "min_perc_match": int (Minimum percent_match)
    """

    cp = clps_check_cfg(cfg_d)
    dfs_dir = os.path.join(op_lib_dir, cfg_d["d"]["steps2dirs"]["5"])
    bc_fp = os.path.join(dfs_dir, cfg_d["d"]["fns"]["5"]["bc_loc_df"])
    # Note that gene count and locus tags aren't in use
    bc_df = import_barcodes_file(bc_fp)

    (
        best_mappings_df,
        best_mappings_d,
        ls_bc_multi_mappings,
        multi_contig_bcs,
        failed_bcs,
        log_list,
    ) = bc_df_collapse_steps(bc_df, cp, dfs_dir, cfg_d, lib_name)

    export_collapsed_bc_files(
        best_mappings_df,
        ls_bc_multi_mappings,
        multi_contig_bcs,
        failed_bcs,
        dfs_dir,
        op_lib_dir,
        cfg_d,
        lib_name,
    )


if __name__ == "__main__":
    _, op_lib_dir, cfg_fp = sys.argv
    lib_name = op_lib_dir.split("/")[-1]
    with open(cfg_fp) as f:
        cfg_d = json.loads(f.read())
    main_collapse_barcodes(op_lib_dir, lib_name, cfg_d)
