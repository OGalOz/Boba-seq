import pandas as pd
import numpy as np
import sys
from typing import List, TypeVar, Dict


T = TypeVar("T")


def parse_paf_file(paf_fp, cfg_d) -> pd.DataFrame:
    """
    .paf format, minimap2:
    Col Type    Description
    1   string  Query sequence name
    2   int Query sequence length
    3   int Query start (0-based; BED-like; closed)
    4   int Query end (0-based; BED-like; open)
    5   char    Relative strand: "+" or "-"
    6   string  Target sequence name
    7   int Target sequence length
    8   int Target start on original strand (0-based)
    9   int Target end on original strand (0-based)
    10  int Number of residue matches
    11  int Alignment block length
    12  int Mapping quality (0-255; 255 for missing)

    Note, other columns:
    tp  A   Type of aln: P/primary, S/secondary and I,i/inversion
    cm  i   Number of minimizers on the chain
    s1  i   Chaining score
    s2  i   Chaining score of the best secondary chain
    NM  i   Total number of mismatches and gaps in the alignment
    AS  i   DP alignment score
    ms  i   DP score of the max scoring segment in the alignment
    nn  i   Number of ambiguous bases in the alignment
    """
    dtype: Dict[str, str] = {
        "read_name": str,
        "qlen": int,
        "qstart": int,
        "qend": int,
        "strand": str,
        "target": str,
        "tlen": int,
        "tstart": int,
        "tend": int,
        "matches": int,
        "aln_len": int,
        "qual": np.int16,
        "v1": str,
        "v2": str,
        "v3": str,
        "v4": str,
        "v5": str,
        "v6": str,
    }
    """
        "aln_type": str,
        "minimizer": str,
        "chain_score": float,
        "chain_s_2nd": str,
        "per_base_div": str,
        "rep_seeds": str
        "tp":str, 
        "cm":str,
        "s1":str,
        "s2":str,
        "NM":str,
        "AS":str,
        "ms":str,
        "nn":str
    """
    col_names: List[str] = list(dtype.keys())
    df = pd.read_table(paf_fp, names=col_names, dtype=dtype)
    for i in range(1, 7):
        del df["v" + str(i)]

    num_reads_searched: int = len(df["read_name"].unique())
    not_mapped_uniq: int = len(df[df["strand"] == "*"]["read_name"].unique())
    mapped_uniq: int = len(df[df["strand"] != "*"]["read_name"].unique())
    init_tot_hits = df[df["strand"] != "*"].shape[0]

    log_d: Dict[str, T] = {
        "reads_in_paf": num_reads_searched,
        "reads_not_mapped": not_mapped_uniq,
        "reads_mapped": mapped_uniq,
    }
    # Exclude reads with quality less than X
    df = df[df["qual"] >= cfg_d["minimap_qual_min"]]

    # Compute stats for Logs/
    hits_high_qual = df.shape[0]
    reads_high_qual = df["read_name"].unique().shape[0]
    perc_high_qual = round(hits_high_qual / init_tot_hits * 100, 2)
    perc_read_high_qual = round(reads_high_qual / mapped_uniq * 100, 2)

    reads_one_hit = df.drop_duplicates(subset="read_name", keep=False).shape[0]
    reads_more_one_hit = df[df.duplicated(subset="read_name", keep="first")].shape[0]
    perc_mapped_with_one_hit = round(reads_one_hit / reads_high_qual * 100, 2)
    perc_mapped_with_more = round(reads_more_one_hit / reads_high_qual * 100, 2)

    # a small % of reads have multiple hits, ignore these and write to file
    df_multi_hits = df[df.duplicated(subset="read_name", keep=False)]
    df_filtered = df.drop_duplicates(subset="read_name", keep=False, ignore_index=True)

    log_d["hits_high_qual"] = hits_high_qual
    log_d["reads_high_qual"] = reads_high_qual
    log_d["reads_high_qual/reads_mapped, perc"] = perc_read_high_qual
    log_d["perc_reads_w_one_hit"] = perc_mapped_with_one_hit
    log_d["perc_reads_w_more_hits (excluded)"] = perc_mapped_with_more

    # nMatches / read_name Length
    df_filtered.loc[:, "perc_match_cov"] = df_filtered["matches"] / df_filtered["qlen"]
    # nMatches / Alignment Length
    df_filtered.loc[:, "perc_match"] = df_filtered["matches"] / df_filtered["aln_len"]

    return df_filtered, df_multi_hits, log_d


def get_read_name_to_info_dict(paf_df: pd.DataFrame) -> Dict[str, List[T]]:
    """
    Output Dictionary maps read name to other info from paf dataframe.
    This allows us to match barcode to other info more quickly.
    At this point there should be no duplicate read_name values (removed in parse_paf_file)
    """
    info_dict: Dict[str, int] = {}
    cols = list(paf_df.columns)
    cols.remove("read_name")
    N = len(cols)
    M = paf_df.shape[0]
    L = 10**5

    for ix, row in paf_df.iterrows():
        #         if ix % L == 0:
        #             print(f"finished making {ix}/{M} read2ix")
        if row["read_name"] not in info_dict:
            info_dict[row["read_name"]] = ix
        else:
            print(f"WARNING: {row['read_name']} appears twice in filtered paf.")
    return info_dict


if __name__ == "__main__":
    _, paf_fp, op_path = sys.argv
    df = parse_paf_file(paf_fp)
    read2info_list: Dict = get_read_name_to_info_dict(df)
