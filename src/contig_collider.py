from typing import List, Tuple, Optional, Dict

"""

One central entry point is 'special_match_contig_names'.

"""


def match_contig_names(A_df, B_df, debug=False):

    a_ctgs: List[str] = sorted(list(A_df["contig"].unique()))
    b_ctgs: List[str] = sorted(list(B_df["contig"].unique()))
    a_ctgs_new, b_ctgs_new = collider(a_ctgs, b_ctgs, presorted=True)

    print("Now renaming the contig values")
    for i in range(len(a_ctgs_new)):
        A_df["contig"].replace(a_ctgs[i], a_ctgs_new[i], inplace=True)
    for i in range(len(b_ctgs_new)):
        B_df["contig"].replace(b_ctgs[i], b_ctgs_new[i], inplace=True)

    return A_df, B_df


# Called from step6.py and validate.py
def special_match_contig_names(a_ctgs: List[str], b_ctgs: List[str]) -> Dict[str, str]:
    """
    This returns a dictionary mapping original 'a' contig name to the
    name with the suffix, if there is a simple suffix. If there is a match
    between the two set of names, then it just returns the original names.
    """
    a_ctgs = get_contig_names_from_strs(a_ctgs)
    a_ctgs = sorted(a_ctgs)
    orig_a_ctgs = a_ctgs[:]
    b_ctgs = sorted(b_ctgs)
    try:
        res = collider(a_ctgs, b_ctgs, presorted=True)
        a_ctgs_new: List[str] = res[0]
        b_ctgs_new: List[str] = res[1]
    except Exception:
        raise Exception(
            "Failed to match contigs from genome fna and Genbank File - "
            + "please rename contigs to be the same."
        )
    op_d = {}
    for i in range(len(a_ctgs_new)):
        op_d[orig_a_ctgs[i]] = a_ctgs_new[i]
    return op_d


def get_contig_names_from_strs(ctgs: List[str]) -> List[str]:
    return [x.split(" ")[0] for x in ctgs]


def collider(
    ctgsA: List[str], ctgsB: List[str], presorted=False
) -> Tuple[List[str], List[str]]:
    """
    Central part of the program
    """
    if not presorted:
        ctgsA = sorted(ctgsA)
        ctgsB = sorted(ctgsB)

    total_miss_bool = check_for_total_miss(ctgsA, ctgsB)

    if total_miss_bool:
        # There is no overlap between them
        ctgsA, ctgsB = attempt_reconcile(ctgsA, ctgsB)

    return ctgsA, ctgsB


def check_for_total_miss(ctgsA, ctgsB) -> bool:
    """
    This checks if there is zero overlap between the
    contig names.
    """

    cA_set = set(ctgsA)
    cB_set = set(ctgsB)
    intersection_set = cA_set.intersection(cB_set)
    if len(intersection_set) > 0:
        # Some overlap
        return False
    else:
        # No overlap
        return True


def attempt_reconcile(ctgsA, ctgsB):
    """
    Note this function is only called if there is
    no overlap between ctgsA and ctgsB.
    Also, ctgsA and ctgsB are sorted.
    Note that length ( new_ctgsA ) = length (ctgsA), the difference would be
    that the new sets have the same suffix

    """

    prfx: Tuple[Optional[str], Optional[str]] = check_suffix(ctgsA, ctgsB)
    if prfx[0] is None:
        raise RuntimeError("Contigs don't match and have no matching suffix")
    else:
        new_ctgsA, new_ctgsB = fix_using_suffix(ctgsA, ctgsB, prfx[0], prfx[1])
    return new_ctgsA, new_ctgsB


def fix_using_suffix(ctgsA, ctgsB, suffix: str, large_tag: str):
    if large_tag == "A":
        ctgsA = [x + suffix for x in ctgsA]
    elif large_tag == "B":
        ctgsB = [x + suffix for x in ctgsB]
    else:
        raise RuntimeError(
            f"Incorrectly running function, large_tag should be 'A' or 'B', instead {large_tag}"
        )

    return ctgsA, ctgsB


def check_suffix(ctgsA, ctgsB):
    first = ctgsA[0]
    for i in range(len(ctgsB)):
        if ctgsB[i].startswith(first):
            extra_tag = ctgsB[i][len(first) :]
            return (extra_tag, "A")

    first = ctgsB[0]
    for i in range(len(ctgsA)):
        if ctgsA[i].startswith(first):
            extra_tag = ctgsA[i][len(first) :]
            return (extra_tag, "B")

    return (None, None)
