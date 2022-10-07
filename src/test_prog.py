"""
Basic testing file 

Write tests for all imported functions from stepX files.

"""

import os
import sys
import json
import pandas as pd
import random
import string
from typing import List, Dict, Tuple, Any, Union, Optional
import unittest

from step1 import (
    run_usearch_search_pcr2_command,
    run_usearch_search_oligodb_command,
    concat_files_to_output,
    extract_file_name,
        )

from step5 import (
    gene_fragment_overlap,
    combine_barcodes_with_other_data,
    create_read_to_barcode_dict,
    reshape_dataframe,
)
from collapse_bcs import (
    get_bc_to_locs,
    get_best_mapping_for_barcodes,
    compute_best_mappings,
    split_locs_by_ranges,
    split_by_range_bins,
    split_by_contigs,
    compute_best_mappings_by_split_range,
    compute_min_start_and_end,
)




class TestStringMethods(unittest.TestCase):
    
    def test_reshape_dataframe(self, bc_df_fp: Optional[str] = None):
        # Set bc_df_fp here:

        if bc_df_fp is None:
            raise Exception("bc_df_fp should be set above...")
        bc_df: pd.DataFrame = pd.read_table(bc_df_fp, sep=",")
        post_reshape_bc_df = reshape_dataframe(bc_df)
        self.assertTrue(True)
    
    
    def test_gene_fragment_overlap(self, inp_genes_df_fp: Optional[str] = None, inp_bc_df_fp: Optional[str] = None):
        """
        inp_genes_df needs cols:
            "pos_from"
            "pos_to"
            "contig"
            "strand"
    
        inp_bc_df_fp needs cols:
            "tstart"
            "tend"
            "contig"
            "strand"
    
        """
        # inp_genes_df_fp = ""
        # inp_bc_df_fp = ""
        if inp_genes_df_fp is None or inp_bc_df_fp is None:
            raise Exception("You must set both inp_genes_df_fp and inp_bc_df_fp")

        genes_df = pd.read_table(inp_genes_df_fp)
        bc_df = pd.read_table(inp_bc_df_fp)
        gene_fragment_overlap(genes_df, bc_df, dbg=True)
        self.assertTrue(True)
    


    def test_search_pcr2(self, inp_fq: Optional[str] = None, op_json_fp: Optional[str] = None) -> None:
        """
        # Writes a JSON file in loc indicated
        # cfg_d must have ( for example ):
        #   "fwd": "GTTCTTATCTTTGCAGTCTC",
        #   "rev": "GAGATTTACGCTTTGGTAAAAGTTGG",
        #   "maxdiffs": 2,
        #   "minamp": 300,
        #   "maxamp": 15000
        """
        # Set inp_fq here:

        # Set op_json_fp here:


        cfg_d = {
            "fwd": "GTTCTTATCTTTGCAGTCTC",
            "rev": "GAGATTTACGCTTTGGTAAAAGTTGG",
            "maxdiffs": 2,
            "minamp": 300,
            "maxamp": 15000,
        }
    
        fwd = cfg_d["fwd"]
        rev = cfg_d["rev"]
        maxdiffs = cfg_d["maxdiffs"]
        minamp = cfg_d["minamp"]
        maxamp = cfg_d["maxamp"]
        tuple_list = []
        seq_num = 1
        if inp_fq is None:
            raise Exception("You need to set inp_fq within this function.")

        with open(inp_fq, "r") as f:
            rn = f.readline()
            while rn != "":
                seq = f.readline()
                l = len(seq)
                if l > minamp and l < maxamp:
                    res = spcr2_check_seq(seq, fwd, rev)
                    if res is not None:
                        tuple_list.append((seq_num, res[0], res[1]))
                sep = f.readline()
                qual = f.readline()
                seq_num += 1
    
                if seq_num % 5000 == 0:
                    print(
                        "seq n",
                        seq_num,
                        "n found",
                        len(tuple_list),
                        "\nlast couple\n",
                        tuple_list[-10:],
                    )
                rn = f.readline()
    
        op_d = {"total": len(tuple_list), "tup_list": tuple_list}
        if op_json_fp is not None:
            with open(op_json_fp, "w") as g:
                g.write(json.dumps(op_d))
            print("wrote output to " + op_json_fp)
        self.assertTrue(True)
    

# Function which creates a sample barcodes dataframe. 
def generate_bc_df(
    op_fp: Optional[str] = None,
    start_range: int = 0,
    end_range: int = 200,
    nSamples: int = 30,
    ncontigs: int = 3,
    min_frag_len: int = 10,
    max_frag_len: int = 40,
):
    """
    We create a sample barcode df to test gene fragment overlap
    """
    start_points: List[int] = [
        random.randint(start_range, end_range) for i in range(nSamples)
    ]
    start_points = sorted(start_points)
    end_points = [x + random.randint(min_frag_len, max_frag_len) for x in start_points]
    if ncontigs >= 20:
        raise Exception("For testing expecting less than 20 contigs")
    uppercase_chain = string.ascii_uppercase
    contig_options: List[str] = [uppercase_chain[i] for i in range(ncontigs)]
    contig_list = [random.choice(contig_options) for i in range(nSamples)]
    strand_list = [random.choice(["+", "-"]) for i in range(nSamples)]

    d = {
        "tstart": start_points,
        "tend": end_points,
        "contig": contig_list,
        "strand": strand_list,
    }
    op_df = pd.DataFrame.from_dict(d)
    if op_fp:
        op_df.to_csv(op_fp, index=False, sep="\t")
        print("Wrote out sample BC df to " + op_fp)
    else:
        print("output dataframe", op_df)


# Function which creates a sample genes dataframe. 
def generate_genes_df(
    op_fp: str,
    start_range: int = 0,
    end_range: int = 200,
    nSamples: int = 10,
    ncontigs: int = 3,
    min_frag_len: int = 7,
    max_frag_len: int = 20,
):
    """
    We create a sample genes df to test gene fragment overlap
    """
    start_points: List[int] = [
        random.randint(start_range, end_range) for i in range(nSamples)
    ]
    start_points = sorted(start_points)
    end_points = [x + random.randint(min_frag_len, max_frag_len) for x in start_points]
    if ncontigs >= 20:
        raise Exception("For testing expecting less than 20 contigs")
    if nSamples >= 20:
        raise Exception("For testing expecting less than 20 genes")
    uppercase_chain = string.ascii_uppercase
    contig_options: List[str] = [uppercase_chain[i] for i in range(ncontigs)]
    contig_list = [random.choice(contig_options) for i in range(nSamples)]
    strand_list = [random.choice(["+", "-"]) for i in range(nSamples)]
    lc_chain = string.ascii_lowercase
    locus_tags: List[str] = [lc_chain[i] for i in range(nSamples)]

    d = {
        "pos_from": start_points,
        "pos_to": end_points,
        "contig": contig_list,
        "strand": strand_list,
        "locus_tag": locus_tags,
    }
    op_df = pd.DataFrame.from_dict(d)
    op_df.to_csv(op_fp, index=False, sep="\t")
    print("Wrote out sample genes df to " + op_fp)
    
    
    
# Util test function to check sequence for substrings.    
def spcr2_check_seq(seq: str, fwd: str, rev: str) -> Optional[Tuple[int, int]]:
    x = seq.find(fwd)
    y = seq.find(rev)
    if x == -1 or y == -1:
        return None
    else:
        return (x, y)


def main():
    args = sys.argv
    inp_fq = args[-2]
    op_json_fp = args[-1]
    search_pcr2_util(inp_fq, op_json_fp, cfg_d)


def test_split_by_range_bins():

    test_case = {
        (12, 20, "contigA"): 5,
        (12, 19, "contigA"): 2,
        (13, 22, "contigA"): 1,
        (13, 33, "contigA"): 2,
        (25, 42, "contigA"): 3,
        (45, 55, "contigA"): 4,
        (120, 150, "contigA"): 13,
    }
    result = split_by_range_bins(test_case)
    print("result\n", result)


def run_tests_by_flag(args: List[str], flag: int):
    """
    Valid flag options below
    """
    if flag == 0:
        unittest.main()
    if flag == 1:
        generate_bc_df("tests/outputs/UNIT_TEST_BC_DF.tsv")
        bc_df_fp = args[0]
    elif flag == 2:
        generate_genes_df("tests/outputs/UNIT_TEST_GENES_DF.tsv")
    else:
        raise Exception("Could not recognize flag number:", flag)

    return None


"""
python -m unittest -v test_module
"""
if __name__ == "__main__":
    unittest.main()

