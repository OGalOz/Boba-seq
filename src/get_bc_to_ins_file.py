"""
Generates a file with a barcode, then the related insert sequences

bc:
in:
in:
in:
.
.
.
bc:
in:
in:
in:
etc

"""

import os
import sys
import logging
import json


def get_bc_to_ins_file(bc_seq_fa, ins_fq_fp, op_dir):

    bc2read = get_bc_to_read_dict(bc_seq_fa)
    read2seq = get_read_to_seq_dict_from_fq(ins_fq_fp)
    # output_tmp_dicts((op_dir, bc2read, read2seq)
    op_fp = os.path.join(op_dir, "bc_to_seqs.txt")
    op_count_fp = os.path.join(op_dir, "bc2counts.tsv")
    create_bc_to_seqs_file(bc2read, read2seq, op_fp, op_count_fp)
    return None


def get_bc_to_seqs_dict_from_file(inp_fp, print_iter=10**4):
    """
    Args:
        inp_fp: Path to file that looks like
            bc:AC...
            in:CCA...
            in:TG...
        Where there is one barcode line and many following 'insert'
        lines
    """
    FH = open(inp_fp, "r")
    bc_to_ins_d = {}
    line_num = 1
    bc_line = FH.readline().rstrip()
    bc_num = 1
    while bc_line != "":
        if bc_line[0:3] != "bc:":
            raise Exception(f"Expecting barcode at line # {line_num}")
        else:
            bc = bc_line[3:]
        inserts = set()
        ins_line = FH.readline().rstrip()
        while ins_line[:3] == "in:":
            line_num += 1
            insert = ins_line[3:]
            inserts.add(insert)
            ins_line = FH.readline().rstrip()
        bc_to_ins_d[bc] = inserts
        bc_line = ins_line
        line_num += 1
        bc_num += 1

    FH.close()
    return bc_to_ins_d


def create_bc_to_seqs_file(bc2read, read2seq, op_fp, op_count_fp):
    op_FH = open(op_fp, "w")
    op_count_FH = open(op_count_fp, "w")
    for bc in bc2read.keys():
        op_FH.write("bc:" + bc + "\n")
        read_ids = bc2read[bc]
        nR = len(read_ids)
        op_count_FH.write(f"{bc}\t{nR}\n")
        for read_id in read_ids:
            op_FH.write("in:" + read2seq[read_id] + "\n")
    op_FH.close()
    op_count_FH.close()
    print("Wrote out to " + op_fp + " and " + op_count_fp)


def output_tmp_dicts(op_dir, bc2read, read2seq):
    tmp_bc_json = os.path.join(op_dir, "tmp_bc.json")
    tmp_ins_json = os.path.join(op_dir, "tmp_ins.json")
    for x in [[bc2read, tmp_bc_json], [read2seq, tmp_ins_json]]:
        with open(x[1], "w") as g:
            g.write(json.dumps(x[0], indent=2))
        print("Wrote to " + x[1])


def get_bc_to_read_dict(bc_fasta_fp):
    """
    Input is file containing read name then barcode, one line each.
    """
    bc2read = {}

    FH = open(bc_fasta_fp, "r")
    read_line = FH.readline()
    line_num = 0
    while read_line != "":
        line_num += 1
        if read_line[0] != ">":
            raise Exception(f"Expecting '>' at first pos, line no. {line_num}")
        read_name = read_line.rstrip()[1:]
        BC_seq = FH.readline().rstrip()
        if len(BC_seq) != 20:
            raise Exception(
                f"Expecting bc to be length 20, line no. {line_num}, "
                + f"instead length: {len(BC_seq)}. "
                + BC_seq
            )
        if BC_seq in bc2read:
            # print("Found duplicate barcode: " + BC_seq + \
            #        f". Line no. {line_num}")
            bc2read[BC_seq].append(read_name)
        else:
            bc2read[BC_seq] = [read_name]
        read_line = FH.readline()

    print(f"Num of unique barcodes: {len(bc2read.keys())}")
    FH.close()

    return bc2read


def get_read_to_seq_dict_from_fq(fq_fp):
    """
    Gets dictionary
        read_name -> list<sequence>
    """
    read2seq = {}
    FH = open(fq_fp, "r")
    read_line = FH.readline()
    line_num = 1
    while read_line != "":
        if read_line[0] != "@":
            raise Exception(
                f"Expecting '@' at first pos, " + f"line no. {line_num}, file {fq_fp}"
            )
        read_name = read_line.rstrip()[1:]
        seq_line = FH.readline().rstrip()
        spacer_line = FH.readline().rstrip()
        if spacer_line != "+":
            raise Exception(
                f"Expecting '+' at " + f"line no. {line_num + 2}, file {fq_fp}"
            )
        qual_line = FH.readline().rstrip()
        if len(seq_line) != len(qual_line):
            raise Exception(
                "Expecting length of seq and qual to be "
                + f"the same at lines {line_num + 1}, {line_num + 3}"
            )
        if read_name in read2seq:
            raise Exception(
                f"Found duplicate read name in fq file: {fq_fp},"
                + f" line no. {line_num}"
            )
        else:
            read2seq[read_name] = seq_line
        line_num += 4
        read_line = FH.readline()
    FH.close()

    return read2seq


def get_read_to_seq_dict_from_fa(fa_fp):
    """
    Gets dictionary
        read_name -> sequence
    """

    read2seq = {}
    FH = open(fa_fp, "r")
    read_line = FH.readline()
    line_num = 0
    print("Parsing fasta file into dictionary, ", fa_fp)
    while read_line != "":
        line_num += 1
        if read_line[0] != ">":
            raise Exception(
                f"Expecting '>' at first pos, " + f"line no. {line_num}, file {fa_fp}"
            )
        read_name = read_line.rstrip()[1:]
        seq = ""
        seq_line = FH.readline()
        while seq_line != "" and seq_line[0] != ">":
            line_num += 1
            seq += seq_line.rstrip()
            seq_line = FH.readline()
        if seq_line == "":
            read_line = ""
        elif seq_line[0] == ">":
            read_line = seq_line
        read2seq[read_name] = seq

    print(f"Num of unique inserts: {len(read2seq.keys())}")

    FH.close()

    return read2seq


def main():
    args = sys.argv
    help_str = "python3 get_bc_to_ins_file.py bc.fasta ins.fq op_dir 1"

    if args[-1] != "1":
        print("Incorrect Args:")
        print(help_str)
    else:
        fn, bc_fp, ins_fp, op_dir, indic = args
        get_bc_to_ins_file(bc_fp, ins_fp, op_dir)
    return None


if __name__ == "__main__":
    main()
