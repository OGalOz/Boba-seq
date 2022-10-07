def get_read_to_seq_dict_from_fa(fa_fp):
    """
    Gets dictionary
        read_name -> sequence
    """

    read2seq = {}
    FH = open(fa_fp, "r")
    read_line = FH.readline()
    line_num = 0
    while read_line != "":
        line_num += 1
        if read_line[0] != ">":
            raise Exception(
                f"Expecting '>' at first pos, " + f"line no. {line_num}, file {fa_fp}"
            )
        read_name = read_line.rstrip()[1:]
        seq = ""
        seq_line = FH.readline()
        while seq_line != "":
            line_num += 1
            seq += seq_line.rstrip()
            seq_line = FH.readline()
            if seq_line[0] == ">":
                read_line = seq_line
                seq_line = ""
        read2seq[read_name] = seq

    print(f"Num of unique inserts: {len(read2seq.keys())}")

    FH.close()

    return read2seq
