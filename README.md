# Bacterial Overexpression Barcoded Shotgun Expression Library Sequencing (Boba-seq)


## Overview
The code in this repository is used to map barcodes to insert sequences and then
to their genomic location (if a reference assembly is provided). 
It was designed to map barcoded expression libraries.

## Dependencies
This program runs on Linux using python3 and uses three programs within it: usearch, vsearch, and minimap2.
You can download them at the following links:
usearch: https://www.drive5.com/usearch/download.html
vsearch: https://github.com/torognes/vsearch
minimap2: https://github.com/lh3/minimap2

* It also requires you to have pandas and matplotlib installed for your python3


## Input Files:
#### Reads FASTQ
This is the FASTQ file containing long read amplicons
#### Genome assembly FNA
This is the file containing the reference genome sequence
#### Genome assembly GFF
This file contains the annotations of the reference assembly - note that contig
names should correspond between this and the genome fna file. 
Example:
```
NZ_JH724188.1	RefSeq	region	1	2641955	.	+	.	ID=NZ_JH724188.1:1..2641955;Dbxref=HMP:1079,taxon:997880;gbkey=Src;genome=genomic;isolation-source=biological product [ENVO:02000043];mol_type=genomic DNA;nat-host=Homo sapiens;strain=CL05T00C42
NZ_JH724188.1	RefSeq	gene	207	698	.	-	.	ID=gene-HMPREF1079_RS21295;Name=HMPREF1079_RS21295;gbkey=Gene;gene_biotype=protein_coding;locus_tag=HMPREF1079_RS21295;old_locus_tag=HMPREF1079_00001
.
.
.
```
#### Oligos FASTA
This file contains the 4 oligos flanking barcode and insert fragment regions
(5' of barcode, 3' of barcode, 5' of insert, 3' of insert).
Could look like this:
```
>5insert
AAAGACACCATGCAC
>3insert
GTGAGCCTCGGTACC
>5bc
GACCTGCAGCGTACG
>3bc
AGAGACCTCGTGGAC
```
#### Config Files
There is one central config file, a JSON file that should be modelled after quik_test.json,
and another configurating JSON file (Default Config) which can be left unmodified, but can have its values changed if
there is a need for new file names. It's up to you to name your config files, but you
must point the central config file to the auxiliary JSON file within the key "default_cfg_path".
Details described below.


## Central Config Description
The central configuration file contains many keys and details.
The following is a description for how to enter each key:
```
    "lib_names": e.g. ["Btheta", "Buniformis"], a list of strings, each being a library name
    "lib_genome_dir": e.g. "$HOME/data/ref_genomes"; a string with a path to a directory containing
    the files listed in the keys "lib_genome_filenames" and "lib_genome_gffs".
    "lib_genome_filenames": e.g. ["Btheta.fasta", "Buniformis.fasta"]; a list of strings representing
    genome fna files within the directory from the key "lib_genome_dir", and the order must be corresponding
    to the lib names within the key "lib_names"

    "lib_genome_gffs": e.g. ["Btheta.gff", "Buniformis.gff"]; a list of strings representing
    genome gff files within the directory from the key "lib_genome_dir", and the order must be corresponding
    to the lib names within the key "lib_names". (GFF files are a table with genes and locations)

    "minimap2_exec_path": e.g. "$HOME/bin/minimap2-2.24_x64-linux/minimap2"; path to minimap 2 executable file.
    "vsearch_exec_path": e.g. "$HOME/.conda/envs/main_env/bin/vsearch"; path to vsearch executable file.
    "metagenome_bool": e.g. false, A bool describing whether this is a metagenomic library or not.

    "minimap_qual_min": e.g. 60; an int representing minimum allowed minimap quality
    "default_cfg_path":  e.g. "$HOME/tests/json/default_cfg.json"; path to default config file (explained below)

    # Below params point to objects. The following are parameters for choosing the best barcode 
    "collapse_params": {
        "min_perc_cov": 0.9, 
        "min_perc_match": 0.9,
        "min_BC_support_ratio": 3,
        "max_frag_diff_to_len": 0.25
    },
    # Below params point to objects. The following are parameters for choosing the best barcode 
    "primer_info": {
        "oligo_db_fp": "oligos_blunt-end.fa",
        "flanking_names": {
            "BC_flanking": {"5": "5bc", 
                            "3": "3bc"},
            "insert_flanking": {"5": "5insert", 
                                "3": "3insert"}
        }
    },

    # Below params define whether you are running a second time, and include
    # path to usearch executable. If you're running on big input FASTQ files
    # and want to save time, replace the value for the key 'files_pre_split' to
    # false.
    "step_1": {
        "files_pre_split": false,
        "usearch_exec_path": "$HOME/bin/usearch",
        "search_pcr2": {
            "fwd": "GTTCTTATCTTTGCAGTCTC",
            "rev": "GAGATTTACGCTTTGGTAAAAGTTGG",
            "maxdiffs": 2,
            "minamp": 300,
            "maxamp": 15000
        },
        "search_oligodb": {
            "maxdiffs": 1,
            "minamp": 300,
            "maxamp": 15000
        },
        "remove_non_concatenated_oligo_ops": false
    },
    "step_2": {
        "oligo_id_cutoff": 90
    },
    "step_3": {
        "max_expected_error": 10
    }
}

```
### Default Config File
The Default Config File ($HOME/tests/json/default_cfg.json)
contains the names of the output folders and files that will be created under
your main output file


## How to Run
The entry point is 
```
src/run_steps.py
```
This file us usually run with the following arguments:
```
python3 src/run_steps.py central_config.json input_dir output_dir {start_point}
```
Where {start_point} is an integer between 1 and 6 - where 1 means start from the beginning,
6 means run from the last step (Plot and stats generation), and every number between starts at that point.
The program normally runs from step 1 through step 6






'''python3 src/run_steps.py quik_test.cfg sample_data/small_test tmp_small 1'''


1. Place all FASTQ library files (inputs) in their own directory 
    'inp_dir' in instruction #4
2. Name a temporary directory for the work to be
    done in, for example tmp_Jan1, this is the 'op_dir' in
    instruction #4. (This directory could either already exist 
    or not, but if it exists then everything
    will be overwritten)
3. Create a config file to have the key lib_names point to
    a list of the prefixes from the FASTQ library files,
    e.g. the 'library names'. Make it match the cfg
    'sample_cfg.json' in terms of keys and structure.
    There are 3 forbidden lib names:
    'concat', 'discard' and 'Logs'.
    This will be 'cfg_fp' in instruction #4.
    If you are running this on metagenomic data, set the
    flag "metagenome_bool" to true, otherwise keep it
    as False.

4. From the central directory, run
    python3 src/run_steps.py cfg_fp inp_dir op_dir {N}
    
    Where {N} is the step number you want to start
    on. To start at the beginning, you would set N to 1,
    for example:

    python3 src/run_steps.py cfg_fp inp_dir op_dir 1

    Suppose you already ran up to step 3 and want to continue
    from step 4, then run 

    python3 src/run_steps.py cfg_fp inp_dir op_dir 4

    where inp_dir and op_dir are the same as before.

Note: If running a second time, use the same directory as before
      as an input (From #1), but this time the directory should
      include a directory within it called "split_files", and that
      should be the one you point to - in this case, set
      the flag (within the config file "step_1" - 
      "file_pre_split" to true

5. If you are NOT running on metagenomic data, then there are a few
    parts to the config file that are very important to get right,
    and the program can't check for you that they are correct, 
    so they might produce incorrect data if not set properly.
    These specifically are the keys you'll need to set correctly:
        "lib_names", "lib_genome_dir", "lib_genome_filenames",
        "lib_genome_gffs"
    "lib_names" is an array e.g. [lib1, lib2, ... lib N],
    lib_genome_dir is the path to a directory that contains all of 
        the files listed in the following keys
    lib_genome_filenames is an array e.g. [fp1, fp2, ..., fpN]
    with the paths to the FNA files (Fasta Nucleotide), where fp1
    is associated with lib1, etc.
    lib_genome_gffs is an array e.g. [fp1, fp2, ..., fpN] with
    the paths to the genome's GFF files, similar to above with the
    library being related to the library in the first array.
    Make sure all input files are readable.



To run the Genome Viewer:
```
    # first import:
    from GenomeViewerNoConditions import GenomeViewerNoconditions
    # then instantiate a version of the class using the genes_count_df, and BC_loc_df
    # If you want to output PDFs, then set the pdf_op_dir argument to a directory path
    # If you want to show plots within the shell, set show_plot arg to True
    dbv = GenomeViewerNoConditions('genes_count_df.tsv', 'BC_loc_df.tsv', show_plot=True, pdf_op_dir="my_pdfs")
    # In order to see a specific gene, use the following command, where you can set gene name/locus_tag
    dbv.show_gene(locus_tag="HMPREF1079_RS10730")
    # The locus_tag is the more unique identifier, so when possible, just use the locus_tag
    # Otherwise, you can use the gene:
    dbv.show_gene(name="folE")
    # To go to the next gene in the genes_count_df, use the following command:
    dbv.show_next_gene()
    # It might be worthwhile to separate the genes_count_df files by contig, and that way
    # when you run show_next_gene, it sticks to a certain contig, otherwise it may shuffle
    # between contigs.
```
