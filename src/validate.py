import os, logging, json
from typing import List, Dict
from contig_collider import special_match_contig_names
from import_gff import DubSeq_import_gff
from get_bc_to_ins_file import get_read_to_seq_dict_from_fa


def load_entire_cfg(cfg_fp):
    with open(cfg_fp, "r") as f:
        cfg_d = json.loads(f.read())
    cfg_d = verify_cfg_d(cfg_d)
    return cfg_d


def verify_cfg_d(cfg_d) -> Dict:
    # Checking the first level of keys
    # Adds default config file to cfg_d
    for x in [
        "lib_names",
        "metagenome_bool",
        "vsearch_exec_path",
        "primer_info",
        "step_1",
        "step_2",
        "step_3",
        "default_cfg_path",
    ]:
        if x not in cfg_d:
            raise Exception(f"Expecting key '{x}' in input config.")

    if not isinstance(cfg_d["metagenome_bool"], bool):
        raise TypeError(
            "'metagenome_bool' key in cfg should be boolean - True or False."
        )
    if not cfg_d["metagenome_bool"]:
        for x in [
            "lib_genome_dir",
            "lib_genome_filenames",
            "lib_genome_gffs",
            "minimap2_exec_path",
            "minimap_qual_min",
        ]:
            if x not in cfg_d:
                raise KeyError(
                    f"If not a metagenome, expecting key {x} in input config."
                )

    # Checking each of the first level of keys
    validate_primer_info(cfg_d["primer_info"])
    validate_step_1(cfg_d["step_1"])
    validate_step_2(cfg_d["step_2"])
    validate_step_3(cfg_d["step_3"])
    validate_lib_names(cfg_d["lib_names"])
    validate_vsearch_exec_path(cfg_d["vsearch_exec_path"])
    if not cfg_d["metagenome_bool"]:
        validate_minimap_and_genomes(cfg_d)
        validate_collapse_params(cfg_d)
    validate_and_add_default_cfg(cfg_d)
    validate_matching_contig_names(cfg_d)
    return cfg_d


def validate_matching_contig_names(cfg_d):
    """
    Desc:
        We check if there's a match between contig names
        for the genome FNA files and the genome GFF files.
        If they don't match in any way, the run will fail.
    """
    files_dir = cfg_d["lib_genome_dir"]
    for i in range(len(cfg_d["lib_names"])):
        lib_name: str = cfg_d["lib_names"][i]
        genome_fna_fp: str = os.path.join(files_dir, cfg_d["lib_genome_filenames"][i])
        gff_fp: str = os.path.join(files_dir, cfg_d["lib_genome_gffs"][i])

        mini_validate_matching_contig_names(lib_name, genome_fna_fp, gff_fp)

    return None


def mini_validate_matching_contig_names(lib_name, genome_fna_fp, gff_fp):

    contig2seq = get_read_to_seq_dict_from_fa(genome_fna_fp)

    fna_contig_names = []
    for ctg in contig2seq:
        new_ctg = ctg.split(" ")[0]
        if new_ctg != ctg:
            print(
                f"WARNING: contig name `{ctg}` from Genome FNA "
                + f"file `{genome_fna_fp}` contains a space."
            )
        fna_contig_names.append(new_ctg)

    gff_df = DubSeq_import_gff(gff_fp)
    gff_contig_names = list(gff_df["contig"].unique())

    if len(gff_contig_names) > len(fna_contig_names):
        raise Exception(
            "Number of contig names in GFF cannot be greater than"
            + " the number of contig names in the Genome FNA file. "
            + "That would imply there are genes in contigs that "
            + "don't exist."
        )

    new_gff_contig_names = []
    for ctg in gff_contig_names:
        new_ctg = ctg.split(" ")[0]
        if new_ctg != ctg:
            print(
                f"WARNING: contig name `{ctg}` from Genome GFF "
                + f"file `{gff_fp}` contains a space."
            )
        new_gff_contig_names.append(new_ctg)

    special_match_contig_names(gff_contig_names, new_gff_contig_names)

    # Succesfully validated contig name matching
    return None


def validate_and_add_default_cfg(cfg_d):
    # ((Note that dicts aren't copied but passed in.))
    path2dflt = cfg_d["default_cfg_path"]
    if not os.path.exists(path2dflt):
        raise ValueError(f"Default config path not found at {path2dflt}")
    with open(path2dflt, "r") as f:
        df_cfg = json.loads(f.read())
    expected_keys = ["steps2dirs", "fns"]
    for x in expected_keys:
        if x not in df_cfg:
            raise KeyError(f"Key {x} missing from default config")
    s2d = df_cfg["steps2dirs"]
    for x in [str(i) for i in range(1, 7)]:
        if x not in s2d:
            raise KeyError(f"Key {x} missing from default config - steps2dirs")
        assert isinstance(s2d[x], str), f"Expecting string for value with key {x}"

    fnd = df_cfg["fns"]
    for x in [str(i) for i in range(1, 6)] + ["P"]:
        if x not in fnd:
            raise KeyError(f"Key {x} missing from default config - 'fns'")
        assert isinstance(
            fnd[x], dict
        ), f"Expecting dict for value with key {x} in 'fns' default."

    cfg_d["d"] = df_cfg


def validate_collapse_params(cfg_d):
    # We return the sub-dictionary for future tests
    if "collapse_params" not in cfg_d:
        raise Exception("config must include key 'collapse_params'")
    cp: Dict = cfg_d["collapse_params"]
    if not isinstance(cp, dict):
        raise TypeError("During collapse step, 'collapse_params' must be a dict.")
    for x in [
        "min_perc_cov",
        "min_perc_match",
        "min_BC_support_ratio",
        "max_frag_diff_to_len",
    ]:
        if x not in cp:
            raise KeyError("'collapse_params' missing param " + x)
        elif not (isinstance(cp[x], float) or isinstance(cp[x], int)):
            raise TypeError(f"Key {x} must be float or int: {cp[x]}")
    return cp


def validate_minimap_and_genomes(cfg_d):
    if "minimap2_exec_path" not in cfg_d:
        raise Exception("Please add 'minimap2_exec_path' to cfg")
    minimap_exec = cfg_d["minimap2_exec_path"]
    if not isinstance(minimap_exec, str):
        raise Exception("minimap2 exec path should be string.")
    if not os.path.exists(minimap_exec):
        raise Exception("minimap2 exec path not found in system.")
    v_exc(minimap_exec)
    if "minimap_qual_min" not in cfg_d:
        raise Exception("'minimap_qual_min' (float or int) missing from config file.")
    if not isinstance(cfg_d["minimap_qual_min"], float) and not isinstance(
        cfg_d["minimap_qual_min"], int
    ):
        raise Exception("'minimap_qual_min' must be float or int.")
    lib_genome_dir = cfg_d["lib_genome_dir"]
    if not isinstance(lib_genome_dir, str):
        raise Exception("genome dir path should be string.")
    if not os.path.exists(lib_genome_dir):
        raise Exception("genome dir path not found in system.")

    # These are FNAs
    genome_filenames = cfg_d["lib_genome_filenames"]
    for gnm_fn in genome_filenames:
        if not isinstance(gnm_fn, str):
            raise Exception("genome filename should be string.")
        fp = os.path.join(lib_genome_dir, gnm_fn)
        if not os.path.exists(fp):
            raise Exception("genome file path not found at " + fp)
        v_read(fp)

    if len(cfg_d["lib_names"]) != len(genome_filenames):
        raise Exception(
            "Expecting the same number of library names as"
            + " genome filenames (one for one)."
        )

    corresponding_str = ""
    for i in range(len(genome_filenames)):
        corresponding_str += f"{cfg_d['lib_names'][i]} <-> {genome_filenames[i]}\n"
    print("These are the corresponding libs to genomes:\n " + corresponding_str)

    # These are genome gffs
    if "lib_genome_gffs" not in cfg_d:
        raise KeyError("'lib_genome_gffs' missing from cfg_d")
    else:
        for x in cfg_d["lib_genome_gffs"]:
            fp = os.path.join(cfg_d["lib_genome_dir"], x)
            if not os.path.exists(fp):
                raise RuntimeError(f"Path to gff at {fp} not found.")
            v_read(fp)


def validate_vsearch_exec_path(inp_str):
    if not isinstance(inp_str, str):
        raise Exception("vsearch exec path must be a string")
    if not os.path.exists(inp_str):
        raise Exception("Could not find path to vsearch executable in OS.")
    # Check if executable permissions:
    v_exc(inp_str)


def validate_lib_names(inp_list):
    if not isinstance(inp_list, list):
        raise Exception(
            "Expecting 'lib_names' key to point to a list, "
            + "instead "
            + str(type(inp_list))
        )
    forbidden_lib_names = ["Logs", "concat", "discard"]
    for lib_name in inp_list:
        if not isinstance(lib_name, str):
            raise Exception("Each component of lib_names should be a string.")

        if lib_name in forbidden_lib_names:
            raise Exception(
                f"lib_name {lib_name} part of "
                + "forbidden_lib_names: "
                + ", ".join(forbidden_lib_names)
            )


def verify_op_dir(op_dir):
    if not os.path.isdir(op_dir):
        raise Exception("Output directory not a directory as expected.")


def validate_primer_info(inp_d):
    test_key = "oligo_db_fp"
    if test_key not in inp_d:
        raise Exception(
            "You must include path to oligos fa file "
            + f"as key '{test_key}' in primer_info section."
        )
    else:
        fp = inp_d[test_key]
        if not isinstance(fp, str):
            raise Exception(
                "Path to oligos fa file must be string, " + "instead " + str(type(fp))
            )
        if not os.path.exists(fp):
            raise Exception(f"Cannot find oligos fa file (key '{test_key}') at " + fp)

    if "flanking_names" not in inp_d:
        raise Exception(
            "You must include labels of flanking primers"
            + " under key 'flanking_names' in primer_info section."
        )
    else:
        fl_d = inp_d["flanking_names"]
        if not isinstance(fl_d, dict):
            raise Exception(
                "Type of object mapped from key 'flanking_names'"
                + " must be a dictionary, instead "
                + str(type(fl_d))
            )
        else:
            for x in ["BC_flanking", "insert_flanking"]:
                if x not in fl_d:
                    raise Exception(x + " must be a key in 'flanking_names'")
                crt = fl_d[x]
                for y in ["5", "3"]:
                    if y not in crt:
                        raise Exception(y + " (string) must be a key in " + x)


def validate_step_1(step_1_d):
    for x in ["usearch_exec_path", "search_pcr2", "remove_non_concatenated_oligo_ops"]:
        if x not in step_1_d:
            raise Exception(f"Expecting key {x} in step_1 config.")

    if not os.path.isfile(step_1_d["usearch_exec_path"]):
        raise Exception(
            "Full path to usearch executable file not"
            f" showing up as file at {step_1_d['usearch_exec_path']}"
        )
    v_exc(step_1_d["usearch_exec_path"])

    if not isinstance(step_1_d["search_pcr2"], dict):
        raise Exception(
            "Expecting dict as search_pcr2 key, instead got "
            + str(type(step_1_d["search_pcr2"]))
        )

    for param in ["remove_non_concatenated_oligo_ops"]:
        if not isinstance(step_1_d[param], bool):
            raise Exception(f"Expecting parameter {param} to be boolean.")


def validate_step_2(inp_d):
    if "oligo_id_cutoff" not in inp_d:
        raise Exception(
            "You must include integer id cutoff (1-99) "
            + "as key 'oligo_id_cutoff' in step 2 section."
        )
    else:
        idc = inp_d["oligo_id_cutoff"]
        if not isinstance(idc, int):
            raise Exception(
                "'oligo_id_cutoff' key should point to an integer,"
                + f" instead {type(idc)}, {idc}"
            )
        else:
            if idc < 0 or idc > 100:
                raise Exception(
                    "oligo_id_cutoff should be between 0 and 100, "
                    + f"instead got {idc}."
                )


def validate_step_3(inp_d):
    if "max_expected_error" not in inp_d:
        raise Exception(
            "You must include integer max_expected_error in step 3 section."
        )
    else:
        idc = inp_d["max_expected_error"]
        if not isinstance(idc, int):
            raise Exception(
                "'max_expected_error' key should point to an integer,"
                + f" instead {type(idc)}, {idc}"
            )
        else:
            if idc < 0:
                raise Exception(
                    "max_expected_error should be a positive integer, "
                    + f"instead got {idc}."
                )


def v_exc(path_to_executable: str):
    executable: bool = os.access(path_to_executable, os.X_OK)
    if not executable:
        raise PermissionError(
            "You do not have permission to execute: " + path_to_executable
        )


def v_read(path_to_readable: str):
    readable: bool = os.access(path_to_readable, os.R_OK)
    if not readable:
        raise PermissionError("You do not have permission to read: " + path_to_readable)
