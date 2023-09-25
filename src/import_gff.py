import sys
import pandas as pd

# Modified from https://github.com/novichkov-lab/DubSeq

def import_gff(gff_path) -> pd.DataFrame:
    # First read all features that has "Parent" property and hash them
#     features, id2features = import_gff_parent(gff_path)

    # First read all features that has "Parent" property and hash them
    id2features = {}
    features=[]
    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or line.startswith("\n"):
                continue                
            vals = line.split("\t")
            # ignore all features that are not gene nor CDS
            if vals[2] != "CDS" and vals[2] != "gene":
                continue
            f_contig = vals[0]
            f_contig = f_contig.split(" ")[0]
            f_pos_from = int(vals[3])
            f_pos_to = int(vals[4])
            f_strand = vals[6]
            f_description = vals[8].strip()

            f_parent = None
            f_name = ""
            f_product = ""
            f_note = ""
            f_pseudo = False
            ## delimiter is " " instead of ";" at times
            if len(f_description.split(";")) > 1:
                for dval in f_description.split(";"):
                    if dval.startswith("Parent="):
                        f_parent = dval[len("Parent=") :].strip()
                    elif dval.startswith("gene="):
                        f_name = dval[len("gene=") :].strip()
                    elif dval.startswith("product="):
                        f_product = dval[len("product=") :].strip()
                    elif dval.startswith("Note="):
                        f_note = dval[len("Note=") :].strip()
                    elif "pseudo=true" in dval:
                        f_pseudo = True
            else:
                print("WARNING: check gff file, unexpected delimiter used for gene annotations.")
            
            if f_parent:
                features = id2features.get(f_parent)
                if not features:
                    features = []
                    id2features[f_parent] = features
                features.append(
                    {
                        "gene_type": vals[2],
                        "gene_name": f_name,
                        "contig": f_contig,
                        "pos_from": f_pos_from,
                        "pos_to": f_pos_to,
                        "strand": f_strand,
                        "pseudo": f_pseudo,
                        "product": f_product,
                    }
                )

    # Now read all "gene" features and collect of children
    genes = []
    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or line.startswith("\n"):
                continue
            vals = line.split("\t")
            if vals[2] == "gene":
                gene_contig = vals[0]
                gene_pos_from = int(vals[3])
                gene_pos_to = int(vals[4])
                gene_strand = vals[6]
                gene_description = vals[8].strip()

                gene_locus_tag = None
                gene_id = None
                if len(vals[8].split(";")) > 1:
                    for term in vals[8].split(";"):
                        (key, value) = term.split("=")
                        if key == "locus_tag":
                            gene_locus_tag = value.strip()
                        elif key == "ID":
                            gene_id = value.strip()
                else:
                    print("WARNING: check gff file, unexpected delimiter used for gene annotations.")

                if not gene_id:
                    continue

                features = id2features.get(gene_id)
                if not features:
                    continue

                # build features related to this gene and locations are correct
                gene_features = []
                for ft in features:
                    if ft["contig"] != gene_contig:
                        continue
                    if ft["strand"] != gene_strand:
                        continue
                    if ft["pos_from"] < gene_pos_from:
                        continue
                    if ft["pos_to"] > gene_pos_to:
                        continue
                    gene_features.append(ft)

                if len(gene_features) == 0:
                    continue

                # if there are more than one feature, check that the type of feature is the same
                gene_types = {}
                for f in gene_features:
                    gene_types[f["gene_type"]] = 1
                if len(gene_types) > 1:
                    raise ValueError(
                        "More than one gene type for a given gene: " + gene_id
                    )

                f = gene_features[0]
                genes.append(
                    {
                        "gene_type": f["gene_type"],
                        "gene_name": f["gene_name"],
                        "locus_tag": gene_locus_tag,
                        "contig": f["contig"],
                        "pos_from": f["pos_from"],
                        "pos_to": f["pos_to"],
                        "strand": f["strand"],
                        "pseudo": f["pseudo"],
                        "product": f["product"],
                    }
                )

    if len(genes)==0:
        return "import_gff failed to import GFF."
    else:
        df_genes=pd.DataFrame(genes)[
            [
                "gene_type",
                "gene_name",
                "locus_tag",
                "contig",
                "pos_from",
                "pos_to",
                "strand",
                "pseudo",
                "product",
            ]
        ]
        df_genes.sort_values(by=["contig","pos_from"], inplace=True)
        return df_genes

def import_gff_alt(gff_path) -> pd.DataFrame:
    # First read all features that has "Parent" property and hash them
#     features, id2features = import_gff_parent(gff_path)
# for gff with CDS, no parent feature, and product is part of gene_description
    # If GFF doesn't have "gene" features and uses "CDS" instead
    CDS = []
    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or line.startswith("\n"):
                continue                
            vals = line.split("\t")
            if vals[2] == "CDS":
                f_contig = vals[0]
                f_contig = f_contig.split(" ")[0]
                f_pos_from = int(vals[3])
                f_pos_to = int(vals[4])
                f_strand = vals[6]
                f_description = vals[8].strip()

                f_name = ""
                f_product = ""
                f_note = ""
                f_locus_tag = ""
                f_pseudo = False
                if len(f_description.split(";")) > 1:
                    for dval in f_description.split(";"):
                        if dval.startswith("Name="): # name of gene, ex. bla
                            f_name = dval[len("Name=") :].strip()
                        elif dval.startswith("ID="): # also different
                            f_locus_tag = dval[len("ID=") :].strip()
                        elif dval.startswith("locus_tag="): # also different
                            f_locus_tag = dval[len("locus_tag=") :].strip()
                        elif dval.startswith("product="):
                            f_product = dval[len("product=") :].strip()
                        elif "pseudo=true" in dval:
                            f_pseudo = True
                else:
                    print("WARNING: check gff file, unexpected delimiter used for gene annotations.")

                    ### YH: add compilation of features here??###
                    
                CDS.append(
                    {
                        "gene_type": vals[2],
                        "gene_name": f_name,
                        "locus_tag": f_locus_tag,
                        "contig": f_contig,
                        "pos_from": f_pos_from,
                        "pos_to": f_pos_to,
                        "strand": f_strand,
                        "pseudo": f_pseudo,
                        "product": f_product,
                    }
                )

    if len(CDS)==0:
        return "import_gff_alt failed to import GFF."
    else:
        df_CDS=pd.DataFrame(CDS)[
            [
                "gene_type",
                "gene_name",
                "locus_tag",
                "contig",
                "pos_from",
                "pos_to",
                "strand",
                "pseudo",
                "product",
            ]
        ]
        df_CDS.sort_values(by=["contig","pos_from"], inplace=True)
        return df_CDS

if __name__ == "__main__":
    _, gff_fp = sys.argv
    import_gff(gff_fp)
#     import_gff_alt(gff_fp)
