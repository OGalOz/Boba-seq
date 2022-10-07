import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import contig_collider

class GenomeViewerNoConditions:
    def __init__(
        self, genes_count_df_fp, BC_loc_df_fp, show_plot=False, pdf_op_dir=None
    ):

        self.show_plot = show_plot
        self.pdf_op_dir = pdf_op_dir
        dtypes_d = {"gene_name": str, "locus_tag": str, "gene_type": str, "contig": str}
        self.genes_count_df = pd.read_table(genes_count_df_fp, dtype=dtypes_d)
        self.genes_count_df["gene_index"] = np.arange(self.genes_count_df.shape[0])
        self.BC_loc_df = pd.read_table(BC_loc_df_fp, dtype=dtypes_d)
        self.BC_loc_df = self.BC_loc_df.drop_duplicates(subset="bc")
        self.BC_loc_df = self.BC_loc_df.rename(
            columns={"tstart": "pos_from", "tend": "pos_to", "target": "contig"}
        )
        self.split_by_contigs()
        self.check_input_tables()

        self.__window_size = 10**4
        self.__gene_y = 18
        self.__gene_x_offset = 200

        self.__cur_gene_index = 0
        self.__cur_contig = None
        self.__score_type = "score_cnnls"
        self.__fr_covered_color = "#00FF00"
        self.__fr_non_covered_color = "#AAAAAA"
        self.__cur_gene_color = "#FF0000"
        self.__gene_color = "#000000"
        self.__gene_score_color = "#FF0000"


    def split_by_contigs(self) -> None:
        """
        Desc:
            We take the two dataframes:
                1. self.BC_loc_df
                2. self.genes_count_df
            and we split them up by contigs.
            First we get the contig names for both.
            Then we check if the contig names match or if they are close.
        """

        gc_df, bc_df = contig_collider.match_contig_names(
            self.genes_count_df, self.BC_loc_df, debug=True
        )
        self.genes_count_df, self.BC_loc_df = gc_df, bc_df
        g_ctgs_new = sorted(self.genes_count_df["contig"].unique())
        bc_ctgs_new = sorted(self.BC_loc_df["contig"].unique())

        g_ctg_d: Dict[str, pd.DataFrame] = {}
        bc_ctg_d: Dict[str, pd.DataFrame] = {}
        for ctg in g_ctgs_new:
            g_ctg_d[ctg] = self.genes_count_df[self.genes_count_df["contig"] == ctg]
            g_ctg_d[ctg]["gene_index"] = np.arange(g_ctg_d[ctg].shape[0])
        for ctg in bc_ctgs_new:
            bc_ctg_d[ctg] = self.BC_loc_df[self.BC_loc_df["contig"] == ctg]
            # bc_ctg_d[ctg]['bc_index'] = np.arange(bc_ctg_d[ctg].shape[0])

        self.g_ctg_d = g_ctg_d
        self.bc_ctg_d = bc_ctg_d
        return None

    def check_input_tables(self):
        # TD
        for x in [
            "gene_index",
            "fragment_count",
            "gene_name",
            "locus_tag",
            "gene_type",
            "contig",
            "pos_from",
            "pos_to",
            "strand",
        ]:
            if x not in self.genes_count_df.columns:
                raise Exception(f"Column {x} missing from genes_count_df")
        for x in [
            "bc",
            "pos_from",
            "pos_to",
            "index",
            "read_name",
            "qlen",
            "qstart",
            "qend",
            "strand",
            "contig",
            "tlen",
            "matches",
            "aln_len",
            "qual",
            "perc_match_cov",
            "perc_match",
            "gene_count",
            "locus_tags",
        ]:
            if x not in self.BC_loc_df.columns:
                raise Exception(f"Column {x} missing from BC_count_df")

        print("Dataframes passed input test")

    def set_gene(
        self, locus_tag=None, index=None, contig=None, name=None
    ) -> pd.DataFrame:
        """
        Currently only allowing access by locus tag OR contig & index
        """
        input_str = f"locus_tag: `{locus_tag}`, index: `{index}`,"
        input_str += f" contig: `{contig}`, name `{name}`"

        genes: pd.DataFrame = None
        ctg_newly_found = False
        if locus_tag is not None:
            genes = self.get_appropriate_genes(locus_tag=locus_tag)
            if genes.shape[0] > 0:
                print("Found matching locus_tag.")
                if contig is None:
                    ctg_newly_found = True
                    contig = genes.iloc[0]["contig"]
                genes = self.get_appropriate_genes(locus_tag=locus_tag, contig=contig)

        if genes is None:
            if index is None and contig is None:
                raise Exception("Could not find appropriate genes for " + input_str)
            else:
                if (index is None and contig is not None) or (
                    index is not None and contig is None
                ):
                    raise Exception(
                        "If index is set, contig must be set, and vice versa."
                    )
                if contig not in self.g_ctg_d:
                    raise Exception(
                        f"Contig {contig} not found in genes contig dataframe."
                    )
                cur_contig_df = self.g_ctg_d[contig]
                if index == cur_contig_df.shape[0]:
                    print("index reached max alloted, resetting to 0")
                    index = 0
                genes = cur_contig_df.iloc[index]
                self.__cur_contig = contig
                self.__cur_gene_index = index
        elif genes.shape[0] > 0:
            self.__cur_contig = genes.iloc[0]["contig"]
            self.__cur_gene_index = genes.iloc[0]["gene_index"]
        else:
            raise RuntimeError(
                f"No genes that match this description were found: {input_str}"
            )


        return genes

    def __filter_range(self, d, pos_from, pos_to) -> pd.DataFrame:
        """
        d is a dataframe.
        pos_from and to are integers
        """
        if pos_from is not None:
            d = d[d.pos_from >= pos_from]
        if pos_to is not None:
            d = d[d.pos_to <= pos_to]
        return d

    def window(self):
        cur_gene = self.current_gene()
        gene_center = cur_gene.pos_from + (cur_gene.pos_to - cur_gene.pos_from) / 2
        # window size is expected to be 10K
        window_from = gene_center - self.__window_size / 2
        window_to = gene_center + self.__window_size / 2
        return (window_from, window_to)

    def set_window_size(self, window_size):
        self.__window_size = window_size

    def zoom_in(self):
        self.__window_size /= 1.2

    def zoom_out(self):
        self.__window_size *= 1.2

    def current_gene(self) -> pd.Series:
        crt_df = self.g_ctg_d[self.__cur_contig]
        single_row_df = crt_df[crt_df["gene_index"] == self.__cur_gene_index]
        res = single_row_df.iloc[0]
        return res

    def next_gene(self):
        self.set_gene(contig=self.__cur_contig, index=self.__cur_gene_index + 1)

    def prev_gene(self):
        self.set_gene(contig=self.__cur_contig, index=self.__cur_gene_index - 1)

    def get_appropriate_genes(
        self, contig=None, locus_tag=None, name=None, pos_from=None, pos_to=None
    ) -> pd.DataFrame:
        #
        """
        gene df requires the following:
            'gene_name', 'locus_tag', 'pos_from', 'pos_to'

        """
        d = self.genes_count_df
        if contig is not None:
            if contig not in self.g_ctg_d:
                raise KeyError(f"Contig {contig} missing from genes dataframe")
            d = self.g_ctg_d[contig]
        if locus_tag is not None:
            d = d[d["locus_tag"] == locus_tag]
        if name is not None:
            d = d[d["gene_name"] == name]

        # below __filter_range might not change anything if pos_from/to are None
        return self.__filter_range(d, pos_from, pos_to)

    def conditions(self, name=None):
        d = self.barseq_layout
        if name is not None:
            d = d[d["name"].str.find(name) != -1]
        return d

    def get_fragments(self, contig=None, pos_from=None, pos_to=None) -> pd.DataFrame:
        if contig is not None:
            if contig not in self.bc_ctg_d:
                raise KeyError(f"Contig {contig} missing from BC dataframe")
            d = self.bc_ctg_d[contig]
        return self.__filter_range(d, pos_from, pos_to)

    @property
    def barseq_layout(self):
        return self.__barseq_layout_df

    def show_next_gene(self):
        self.next_gene()
        self.show()

    def show_prev_gene(self):
        self.prev_gene()
        self.show()

    def show_gene(self, locus_tag=None, name=None):
        self.set_gene(locus_tag=locus_tag, name=name)
        self.show()

    def show_zoom_in(self):
        self.zoom_in()
        self.show()

    def show_zoom_out(self):
        self.zoom_out()
        self.show()

    def show(self):
        cur_gene: pd.Series = self.current_gene()
        print("current gene", cur_gene)

        cur_contig: str = self.__cur_contig

        print("current contig", cur_contig)

        (window_from, window_to) = self.window()

        print("Current window:", (window_from, window_to))

        genes: pd.DataFrame = self.get_appropriate_genes(
            contig=cur_contig, pos_from=window_from, pos_to=window_to
        )

        fragments_df: pd.DataFrame = self.get_fragments(
            contig=cur_contig, pos_from=window_from, pos_to=window_to
        )

        fragments_df = fragments_df.reset_index()

        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_subplot(111)
        if not pd.isna(cur_gene["locus_tag"]):
            ax.set_title(
                "Locus Tag: " + cur_gene["locus_tag"] + ". Contig: " + cur_contig,
                fontsize=15,
            )
        else:
            ax.set_title("Locus Tag: NA. Contig: " + cur_contig, fontsize=15)
        ax.grid(True)

        # Do genes
        # We stagger vertical location of text
        location_marker: int = 0
        for ix, gene in genes.iterrows():
            color = (
                self.__cur_gene_color
                if gene["gene_index"] == cur_gene["gene_index"]
                else self.__gene_color
            )

            if location_marker == 4:
                location_marker = 0
            arrowstyle = "->" if gene.strand == "+" else "<-"

            ax.annotate(
                gene["locus_tag"],
                xy=(gene.pos_from, 14),
                xytext=(gene.pos_from, 14 - (1.8 * location_marker)),
                fontsize=8,
                color=color,
            )

            ax.annotate(
                "",
                xy=(gene.pos_to, 18),
                xytext=(gene.pos_from, 18),
                fontsize=20,
                arrowprops=dict(arrowstyle=arrowstyle, color=color),
            )
            location_marker += 1

        # Do fragments
        for ix, fragment in fragments_df.iterrows():
            color = (
                self.__fr_covered_color
                if fragment.pos_from <= cur_gene.pos_from
                and fragment.pos_to >= cur_gene.pos_to
                else self.__fr_non_covered_color
            )

            ax.annotate(
                "",
                xy=(fragment.pos_to, ix + 20),
                xytext=(fragment.pos_from, ix + 20),
                fontsize=20,
                arrowprops=dict(arrowstyle="-", color=color),
            )

        x_min = genes.pos_from.min() - self.__gene_x_offset
        x_max = genes.pos_to.max() + self.__gene_x_offset
        y_min = fragments_df.shape[0] + 20
        y_max = 0

        ax.axis([x_min, x_max, y_min, y_max])
        plt.ylabel("fragment")
        if self.pdf_op_dir is not None:
            if os.path.exists(self.pdf_op_dir):
                print(cur_gene)
                plot_fn = cur_gene["contig"] + "_" + cur_gene["locus_tag"] + "_plot.pdf"
                op_figure_path = os.path.join(self.pdf_op_dir, plot_fn)
                plt.savefig(op_figure_path)
                print("Wrote PDF file at " + op_figure_path)
        if self.show_plot:
            plt.show()


def main():
    args = sys.argv
    gene_count_df_fp, BC_loc_df, pdf_op_dir = args[1:]
    dbv = GenomeViewerNoConditions(gene_count_df_fp, BC_loc_df, pdf_op_dir=pdf_op_dir)
    dbv.show_gene(locus_tag="HMPREF1079_RS05170")

    ## 1
    dbv.show_next_gene()


if __name__ == "__main__":
    main()
