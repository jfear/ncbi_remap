"""Figure 7: Selection of genes with similar expression.
(A) Schematic showing how selection works.
(B) UMAP of gene space
(C) Table of top 5 genes.
"""
import os

import svgutils.compose as svg

from ncbi_remap.plotting.size_conversion import FigSize


def main():
    panel_labels_kwargs = dict(
        size=snakemake.params.label_size,
        weight=snakemake.params.label_weight,
        font=snakemake.params.label_font,
    )

    figsize = FigSize(snakemake.params.figsize)

    svg.Figure(
        figsize.width_cm,
        figsize.height_cm,
        # Schematic showing how selection works.
        svg.Panel(
            svg.SVG(snakemake.input.drawing).move(10, 10),
            svg.Text("A", 5, 10, **panel_labels_kwargs),
        ),
        # UMAP of gene space
        svg.Panel(
            svg.SVG(snakemake.input.umap).move(10, 10),
            svg.Text("B", 5, 10, **panel_labels_kwargs),
        ),
        # Table of top 5 genes.
        svg.Panel(
            svg.SVG(snakemake.input.tbl).move(10, 10),
            svg.Text("C", 5, 10, **panel_labels_kwargs),
        ),
    ).save(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                drawing="../../data/drawings/overview_schematic.svg",
                distribution="../../output/paper-wf/figure_panels/sample_submission.svg",
            )
        )

    main()
