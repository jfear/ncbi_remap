"""Figure 5: Comparing newly generated data.
(A) Schematic showing how comparison works.
(B) UMAP of comparison.
(C) Table of top 5 by distance.
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
        # Schematic showing how comparison works.
        svg.Panel(
            svg.SVG(snakemake.input.drawing).move(10, 10),
            svg.Text("A", 5, 10, **panel_labels_kwargs),
        ),
        # UMAP of comparison.
        svg.Panel(
            svg.SVG(snakemake.input.umap).move(10, 10),
            svg.Text("B", 5, 10, **panel_labels_kwargs),
        ),
        # Table of top 5 by distance.
        svg.Panel(
            svg.SVG(snakemake.input.table).move(10, 10),
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
