"""Figure 4: Genome annotation with aggregated data.
(A) Annotation of low expression genes.
(B) Tissue specific annotation.
(C) Individual tracks.
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
        # Annotation of low expression genes.
        svg.Panel(
            svg.SVG(snakemake.input.drawing).move(10, 10),
            svg.Text("A", 5, 10, **panel_labels_kwargs),
        ),
        # Tissue specific annotation.
        svg.SVG(snakemake.input.tissue).move(10, 10),
        svg.Text("B", 5, 10, **panel_labels_kwargs),
        # Individual tracks.
        svg.SVG(snakemake.input.individual).move(10, 10),
        svg.Text("C", 5, 10, **panel_labels_kwargs),
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
