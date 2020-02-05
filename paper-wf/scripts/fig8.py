"""Figure 8: Data access.
(A) GEO is main source.
(B) FlyBase has aggregated tracks
(C) DroSRA Tool
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
        # GEO is main source.
        svg.Panel(
            svg.SVG(snakemake.input.geo).move(10, 10),
            svg.Text("A", 5, 10, **panel_labels_kwargs),
        ),
        # FlyBase has aggregated tracks
        svg.Panel(
            svg.SVG(snakemake.input.flybase).move(10, 10),
            svg.Text("B", 5, 10, **panel_labels_kwargs),
        ),
        # DroSRA Tool
        svg.Panel(
            svg.SVG(snakemake.input.drosra).move(10, 10),
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
