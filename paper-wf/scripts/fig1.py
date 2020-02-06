import os

import svgutils.compose as svg

from ncbi_remap.plotting.size_conversion import FigSize


def main():
    panel_labels_kwargs = dict(
        size=snakemake.params.label_size,
        weight=snakemake.params.label_weight,
        font=snakemake.params.label_font
    )

    figsize = FigSize(snakemake.params.figsize)

    svg.Figure(figsize.width_cm, figsize.height_cm, 
        svg.Panel(
            svg.SVG(snakemake.input.drawing).scale(2.95).move(0, 0),
            svg.Text("A", 0, 10, **panel_labels_kwargs),
            svg.Text("B", 120, 17, **panel_labels_kwargs),
            svg.Text("C", 120, 75, **panel_labels_kwargs),
            svg.Text("D", 120, 135, **panel_labels_kwargs),
            svg.Grid(10, 10, 2)
        ),
        svg.Panel(
            svg.SVG(snakemake.input.distribution).move(10, 0),
            svg.Text("E", 0, 10, **panel_labels_kwargs),
        ).move(0, 200),
    ).save(snakemake.output[0])




if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(input=dict(
            drawing="../../data/drawings/overview_schematic.svg",
            distribution="../../output/paper-wf/figure_panels/sample_submission.svg",
        ))

    main()
