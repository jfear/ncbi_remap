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
            svg.SVG(snakemake.input.drawing).scale(3.5).move(10, 10),
            svg.Text("A", 5, 10, **panel_labels_kwargs),
            svg.Text("B", 154, 27, **panel_labels_kwargs),
            svg.Text("C", 154, 97, **panel_labels_kwargs),
            svg.Text("D", 154, 167, **panel_labels_kwargs),
        ),
        svg.Panel(
            svg.SVG(snakemake.input.distribution).scale(1.3).move(10, 10),
            svg.Text("E", 5, 10, **panel_labels_kwargs),
        ).move(410, 0),
    ).save(snakemake.output[0])




if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(input=dict(
            drawing="../../data/drawings/overview_schematic.svg",
            distribution="../../output/paper-wf/figure_panels/sample_submission.svg",
        ))

    main()
