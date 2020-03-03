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
            svg.SVG(snakemake.input.distribution).scale(1.18).move(18, 0),
            svg.Text("A", 0, 10, **panel_labels_kwargs),
        ),
        svg.Panel(
            svg.SVG(snakemake.input.drawing).scale(2.80).move(0, 0),
            svg.Text("B", 0, 10, **panel_labels_kwargs),
            svg.Text("C", 123, 22, **panel_labels_kwargs),
            svg.Text("D", 123, 80, **panel_labels_kwargs),
            svg.Text("E", 123, 136, **panel_labels_kwargs),
            svg.Text("F", 123, 193, **panel_labels_kwargs),
        ).move(0, 155),
    ).save(snakemake.output[0])




if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(input=dict(
            drawing="../../data/drawings/overview_schematic.svg",
            distribution="../../output/paper-wf/figure_panels/sample_submission.svg",
        ))

    main()
