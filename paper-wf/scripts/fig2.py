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
        # UMAP of Strategy (RNA-Seq and ChIP-Seq)
        svg.Panel(
            svg.SVG(snakemake.input.umap_library_strategy_panel).move(10, 10),
            svg.Text("A", 5, 10, **panel_labels_kwargs),
        ),
        # Outlier detection score of library strategy and selection
        svg.Panel(
            svg.SVG(snakemake.input.isolation_score_panel).move(10, 10),
            svg.Text("B", 5, 10, **panel_labels_kwargs),
            svg.Text("C", 260, 10, **panel_labels_kwargs),
        ).move(250, 0),
        # Schematic of the library strategy random forest
        svg.Panel(
            svg.SVG(snakemake.input.drawing).move(10, 10),
            svg.Text("D", 5, 10, **panel_labels_kwargs),
        ).move(0, 190),
        # Feature importances from the random forest
        svg.Panel(
            svg.SVG(snakemake.input.importance).move(10, 10),
            svg.Text("E", 35, 10, **panel_labels_kwargs),
        ).move(480, 190),
        # Panel of features
        svg.Panel(),
        # UMAP of RNA-Seq with library selection
        svg.Panel(),
        # Outlier detection score of updated library strategy and selection
        svg.Panel(),
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
