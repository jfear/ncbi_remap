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
            svg.SVG(snakemake.input.umap_library_strategy_panel).move(10, 0),
            svg.Text("A", 0, 10, **panel_labels_kwargs),
        ),
        # Outlier detection score of library strategy and selection
        svg.Panel(
            svg.SVG(snakemake.input.isolation_score_panel).move(10, 0),
            svg.Text("B", 0, 10, **panel_labels_kwargs),
            svg.Text("C", 170, 10, **panel_labels_kwargs),
        ).move(175, 0),
        # Schematic of the library strategy random forest
        svg.Panel(
            svg.SVG("../data/drawings/large_16_9.svg").scale(3.78 * 1.4).move(10, 0),
            svg.Text("D", 0, 10, **panel_labels_kwargs),
            # svg.SVG("../data/drawings/small_16_9.svg").move(50, 10),
        ).move(0, 130),
        # Feature importances from the random forest
        svg.Panel(
            svg.SVG(snakemake.input.importance).move(10, 0),
            svg.Text("E", 100, 10, **panel_labels_kwargs),
            svg.Text("F", 180, 10, **panel_labels_kwargs),
        ).move(300, 130),
        # Panel of features
        svg.Panel(
            svg.SVG(snakemake.input.top_features).move(10, 0),
            svg.Text("G", 0, 10, **panel_labels_kwargs),
        ).move(0, 310),
        # UMAP of RNA-Seq with library selection
        svg.Panel(
            svg.SVG("../data/drawings/small_4_3.svg").scale(3.78 * .65).move(10, 0),
            svg.Text("H", 0, 10, **panel_labels_kwargs),

        ).move(330, 310),
        # Outlier detection score of updated library strategy and selection
        svg.Panel(
            svg.SVG("../data/drawings/small_16_9.svg").scale(3.78 * .8).move(10, 0),
            svg.Text("I", 0, 10, **panel_labels_kwargs),
        ).move(330, 390),
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
