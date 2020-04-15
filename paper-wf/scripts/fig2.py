import sys
import svgutils.compose as svg

sys.path.insert(0, "../src")
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
        
        # Schematic
        svg.Panel(
            svg.Text("A", 0, 10, **panel_labels_kwargs),
            # svg.SVG(snakemake.input.drawing).move(10, 0),
        ),

        # UMAP
        svg.Panel(
            svg.Text("B", 0, 10, **panel_labels_kwargs),
            svg.SVG(snakemake.input.umap).scale(.6).move(10, 0),
        ).move(315, 0),

        # Feature Importance
        svg.Panel(
            svg.Text("C", 0, 10, **panel_labels_kwargs),
            svg.SVG(snakemake.input.feature_importance).scale(.6).move(10, 0),
        ).move(0, 345),

        # SHAP Feature Effects
        svg.Panel(
            svg.Text("D", 0, 10, **panel_labels_kwargs),
            # svg.SVG(snakemake.input.shap_interaction_panel).scale(.5).move(10, 0),
        ).move(185, 345),

        # Feature Interaction Panels
        svg.Panel(
            svg.Text("E", 0, 10, **panel_labels_kwargs),
            # svg.SVG(snakemake.input.feature_interaction_panel).scale(.65).move(10, 0),
        ).move(375, 345),

        # svg.Grid(20, 20, 4)

    ).save(snakemake.output[0])


if __name__ == "__main__":
    main()
