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

        # UMAPs
        svg.Panel(
            svg.Text("B", 0, 10, **panel_labels_kwargs),

            # UMAP RNA-Seq Outliers
            svg.SVG(snakemake.input.umap_rnaseq_outliers).scale(.6).move(10, 0),

            # UMAP of Strategy (RNA-Seq, EST, WGS, and ChIP-Seq)
            svg.Panel(
                svg.SVG(snakemake.input.umap_library_strategy),
                svg.Line([(20, 10), (20, 330)], width=1),   # Left Line
                svg.Line([(450, 10), (450, 330)], width=1),  # Right Line
                svg.Line([(20, 10), (450, 10)], width=1),  # Top Line
                svg.Line([(20, 330), (450, 330)], width=1),  # bottom Line
            ).scale(.2).move(28, 125),
        ).move(320, 0),

        # Feature Importance
        svg.Panel(
            svg.Text("C", 0, 10, **panel_labels_kwargs),
            svg.SVG(snakemake.input.feature_importance).scale(.6).move(10, 0),
        ).move(0, 210),

        # SHAP Feature Effects
        svg.Panel(
            svg.Text("D", 0, 10, **panel_labels_kwargs),
            svg.SVG(snakemake.input.shap_interaction_panel).scale(.5).move(10, 0),
        ).move(185, 210),

        # Feature Interaction Panels
        svg.Panel(
            svg.Text("E", 0, 10, **panel_labels_kwargs),
            svg.SVG(snakemake.input.feature_interaction_panel).scale(.65).move(10, 0),
        ).move(375, 210),

        # svg.Grid(10, 10, 2)

    ).save(snakemake.output[0])


if __name__ == "__main__":
    main()
