import sys
import svgutils.compose as svg

sys.path.insert(0, "../src")
from ncbi_remap.plotting.size_conversion import FigSize


def main():
    figsize = FigSize(snakemake.params.figsize)

    svg.Figure(
        figsize.width_cm,
        figsize.height_cm,
        
        svg.Panel(
            svg.SVG(snakemake.input[0]).move(0, 0),
            svg.SVG(snakemake.input[1]).move(430, 0),
            svg.SVG(snakemake.input[2]).move(0, 360),
            svg.SVG(snakemake.input[3]).move(430, 360),
            svg.SVG(snakemake.input[4]).move(0, 720),
            svg.SVG(snakemake.input[5]).move(430, 720),
        ).scale(.4)
    ).save(snakemake.output[0])


if __name__ == "__main__":
    main()
