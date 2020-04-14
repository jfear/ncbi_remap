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
            svg.SVG(snakemake.input[1]).move(420, 0),
            svg.SVG(snakemake.input[2]).move(0, 420),
            svg.SVG(snakemake.input[3]).move(420, 420),
        ).scale(.4)
    ).save(snakemake.output[0])


if __name__ == "__main__":
    main()
