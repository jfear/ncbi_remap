"""Figure 1."""

from svgutils.compose import Figure, Panel, SVG, Grid
from common import label

Figure(
#     "8.5cm", "8.5cm",
#     "11.6cm", "11.6cm",
#     "17.8cm", "11.6cm",
    720, 720,
    Panel(
        SVG(snakemake.input.pie),
        label('A')
    ),
    Panel(
        SVG(snakemake.input.mapping),
        label('B')
    ).move(0, 142),
    Panel(
        SVG(snakemake.input.libsize),
        label('C')
    ).move(0, 284),
    Panel(
        SVG(snakemake.input.libsize),
        label('C')
    ).move(0, 284),
    Panel(
        SVG(snakemake.input.dups),
        label('D')
    ).move(0, 426),
#     Grid(20, 20)
).save(snakemake.output[0])

# 1-column "8.5cm", "8.5cm",
# 1.5-column "11.6cm", "11.6cm",
# 2-column "17.8cm", "17.8cm",
