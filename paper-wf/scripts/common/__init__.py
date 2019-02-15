"""Shared bits of code."""

from pathlib import Path
from pickle import load

import matplotlib as mpl
import seaborn as sns
from svgutils.compose import Text

from ncbi_remap.config import config, CONFIG_DIR
from ncbi_remap.plotting import add_styles


def label(text):
    return Text(text, 5, 10, size=12, weight="bold", font="sans-serif")


# Setup general plotting info
add_styles(Path(CONFIG_DIR, 'stylelib'))
mpl.style.use(['common', 'paper'])

# make mappers
try:
    with Path('../output/paper-wf/fbgn2symbol.pkl').open('rb') as fh:
        fbgn2symbol = load(fh)

    with Path('../output/paper-wf/symbol2fbgn.pkl').open('rb') as fh:
        symbol2fbgn = load(fh)
except FileNotFoundError:
    print('Please make sure to run the rules fbgn2symbol and symbol2fbgn first.')
    raise
