# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.0.3
#   kernelspec:
#     display_name: Python [conda env:ncbi_remap]
#     language: python
#     name: conda-env-ncbi_remap-py
# ---

# %%
import os
import sys
from pathlib import Path

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import gffutils
from gffutils.pybedtools_integration import to_bedtool
import pybedtools

# Project level imports
from ncbi_remap.notebook import Nb

# %%
db = gffutils.FeatureDB('../references/dm6/r6-11/gtf/dm6_r6-11.gtf.db')

# %%
genes = to_bedtool(db.features_of_type('gene', order_by='start')).sort().saveas()

# %%
genes[0]

# %%

# %%

# %%

# %%

# %%

# %%
pb = pybedtools.BedTool('../../dmel_pacbio/output/pacbio-wf/w1118_testi1/w1118_testi1.collapsed.gff').sort().saveas()

# %%
feature = pb[0]

# %%
near = mRNA.closest(pybedtools.BedTool([feature]).saveas())

# %%
for i in near:
    print(i)

# %%
