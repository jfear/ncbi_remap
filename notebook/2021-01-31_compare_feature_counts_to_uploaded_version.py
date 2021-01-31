# %%
from pathlib import Path
import pandas as pd
from pandas.testing import assert_frame_equal

# %%
srxs = Path("../output/rnaseq-counts-wf/test.txt").read_text().strip().splitlines()

# %%
def fix_name(x):
    if "bam" in x:
        return Path(x).stem
    return x

# %%
bad = []
for srx in srxs:
    old = next(Path("/mnt/raid0/Globus/GSE117217_RAW").glob(f"*{srx}.bam.counts.txt.gz"))
    new = Path(f"../output/rnaseq-counts-wf/gene_counts/{srx}.bam.counts.txt.gz")

    try:
        assert_frame_equal(
            pd.read_csv(old, comment="#", sep="\t").rename(fix_name, axis=1),
            pd.read_csv(new, comment="#", sep="\t").rename(fix_name, axis=1)
        )
        print("Same", srx, old, new)
    except AssertionError:
        print("Different", srx, old, new)
        bad.append((old, new))
    except FileNotFoundError:
        print("Missing", srx, old, new)

# %%
bad = []
for srx in srxs:
    old = next(Path("/mnt/raid0/Globus/GSE117217_RAW").glob(f"*{srx}.bam.counts.jcounts.txt.gz"))
    new = Path(f"../output/rnaseq-counts-wf/gene_counts/{srx}.bam.counts.jcounts.txt.gz")

    try:
        assert_frame_equal(
            pd.read_csv(old, comment="#", sep="\t").rename(fix_name, axis=1),
            pd.read_csv(new, comment="#", sep="\t").rename(fix_name, axis=1)
        )
        print("Same", srx, old, new)
    except AssertionError:
        print("Different", srx, old, new)
        bad.append((old, new))
    except FileNotFoundError:
        print("Missing", srx, old, new)


# %%
bad = []
for srx in srxs:
    old = next(Path("/mnt/raid0/Globus/GSE117217_RAW").glob(f"*{srx}.bam.intergenic.counts.txt.gz"))
    new = Path(f"../output/rnaseq-counts-wf/intergenic_counts/{srx}.bam.intergenic.counts.txt.gz")

    try:
        assert_frame_equal(
            pd.read_csv(old, comment="#", sep="\t").rename(fix_name, axis=1),
            pd.read_csv(new, comment="#", sep="\t").rename(fix_name, axis=1)
        )
        print("Same", srx, old, new)
    except AssertionError:
        print("Different", srx, old, new)
        bad.append((old, new))
    except FileNotFoundError:
        print("Missing", srx, old, new)

# %%
bad = []
for srx in srxs:
    old = next(Path("/mnt/raid0/Globus/GSE117217_RAW").glob(f"*{srx}.bam.intergenic.counts.jcounts.txt.gz"))
    new = Path(f"../output/rnaseq-counts-wf/intergenic_counts/{srx}.bam.intergenic.counts.jcounts.txt.gz")

    try:
        assert_frame_equal(
            pd.read_csv(old, comment="#", sep="\t").rename(fix_name, axis=1),
            pd.read_csv(new, comment="#", sep="\t").rename(fix_name, axis=1)
        )
        print("Same", srx, old, new)
    except AssertionError:
        print("Different", srx, old, new)
        bad.append((old, new))
    except FileNotFoundError:
        print("Missing", srx, old, new)
# %%
