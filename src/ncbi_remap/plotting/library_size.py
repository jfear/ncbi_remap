"""Distribution sample submssion by year."""
from typing import Union, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs


PLOT_DEFAULTS = dict()
JOINT_AXES_DEFAULTS = dict(ylabel="Duplication Rate (%)", xlabel="Library Size (# Reads)")
LIBSIZE_AXES_DEFAULTS = dict()
DUPLICATION_AXES_DEFAULTS = dict()


def get_data(file_name: str) -> Tuple[pd.Series, pd.Series]:
    return (
        pd.read_table(
            file_name,
            usecols=[
                "srx",
                "date_created",
                "layout",
                "libsize_R1",
                "libsize_R2",
                "pct_duplication",
            ],
        )
        .assign(date_created=lambda x: pd.to_datetime(x.date_created))
        .assign(year=lambda x: x.date_created.dt.year)
        .assign(libsize=lambda x: x.apply(get_libsize, axis=1))
        .assign(log_libsize=lambda x: np.log10(x.libsize))
        .sort_values("libsize", ascending=True)
    )


def get_libsize(x: pd.Series) -> float:
    if x.layout in ["SE", "keep_R1"]:
        return x.libsize_R1
    elif x.layout == "keep_R2":
        return x.libsize_R2
    elif x.layout == "PE":
        try:
            return max(x.libsize_R1, x.libsize_R2)
        except Exception as e:
            print(x)
            raise e


class Plot(NcbiPlotter):
    def __init__(
        self,
        file_name: str,
        plot_kwargs: Union[None, dict] = None,
        joint_ax_kwargs: Union[None, dict] = None,
        libsize_ax_kwargs: Union[None, dict] = None,
        duplication_ax_kwargs: Union[None, dict] = None,
        ax: plt.Axes = None,
    ):
        """
        Example
        ----------
        >>> from ncbi_remap.plotting.library_size import Plot
        >>> plot = Plot("../../../output/paper-wf/srx_prealn.tsv")
        """
        self.update_figsize()
        self.ax_joint = None
        self.ax_marg_libsize = None
        self.ax_marg_duplication = None

        self.plot_kwargs = update_kwargs(PLOT_DEFAULTS, plot_kwargs)

        self.joint_axes_kwargs = update_kwargs(JOINT_AXES_DEFAULTS, joint_ax_kwargs)
        self.libsize_axes_kwargs = update_kwargs(LIBSIZE_AXES_DEFAULTS, libsize_ax_kwargs)
        self.duplication_axes_kwargs = update_kwargs(DUPLICATION_AXES_DEFAULTS, duplication_ax_kwargs)

        self.data = get_data(file_name)
        self.plot()
        self.tweak()

    def plot(self):
        g = sns.jointplot(
            x="log_libsize",
            y="pct_duplication",
            data=self.data,
            kind="kde",
            **self.plot_kwargs,
        ) # type: sns.JointGrid

        self.ax_joint = g.ax_joint
        self.ax_marg_libsize = g.ax_marg_x
        self.ax_marg_duplication = g.ax_marg_y

    def tweak(self):
        self.ax_joint.set(**self.joint_axes_kwargs)
        self.ax_marg_libsize.set(**self.libsize_axes_kwargs)
        self.ax_marg_duplication.set(**self.duplication_axes_kwargs)
        self.ax_joint.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, pos: f"$10^{x:.0f}$"))
