"""Set of helper functions for the notebook."""
import os
from pathlib import Path
from datetime import datetime
from subprocess import check_output

import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns
from IPython import get_ipython

from .config import PROJECT_DIR, CONFIG_DIR, REFERENCES_DIR, config
from .plotting import add_styles


class Nb(object):
    def __init__(self, nb_name=None, project_dir=None, subproject_dir=None,
                 config_dir=None, ref_dir=None, fig_dir=None, table_dir=None,
                 formats=None, styles=None, styles_wide=None, styles_full=None,
                 watermark=None, **kwargs):
        """Helper method for working consistently in notebook.

        Stores a set a bunch of useful attributes. Turns on a bunch of commonly
        used notebook magics. If matplotlib stylelib exists in the config_dir
        then it imports user defined styles.

        Parameters
        ----------
        nb_name : str
            Name of the current notebook.
        project_dir : str
            Name of the project directory.
        config_dir : str
            Name of the config directory.
        ref_dir : str
            Name of the references directory.
        subproject_dir : str
            Name of the subproject directory for placing output.
        fig_dir : str
            Name of the figures directory.
        table_dir : str
            Name of the tables directory.
        formats : str or list
            Default list of formats to use for plotting. For example 'png' or
            ['png', 'svg'].
        styles : str or list
            Default list of matplotlib.style.library to use for plotting. For
            example 'seaborn-notebook' or ['seaborn-notebook',
            'seaborn-paper'].
        watermark : bool
            If true turn on watermarking.
        **kwargs
            Additional arguments that are stored as attributes

        Attributes
        ----------
        nb_name : str
            Name of the current notebook.
        project_dir : str
            Name of the project directory.
        subproject_dir : str
            Directory to save outputs from this subproject.
        config_dir : str
            Name of the config directory.
        ref_dir : str
            Name of the references directory.
        fig_dir : str
            Name of the figures directory.
        table_dir : str
            Name of the tables directory.
        formats : str or list
            Default list of formats to use for plotting. For example 'png' or
            ['png', 'svg'].
        styles : str or list
            Default list of matplotlib.style.library to use for plotting. For
            example 'seaborn-notebook' or ['seaborn-notebook',
            'seaborn-paper'].
        styles_wide : str or list
            Default list of matplotlib.style.library to use for plotting wide
            (two column) images. For example 'seaborn-notebook' or
            ['seaborn-notebook', 'seaborn-paper'].
        styles_full : str or list
            Default list of matplotlib.style.library to use for plotting wide
            (two column) images. For example 'seaborn-notebook' or
            ['seaborn-notebook', 'seaborn-paper'].
        date : str
            Current date, generated upon creation.
        conda_env : str
            Name of the current conda environment location.
        fasta : str
            Path to fasta file.
        chromsizes : str
            Path to chromsizes file.
        gtf : str
            Path to gtf file.
        gtf_db : str
            Path to gtf_db file.
        annot : str
            Path to annot file.
        syn : str
            Path to syn file.

        """
        self.nb_name = nb_name
        self.project_dir = project_dir
        self.config_dir = config_dir
        self.ref_dir = ref_dir
        self.fig_dir = fig_dir
        self.table_dir = table_dir
        self.formats = formats
        self.styles = styles
        self.styles_wide = styles_wide
        self.styles_full = styles_full
        self.date = datetime.now().strftime("%Y-%m-%d")
        self.conda_env = self.get_conda()

        # Add useful reference paths
        assembly = kwargs['assembly']
        tag = kwargs['tag']
        self.fasta = os.path.join(self.ref_dir, assembly, tag, 'fasta',
                                  f'{assembly}_{tag}.fasta')

        self.chromsizes = os.path.join(self.ref_dir, assembly, tag, 'fasta',
                                       f'{assembly}_{tag}.chromsizes')

        self.gtf = os.path.join(self.ref_dir, assembly, tag, 'gtf',
                                f'{assembly}_{tag}.gtf')

        self.gtf_db = os.path.join(self.ref_dir, assembly, tag, 'gtf',
                                   f'{assembly}_{tag}.gtf.db')

        self.annot = os.path.join(self.ref_dir, assembly, tag,
                                  'fb_annotation',
                                  f'{assembly}_{tag}.fb_annotation')

        self.syn = os.path.join(self.ref_dir, assembly, tag,
                                'fb_synonym',
                                f'{assembly}_{tag}.fb_synonym')

        # Add useful mappers
        _annot = pd.read_csv(self.annot, sep='\t', index_col=1)
        self.fbgn2symbol = _annot['gene_symbol'].to_dict()
        self.symbol2fbgn = {v: k for k, v in self.fbgn2symbol.items()}

        try:
            self.fbgn2chrom = pd.read_csv(
                os.path.join(self.project_dir, 'output/fbgn2chrom.tsv'),
                sep='\t', index_col=0)
        except Exception:
            print(
                'Please check output/fbgn2chrom.tsv. '
                'If it does not exist, run bin/fbgn2chrom.py'
            )

        # Add Colors
        colors = kwargs.pop('colors')
        self.colors = colors['cycle']

        self.color_chrom = colors['chrom']

        self.color_female = colors['female']
        self.color_male = colors['male']
        self.colors_sex = sns.color_palette([self.color_female,
                                             self.color_male])

        self.color_c1 = colors['c1']
        self.color_c2 = colors['c2']
        self.color_c3 = colors['c3']

        # Add any key word args
        self._config_attrs = kwargs.keys()
        for k, v in kwargs.items():
            setattr(self, k, v)

        # turn on magics
        self._start_magics(watermark=watermark)

        # Set up plotting
        self._setup_plotting()

        # Turn off scientific notation
        np.set_printoptions(precision=5, suppress=True)

    def _start_magics(self, watermark=None):
        """Start up the notebook magics I commonly use."""
        mgc = get_ipython().magic

        ## Activate the autoreload extension for easy reloading of external packages
        mgc('reload_ext autoreload')
        mgc('autoreload 2')

        ## Trun on the water mark
        if watermark:
            mgc('reload_ext watermark')
            mgc('watermark -u -d -g')

        ## Plot inline
        mgc('matplotlib inline')

    def _setup_plotting(self):
        styles = os.path.join(self.config_dir, 'stylelib')
        if os.path.exists(styles):
            add_styles(styles)

        mpl.style.use(['common', 'notebook'])
        sns.set_palette(self.colors)

    def get_conda(self):
        conda_info = check_output(['conda', 'info']).decode('utf-8')
        for x in conda_info.split('\n'):
            if 'envs directories' in x:
                return x.split(':')[1].strip()

    @classmethod
    def setup_notebook(cls, nb_name=None, subproject_dir=None, watermark=True,
                       **kwargs):
        """Helper function to consistently setup notebooks.

        Functions detects working folder and sets up a ncbi_remap.notebook.Nb
        object with sane defaults.

        Parameters
        ----------
        nb_name : str
            Name of the current notebook.
        subproject_dir : str
            Directory to save outputs from this subproject.
        watermark : bool
            If truen then output watermark information.
        kwargs
            Additional arguments to pass to Nb.

        """
        # set defaults
        defaults = {
            'nb_name': nb_name,
            'project_dir': PROJECT_DIR,
            'subproject_dir': subproject_dir,
            'config_dir': CONFIG_DIR,
            'ref_dir': REFERENCES_DIR,
            'fig_dir': './figures',
            'table_dir': './tables',
            'formats': ['png', 'pdf', 'svg'],
            'styles': ['notebook', 'paper'],
            'styles_wide': ['notebook-wide', 'paper-wide'],
            'styles_full': ['notebook-full', 'paper-full'],
            'watermark': watermark
        }

        defaults.update(kwargs)

        # Import external config
        defaults.update(config)

        return cls(**defaults)

    def fig_name(self, fname):
        if self.nb_name is not None:
            fname = '_'.join([self.nb_name, fname])

        return os.path.join(self.fig_dir, fname)

    def table_name(self, fname):
        if self.nb_name is not None:
            fname = '_'.join([self.nb_name, fname])

        return os.path.join(self.table_dir, fname)

    def __repr__(self):
        return str(self)

    def __str__(self):
        keys = ['nb_name', 'project_dir', 'config_dir', 'fig_dir', 'table_dir',
                'formats', 'styles', 'styles_wide', 'styles_full', 'date']
        keys.extend(self._config_attrs)
        res = []
        for key in keys:
            value = self.__dict__[key]
            if value is None:
                value = 'None'
            res.append('{}:\t{}'.format(key, value))

        return '\n'.join(res)
