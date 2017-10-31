"""Set of helper functions for the notebook."""
import os
from datetime import datetime
from yaml import load

import numpy as np
import seaborn as sns
from IPython import get_ipython

from .plotting import add_styles


class Nb(object):
    def __init__(self, nb_name=None, project_dir=None, config_dir=None, fig_dir=None, table_dir=None,
                 cache=None, formats=None, styles=None, watermark=None, **kwargs):
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
        fig_dir : str
            Name of the figures directory.
        table_dir : str
            Name of the tables directory.
        cache : str
            Name of the cache directory.
        formats : str or list
            Default list of formats to use for plotting. For example 'png' or ['png', 'svg'].
        styles : str or list
            Default list of matplotlib.style.library to use for plotting. For example 'seaborn-notebook' or ['seaborn-notebook', 'seaborn-paper'].
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
        config_dir : str
            Name of the config directory.
        fig_dir : str
            Name of the figures directory.
        table_dir : str
            Name of the tables directory.
        cache : str
            Name of the cache directory.
        formats : str or list
            Default list of formats to use for plotting. For example 'png' or ['png', 'svg'].
        styles : str or list
            Default list of matplotlib.style.library to use for plotting. For example 'seaborn-notebook' or ['seaborn-notebook', 'seaborn-paper'].
        date : str
            Current date, generated upon creation.

        """
        self.nb_name = nb_name
        self.project_dir = project_dir
        self.config_dir = config_dir
        self.fig_dir = fig_dir
        self.table_dir = table_dir
        self.cache = cache
        self.formats = formats
        self.styles = styles
        self.date = datetime.now().strftime("%Y-%m-%d")

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

        sns.set_context('notebook')

    @classmethod
    def setup_notebook(cls, nb_name=None, config_name='common.yml', watermark=True):
        """Helper function to consistently setup notebooks.

        Functions detects working folder and sets up a ncbi_remap.notebook.Nb
        object with sane defaults.

        Parameters
        ----------
        nb_name : str
            Name of the current notebook.
        config_name : str
            Name of the a YAML config file.
        watermark : bool
            If truen then output watermark information.

        """
        if nb_name is None:
            cache_dir = 'cache'
        else:
            cache_dir = os.path.join('cache', nb_name)

        # Figure out current project and config folder
        prj = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../'))
        cfg = os.path.join(prj, 'config')

        # set defaults
        defaults = {
            'nb_name': nb_name,
            'project_dir': prj,
            'config_dir': cfg,
            'fig_dir': os.path.join(prj, 'output/figures'),
            'table_dir': os.path.join(prj, 'output/tables'),
            'cache': os.path.join(prj, 'output', cache_dir),
            'formats': ['png', 'pdf'],
            'styles': ['notebook', 'paper', 'talk', 'poster'],
            'watermark': watermark
        }


        # Import external config
        fname = os.path.join(cfg, 'common.yml')
        with open(fname) as fh:
            config = load(fh)
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
                'cache', 'formats', 'styles', 'date']
        keys.extend(self._config_attrs)
        res = []
        for key in keys:
            value = self.__dict__[key]
            if value is None:
                value = 'None'
            res.append('{}:\t{}'.format(key, value))

        return '\n'.join(res)
