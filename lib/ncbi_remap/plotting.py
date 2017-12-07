import os
from functools import wraps

import matplotlib as mpl
import matplotlib.pyplot as plt

def add_styles(dirname):
    mpl.style.core.USER_LIBRARY_PATHS.append(dirname)
    mpl.style.core.update_user_library(mpl.style.library)


def make_figs(fname=None, styles=None, formats=None):
    def _plot_all(func, styles=styles, formats=formats):
        @wraps(func)
        def wrapper(*args, styles=styles, formats=formats, **kwargs):
            def plot_style(style, formats):
                if isinstance(formats, str):
                    formats = [formats]
                elif formats is None:
                    formats = ['png', 'eps']

                with plt.style.context(style):
                    func(*args,  **kwargs)
                    plt.tight_layout()
                    if (not 'notebook' in style) & (fname is not None):
                        fn = fname + '_' + style
                        for f in formats:
                            plt.savefig('{}.{}'.format(fn, f))
                        plt.close()

            if isinstance(styles, str):
                styles = [styles]
            elif styles is None:
                styles = ['notebook']

            for style in styles:
                plot_style(style, formats)

        return wrapper
    return _plot_all
