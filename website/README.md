# Tools for creating http://jfear.github.io/larval_gonad/

## Building the Website

Clone the repository & make sure submodules are included

```
$ git clone https://github.com/jfear/larval_gonad.git
$ cd larval_gonad
$ git submodule update --init --recursive
$ cd website
```

Install the required packages:

```
$ conda create -n larval_gonad --file ../requirements.txt
$ source activate larval_gonad
$ conda install --file ../dev_requirements.txt
```

Update notebooks with navigation information.
```
$ make add_nav
```

Copy the notebook content to the right location (this script also modifies some links for the HTML):
```
$ make copy_notebook
```

Build the html and serve locally:

```
$ make html
$ make serve
$ open http://localhost:8000
```

Deploy to github pages
```
$ make publish-to-github
```

## References

The website is generated using the [Pelican](http://docs.getpelican.com/)
static site generator and [Jupyter Notebooks](http://jupyter.org/). Much of the code for generating this site is adapted
from: [Python Data Science Handbook](https://github.com/jakevdp/PythonDataScienceHandbook) and
[Jake's Blog](https://github.com/jakevdp/jakevdp.github.io-source).
