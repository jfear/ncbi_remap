"""
This script copies all notebooks from the book into the website directory, and
creates pages which wrap them and link together.
"""
import os
import nbformat
import shutil
from pathlib import Path
from generate_contents import gen_contents


PAGEFILE = """title: {title}
url:
save_as: {htmlfile}
Template: {template}

{{% notebook notebooks/{notebook_file} cells[{cells}] %}}
"""


def abspath_from_here(*args):
    here = os.path.dirname(__file__)
    path = os.path.join(here, *args)
    return os.path.abspath(path)


# TODO: Double check paths
NB_SOURCE_DIR = abspath_from_here('..', '..', 'docs')
NB_DEST_DIR = abspath_from_here('..', 'content', 'notebooks')
PAGE_DEST_DIR = abspath_from_here('..', 'content', 'pages')

Path(NB_DEST_DIR).mkdir(parents=True, exist_ok=True)
Path(PAGE_DEST_DIR).mkdir(parents=True, exist_ok=True)


def copy_notebooks():
    nblist = sorted(nb for nb in os.listdir(NB_SOURCE_DIR)
                    if nb.endswith('.ipynb'))
    name_map = {nb: nb.rsplit('.', 1)[0].lower() + '.html'
                for nb in nblist}

    # TODO: Double check paths
    figsource = abspath_from_here('..', '..', 'docs', 'figures')
    figdest = abspath_from_here('..', 'content', 'figures')
    figurelist = []
    figure_map = {}

    if os.path.exists(figdest):
        shutil.rmtree(figdest)

    if os.path.exists(figsource):
        shutil.copytree(figsource, figdest)

        # TODO: Double check paths and EDIT ProjectName
        figurelist = os.listdir(abspath_from_here('..', 'content', 'figures'))
        figure_map = {
            os.path.join('figures', fig):
                os.path.join('/ncbi_remap/figures', fig) for fig in figurelist
        }

    for nb in nblist:
        base, ext = os.path.splitext(nb)
        print('-', nb)

        content = nbformat.read(os.path.join(NB_SOURCE_DIR, nb), as_version=4)

        if nb == 'Index.ipynb':
            cells = '1:'
            template = 'page'
            # TODO: Edit Index title
            title = 'NCBI Remap'
            content.cells[3].source = '\n'.join(gen_contents())
            content.cells[3].cell_type = 'markdown'
        else:
            cells = '1:'
            template = 'page'
            title = content.cells[0].source
            if not title.startswith('#') or len(title.splitlines()) > 1:
                raise ValueError('title not found in third cell')
            title = title.lstrip('#').strip()

        # Replace internal URLs and figure links in notebook
        for cell in content.cells:
            if cell.cell_type == 'markdown':
                for nbname, htmlname in name_map.items():
                    if nbname in cell.source:
                        cell.source = cell.source.replace(nbname, htmlname)
                for figname, newfigname in figure_map.items():
                    if figname in cell.source:
                        cell.source = cell.source.replace(figname, newfigname)
            elif cell.cell_type == 'code':
                cell.source = "#ignore\n" + cell.source

        nbformat.write(content, os.path.join(NB_DEST_DIR, nb))

        pagefile = os.path.join(PAGE_DEST_DIR, base + '.md')
        htmlfile = base.lower() + '.html'
        with open(pagefile, 'w') as f:
            f.write(
                PAGEFILE.format(
                    title=title,
                    htmlfile=htmlfile,
                    notebook_file=nb,
                    template=template,
                    cells=cells
                )
            )


if __name__ == '__main__':
    copy_notebooks()

