import nbformat as nbf
from glob import glob

# Collect a list of all notebooks in the content folder
notebooks = glob("0*.ipynb")

text = '# Solution'
replacement = '# Solution goes here'

# Search through each notebook
for ipath in notebooks:
    ntbk = nbf.read(ipath, nbf.NO_CONVERT)

    # if the first cell is a raw cell, remove it
    first = ntbk.cells[0]
    if first['cell_type'] == 'raw':
        del ntbk.cells[0]

    for cell in ntbk.cells:
        # remove tags
        if 'tags' in cell['metadata']:
            cell['metadata']['tags'] = []

        # remove output
        if 'outputs' in cell:
            cell['outputs'] = []

        # remove solutions
        if cell['source'].startswith(text):
            cell['source'] = replacement

    nbf.write(ntbk, ipath)
