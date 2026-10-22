import nbformat as nbf
from glob import glob

# Collect a list of all notebooks in the content folder
notebooks = glob("0*.ipynb")

text = '# Solution'
replacement = '# Solution goes here'

colab = '[Click here to run this notebook on Colab]'
colab_url = 'https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/'

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

        # the Colab link points to this version
        if cell['source'].startswith(colab):
            cell['source'] = f'{colab}({colab_url}{ipath}).'

    nbf.write(ntbk, ipath)
