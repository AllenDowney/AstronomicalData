"""Copy the solution notebooks into jb/ and prepare them for the book.

For each ../soln/NN_name.ipynb, writes NN_name.ipynb here:
- code cells that start with "# Solution" get the hide-cell tag;
- the page description goes in the notebook metadata (see GreenTeaPress plans/seo.md).
"""

import nbformat

DESCRIPTIONS = {
    '01_query': "Query the Gaia archive from Python: use astroquery to list tables, read their metadata, "
                "and write ADQL queries that download stars as an Astropy Table.",
    '02_coords': "Work with units and celestial coordinates in Astropy, transform them to the GD-1 frame, "
                 "and select the stars in a polygon of sky with an ADQL query.",
    '03_motion': "Transform Gaia positions and proper motions to the GD-1 frame, correct for the Sun's "
                 "reflex motion, and select candidate GD-1 stars with pandas.",
    '04_select': "Select GD-1 candidates by proper motion in the Gaia query: a convex hull, "
                 "a bounding box with BETWEEN, and a polygon test with Matplotlib.",
    '05_join': "Join Gaia data with Pan-STARRS photometry in ADQL, using the best-neighbour "
               "cross-match table, and save the candidate stars with pandas.",
    '06_photo': "Select GD-1 stars by photometry: a color-magnitude diagram, a MIST isochrone, "
                "and a polygon around it that picks out stars of the right age and metallicity.",
    '07_plot': "Make a publication-quality figure with Matplotlib: annotations, rcParams, style sheets, "
               "LaTeX labels, and a multi-panel figure like Figure 1 of the paper.",
}


def prepare(name, description):
    nb = nbformat.read(f'../soln/{name}.ipynb', 4)

    for cell in nb.cells:
        if cell.cell_type == 'code' and cell.source.startswith('# Solution'):
            tags = cell.metadata.setdefault('tags', [])
            if 'hide-cell' not in tags:
                tags.append('hide-cell')

    nb.metadata['myst'] = {'html_meta': {'description': description}}
    nb.metadata['og:description'] = description
    nbformat.write(nb, f'{name}.ipynb')
    print(f'{name} ({len(description)} characters)')


for name, description in DESCRIPTIONS.items():
    prepare(name, description)
