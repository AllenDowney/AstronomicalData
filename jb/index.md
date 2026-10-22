---
myst:
  html_meta:
    description: "Python tools for astronomical data in seven Jupyter notebooks: query Gaia with ADQL, transform coordinates, and replicate a study of the GD-1 stream."
og:description: "Python tools for astronomical data in seven Jupyter notebooks: query Gaia with ADQL, transform coordinates, and replicate a study of the GD-1 stream."
---

# Astronomical Data in Python: Notebooks

[Allen Downey](https://allendowney.com)

This is self-study material: the workshop is not offered live now, but the notebooks are maintained and run on Colab.

*Astronomical Data in Python* is an introduction to tools and practices for working with astronomical data. Topics covered include:

* Writing queries that select and download data from a database.

* Using data stored in an Astropy `Table` or Pandas `DataFrame`.

* Working with coordinates and other quantities with units.

* Storing data in various formats.

* Performing database join operations that combine data from multiple tables.

* Visualizing data and preparing publication-quality figures.

As a running example, we will replicate part of the analysis in a 2018 paper, "[Off the beaten path: Gaia reveals GD-1 stars outside of the main stream](https://arxiv.org/abs/1805.00425)" by Adrian M. Price-Whelan and Ana Bonaca.

This material was developed in collaboration with [The Carpentries](https://carpentries.org/) and the Astronomy Curriculum Development Committee, and supported by funding from the American Institute of Physics through the American Astronomical Society. It became the Data Carpentry lesson [*Foundations of Astronomical Data Science*](https://datacarpentry.github.io/astronomy-python/), which The Carpentries retired in 2026. The lesson pages remain online, and they include notes for instructors, but these notebooks are the version that is maintained.

I am grateful for contributions from the members of the committee -- Azalee Bostroem, Rodolfo Montez, and Phil Rosenfield -- and from Erin Becker, Brett Morris and Adrian Price-Whelan.

**Level:** intermediate. You should be familiar with basic Python, including lists, dictionaries, and functions, but not necessarily Astropy or Pandas; and with astronomy at the undergraduate level.

**Length:** seven notebooks, about 20 hours in total.

## The notebooks

The pages on this site show the notebooks with solutions. To do the exercises, run the versions with blanks on Colab:

1. [Queries](01_query): connecting to the Gaia archive and writing ADQL queries. [Exercises on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/01_query.ipynb)
2. [Coordinates and units](02_coords): `Quantity` objects, coordinate frames, and selecting stars in a region of the sky. [Exercises on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/02_coords.ipynb)
3. [Proper motion](03_motion): transforming to the GD-1 frame and selecting stars by proper motion with Pandas. [Exercises on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/03_motion.ipynb)
4. [Transformation and selection](04_select): doing the proper-motion selection as part of the query. [Exercises on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/04_select.ipynb)
5. [Joining tables](05_join): combining Gaia data with Pan-STARRS photometry. [Exercises on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/05_join.ipynb)
6. [Photometry](06_photo): a color-magnitude diagram and a second round of selection. [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/06_photo.ipynb)
7. [Visualization](07_plot): a publication-quality figure with several panels. [Exercises on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/07_plot.ipynb)

The notebooks query the Gaia archive, which limits anonymous users to queries of one minute and 50,000 rows; Notebook 1 explains how to register and log in. The code and the data files are in [the GitHub repository](https://github.com/AllenDowney/AstronomicalData), which explains how to run the notebooks on your own computer.

See also [my other workshops](https://allendowney.github.io/#workshops).

The code is available under the [MIT License](https://opensource.org/license/MIT). The text is available under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).
