# Astronomical Data in Python

[![tests](https://github.com/AllenDowney/AstronomicalData/actions/workflows/tests.yml/badge.svg)](https://github.com/AllenDowney/AstronomicalData/actions/workflows/tests.yml)

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

This is one of [my workshops and tutorials](https://allendowney.github.io/#workshops).

**Level:** intermediate. This material should be accessible to people familiar with basic Python, but not necessarily the libraries we will use, like Astropy or Pandas. If you are familiar with Python lists and dictionaries, and you know how to write a function that takes parameters and returns a value, that should be enough. We assume that you are familiar with astronomy at the undergraduate level, but we will not assume specialized knowledge of the datasets or analysis methods we'll use.

**Length:** seven notebooks, about 20 hours in total; it was first taught as five 4-hour sessions.

## Outline

1. **Queries:** connecting to the Gaia archive, exploring the tables it contains, and writing ADQL queries that download the results as an Astropy `Table`.
2. **Coordinates and units:** `Quantity` objects with units, transforming coordinates between frames, and selecting the stars in a polygonal region of the sky; saving the results in a FITS file.
3. **Proper motion:** transforming coordinates and proper motion to the GD-1 frame, working with a Pandas `DataFrame`, and selecting stars whose proper motion is consistent with GD-1.
4. **Transformation and selection:** doing the proper-motion selection as part of the query, so we can search a bigger region of the sky and download less data.
5. **Joining tables:** using `JOIN` to combine Gaia data with Pan-STARRS photometry.
6. **Photometry:** a color-magnitude diagram, a polygon around the isochrone we expect for GD-1, and a second round of selection.
7. **Visualization:** annotations, customization, style sheets, and a figure with multiple panels, like Figure 1 in the paper.

Each notebook comes in two versions: one with blanks for the exercises, and one with solutions.

| Notebook | Exercises | Solutions |
|---|---|---|
| 1. Queries | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/01_query.ipynb) | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/soln/01_query.ipynb) |
| 2. Coordinates and units | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/02_coords.ipynb) | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/soln/02_coords.ipynb) |
| 3. Proper motion | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/03_motion.ipynb) | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/soln/03_motion.ipynb) |
| 4. Transformation and selection | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/04_select.ipynb) | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/soln/04_select.ipynb) |
| 5. Joining tables | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/05_join.ipynb) | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/soln/05_join.ipynb) |
| 6. Photometry | | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/soln/06_photo.ipynb) |
| 7. Visualization | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/07_plot.ipynb) | [Run on Colab](https://colab.research.google.com/github/AllenDowney/AstronomicalData/blob/main/soln/07_plot.ipynb) |

Notebook 6 has no exercises, so it has only one version. You can also [read the notebooks, with solutions, online](https://allendowney.github.io/AstronomicalData/).

### The Gaia archive

The notebooks query the [Gaia archive](https://gea.esac.esa.int/archive/), which limits what anonymous users can do: a query can run for at most one minute and return at most 50,000 rows ([details here](https://www.cosmos.esa.int/web/gaia-users/archive/collaborate#GaiaArchiveUserQuotas)). When the archive is busy, a query can time out; if that happens, run the cell again. Notebook 1 explains how to register and log in, which raises these limits.

The notebooks save their results, and the later notebooks download saved copies from the `data` folder in this repository, so you can start with any notebook.

## Recordings

There are no recordings. The material was first taught as a Carpentries "community test drive", online, November 16–20, 2020, by Allen Downey and Azalee Bostroem. The tag [`carpentries-2020`](https://github.com/AllenDowney/AstronomicalData/tree/carpentries-2020) is the version of the code from that workshop.

## Running the notebooks

The easiest way to run the notebooks is on Colab: click one of the links above. You need a Google account, but you don't have to install anything. Each notebook installs the libraries Colab doesn't have and downloads the files it needs.

To run them on your own computer, you need Python 3 with Jupyter, NumPy, Pandas, PyTables, SciPy, Matplotlib, Astropy, astroquery, and Gala. You can [download the code in this zip file](https://github.com/AllenDowney/AstronomicalData/archive/main.zip), or clone the repository:

```
git clone https://github.com/AllenDowney/AstronomicalData
cd AstronomicalData
```

Then install the packages with pip:

```
pip install -r requirements.txt
```

Or, if you use conda or mamba, create an environment:

```
mamba env create -f environment.yml
mamba activate AstronomicalData
```

Then run `jupyter lab` and open `test_setup.ipynb`, which checks that you have everything you need.

## Next steps

See [my other workshops](https://allendowney.github.io/#workshops).

## License

The code is available under the [MIT License](LICENSE). The text is available under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).
