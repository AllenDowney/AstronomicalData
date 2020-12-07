# Astronomical Data in Python

*Astronomical Data in Python* is an introduction to tools and practices for working with astronomical data.  Topics covered include:

* Writing queries that select and download data from a database.

* Using data stored in an Astropy `Table` or Pandas `DataFrame`.

* Working with coordinates and other quantities with units.

* Storing data in various formats.

* Performing database join operations that combine data from multiple tables.

* Visualizing data and preparing publication-quality figures.

As a running example, we will replicate part of the analysis in a recent paper, "[Off the beaten path: Gaia reveals GD-1 stars outside of the main stream](https://arxiv.org/abs/1805.00425)" by Adrian M. Price-Whelan and Ana Bonaca.

As the abstract explains, "Using data from the Gaia second data release combined with Pan-STARRS photometry, we present a sample of highly-probable members of the longest cold stream in the Milky Way, GD-1."

GD-1 is a [stellar stream](https://en.wikipedia.org/wiki/List_of_stellar_streams), which is "an association of stars orbiting a galaxy that was once a globular cluster or dwarf galaxy that has now been torn apart and stretched out along its orbit by tidal forces."

[This article in *Science* magazine](https://www.sciencemag.org/news/2018/10/streams-stars-reveal-galaxy-s-violent-history-and-perhaps-its-unseen-dark-matter) explains some of the background, including the process that led to the paper and a discussion of the scientific implications:

* "The streams are particularly useful for ... galactic archaeology --- rewinding the cosmic clock to reconstruct the assembly of the Milky Way."

* "They also are being used as exquisitely sensitive scales to measure the galaxy's mass."

* "... the streams are well-positioned to reveal the presence of dark matter ... because the streams are so fragile, theorists say, collisions with marauding clumps of dark matter could leave telltale scars, potential clues to its nature."

## Data

The datasets we will work with are:
 
* [Gaia](https://en.wikipedia.org/wiki/Gaia_(spacecraft)), which is "a space observatory of the European Space Agency (ESA), launched in 2013 ... designed for astrometry: measuring the positions, distances and motions of stars with unprecedented precision", and

* [Pan-STARRS](https://en.wikipedia.org/wiki/Pan-STARRS), The Panoramic Survey Telescope and Rapid Response System, which is a survey designed to monitor the sky for transient objects, producing a catalog with accurate astronometry and photometry of detected sources.

Both of these datasets are very large, which can make them challenging to work with.  It might not be possible, or practical, to download the entire dataset.
One of the goals of this workshop is to provide tools for working with large datasets.

## Prerequisites

These notebooks are meant for people who are familiar with basic Python, but not necessarily the libraries we will use, like Astropy or Pandas.  If you are familiar with Python lists and dictionaries, and you know how to write a function that takes parameters and returns a value, you know enough Python to get started.

We assume that you have some familiarity with operating systems, like the ability to use a command-line interface.  But we don't assume you have any prior experience with databases.

We assume that you are familiar with astronomy at the undergraduate level, but we will not assume specialized knowledge of the datasets or analysis methods we'll use. 

## Acknowledgements

This material was developed in collaboration with [The Carpentries](https://carpentries.org/) and the Astronomy Curriculum Development Committee, and supported by funding from the American Institute of Physics through the American Astronomical Society.

I am grateful for contributions from the members of the committee -- Azalee Bostroem, Rodolfo Montez, and Phil Rosenfield -- and from Erin Becker, Brett Morris and Adrian Price-Whelan.

This material is also available in the form of [Carpentries lessons](https://datacarpentry.github.io/astronomy-python), but you should be
aware that these versions might diverge in the future.