#!/usr/bin/env python
# coding: utf-8

# # Lesson 2
# 
# This is the second in a series of lessons related to astronomy data.
# 
# As a running example, we are replicating parts of the analysis in a recent paper, "[Off the beaten path: Gaia reveals GD-1 stars outside of the main stream](https://arxiv.org/abs/1805.00425)" by Adrian M. Price-Whelan and Ana Bonaca.
# 
# In the first notebook, we wrote ADQL queries and used them to select and download data from the Gaia server.
# 
# In this notebook, we'll pick up where we left off and write a query to select stars from the region of the sky where we expect GD-1 to be.

# We'll start with an example that does a "cone search"; that is, it selects stars that appear in a circular region of the sky.
# 
# Then, to select stars in the vicinity of GD-1, we'll:
# 
# * Use `Quantity` objects to represent measurements with units.
# 
# * Use the `Gala` library to convert coordinates from one frame to another.
# 
# * Use the ADQL keywords `POLYGON`, `CONTAINS`, and `POINT` to select stars that fall within a polygonal region.
# 
# * Submit a query and download the results.
# 
# * Store the results in a FITS file.
# 
# After completing this lesson, you should be able to
# 
# * Use Python string formatting to compose more complex ADQL queries.
# 
# * Work with coordinates and other quantities that have units.
# 
# * Download the results of a query and store them in a file.

# ## Installing libraries
# 
# If you are running this notebook on Colab, you can run the following cell to install Astroquery and a the other libraries we'll use.
# 
# If you are running this notebook on your own computer, you might have to install these libraries yourself.  
# 
# If you are using this notebook as part of a Carpentries workshop, you should have received setup instructions.
# 
# TODO: Add a link to the instructions.
# 

# In[1]:


# If we're running on Colab, install libraries

import sys
IN_COLAB = 'google.colab' in sys.modules

if IN_COLAB:
    get_ipython().system('pip install astroquery astro-gala pyia')


# ## Selecting a region

# One of the most common ways to restrict a query is to select stars in a particular region of the sky.
# 
# For example, here's a query from the [Gaia archive documentation](https://gea.esac.esa.int/archive-help/adql/examples/index.html) that selects "all the objects ... in a circular region centered at (266.41683, -29.00781) with a search radius of 5 arcmin (0.08333 deg)."

# In[2]:


query = """
SELECT 
TOP 10 source_id
FROM gaiadr2.gaia_source
WHERE 1=CONTAINS(
  POINT(ra, dec),
  CIRCLE(266.41683, -29.00781, 0.08333333))
"""


# This query uses three keywords that are specific to ADQL (not SQL):
# 
# * `POINT`: a location in [ICRS coordinates](https://en.wikipedia.org/wiki/International_Celestial_Reference_System), specified in degrees of right ascension and declination.
# 
# * `CIRCLE`: a circle where the first two values are the coordinates of the center and the third is the radius in degrees.
# 
# * `CONTAINS`: a function that returns `1` if a `POINT` is contained in a shape and `0` otherwise.
# 
# Here is the [documentation of `CONTAINS`](http://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.html#tth_sEc4.2.12).
# 
# A query like this is called a cone search because it selects stars in a cone.
# 
# Here's how we run it.

# In[3]:


from astroquery.gaia import Gaia

job = Gaia.launch_job(query)
result = job.get_results()
result


# **Exercise:** When you are debugging queries like this, you can use `TOP` to limit the size of the results, but then you still don't know how big the results will be.
# 
# An alternative is to use `COUNT`, which asks for the number of rows that would be selected, but it does not return them.
# 
# In the previous query, replace `TOP 10 source_id` with `COUNT(source_id)` and run the query again.  How many stars has Gaia identified in the cone we searched?

# ## Getting GD-1 Data
# 
# From the Price-Whelan and Bonaca paper, we will try to reproduce Figure 1, which includes this representation of stars likely to belong to GD-1:
# 
# <img src="https://github.com/datacarpentry/astronomy-python/raw/gh-pages/fig/gd1-4.png">

# Along the axis of right ascension ($\phi_1$) the figure extends from -100 to 20 degrees.
# 
# Along the axis of declination ($\phi_2$) the figure extends from about -8 to 4 degrees.
# 
# Ideally, we would select all stars from this rectangle, but there are more than 10 million of them, so
# 
# * That would be difficult to work with,
# 
# * As anonymous users, we are limited to 3 million rows in a single query, and
# 
# * While we are developing and testing code, it will be faster to work with a smaller dataset.
# 
# So we'll start by selecting stars in a smaller rectangle, from -55 to -45 degrees right ascension and -8 to 4 degrees of declination.
# 
# But first we let's see how to represent quantities with units like degrees.

# ## Working with coordinates
# 
# Coordinates are physical quantities, which means that they have two parts, a value and a unit.
# 
# For example, the coordinate $30^{\circ}$ has value 30 and its units are degrees.
# 
# Until recently, most scientific computation was done with values only; units were left out of the program altogether, [often with disastrous results](https://en.wikipedia.org/wiki/Mars_Climate_Orbiter#Cause_of_failure).
# 
# Astropy provides tools for including units explicitly in computations, which makes it possible to detect errors before they cause disasters.
# 
# To use Astropy units, we import them like this:

# In[4]:


import astropy.units as u

u


# `u` is an object that contains most common units and all SI units.
# 
# You can use `dir` to list them, but you should also [read the documentation](https://docs.astropy.org/en/stable/units/).

# In[5]:


dir(u)


# To create a quantity, we multiply a value by a unit.

# In[6]:


coord = 30 * u.deg
type(coord)


# The result is a `Quantity` object.
# 
# Jupyter knows how to display `Quantities` like this:

# In[7]:


coord


# ## Selecting a rectangle
# 
# Now we'll select a rectangle from -55 to -45 degrees right ascension and -8 to 4 degrees of declination.
# 
# We'll define variables to contain these limits.

# In[8]:


phi1_min = -55
phi1_max = -45
phi2_min = -8
phi2_max = 4


# To represent a rectangle, we'll use two lists of coordinates and multiply by their units.

# In[9]:


phi1_rect = [phi1_min, phi1_min, phi1_max, phi1_max] * u.deg
phi2_rect = [phi2_min, phi2_max, phi2_max, phi2_min] * u.deg


# `phi1_rect` and `phi2_rect` represent the coordinates of the corners of a rectangle.  
# 
# But they are in "[a Heliocentric spherical coordinate system defined by the orbit of the GD1 stream](https://gala-astro.readthedocs.io/en/latest/_modules/gala/coordinates/gd1.html)"
# 
# In order to use them in a Gaia query, we have to convert them to [International Celestial Reference System](https://en.wikipedia.org/wiki/International_Celestial_Reference_System) (ICRS) coordinates.  We can do that by storing the coordinates in a `GD1Koposov10` object provided by [Gala](https://gala-astro.readthedocs.io/en/latest/coordinates/).

# In[10]:


import gala.coordinates as gc

corners = gc.GD1Koposov10(phi1=phi1_rect, phi2=phi2_rect)
type(corners)


# We can display the result like this:

# In[11]:


corners


# Now we can use `transform_to` to convert to ICRS coordinates.

# In[12]:


import astropy.coordinates as coord

corners_icrs = corners.transform_to(coord.ICRS)
type(corners_icrs)


# The result is an `ICRS` object.

# In[13]:


corners_icrs


# Notice that a rectangle in one coordinate system is not necessarily a rectangle in another.  In this example, the result is a polygon.

# ## Selecting a polygon
# 
# In order to use this polygon as part of an ADQL query, we have to convert it to a string with a comma-separated list of coordinates, as in this example:
# 
# ```
# """
# POLYGON(143.65, 20.98, 
#         134.46, 26.39, 
#         140.58, 34.85, 
#         150.16, 29.01)
# """
# ```

# `corners_icrs` behaves like a list, so we can use a `for` loop to iterate through the points.

# In[14]:


for point in corners_icrs:
    print(point)


# From that, we can select the coordinates `ra` and `dec`:

# In[15]:


for point in corners_icrs:
    print(point.ra, point.dec)


# The results are quantities with units, but if we select the `value` part, we get a dimensionless floating-point number.

# In[16]:


for point in corners_icrs:
    print(point.ra.value, point.dec.value)


# We can use string `format` to convert these numbers to strings.

# In[17]:


point_base = "{point.ra.value}, {point.dec.value}"

t = [point_base.format(point=point)
     for point in corners_icrs]
t


# The result is a list of strings, which we can join into a single string using `join`.

# In[18]:


point_list = ', '.join(t)
point_list


# Notice that we invoke `join` on a string and pass the list as an argument.
# 
# Before we can assemble the query, we need `columns` again (as we saw in the previous notebook).

# In[19]:


columns = 'source_id, ra, dec, pmra, pmdec, parallax, parallax_error, radial_velocity'


# Here's the base for the query, with format specifiers for `columns` and `point_list`.

# In[20]:


query_base = """SELECT {columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON({point_list}))
"""


# And here's the result:

# In[21]:


query = query_base.format(columns=columns, 
                          point_list=point_list)
print(query)


# As always, we should take a minute to proof-read the query before we launch it.
# 
# The result will be bigger than our previous queries, so it will take a little longer.

# In[22]:


job = Gaia.launch_job_async(query)
print(job)


# Here are the results.

# In[23]:


results = job.get_results()
len(results)


# There are more than 100,000 stars in this polygon, but that's a manageable size to work with.

# ## Saving results
# 
# This is the set of stars we'll work with in the next step.  But since we have a substantial dataset now, this is a good time to save it.
# 
# Storing the data in a file means we can shut down this notebook and pick up where we left off without running the previous query again.
# 
# Astropy `Table` objects provide `write`, which writes the table to disk.

# In[24]:


filename = 'gd1_results.fits'
results.write(filename, overwrite=True)


# Because the filename ends with `fits`, the table is written in the [FITS format](https://en.wikipedia.org/wiki/FITS), which preserves the metadata associated with the table.
# 
# If the file already exists, the `overwrite` argument causes it to be overwritten.
# 
# To see how big the file is, we can use `ls` with the `-lh` option, which prints information about the file including its size in human-readable form.

# In[25]:


get_ipython().system('ls -lh gd1_results.fits')


# The file is about 8.6 MB.  If you are using Windows, `ls` might not work; in that case, try:
# 
# ```
# !dir gd1_results.fits
# ```

# ## Summary
# 
# In this notebook, we composed more complex queries to select stars within a polygonal region of the sky.  Then we downloaded the results and saved them in a FITS file.
# 
# In the next notebook, we'll reload the data from this file and replicate the next step in the analysis, using proper motion to identify stars likely to be in GD-1.

# ## Best practices
# 
# * For measurements with units, use `Quantity` objects that represent units explicitly and check for errors.
# 
# * Use the `format` function to compose queries; it is often faster and less error-prone.
# 
# * Develop queries incrementally: start with something simple, test it, and add a little bit at a time.
# 
# * Once you have a query working, save the data in a local file.  If you shut down the notebook and come back to it later, you can reload the file; you don't have to run the query again.
