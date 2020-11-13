#!/usr/bin/env python
# coding: utf-8

# # Chapter 6
# 
# This is the sixth in a series of notebooks related to astronomy data.
# 
# As a continuing example, we will replicate part of the analysis in a recent paper, "[Off the beaten path: Gaia reveals GD-1 stars outside of the main stream](https://arxiv.org/abs/1805.00425)" by Adrian M. Price-Whelan and Ana Bonaca.
# 
# In the previous lesson we downloaded photometry data from Pan-STARRS, which is available from the same server we've been using to get Gaia data. 
# 
# The next step in the analysis is to select candidate stars based on the photometry data.  The following figure from the paper is a color-magnitude diagram for the stars selected based on proper motion:
# 
# <img width="300" src="https://github.com/datacarpentry/astronomy-python/raw/gh-pages/fig/gd1-3.png">
# 
# In red is a theoretical isochrone, showing where we expect the stars in GD-1 to fall based on the metallicity and age of their original globular cluster. 
# 
# By selecting stars in the shaded area, we can further distinguish the main sequence of GD-1 from younger background stars.

# ## Outline
# 
# Here are the steps in this notebook:
# 
# 1. We'll reload the data from the previous notebook and make a color-magnitude diagram.
# 
# 2. Then we'll specify a polygon in the diagram that contains stars with the photometry we expect.
# 
# 3. Then we'll merge the photometry data with the list of candidate stars, storing the result in a Pandas `DataFrame`.
# 
# After completing this lesson, you should be able to
# 
# * Use Matplotlib to specify a `Polygon` and determine which points fall inside it.
# 
# * Use Pandas to merge data from multiple `DataFrames`, much like a database `JOIN` operation.

# ## Installing libraries
# 
# If you are running this notebook on Colab, you can run the following cell to install Astroquery and a the other libraries we'll use.
# 
# If you are running this notebook on your own computer, you might have to install these libraries yourself.  
# 
# If you are using this notebook as part of a Carpentries workshop, you should have received setup instructions.
# 
# TODO: Add a link to the instructions.

# In[1]:


# If we're running on Colab, install libraries

import sys
IN_COLAB = 'google.colab' in sys.modules

if IN_COLAB:
    get_ipython().system('pip install astroquery astro-gala pyia python-wget')


# ## Reload the data
# 
# The following cell downloads the photometry data we created in the previous notebook.

# In[2]:


import os
from wget import download

filename = 'gd1_photo.fits'
filepath = 'https://github.com/AllenDowney/AstronomicalData/raw/main/data/'

if not os.path.exists(filename):
    print(download(filepath+filename))


# Now we can read the data back into an Astropy `Table`.

# In[3]:


from astropy.table import Table

photo_table = Table.read(filename)


# ## Plotting photometry data
# 
# Now that we have photometry data from Pan-STARRS, we can replicate the [color-magnitude diagram](https://en.wikipedia.org/wiki/Galaxy_color%E2%80%93magnitude_diagram) from the original paper:
# 
# <img width="300" src="https://github.com/datacarpentry/astronomy-python/raw/gh-pages/fig/gd1-3.png">
# 
# The y-axis shows the apparent magnitude of each source with the [g filter](https://en.wikipedia.org/wiki/Photometric_system).
# 
# The x-axis shows the difference in apparent magnitude between the g and i filters, which indicates color.
# 
# Stars with lower values of (g-i) are brighter in g-band than in i-band, compared to other stars, which means they are bluer.
# 
# Stars in the lower-left quadrant of this diagram are less bright and less metallic than the others, which means they are [likely to be older](http://spiff.rit.edu/classes/ladder/lectures/ordinary_stars/ordinary.html).
# 
# Since we expect the stars in GD-1 to be older than the background stars, the stars in the lower-left are more likely to be in GD-1.

# In[4]:


import matplotlib.pyplot as plt

def plot_cmd(table):
    """Plot a color magnitude diagram.
    
    table: Table or DataFrame with photometry data
    """
    y = table['g_mean_psf_mag']
    x = table['g_mean_psf_mag'] - table['i_mean_psf_mag']

    plt.plot(x, y, 'ko', markersize=0.3, alpha=0.3)

    plt.xlim([0, 1.5])
    plt.ylim([14, 22])
    plt.gca().invert_yaxis()

    plt.ylabel('$g_0$')
    plt.xlabel('$(g-i)_0$')


# `plot_cmd` uses a new function, `invert_yaxis`, to invert the `y` axis, which is conventional when plotting magnitudes, since lower magnitude indicates higher brightness.
# 
# `invert_yaxis` is a little different from the other functions we've used.  You can't call it like this:
# 
# ```
# plt.invert_yaxis()          # doesn't work
# ```
# 
# You have to call it like this:
# 
# ```
# plt.gca().invert_yaxis()          # works
# ```
# 
# `gca` stands for "get current axis".  It returns an object that represents the axes of the current figure, and that object provides `invert_yaxis`.
# 
# **In case anyone asks:** The most likely reason for this inconsistency in the interface is that `invert_yaxis` is a lesser-used function, so it's not made available at the top level of the interface.

# Here's what the results look like.

# In[5]:


plot_cmd(photo_table)


# Our figure does not look exactly like the one in the paper because we are working with a smaller region of the sky, so we don't have as many stars.  But we can see an overdense region in the lower left that contains stars with the photometry we expect for GD-1.
# 
# The authors of the original paper derive a detailed polygon that defines a boundary between stars that are likely to be in GD-1 or not.
# 
# As a simplification, we'll choose a boundary by eye that seems to contain the overdense region.

# ## Drawing a polygon
# 
# Matplotlib provides a function called `ginput` that lets us click on the figure and make a list of coordinates.
# 
# It's a little tricky to use `ginput` in a Jupyter notebook.  
# Before calling `plt.ginput` we have to tell Matplotlib to use `TkAgg` to draw the figure in a new window.
# 
# When you run the following cell, a figure should appear in a new window.  Click on it 10 times to draw a polygon around the overdense area.  A red cross should appear where you click.

# In[6]:


import matplotlib as mpl

if IN_COLAB:
    coords = None
else:
    mpl.use('TkAgg')
    plot_cmd(photo_table)
    coords = plt.ginput(10)
    mpl.use('agg')


# The argument to `ginput` is the number of times the user has to click on the figure.
# 
# The result from `ginput` is a list of coordinate pairs.

# In[7]:


coords


# If `ginput` doesn't work for you, you could use the following coordinates.

# In[8]:


if coords is None:
    coords = [(0.2, 17.5), 
              (0.2, 19.5), 
              (0.65, 22),
              (0.75, 21),
              (0.4, 19),
              (0.4, 17.5)]


# The next step is to convert the coordinates to a format we can use to plot them, which is a sequence of `x` coordinates and a sequence of `y` coordinates.  The NumPy function `transpose` does what we want. 

# In[9]:


import numpy as np

xs, ys = np.transpose(coords)
xs, ys


# To display the polygon, we'll draw the figure again and use `plt.plot` to draw the polygon.

# In[10]:


plot_cmd(photo_table)
plt.plot(xs, ys);


# If it looks like your polygon does a good job surrounding the overdense area, go on to the next section.  Otherwise you can try again.
# 
# If you want a polygon with more points (or fewer), you can change the argument to `ginput`.
# 
# The polygon does not have to be "closed".  When we use this polygon in the next section, the last and first points will be connected by a straight line.
# 

# ## Which points are in the polygon?
# 
# Matplotlib provides a `Path` object that we can use to check which points fall in the polygon we selected.
# 
# Here's how we make a `Path` using a list of coordinates.

# In[11]:


from matplotlib.path import Path

path = Path(coords)
path


# `Path` provides `contains_points`, which figures out which points are inside the polygon.
# 
# To test it, we'll create a list with two points, one inside the polygon and one outside.

# In[12]:


points = [(0.4, 20), 
          (0.4, 30)]


# Now we can make sure `contains_points` does what we expect.

# In[13]:


inside = path.contains_points(points)
inside


# The result is an array of Boolean values.
# 
# We are almost ready to select stars whose photometry data falls in this polygon.  But first we need to do some data cleaning.

# ## Reloading the data
# 
# Now we need to combine the photometry data with the list of candidate stars we identified in a previous notebook.  The following cell downloads it:
# 
# 

# In[14]:


import os
from wget import download

filename = 'gd1_candidates.hdf5'
filepath = 'https://github.com/AllenDowney/AstronomicalData/raw/main/data/'

if not os.path.exists(filename):
    print(download(filepath+filename))


# In[15]:


import pandas as pd

candidate_df = pd.read_hdf(filename, 'candidate_df')


# `candidate_df` is the Pandas DataFrame that contains the results from Notebook XX, which selects stars likely to be in GD-1 based on proper motion.  It also includes position and proper motion transformed to the ICRS frame.

# ## Merging photometry data
# 
# Before we select stars based on photometry data, we have to solve two problems:
# 
# 1. We only have Pan-STARRS data for some stars in `candidate_df`.
# 
# 2. Even for the stars where we have Pan-STARRS data in `photo_table`, some photometry data is missing.
# 
# We will solve these problems in two step:
# 
# 1. We'll merge the data from `candidate_df` and `photo_table` into a single Pandas `DataFrame`.
# 
# 2. We'll use Pandas functions to deal with missing data.
# 
# `candidate_df` is already a `DataFrame`, but `results` is an Astropy `Table`.  Let's convert it to Pandas:

# In[16]:


photo_df = photo_table.to_pandas()

for colname in photo_df.columns:
    print(colname)


# Now we want to combine `candidate_df` and `photo_df` into a single table, using `source_id` to match up the rows.
# 
# You might recognize this task; it's the same as the JOIN operation in ADQL/SQL.
# 
# Pandas provides a function called `merge` that does what we want.  Here's how we use it.

# In[17]:


merged = pd.merge(candidate_df, 
                  photo_df, 
                  on='source_id', 
                  how='left')
merged.head()


# The first argument is the "left" table, the second argument is the "right" table, and the keyword argument `on='source_id'` specifies a column to use to match up the rows.
# 
# The argument `how='left'` means that the result should have all rows from the left table, even if some of them don't match up with a row in the right table.
# 
# If you are interested in the other options for `how`, you can [read the documentation of `merge`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.merge.html).
# 
# You can also do different types of join in ADQL/SQL; [you can read about that here](https://www.w3schools.com/sql/sql_join.asp).
# 
# The result is a `DataFrame` that contains the same number of rows as `candidate_df`. 

# In[18]:


len(candidate_df), len(photo_df), len(merged)


# And all columns from both tables.

# In[19]:


for colname in merged.columns:
    print(colname)


# **Detail** You might notice that Pandas also provides a function called `join`; it does almost the same thing, but the interface is slightly different.  We think `merge` is a little easier to use, so that's what we chose.  It's also more consistent with JOIN in SQL, so if you learn how to use `pd.merge`, you are also learning how to use SQL JOIN.
# 
# Also, someone might ask why we have to use Pandas to do this join; why didn't we do it in ADQL.  The answer is that we could have done that, but since we already have the data we need, we should probably do the computation locally rather than make another round trip to the Gaia server.

# ## Missing data
# 
# Let's add columns to the merged table for magnitude and color.

# In[20]:


merged['mag'] = merged['g_mean_psf_mag']
merged['color'] = merged['g_mean_psf_mag'] - merged['i_mean_psf_mag']


# These columns contain the special value `NaN` where we are missing data.
# 
# We can use `notnull` to see which rows contain value data, that is, not null values.

# In[21]:


merged['color'].notnull()


# And `sum` to count the number of valid values.

# In[22]:


merged['color'].notnull().sum()


# For scientific purposes, it's not obvious what we should do with candidate stars if we don't have photometry data.  Should we give them the benefit of the doubt or leave them out?
# 
# In part the answer depends on the goal: are we trying to identify more stars that might be in GD-1, or a smaller set of stars that have higher probability?
# 
# In the next section, we'll leave them out, but you can experiment with the alternative.

# ## Selecting based on photometry
# 
# Now let's see how many of these points are inside the polygon we chose.
# 
# We can use a list of column names to select `color` and `mag`.

# In[23]:


points = merged[['color', 'mag']]
points.head()


# The result is a `DataFrame` that can be treated as a sequence of coordinates, so we can pass it to `contains_points`:

# In[24]:


inside = path.contains_points(points)
inside


# The result is a Boolean array.  We can use `sum` to see how many stars fall in the polygon.

# In[25]:


inside.sum()


# Now we can use `inside` as a mask to select stars that fall inside the polygon.

# In[26]:


selected = merged[inside]


# Let's make a color-magnitude plot one more time, highlighting the selected stars with green `x` marks.

# In[27]:


plot_cmd(photo_table)
plt.plot(xs, ys)

plt.plot(selected['color'], selected['mag'], 'gx');


# It looks like the selected stars are, in fact, inside the polygon, which means they have photometry data consistent with GD-1.
# 
# Finally, we can plot the coordinates of the selected stars:

# In[28]:


plt.figure(figsize=(10,2.5))

x = selected['phi1']
y = selected['phi2']

plt.plot(x, y, 'ko', markersize=0.7, alpha=0.9)

plt.xlabel('ra (degree GD1)')
plt.ylabel('dec (degree GD1)')

plt.axis('equal');


# This example includes two new Matplotlib commands:
# 
# * `figure` creates the figure.  In previous examples, we didn't have to use this function; the figure was created automatically.  But when we call it explicitly, we can provide arguments like `figsize`, which sets the size of the figure.
# 
# * `axis` with the parameter `equal` sets up the axes so a unit is the same size along the `x` and `y` axes.
# 
# In an example like this, where `x` and `y` represent coordinates in space, equal axes ensures that the distance between points is represented accurately.   

# ## Write the data
# 
# Let's write the merged DataFrame to a file.

# In[29]:


filename = 'gd1_merged.hdf5'

merged.to_hdf(filename, 'merged')
selected.to_hdf(filename, 'selected')


# In[30]:


get_ipython().system('ls -lh gd1_merged.hdf5')


# If you are using Windows, `ls` might not work; in that case, try:
# 
# ```
# !dir gd1_merged.hdf5
# ```

# ## Save the polygon
# 
# [Reproducibile research](https://en.wikipedia.org/wiki/Reproducibility#Reproducible_research) is "the idea that ... the full computational environment used to produce the results in the paper such as the code, data, etc. can be used to reproduce the results and create new work based on the research."
# 
# This Jupyter notebook is an example of reproducible research because it contains all of the code needed to reproduce the results, including the database queries that download the data and and analysis.
# 
# However, when we used `ginput` to define a polygon by hand, we introduced a non-reproducible element to the analysis.  If someone running this notebook chooses a different polygon, they will get different results.  So it is important to record the polygon we chose as part of the data analysis pipeline.
# 
# Since `coords` is a NumPy array, we can't use `to_hdf` to save it in a file.  But we can convert it to a Pandas `DataFrame` and save that.
# 
# As an alternative, we could use [PyTables](http://www.pytables.org/index.html), which is the library Pandas uses to read and write files.  It is a powerful library, but not easy to use directly.  So let's take advantage of Pandas.

# In[31]:


coords_df = pd.DataFrame(coords)


# In[32]:


filename = 'gd1_polygon.hdf5'
coords_df.to_hdf(filename, 'coords_df')


# We can read it back like this.

# In[33]:


coords2_df = pd.read_hdf(filename, 'coords_df')
coords2 = coords2_df.to_numpy()


# And verify that the data we read back is the same.

# In[34]:


np.all(coords2 == coords)


# ## Summary
# 
# In this notebook, we worked with two datasets: the list of candidate stars from Gaia and the photometry data from Pan-STARRS.
# 
# We drew a color-magnitude diagram and used it to identify stars we think are likely to be in GD-1.
# 
# Then we used a Pandas `merge` operation to combine the data into a single `DataFrame`.

# ## Best practices
# 
# * If you want to perform something like a database `JOIN` operation with data that is in a Pandas `DataFrame`, you can use the `join` or `merge` function.  In many cases, `merge` is easier to use because the arguments are more like SQL.
# 
# * Use Matplotlib options to control the size and aspect ratio of figures to make them easier to interpret.  In this example, we scaled the axes so the size of a degree is equal along both axes.
# 
# * Matplotlib also provides operations for working with points, polygons, and other geometric entities, so it's not just for making figures.
# 
# * Be sure to record every element of the data analysis pipeline that would be needed to replicate the results.

# In[ ]:




