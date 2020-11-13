#!/usr/bin/env python
# coding: utf-8

# # Lesson 1

# ## Introduction
# 
# This workshop is an introduction to tools and practices for working with astronomical data.  Topics covered include:
# 
# * Writing queries that select and download data from a database.
# 
# * Using data stored in an Astropy `Table` or Pandas `DataFrame`.
# 
# * Working with coordinates and other quantities with units.
# 
# * Storing data in various formats.
# 
# * Performing database join operations that combine data from multiple tables.
# 
# * Visualizing data and preparing publication-quality figures.

# As a running example, we will replicate part of the analysis in a recent paper, "[Off the beaten path: Gaia reveals GD-1 stars outside of the main stream](https://arxiv.org/abs/1805.00425)" by Adrian M. Price-Whelan and Ana Bonaca.
# 
# As the abstract explains, "Using data from the Gaia second data release combined with Pan-STARRS photometry, we present a sample of highly-probable members of the longest cold stream in the Milky Way, GD-1."
# 
# GD-1 is a [stellar stream](https://en.wikipedia.org/wiki/List_of_stellar_streams), which is "an association of stars orbiting a galaxy that was once a globular cluster or dwarf galaxy that has now been torn apart and stretched out along its orbit by tidal forces."

# [This article in *Science* magazine](https://www.sciencemag.org/news/2018/10/streams-stars-reveal-galaxy-s-violent-history-and-perhaps-its-unseen-dark-matter) explains some of the background, including the process that led to the paper and an discussion of the scientific implications:
# 
# * "The streams are particularly useful for ... galactic archaeology --- rewinding the cosmic clock to reconstruct the assembly of the Milky Way."
# 
# * "They also are being used as exquisitely sensitive scales to measure the galaxy's mass."
# 
# * "... the streams are well-positioned to reveal the presence of dark matter ... because the streams are so fragile, theorists say, collisions with marauding clumps of dark matter could leave telltale scars, potential clues to its nature."

# ## Prerequisites
# 
# This workshop is meant for people who are familiar with basic Python, but not necessarily the libraries we will use, like Astropy or Pandas.  If you are familiar with Python lists and dictionaries, and you know how to write a function that takes parameters and returns a value, you know enough Python for this workshop.
# 
# We assume that you have some familiarity with operating systems, like the ability to use a command-line interface.  But we don't assume you have any prior experience with databases.
# 
# We assume that you are familiar with astronomy at the undergraduate level, but we will not assume specialized knowledge of the datasets or analysis methods we'll use. 

# ## Data
# 
# The datasets we will work with are:
#  
# * [Gaia](https://en.wikipedia.org/wiki/Gaia_(spacecraft)), which is "a space observatory of the European Space Agency (ESA), launched in 2013 ... designed for astrometry: measuring the positions, distances and motions of stars with unprecedented precision", and
# 
# * [Pan-STARRS](https://en.wikipedia.org/wiki/Pan-STARRS), The Panoramic Survey Telescope and Rapid Response System, which is a survey designed to monitor the sky for transient objects, producing a catalog with accurate astronometry and photometry of detected sources.
# 
# Both of these datasets are very large, which can make them challenging to work with.  It might not be possible, or practical, to download the entire dataset.
# One of the goals of this workshop is to provide tools for working with large datasets.

# ## Lesson 1
# 
# The first lesson demonstrates the steps for selecting and downloading data from the Gaia Database:
# 
# 1. First we'll make a connection to the Gaia server,
# 
# 2. We will explore information about the database and the tables it contains,
# 
# 3. We will write a query and send it to the server, and finally
# 
# 4. We will download the response from the server.
# 
# After completing this lesson, you should be able to
# 
# * Compose a basic query in ADQL.
# 
# * Use queries to explore a database and its tables.
# 
# * Use queries to download data.
# 
# * Develop, test, and debug a query incrementally.

# ## Query Language
# 
# In order to select data from a database, you have to compose a query, which is like a program written in a "query language".
# The query language we'll use is ADQL, which stands for "Astronomical Data Query Language".
# 
# ADQL is a dialect of [SQL](https://en.wikipedia.org/wiki/SQL) (Structured Query Language), which is by far the most commonly used query language.  Almost everything you will learn about ADQL also works in SQL.
# 
# [The reference manual for ADQL is here](http://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.html).
# But you might find it easier to learn from [this ADQL Cookbook](https://www.gaia.ac.uk/data/gaia-data-release-1/adql-cookbook).

# ## Installing libraries
# 
# The library we'll use to get Gaia data is [Astroquery](https://astroquery.readthedocs.io/en/latest/).
# 
# If you are running this notebook on Colab, you can run the following cell to install Astroquery and the other libraries we'll use.
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


# ## Connecting to Gaia
# 
# Astroquery provides `Gaia`, which is an [object that represents a connection to the Gaia database](https://astroquery.readthedocs.io/en/latest/gaia/gaia.html).
# 
# We can connect to the Gaia database like this:

# In[2]:


from astroquery.gaia import Gaia


# #### Optional detail 
# 
# > Running this import statement has the effect of creating a [TAP+](http://www.ivoa.net/documents/TAP/) connection; TAP stands for "Table Access Protocol".  It is a network protocol for sending queries to the database and getting back the results.  We're not sure why it seems to create two connections.

# ## Databases and Tables
# 
# What is a database, anyway?  Most generally, it can be any collection of data, but when we are talking about ADQL or SQL:
# 
# * A database is a collection of one or more named tables.
# 
# * Each table is a 2-D array with one or more named columns of data.
# 
# We can use `Gaia.load_tables` to get the names of the tables in the Gaia database.  With the option `only_names=True`, it loads information about the tables, called the "metadata", not the data itself.

# In[3]:


tables = Gaia.load_tables(only_names=True)


# In[4]:


for table in (tables):
    print(table.get_qualified_name())


# So that's a lot of tables.  The ones we'll use are:
# 
# * `gaiadr2.gaia_source`, which contains Gaia data from [data release 2](https://www.cosmos.esa.int/web/gaia/data-release-2),
# 
# * `gaiadr2.panstarrs1_original_valid`, which contains the photometry data we'll use from PanSTARRS, and
# 
# * `gaiadr2.panstarrs1_best_neighbour`, which we'll use to cross-match each star observed by Gaia with the same star observed by PanSTARRS.
# 
# We can use `load_table` (not `load_tables`) to get the metadata for a single table.  The name of this function is misleading, because it only downloads metadata. 

# In[5]:


meta = Gaia.load_table('gaiadr2.gaia_source')
meta


# Jupyter shows that the result is an object of type `TapTableMeta`, but it does not display the contents.
# 
# To see the metadata, we have to print the object.

# In[6]:


print(meta)


# Notice one gotcha: in the list of table names, this table appears as `gaiadr2.gaiadr2.gaia_source`, but when we load the metadata, we refer to it as `gaiadr2.gaia_source`.
# 
# **Exercise:** Go back and try
# 
# ```
# meta = Gaia.load_table('gaiadr2.gaiadr2.gaia_source')
# ```
# 
# What happens?  Is the error message helpful?  If you had not made this error deliberately, would you have been able to figure it out?

# ## Columns
# 
# The following loop prints the names of the columns in the table.

# In[7]:


for column in meta.columns:
    print(column.name)


# You can probably guess what many of these columns are by looking at the names, but you should resist the temptation to guess.
# To find out what the columns mean, [read the documentation](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html).
# 
# If you want to know what can go wrong when you don't read the documentation, [you might like this article](https://www.vox.com/future-perfect/2019/6/4/18650969/married-women-miserable-fake-paul-dolan-happiness).

# **Exercise:** One of the other tables we'll use is `gaiadr2.gaiadr2.panstarrs1_original_valid`.  Use `load_table` to get the metadata for this table.  How many columns are there and what are their names?
# 
# Hint: Remember the gotcha we mentioned earlier.

# In[8]:


# Solution

meta2 = Gaia.load_table('gaiadr2.panstarrs1_original_valid')
print(meta2)


# In[9]:


# Solution

for column in meta2.columns:
    print(column.name)


# ## Writing queries
# 
# By now you might be wondering how we actually download the data.  With tables this big, you generally don't.  Instead, you use queries to select only the data you want.
# 
# A query is a string written in a query language like SQL; for the Gaia database, the query language is a dialect of SQL called ADQL.
# 
# Here's an example of an ADQL query.

# In[10]:


query1 = """SELECT 
TOP 10
source_id, ref_epoch, ra, dec, parallax 
FROM gaiadr2.gaia_source"""


# **Python note:** We use a [triple-quoted string](https://docs.python.org/3/tutorial/introduction.html#strings) here so we can include line breaks in the query, which makes it easier to read.
# 
# The words in uppercase are ADQL keywords:
# 
# * `SELECT` indicates that we are selecting data (as opposed to adding or modifying data).
# 
# * `TOP` indicates that we only want the first 10 rows of the table, which is useful for testing a query before asking for all of the data.
# 
# * `FROM` specifies which table we want data from.
# 
# The third line is a list of column names, indicating which columns we want.  
# 
# In this example, the keywords are capitalized and the column names are lowercase.  This is a common style, but it is not required.  ADQL and SQL are not case-sensitive.

# To run this query, we use the `Gaia` object, which represents our connection to the Gaia database, and invoke `launch_job`:

# In[11]:


job1 = Gaia.launch_job(query1)
job1


# The result is an object that represents the job running on a Gaia server.
# 
# If you print it, it displays metadata for the forthcoming table.

# In[12]:


print(job1)


# Don't worry about `Results: None`.  That does not actually mean there are no results.
# 
# However, `Phase: COMPLETED` indicates that the job is complete, so we can get the results like this:

# In[13]:


results1 = job1.get_results()
type(results1)


# **Optional detail:**  Why is `table` repeated three times?  The first is the name of the module, the second is the name of the submodule, and the third is the name of the class.  Most of the time we only care about the last one.  It's like the Linnean name for gorilla, which is *Gorilla Gorilla Gorilla*.

# The result is an [Astropy Table](https://docs.astropy.org/en/stable/table/), which is similar to a table in an SQL database except:
# 
# * SQL databases are stored on disk drives, so they are persistent; that is, they "survive" even if you turn off the computer.  An Astropy `Table` is stored in memory; it disappears when you turn off the computer (or shut down this Jupyter notebook).
# 
# * SQL databases are designed to process queries.  An Astropy `Table` can perform some query-like operations, like selecting columns and rows.  But these operations use Python syntax, not SQL.
# 
# Jupyter knows how to display the contents of a `Table`.

# In[14]:


results1


# Each column has a name, units, and a data type.
# 
# For example, the units of `ra` and `dec` are degrees, and their data type is `float64`, which is a 64-bit floating-point number, used to store measurements with a fraction part.
# 
# This information comes from the Gaia database, and has been stored in the Astropy `Table` by Astroquery.

# **Exercise:** Read [the documentation of this table](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html) and choose a column that looks interesting to you.  Add the column name to the query and run it again.  What are the units of the column you selected?  What is its data type?

# ## Asynchronous queries
# 
# `launch_job` asks the server to run the job "synchronously", which normally means it runs immediately.  But synchronous jobs are limited to 2000 rows.  For queries that return more rows, you should run "asynchronously", which mean they might take longer to get started.
# 
# If you are not sure how many rows a query will return, you can use the SQL command `COUNT` to find out how many rows are in the result without actually returning them.  We'll see an example of this later.
# 
# The results of an asynchronous query are stored in a file on the server, so you can start a query and come back later to get the results.
# 
# For anonymous users, files are kept for three days.
# 
# As an example, let's try a query that's similar to `query1`, with two changes:
# 
# * It selects the first 3000 rows, so it is bigger than we should run synchronously.
# 
# * It uses a new keyword, `WHERE`.

# In[15]:


query2 = """SELECT TOP 3000
source_id, ref_epoch, ra, dec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1
"""


# A `WHERE` clause indicates which rows we want; in this case, the query selects only rows "where" `parallax` is less than 1.  This has the effect of selecting stars with relatively low parallax, which are farther away.  We'll use this clause to exclude nearby stars that are unlikely to be part of GD-1.
# 
# `WHERE` is one of the most common clauses in ADQL/SQL, and one of the most useful, because it allows us to select only the rows we need from the database.
# 
# We use `launch_job_async` to submit an asynchronous query.

# In[16]:


job2 = Gaia.launch_job_async(query2)
print(job2)


# And here are the results.

# In[17]:


results2 = job2.get_results()
results2


# You might notice that some values of `parallax` are negative.  As [this FAQ explains](https://www.cosmos.esa.int/web/gaia/archive-tips#negative%20parallax), "Negative parallaxes are caused by errors in the observations."  Negative parallaxes have "no physical meaning," but they can be a "useful diagnostic on the quality of the astrometric solution."
# 
# Later we will see an example where we use `parallax` and `parallax_error` to identify stars where the distance estimate is likely to be inaccurate.

# **Exercise:** The clauses in a query have to be in the right order.  Go back and change the order of the clauses in `query2` and run it again.  
# 
# The query should fail, but notice that you don't get much useful debugging information.  
# 
# For this reason, developing and debugging ADQL queries can be really hard.  A few suggestions that might help:
# 
# * Whenever possible, start with a working query, either an example you find online or a query you have used in the past.
# 
# * Make small changes and test each change before you continue.
# 
# * While you are debugging, use `TOP` to limit the number of rows in the result.  That will make each attempt run faster, which reduces your testing time.  
# 
# * Launching test queries synchronously might make them start faster, too.

# ## Operators
# 
# In a `WHERE` clause, you can use any of the [SQL comparison operators](https://www.w3schools.com/sql/sql_operators.asp):
# 
# * `>`: greater than
# * `<`: less than
# * `>=`: greater than or equal
# * `<=`: less than or equal
# * `=`: equal
# * `!=` or `<>`: not equal
# 
# Most of these are the same as Python, but some are not.  In particular, notice that the equality operator is `=`, not `==`.
# Be careful to keep your Python out of your ADQL!
# 
# You can combine comparisons using the logical operators:
# 
# * AND: true if both comparisons are true
# * OR: true if either or both comparisons are true
# 
# Finally, you can use `NOT` to invert the result of a comparison. 

# **Exercise:** [Read about SQL operators here](https://www.w3schools.com/sql/sql_operators.asp) and then modify the previous query to select rows where `bp_rp` is between `-0.75` and `2`.
# 
# You can [read about this variable here](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html).

# In[18]:


# Solution

# This is what most people will probably do

query = """SELECT TOP 10
source_id, ref_epoch, ra, dec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1 
  AND bp_rp > -0.75 AND bp_rp < 2
"""


# In[19]:


# Solution

# But if someone notices the BETWEEN operator, 
# they might do this

query = """SELECT TOP 10
source_id, ref_epoch, ra, dec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1 
  AND bp_rp BETWEEN -0.75 AND 2
"""


# This [Hertzsprung-Russell diagram](https://sci.esa.int/web/gaia/-/60198-gaia-hertzsprung-russell-diagram) shows the BP-RP color and luminosity of stars in the Gaia catalog.
# 
# Selecting stars with `bp-rp` less than 2 excludes many [class M dwarf stars](https://xkcd.com/2360/), which are low temperature, low luminosity.  A star like that at GD-1's distance would be hard to detect, so if it is detected, it it more likely to be in the foreground.

# ## Cleaning up
# 
# Asynchronous jobs have a `jobid`.

# In[20]:


job1.jobid, job2.jobid


# Which you can use to remove the job from the server.

# In[21]:


Gaia.remove_jobs([job2.jobid])


# If you don't remove it job from the server, it will be removed eventually, so don't feel too bad if you don't clean up after yourself.

# ## Formatting queries
# 
# So far the queries have been string "literals", meaning that the entire string is part of the program.
# But writing queries yourself can be slow, repetitive, and error-prone.
# 
# It is often a good idea to write Python code that assembles a query for you.  One useful tool for that is the [string `format` method](https://www.w3schools.com/python/ref_string_format.asp).
# 
# As an example, we'll divide the previous query into two parts; a list of column names and a "base" for the query that contains everything except the column names.
# 
# Here's the list of columns we'll select.  

# In[22]:


columns = 'source_id, ra, dec, pmra, pmdec, parallax, parallax_error, radial_velocity'


# And here's the base; it's a string that contains at least one format specifier in curly brackets (braces).

# In[23]:


query3_base = """SELECT TOP 10 
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2
"""


# This base query contains one format specifier, `{columns}`, which is a placeholder for the list of column names we will provide.
# 
# To assemble the query, we invoke `format` on the base string and provide a keyword argument that assigns a value to `columns`.

# In[24]:


query3 = query3_base.format(columns=columns)


# The result is a string with line breaks.  If you display it, the line breaks appear as `\n`.

# In[25]:


query3


# But if you print it, the line breaks appear as... line breaks.

# In[26]:


print(query3)


# Notice that the format specifier has been replaced with the value of `columns`.
# 
# Let's run it and see if it works:

# In[27]:


job3 = Gaia.launch_job(query3)
print(job3)


# In[28]:


results3 = job3.get_results()
results3


# Good so far.

# **Exercise:** This query always selects sources with `parallax` less than 1.  But suppose you want to take that upper bound as an input.
# 
# Modify `query3_base` to replace `1` with a format specifier like `{max_parallax}`.  Now, when you call `format`, add a keyword argument that assigns a value to `max_parallax`, and confirm that the format specifier gets replaced with the value you provide.

# In[29]:


# Solution

query4_base = """SELECT TOP 10
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < {max_parallax} AND 
bp_rp BETWEEN -0.75 AND 2
"""


# In[30]:


# Solution

query4 = query4_base.format(columns=columns,
                          max_parallax=0.5)
print(query)


# **Style note:**  You might notice that the variable names in this notebook are numbered, like `query1`, `query2`, etc.  
# 
# The advantage of this style is that it isolates each section of the notebook from the others, so if you go back and run the cells out of order, it's less likely that you will get unexpected interactions.
# 
# A drawback of this style is that it can be a nuisance to update the notebook if you add, remove, or reorder a section.
# 
# What do you think of this choice?  Are there alternatives you prefer?

# ## Summary
# 
# This notebook demonstrates the following steps:
# 
# 1. Making a connection to the Gaia server,
# 
# 2. Exploring information about the database and the tables it contains,
# 
# 3. Writing a query and sending it to the server, and finally
# 
# 4. Downloading the response from the server as an Astropy `Table`.

# ## Best practices
# 
# * If you can't download an entire dataset (or it's not practical) use queries to select the data you need.
# 
# * Read the metadata and the documentation to make sure you understand the tables, their columns, and what they mean.
# 
# * Develop queries incrementally: start with something simple, test it, and add a little bit at a time.
# 
# * Use ADQL features like `TOP` and `COUNT` to test before you run a query that might return a lot of data.
# 
# * If you know your query will return fewer than 3000 rows, you can run it synchronously, which might complete faster (but it doesn't seem to make much difference).  If it might return more than 3000 rows, you should run it asynchronously.
# 
# * ADQL and SQL are not case-sensitive, so you don't have to capitalize the keywords, but you should.
# 
# * ADQL and SQL don't require you to break a query into multiple lines, but you should.
# 

# Jupyter notebooks can be good for developing and testing code, but they have some drawbacks.  In particular, if you run the cells out of order, you might find that variables don't have the values you expect.
# 
# There are a few things you can do to mitigate these problems:
# 
# * Make each section of the notebook self-contained.  Try not to use the same variable name in more than one section.
# 
# * Keep notebooks short.  Look for places where you can break your analysis into phases with one notebook per phase.
