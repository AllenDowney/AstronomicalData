# Chapter 1

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

[This article in *Science* magazine](https://www.sciencemag.org/news/2018/10/streams-stars-reveal-galaxy-s-violent-history-and-perhaps-its-unseen-dark-matter) explains some of the background, including the process that led to the paper and an discussion of the scientific implications:

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

## Outline

The first lesson demonstrates the steps for selecting and downloading data from the Gaia Database:

1. First we'll make a connection to the Gaia server,

2. We will explore information about the database and the tables it contains,

3. We will write a query and send it to the server, and finally

4. We will download the response from the server.

After completing this lesson, you should be able to

* Compose a basic query in ADQL.

* Use queries to explore a database and its tables.

* Use queries to download data.

* Develop, test, and debug a query incrementally.

## Query Language

In order to select data from a database, you have to compose a query, which is like a program written in a "query language".
The query language we'll use is ADQL, which stands for "Astronomical Data Query Language".

ADQL is a dialect of [SQL](https://en.wikipedia.org/wiki/SQL) (Structured Query Language), which is by far the most commonly used query language.  Almost everything you will learn about ADQL also works in SQL.

[The reference manual for ADQL is here](http://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.html).
But you might find it easier to learn from [this ADQL Cookbook](https://www.gaia.ac.uk/data/gaia-data-release-1/adql-cookbook).

## Installing libraries

The library we'll use to get Gaia data is [Astroquery](https://astroquery.readthedocs.io/en/latest/).

If you are running this notebook on Colab, you can run the following cell to install Astroquery and the other libraries we'll use.

If you are running this notebook on your own computer, you might have to install these libraries yourself.  

If you are using this notebook as part of a Carpentries workshop, you should have received setup instructions.

TODO: Add a link to the instructions.



```python
# If we're running on Colab, install libraries

import sys
IN_COLAB = 'google.colab' in sys.modules

if IN_COLAB:
    !pip install astroquery astro-gala pyia
```

## Connecting to Gaia

Astroquery provides `Gaia`, which is an [object that represents a connection to the Gaia database](https://astroquery.readthedocs.io/en/latest/gaia/gaia.html).

We can connect to the Gaia database like this:


```python
from astroquery.gaia import Gaia
```

    Created TAP+ (v1.2.1) - Connection:
    	Host: gea.esac.esa.int
    	Use HTTPS: True
    	Port: 443
    	SSL Port: 443
    Created TAP+ (v1.2.1) - Connection:
    	Host: geadata.esac.esa.int
    	Use HTTPS: True
    	Port: 443
    	SSL Port: 443


Running this import statement has the effect of creating a [TAP+](http://www.ivoa.net/documents/TAP/) connection; TAP stands for "Table Access Protocol".  It is a network protocol for sending queries to the database and getting back the results.  We're not sure why it seems to create two connections.

## Databases and Tables

What is a database, anyway?  Most generally, it can be any collection of data, but when we are talking about ADQL or SQL:

* A database is a collection of one or more named tables.

* Each table is a 2-D array with one or more named columns of data.

We can use `Gaia.load_tables` to get the names of the tables in the Gaia database.  With the option `only_names=True`, it loads information about the tables, called the "metadata", not the data itself.


```python
tables = Gaia.load_tables(only_names=True)
```

    INFO: Retrieving tables... [astroquery.utils.tap.core]
    INFO: Parsing tables... [astroquery.utils.tap.core]
    INFO: Done. [astroquery.utils.tap.core]



```python
for table in (tables):
    print(table.get_qualified_name())
```

    external.external.apassdr9
    external.external.gaiadr2_geometric_distance
    external.external.galex_ais
    external.external.ravedr5_com
    external.external.ravedr5_dr5
    external.external.ravedr5_gra
    external.external.ravedr5_on
    external.external.sdssdr13_photoprimary
    external.external.skymapperdr1_master
    external.external.tmass_xsc
    public.public.hipparcos
    public.public.hipparcos_newreduction
    public.public.hubble_sc
    public.public.igsl_source
    public.public.igsl_source_catalog_ids
    public.public.tycho2
    public.public.dual
    tap_config.tap_config.coord_sys
    tap_config.tap_config.properties
    tap_schema.tap_schema.columns
    tap_schema.tap_schema.key_columns
    tap_schema.tap_schema.keys
    tap_schema.tap_schema.schemas
    tap_schema.tap_schema.tables
    gaiadr1.gaiadr1.aux_qso_icrf2_match
    gaiadr1.gaiadr1.ext_phot_zero_point
    gaiadr1.gaiadr1.allwise_best_neighbour
    gaiadr1.gaiadr1.allwise_neighbourhood
    gaiadr1.gaiadr1.gsc23_best_neighbour
    gaiadr1.gaiadr1.gsc23_neighbourhood
    gaiadr1.gaiadr1.ppmxl_best_neighbour
    gaiadr1.gaiadr1.ppmxl_neighbourhood
    gaiadr1.gaiadr1.sdss_dr9_best_neighbour
    gaiadr1.gaiadr1.sdss_dr9_neighbourhood
    gaiadr1.gaiadr1.tmass_best_neighbour
    gaiadr1.gaiadr1.tmass_neighbourhood
    gaiadr1.gaiadr1.ucac4_best_neighbour
    gaiadr1.gaiadr1.ucac4_neighbourhood
    gaiadr1.gaiadr1.urat1_best_neighbour
    gaiadr1.gaiadr1.urat1_neighbourhood
    gaiadr1.gaiadr1.cepheid
    gaiadr1.gaiadr1.phot_variable_time_series_gfov
    gaiadr1.gaiadr1.phot_variable_time_series_gfov_statistical_parameters
    gaiadr1.gaiadr1.rrlyrae
    gaiadr1.gaiadr1.variable_summary
    gaiadr1.gaiadr1.allwise_original_valid
    gaiadr1.gaiadr1.gsc23_original_valid
    gaiadr1.gaiadr1.ppmxl_original_valid
    gaiadr1.gaiadr1.sdssdr9_original_valid
    gaiadr1.gaiadr1.tmass_original_valid
    gaiadr1.gaiadr1.ucac4_original_valid
    gaiadr1.gaiadr1.urat1_original_valid
    gaiadr1.gaiadr1.gaia_source
    gaiadr1.gaiadr1.tgas_source
    gaiadr2.gaiadr2.aux_allwise_agn_gdr2_cross_id
    gaiadr2.gaiadr2.aux_iers_gdr2_cross_id
    gaiadr2.gaiadr2.aux_sso_orbit_residuals
    gaiadr2.gaiadr2.aux_sso_orbits
    gaiadr2.gaiadr2.dr1_neighbourhood
    gaiadr2.gaiadr2.allwise_best_neighbour
    gaiadr2.gaiadr2.allwise_neighbourhood
    gaiadr2.gaiadr2.apassdr9_best_neighbour
    gaiadr2.gaiadr2.apassdr9_neighbourhood
    gaiadr2.gaiadr2.gsc23_best_neighbour
    gaiadr2.gaiadr2.gsc23_neighbourhood
    gaiadr2.gaiadr2.hipparcos2_best_neighbour
    gaiadr2.gaiadr2.hipparcos2_neighbourhood
    gaiadr2.gaiadr2.panstarrs1_best_neighbour
    gaiadr2.gaiadr2.panstarrs1_neighbourhood
    gaiadr2.gaiadr2.ppmxl_best_neighbour
    gaiadr2.gaiadr2.ppmxl_neighbourhood
    gaiadr2.gaiadr2.ravedr5_best_neighbour
    gaiadr2.gaiadr2.ravedr5_neighbourhood
    gaiadr2.gaiadr2.sdssdr9_best_neighbour
    gaiadr2.gaiadr2.sdssdr9_neighbourhood
    gaiadr2.gaiadr2.tmass_best_neighbour
    gaiadr2.gaiadr2.tmass_neighbourhood
    gaiadr2.gaiadr2.tycho2_best_neighbour
    gaiadr2.gaiadr2.tycho2_neighbourhood
    gaiadr2.gaiadr2.urat1_best_neighbour
    gaiadr2.gaiadr2.urat1_neighbourhood
    gaiadr2.gaiadr2.sso_observation
    gaiadr2.gaiadr2.sso_source
    gaiadr2.gaiadr2.vari_cepheid
    gaiadr2.gaiadr2.vari_classifier_class_definition
    gaiadr2.gaiadr2.vari_classifier_definition
    gaiadr2.gaiadr2.vari_classifier_result
    gaiadr2.gaiadr2.vari_long_period_variable
    gaiadr2.gaiadr2.vari_rotation_modulation
    gaiadr2.gaiadr2.vari_rrlyrae
    gaiadr2.gaiadr2.vari_short_timescale
    gaiadr2.gaiadr2.vari_time_series_statistics
    gaiadr2.gaiadr2.panstarrs1_original_valid
    gaiadr2.gaiadr2.gaia_source
    gaiadr2.gaiadr2.ruwe


So that's a lot of tables.  The ones we'll use are:

* `gaiadr2.gaia_source`, which contains Gaia data from [data release 2](https://www.cosmos.esa.int/web/gaia/data-release-2),

* `gaiadr2.panstarrs1_original_valid`, which contains the photometry data we'll use from PanSTARRS, and

* `gaiadr2.panstarrs1_best_neighbour`, which we'll use to cross-match each star observed by Gaia with the same star observed by PanSTARRS.

We can use `load_table` (not `load_tables`) to get the metadata for a single table.  The name of this function is misleading, because it only downloads metadata. 


```python
meta = Gaia.load_table('gaiadr2.gaia_source')
meta
```

    Retrieving table 'gaiadr2.gaia_source'
    Parsing table 'gaiadr2.gaia_source'...
    Done.





    <astroquery.utils.tap.model.taptable.TapTableMeta at 0x7f922376e0a0>



Jupyter shows that the result is an object of type `TapTableMeta`, but it does not display the contents.

To see the metadata, we have to print the object.


```python
print(meta)
```

    TAP Table name: gaiadr2.gaiadr2.gaia_source
    Description: This table has an entry for every Gaia observed source as listed in the
    Main Database accumulating catalogue version from which the catalogue
    release has been generated. It contains the basic source parameters,
    that is only final data (no epoch data) and no spectra (neither final
    nor epoch).
    Num. columns: 96


Notice one gotcha: in the list of table names, this table appears as `gaiadr2.gaiadr2.gaia_source`, but when we load the metadata, we refer to it as `gaiadr2.gaia_source`.

**Exercise:** Go back and try

```
meta = Gaia.load_table('gaiadr2.gaiadr2.gaia_source')
```

What happens?  Is the error message helpful?  If you had not made this error deliberately, would you have been able to figure it out?

## Columns

The following loop prints the names of the columns in the table.


```python
for column in meta.columns:
    print(column.name)
```

    solution_id
    designation
    source_id
    random_index
    ref_epoch
    ra
    ra_error
    dec
    dec_error
    parallax
    parallax_error
    parallax_over_error
    pmra
    pmra_error
    pmdec
    pmdec_error
    ra_dec_corr
    ra_parallax_corr
    ra_pmra_corr
    ra_pmdec_corr
    dec_parallax_corr
    dec_pmra_corr
    dec_pmdec_corr
    parallax_pmra_corr
    parallax_pmdec_corr
    pmra_pmdec_corr
    astrometric_n_obs_al
    astrometric_n_obs_ac
    astrometric_n_good_obs_al
    astrometric_n_bad_obs_al
    astrometric_gof_al
    astrometric_chi2_al
    astrometric_excess_noise
    astrometric_excess_noise_sig
    astrometric_params_solved
    astrometric_primary_flag
    astrometric_weight_al
    astrometric_pseudo_colour
    astrometric_pseudo_colour_error
    mean_varpi_factor_al
    astrometric_matched_observations
    visibility_periods_used
    astrometric_sigma5d_max
    frame_rotator_object_type
    matched_observations
    duplicated_source
    phot_g_n_obs
    phot_g_mean_flux
    phot_g_mean_flux_error
    phot_g_mean_flux_over_error
    phot_g_mean_mag
    phot_bp_n_obs
    phot_bp_mean_flux
    phot_bp_mean_flux_error
    phot_bp_mean_flux_over_error
    phot_bp_mean_mag
    phot_rp_n_obs
    phot_rp_mean_flux
    phot_rp_mean_flux_error
    phot_rp_mean_flux_over_error
    phot_rp_mean_mag
    phot_bp_rp_excess_factor
    phot_proc_mode
    bp_rp
    bp_g
    g_rp
    radial_velocity
    radial_velocity_error
    rv_nb_transits
    rv_template_teff
    rv_template_logg
    rv_template_fe_h
    phot_variable_flag
    l
    b
    ecl_lon
    ecl_lat
    priam_flags
    teff_val
    teff_percentile_lower
    teff_percentile_upper
    a_g_val
    a_g_percentile_lower
    a_g_percentile_upper
    e_bp_min_rp_val
    e_bp_min_rp_percentile_lower
    e_bp_min_rp_percentile_upper
    flame_flags
    radius_val
    radius_percentile_lower
    radius_percentile_upper
    lum_val
    lum_percentile_lower
    lum_percentile_upper
    datalink_url
    epoch_photometry_url


You can probably guess what many of these columns are by looking at the names, but you should resist the temptation to guess.
To find out what the columns mean, [read the documentation](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html).

If you want to know what can go wrong when you don't read the documentation, [you might like this article](https://www.vox.com/future-perfect/2019/6/4/18650969/married-women-miserable-fake-paul-dolan-happiness).

**Exercise:** One of the other tables we'll use is `gaiadr2.gaiadr2.panstarrs1_original_valid`.  Use `load_table` to get the metadata for this table.  How many columns are there and what are their names?

Hint: Remember the gotcha we mentioned earlier.


```python
# Solution

meta2 = Gaia.load_table('gaiadr2.panstarrs1_original_valid')
print(meta2)
```

    Retrieving table 'gaiadr2.panstarrs1_original_valid'
    Parsing table 'gaiadr2.panstarrs1_original_valid'...
    Done.
    TAP Table name: gaiadr2.gaiadr2.panstarrs1_original_valid
    Description: The Panoramic Survey Telescope and Rapid Response System (Pan-STARRS) is
    a system for wide-field astronomical imaging developed and operated by
    the Institute for Astronomy at the University of Hawaii. Pan-STARRS1
    (PS1) is the first part of Pan-STARRS to be completed and is the basis
    for Data Release 1 (DR1). The PS1 survey used a 1.8 meter telescope and
    its 1.4 Gigapixel camera to image the sky in five broadband filters (g,
    r, i, z, y).
    
    The current table contains a filtered subsample of the 10 723 304 629
    entries listed in the original ObjectThin table.
    We used only ObjectThin and MeanObject tables to extract
    panstarrs1OriginalValid table, this means that objects detected only in
    stack images are not included here. The main reason for us to avoid the
    use of objects detected in stack images is that their astrometry is not
    as good as the mean objects astrometry: “The stack positions (raStack,
    decStack) have considerably larger systematic astrometric errors than
    the mean epoch positions (raMean, decMean).” The astrometry for the
    MeanObject positions uses Gaia DR1 as a reference catalog, while the
    stack positions use 2MASS as a reference catalog.
    
    In details, we filtered out all objects where:
    
    -   nDetections = 1
    
    -   no good quality data in Pan-STARRS, objInfoFlag 33554432 not set
    
    -   mean astrometry could not be measured, objInfoFlag 524288 set
    
    -   stack position used for mean astrometry, objInfoFlag 1048576 set
    
    -   error on all magnitudes equal to 0 or to -999;
    
    -   all magnitudes set to -999;
    
    -   error on RA or DEC greater than 1 arcsec.
    
    The number of objects in panstarrs1OriginalValid is 2 264 263 282.
    
    The panstarrs1OriginalValid table contains only a subset of the columns
    available in the combined ObjectThin and MeanObject tables. A
    description of the original ObjectThin and MeanObjects tables can be
    found at:
    https://outerspace.stsci.edu/display/PANSTARRS/PS1+Database+object+and+detection+tables
    
    Download:
    http://mastweb.stsci.edu/ps1casjobs/home.aspx
    Documentation:
    https://outerspace.stsci.edu/display/PANSTARRS
    http://pswww.ifa.hawaii.edu/pswww/
    References:
    The Pan-STARRS1 Surveys, Chambers, K.C., et al. 2016, arXiv:1612.05560
    Pan-STARRS Data Processing System, Magnier, E. A., et al. 2016,
    arXiv:1612.05240
    Pan-STARRS Pixel Processing: Detrending, Warping, Stacking, Waters, C.
    Z., et al. 2016, arXiv:1612.05245
    Pan-STARRS Pixel Analysis: Source Detection and Characterization,
    Magnier, E. A., et al. 2016, arXiv:1612.05244
    Pan-STARRS Photometric and Astrometric Calibration, Magnier, E. A., et
    al. 2016, arXiv:1612.05242
    The Pan-STARRS1 Database and Data Products, Flewelling, H. A., et al.
    2016, arXiv:1612.05243
    
    Catalogue curator:
    SSDC - ASI Space Science Data Center
    https://www.ssdc.asi.it/
    Num. columns: 26



```python
# Solution

for column in meta2.columns:
    print(column.name)
```

    obj_name
    obj_id
    ra
    dec
    ra_error
    dec_error
    epoch_mean
    g_mean_psf_mag
    g_mean_psf_mag_error
    g_flags
    r_mean_psf_mag
    r_mean_psf_mag_error
    r_flags
    i_mean_psf_mag
    i_mean_psf_mag_error
    i_flags
    z_mean_psf_mag
    z_mean_psf_mag_error
    z_flags
    y_mean_psf_mag
    y_mean_psf_mag_error
    y_flags
    n_detections
    zone_id
    obj_info_flag
    quality_flag


## Writing queries

By now you might be wondering how we actually download the data.  With tables this big, you generally don't.  Instead, you use queries to select only the data you want.

A query is a string written in a query language like SQL; for the Gaia database, the query language is a dialect of SQL called ADQL.

Here's an example of an ADQL query.


```python
query1 = """SELECT 
TOP 10
source_id, ref_epoch, ra, dec, parallax 
FROM gaiadr2.gaia_source"""
```

**Python note:** We use a [triple-quoted string](https://docs.python.org/3/tutorial/introduction.html#strings) here so we can include line breaks in the query, which makes it easier to read.

The words in uppercase are ADQL keywords:

* `SELECT` indicates that we are selecting data (as opposed to adding or modifying data).

* `TOP` indicates that we only want the first 10 rows of the table, which is useful for testing a query before asking for all of the data.

* `FROM` specifies which table we want data from.

The third line is a list of column names, indicating which columns we want.  

In this example, the keywords are capitalized and the column names are lowercase.  This is a common style, but it is not required.  ADQL and SQL are not case-sensitive.

To run this query, we use the `Gaia` object, which represents our connection to the Gaia database, and invoke `launch_job`:


```python
job1 = Gaia.launch_job(query1)
job1
```




    <astroquery.utils.tap.model.job.Job at 0x7f9222e9cb20>



The result is an object that represents the job running on a Gaia server.

If you print it, it displays metadata for the forthcoming table.


```python
print(job1)
```

    <Table length=10>
       name    dtype  unit                            description                            
    --------- ------- ---- ------------------------------------------------------------------
    source_id   int64      Unique source identifier (unique within a particular Data Release)
    ref_epoch float64   yr                                                    Reference epoch
           ra float64  deg                                                    Right ascension
          dec float64  deg                                                        Declination
     parallax float64  mas                                                           Parallax
    Jobid: None
    Phase: COMPLETED
    Owner: None
    Output file: sync_20201005090721.xml.gz
    Results: None


Don't worry about `Results: None`.  That does not actually mean there are no results.

However, `Phase: COMPLETED` indicates that the job is complete, so we can get the results like this:


```python
results1 = job1.get_results()
type(results1)
```




    astropy.table.table.Table



**Optional detail:**  Why is `table` repeated three times?  The first is the name of the module, the second is the name of the submodule, and the third is the name of the class.  Most of the time we only care about the last one.  It's like the Linnean name for gorilla, which is *Gorilla Gorilla Gorilla*.

The result is an [Astropy Table](https://docs.astropy.org/en/stable/table/), which is similar to a table in an SQL database except:

* SQL databases are stored on disk drives, so they are persistent; that is, they "survive" even if you turn off the computer.  An Astropy `Table` is stored in memory; it disappears when you turn off the computer (or shut down this Jupyter notebook).

* SQL databases are designed to process queries.  An Astropy `Table` can perform some query-like operations, like selecting columns and rows.  But these operations use Python syntax, not SQL.

Jupyter knows how to display the contents of a `Table`.


```python
results1
```




<i>Table length=10</i>
<table id="table140265627585264" class="table-striped table-bordered table-condensed">
<thead><tr><th>source_id</th><th>ref_epoch</th><th>ra</th><th>dec</th><th>parallax</th></tr></thead>
<thead><tr><th></th><th>yr</th><th>deg</th><th>deg</th><th>mas</th></tr></thead>
<thead><tr><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
<tr><td>4530738361793769600</td><td>2015.5</td><td>281.56725362448725</td><td>20.40682117430378</td><td>0.9785380604519425</td></tr>
<tr><td>4530752651135081216</td><td>2015.5</td><td>281.0861565355257</td><td>20.523350496351846</td><td>0.2674800612552977</td></tr>
<tr><td>4530743343951405568</td><td>2015.5</td><td>281.37114418299177</td><td>20.474147574053124</td><td>-0.43911323550176806</td></tr>
<tr><td>4530755060627162368</td><td>2015.5</td><td>281.2676236268299</td><td>20.558523922346158</td><td>1.1422630184554958</td></tr>
<tr><td>4530746844341315968</td><td>2015.5</td><td>281.1370431749541</td><td>20.377852388898184</td><td>1.0092247424630945</td></tr>
<tr><td>4530768456615026432</td><td>2015.5</td><td>281.8720921436347</td><td>20.31829694530366</td><td>-0.06900136127674149</td></tr>
<tr><td>4530763513119137280</td><td>2015.5</td><td>281.9211808864116</td><td>20.20956829578524</td><td>0.1266016679823622</td></tr>
<tr><td>4530736364618539264</td><td>2015.5</td><td>281.4913475613274</td><td>20.346579041327693</td><td>0.3894019486060072</td></tr>
<tr><td>4530735952305177728</td><td>2015.5</td><td>281.4085549165704</td><td>20.311030903719928</td><td>0.2041189982608354</td></tr>
<tr><td>4530751281056022656</td><td>2015.5</td><td>281.0585328377638</td><td>20.460309556214753</td><td>0.10294642821734962</td></tr>
</table>



Each column has a name, units, and a data type.

For example, the units of `ra` and `dec` are degrees, and their data type is `float64`, which is a 64-bit floating-point number, used to store measurements with a fraction part.

This information comes from the Gaia database, and has been stored in the Astropy `Table` by Astroquery.

**Exercise:** Read [the documentation of this table](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html) and choose a column that looks interesting to you.  Add the column name to the query and run it again.  What are the units of the column you selected?  What is its data type?

## Asynchronous queries

`launch_job` asks the server to run the job "synchronously", which normally means it runs immediately.  But synchronous jobs are limited to 2000 rows.  For queries that return more rows, you should run "asynchronously", which mean they might take longer to get started.

If you are not sure how many rows a query will return, you can use the SQL command `COUNT` to find out how many rows are in the result without actually returning them.  We'll see an example of this later.

The results of an asynchronous query are stored in a file on the server, so you can start a query and come back later to get the results.

For anonymous users, files are kept for three days.

As an example, let's try a query that's similar to `query1`, with two changes:

* It selects the first 3000 rows, so it is bigger than we should run synchronously.

* It uses a new keyword, `WHERE`.


```python
query2 = """SELECT TOP 3000
source_id, ref_epoch, ra, dec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1
"""
```

A `WHERE` clause indicates which rows we want; in this case, the query selects only rows "where" `parallax` is less than 1.  This has the effect of selecting stars with relatively low parallax, which are farther away.  We'll use this clause to exclude nearby stars that are unlikely to be part of GD-1.

`WHERE` is one of the most common clauses in ADQL/SQL, and one of the most useful, because it allows us to select only the rows we need from the database.

We use `launch_job_async` to submit an asynchronous query.


```python
job2 = Gaia.launch_job_async(query2)
print(job2)
```

    INFO: Query finished. [astroquery.utils.tap.core]
    <Table length=3000>
       name    dtype  unit                            description                            
    --------- ------- ---- ------------------------------------------------------------------
    source_id   int64      Unique source identifier (unique within a particular Data Release)
    ref_epoch float64   yr                                                    Reference epoch
           ra float64  deg                                                    Right ascension
          dec float64  deg                                                        Declination
     parallax float64  mas                                                           Parallax
    Jobid: 1601903242219O
    Phase: COMPLETED
    Owner: None
    Output file: async_20201005090722.vot
    Results: None


And here are the results.


```python
results2 = job2.get_results()
results2
```




<i>Table length=3000</i>
<table id="table140265625141056" class="table-striped table-bordered table-condensed">
<thead><tr><th>source_id</th><th>ref_epoch</th><th>ra</th><th>dec</th><th>parallax</th></tr></thead>
<thead><tr><th></th><th>yr</th><th>deg</th><th>deg</th><th>mas</th></tr></thead>
<thead><tr><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
<tr><td>4530738361793769600</td><td>2015.5</td><td>281.56725362448725</td><td>20.40682117430378</td><td>0.9785380604519425</td></tr>
<tr><td>4530752651135081216</td><td>2015.5</td><td>281.0861565355257</td><td>20.523350496351846</td><td>0.2674800612552977</td></tr>
<tr><td>4530743343951405568</td><td>2015.5</td><td>281.37114418299177</td><td>20.474147574053124</td><td>-0.43911323550176806</td></tr>
<tr><td>4530768456615026432</td><td>2015.5</td><td>281.8720921436347</td><td>20.31829694530366</td><td>-0.06900136127674149</td></tr>
<tr><td>4530763513119137280</td><td>2015.5</td><td>281.9211808864116</td><td>20.20956829578524</td><td>0.1266016679823622</td></tr>
<tr><td>4530736364618539264</td><td>2015.5</td><td>281.4913475613274</td><td>20.346579041327693</td><td>0.3894019486060072</td></tr>
<tr><td>4530735952305177728</td><td>2015.5</td><td>281.4085549165704</td><td>20.311030903719928</td><td>0.2041189982608354</td></tr>
<tr><td>4530751281056022656</td><td>2015.5</td><td>281.0585328377638</td><td>20.460309556214753</td><td>0.10294642821734962</td></tr>
<tr><td>4530740938774409344</td><td>2015.5</td><td>281.3762569536416</td><td>20.436140058941206</td><td>0.9242670062090182</td></tr>
<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
<tr><td>4467710915011802624</td><td>2015.5</td><td>269.9680969307347</td><td>1.1429085038160882</td><td>0.42361471245557913</td></tr>
<tr><td>4467706551328679552</td><td>2015.5</td><td>270.033164589881</td><td>1.0565747323689927</td><td>0.922888231734588</td></tr>
<tr><td>4467712255037300096</td><td>2015.5</td><td>270.7724717923047</td><td>0.6581664892880896</td><td>-2.669179465293931</td></tr>
<tr><td>4467735001181761792</td><td>2015.5</td><td>270.3628606248308</td><td>0.8947079323599124</td><td>0.6117399163086398</td></tr>
<tr><td>4467737101421916672</td><td>2015.5</td><td>270.5110834661444</td><td>0.9806225910160181</td><td>-0.39818224846127004</td></tr>
<tr><td>4467707547757327488</td><td>2015.5</td><td>269.88746280594927</td><td>1.0212759940136962</td><td>0.7741412301054209</td></tr>
<tr><td>4467732772094573056</td><td>2015.5</td><td>270.55997182760126</td><td>0.9037072088489417</td><td>-1.7920417800164183</td></tr>
<tr><td>4467732355491087744</td><td>2015.5</td><td>270.6730790702491</td><td>0.9197224705139885</td><td>-0.3464446494840354</td></tr>
<tr><td>4467717099766944512</td><td>2015.5</td><td>270.57667173120825</td><td>0.726277659009568</td><td>0.05443955111134051</td></tr>
<tr><td>4467719058265781248</td><td>2015.5</td><td>270.7248052971514</td><td>0.8205551921782785</td><td>0.3733943917490343</td></tr>
</table>



You might notice that some values of `parallax` are negative.  As [this FAQ explains](https://www.cosmos.esa.int/web/gaia/archive-tips#negative%20parallax), "Negative parallaxes are caused by errors in the observations."  Negative parallaxes have "no physical meaning," but they can be a "useful diagnostic on the quality of the astrometric solution."

Later we will see an example where we use `parallax` and `parallax_error` to identify stars where the distance estimate is likely to be inaccurate.

**Exercise:** The clauses in a query have to be in the right order.  Go back and change the order of the clauses in `query2` and run it again.  

The query should fail, but notice that you don't get much useful debugging information.  

For this reason, developing and debugging ADQL queries can be really hard.  A few suggestions that might help:

* Whenever possible, start with a working query, either an example you find online or a query you have used in the past.

* Make small changes and test each change before you continue.

* While you are debugging, use `TOP` to limit the number of rows in the result.  That will make each attempt run faster, which reduces your testing time.  

* Launching test queries synchronously might make them start faster, too.

## Operators

In a `WHERE` clause, you can use any of the [SQL comparison operators](https://www.w3schools.com/sql/sql_operators.asp); here are the most common ones:

| Symbol | Operation
|--------| :---
| `>` | greater than
| `<` | less than
| `>=` | greater than or equal
| `<=` | less than or equal
| `=` | equal
| `!=` or `<>` | not equal

Most of these are the same as Python, but some are not.  In particular, notice that the equality operator is `=`, not `==`.
Be careful to keep your Python out of your ADQL!

You can combine comparisons using the logical operators:

* AND: true if both comparisons are true
* OR: true if either or both comparisons are true

Finally, you can use `NOT` to invert the result of a comparison. 

**Exercise:** [Read about SQL operators here](https://www.w3schools.com/sql/sql_operators.asp) and then modify the previous query to select rows where `bp_rp` is between `-0.75` and `2`.

You can [read about this variable here](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html).


```python
# Solution

# This is what most people will probably do

query = """SELECT TOP 10
source_id, ref_epoch, ra, dec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1 
  AND bp_rp > -0.75 AND bp_rp < 2
"""
```


```python
# Solution

# But if someone notices the BETWEEN operator, 
# they might do this

query = """SELECT TOP 10
source_id, ref_epoch, ra, dec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1 
  AND bp_rp BETWEEN -0.75 AND 2
"""
```

This [Hertzsprung-Russell diagram](https://sci.esa.int/web/gaia/-/60198-gaia-hertzsprung-russell-diagram) shows the BP-RP color and luminosity of stars in the Gaia catalog.

Selecting stars with `bp-rp` less than 2 excludes many [class M dwarf stars](https://xkcd.com/2360/), which are low temperature, low luminosity.  A star like that at GD-1's distance would be hard to detect, so if it is detected, it it more likely to be in the foreground.

## Cleaning up

Asynchronous jobs have a `jobid`.


```python
job1.jobid, job2.jobid
```




    (None, '1601903242219O')



Which you can use to remove the job from the server.


```python
Gaia.remove_jobs([job2.jobid])
```

    Removed jobs: '['1601903242219O']'.


If you don't remove it job from the server, it will be removed eventually, so don't feel too bad if you don't clean up after yourself.

## Formatting queries

So far the queries have been string "literals", meaning that the entire string is part of the program.
But writing queries yourself can be slow, repetitive, and error-prone.

It is often a good idea to write Python code that assembles a query for you.  One useful tool for that is the [string `format` method](https://www.w3schools.com/python/ref_string_format.asp).

As an example, we'll divide the previous query into two parts; a list of column names and a "base" for the query that contains everything except the column names.

Here's the list of columns we'll select.  


```python
columns = 'source_id, ra, dec, pmra, pmdec, parallax, parallax_error, radial_velocity'
```

And here's the base; it's a string that contains at least one format specifier in curly brackets (braces).


```python
query3_base = """SELECT TOP 10 
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2
"""
```

This base query contains one format specifier, `{columns}`, which is a placeholder for the list of column names we will provide.

To assemble the query, we invoke `format` on the base string and provide a keyword argument that assigns a value to `columns`.


```python
query3 = query3_base.format(columns=columns)
```

The result is a string with line breaks.  If you display it, the line breaks appear as `\n`.


```python
query3
```




    'SELECT TOP 10 \nsource_id, ra, dec, pmra, pmdec, parallax, parallax_error, radial_velocity\nFROM gaiadr2.gaia_source\nWHERE parallax < 1\n  AND bp_rp BETWEEN -0.75 AND 2\n'



But if you print it, the line breaks appear as... line breaks.


```python
print(query3)
```

    SELECT TOP 10 
    source_id, ra, dec, pmra, pmdec, parallax, parallax_error, radial_velocity
    FROM gaiadr2.gaia_source
    WHERE parallax < 1
      AND bp_rp BETWEEN -0.75 AND 2
    


Notice that the format specifier has been replaced with the value of `columns`.

Let's run it and see if it works:


```python
job3 = Gaia.launch_job(query3)
print(job3)
```

    <Table length=10>
          name       dtype    unit                              description                             n_bad
    --------------- ------- -------- ------------------------------------------------------------------ -----
          source_id   int64          Unique source identifier (unique within a particular Data Release)     0
                 ra float64      deg                                                    Right ascension     0
                dec float64      deg                                                        Declination     0
               pmra float64 mas / yr                         Proper motion in right ascension direction     0
              pmdec float64 mas / yr                             Proper motion in declination direction     0
           parallax float64      mas                                                           Parallax     0
     parallax_error float64      mas                                         Standard error of parallax     0
    radial_velocity float64   km / s                                                    Radial velocity    10
    Jobid: None
    Phase: COMPLETED
    Owner: None
    Output file: sync_20201005090726.xml.gz
    Results: None



```python
results3 = job3.get_results()
results3
```




<i>Table length=10</i>
<table id="table140265627700432" class="table-striped table-bordered table-condensed">
<thead><tr><th>source_id</th><th>ra</th><th>dec</th><th>pmra</th><th>pmdec</th><th>parallax</th><th>parallax_error</th><th>radial_velocity</th></tr></thead>
<thead><tr><th></th><th>deg</th><th>deg</th><th>mas / yr</th><th>mas / yr</th><th>mas</th><th>mas</th><th>km / s</th></tr></thead>
<thead><tr><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
<tr><td>4467710915011802624</td><td>269.9680969307347</td><td>1.1429085038160882</td><td>2.0233280236600626</td><td>-2.5692427875510266</td><td>0.42361471245557913</td><td>0.470352406647465</td><td>--</td></tr>
<tr><td>4467706551328679552</td><td>270.033164589881</td><td>1.0565747323689927</td><td>-3.414829591355289</td><td>-3.8437215857495737</td><td>0.922888231734588</td><td>0.927008559859825</td><td>--</td></tr>
<tr><td>4467712255037300096</td><td>270.7724717923047</td><td>0.6581664892880896</td><td>-3.5620173752896025</td><td>-6.595792323153987</td><td>-2.669179465293931</td><td>0.9719742773203504</td><td>--</td></tr>
<tr><td>4467735001181761792</td><td>270.3628606248308</td><td>0.8947079323599124</td><td>2.1307079926489205</td><td>0.8826727710910712</td><td>0.6117399163086398</td><td>0.509812721702093</td><td>--</td></tr>
<tr><td>4467737101421916672</td><td>270.5110834661444</td><td>0.9806225910160181</td><td>0.17532366511560785</td><td>-5.113270239706202</td><td>-0.39818224846127004</td><td>0.7549581886719651</td><td>--</td></tr>
<tr><td>4467707547757327488</td><td>269.88746280594927</td><td>1.0212759940136962</td><td>-2.6382230817672987</td><td>-3.707776532049287</td><td>0.7741412301054209</td><td>0.3022057897812064</td><td>--</td></tr>
<tr><td>4467732355491087744</td><td>270.6730790702491</td><td>0.9197224705139885</td><td>-2.2735991502653037</td><td>-11.864952855984358</td><td>-0.3464446494840354</td><td>0.4937921513912002</td><td>--</td></tr>
<tr><td>4467717099766944512</td><td>270.57667173120825</td><td>0.726277659009568</td><td>-3.4598362614808367</td><td>-4.601426893365921</td><td>0.05443955111134051</td><td>0.8867339293525688</td><td>--</td></tr>
<tr><td>4467719058265781248</td><td>270.7248052971514</td><td>0.8205551921782785</td><td>-3.255079498426542</td><td>-9.249285069111085</td><td>0.3733943917490343</td><td>0.390952370410666</td><td>--</td></tr>
<tr><td>4467722326741572352</td><td>270.87431291888504</td><td>0.8595565975869158</td><td>0.10696398351859826</td><td>1.2035993780158853</td><td>-0.11850943432864373</td><td>0.1660452431882023</td><td>--</td></tr>
</table>



Good so far.

**Exercise:** This query always selects sources with `parallax` less than 1.  But suppose you want to take that upper bound as an input.

Modify `query3_base` to replace `1` with a format specifier like `{max_parallax}`.  Now, when you call `format`, add a keyword argument that assigns a value to `max_parallax`, and confirm that the format specifier gets replaced with the value you provide.


```python
# Solution

query4_base = """SELECT TOP 10
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < {max_parallax} AND 
bp_rp BETWEEN -0.75 AND 2
"""
```


```python
# Solution

query4 = query4_base.format(columns=columns,
                          max_parallax=0.5)
print(query)
```

    SELECT TOP 10
    source_id, ra, dec, pmra, pmdec, parallax, parallax_error, radial_velocity
    FROM gaiadr2.gaia_source
    WHERE parallax < 0.5 AND 
    bp_rp BETWEEN -0.75 AND 2
    


**Style note:**  You might notice that the variable names in this notebook are numbered, like `query1`, `query2`, etc.  

The advantage of this style is that it isolates each section of the notebook from the others, so if you go back and run the cells out of order, it's less likely that you will get unexpected interactions.

A drawback of this style is that it can be a nuisance to update the notebook if you add, remove, or reorder a section.

What do you think of this choice?  Are there alternatives you prefer?

## Summary

This notebook demonstrates the following steps:

1. Making a connection to the Gaia server,

2. Exploring information about the database and the tables it contains,

3. Writing a query and sending it to the server, and finally

4. Downloading the response from the server as an Astropy `Table`.

## Best practices

* If you can't download an entire dataset (or it's not practical) use queries to select the data you need.

* Read the metadata and the documentation to make sure you understand the tables, their columns, and what they mean.

* Develop queries incrementally: start with something simple, test it, and add a little bit at a time.

* Use ADQL features like `TOP` and `COUNT` to test before you run a query that might return a lot of data.

* If you know your query will return fewer than 3000 rows, you can run it synchronously, which might complete faster (but it doesn't seem to make much difference).  If it might return more than 3000 rows, you should run it asynchronously.

* ADQL and SQL are not case-sensitive, so you don't have to capitalize the keywords, but you should.

* ADQL and SQL don't require you to break a query into multiple lines, but you should.


Jupyter notebooks can be good for developing and testing code, but they have some drawbacks.  In particular, if you run the cells out of order, you might find that variables don't have the values you expect.

There are a few things you can do to mitigate these problems:

* Make each section of the notebook self-contained.  Try not to use the same variable name in more than one section.

* Keep notebooks short.  Look for places where you can break your analysis into phases with one notebook per phase.
