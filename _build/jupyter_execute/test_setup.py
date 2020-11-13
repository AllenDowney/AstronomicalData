#!/usr/bin/env python
# coding: utf-8

# # Astronomical Data in Python
# 
# This notebook imports the libraries we need for the workshop.
# 
# If any of them are missing, you'll get an error message.
# 
# If you don't get any error messages, you are all set.

# In[7]:


from wget import download


# In[1]:


import pandas as pd
import numpy as np


# In[2]:


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import Polygon


# In[3]:


import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table


# In[4]:


import gala.coordinates as gc
from pyia import GaiaData


# In[5]:


# Note: running this import statement opens a connection
# to a Gaia server, so it will fail if you are not connected
# to the internet.

from astroquery.gaia import Gaia


# During the workshop, we might put some code on Slack and ask you to cut and paste it into the notebook.
# 
# If you are on a Mac, you might encounter a problem: 

# In[ ]:




