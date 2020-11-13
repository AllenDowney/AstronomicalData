#!/usr/bin/env python
# coding: utf-8

# # The Notebook of Last Resort

# If you are not able to get everything installed that we need for the workshop, you have the option of running this notebook on Colab.
# 
# Before you get started, you probably want to press the Save button!

# In[1]:


# If we're running on Colab, install libraries

import sys
IN_COLAB = 'google.colab' in sys.modules

if IN_COLAB:
    get_ipython().system('pip install astroquery astro-gala pyia')


# That should be everything you need.  Now you can type code and run it in the following cells.
