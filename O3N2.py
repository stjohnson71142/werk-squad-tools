
# coding: utf-8

# In[2]:


import astropy.units as u
import astropy.constants as ac
from astropy.cosmology import FlatLambdaCDM
import numpy as np
import pandas as pd
import os
from os.path import expanduser
from astropy.table import Table
from glob import glob


# In[37]:


def O3N2(gl):
    """ Calculate 'O3N2' from available ‘O3HB' and 'N2'
    
    Parameters
    ----------
    gl : object
        The information of galaxies
    
    Returns
    -------
    array 
        A NumPy array of available 'O3N2' values
    """
    # Define the calculation for 'O3N2'
    def eq(o3hb, n2):
        value = np.log10(o3hb/n2)
        return value
    # Create a new array of float zeros
    array = np.zeros(782)
    # Calculate 'O3N2' for each target 
    for i in gl.index:
        if (gl['O3HB'][i].dtype == np.float) and (gl['N2'][i].dtype == np.float):
            value = eq(gl['O3HB'][i],gl['N2'][i])
            array[i-2] = np.float(value)
        else:
            array[i-2] = np.nan # for not available 'O3N2'
    # Return
    return array


# In[39]:


def Z_O3N2L(gl):
    """ Calculate the metallicities from 'O3N2' based on the Pettini and Pagel equation 1
    
    Parameters
    ----------
    gl : object
        The information of galaxies 
    
    Returns
    -------
    array 
        A NumPy array of available metallicities from 'O3N2'
    """
    # Create a new array of zeros
    array = np.zeros(782)
    # Calculate the metallicities based on PP eq.1 if 'O3N2' < 1.9
    for i in gl.index:
        if (gl['O3N2'][i].dtype == np.float) and (gl['O3N2'][i] < 1.9):
            metal = 8.73-(0.32*gl['O3N2'][i])
            array[i-2] = np.float(metal)
        else:
            array[i-2] = np.nan # for not available metallicity
    # Return
    return array

