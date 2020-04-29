import numpy as np
import astropy.units as u

def Ha_luminosity(distance, Ha_flux):
    """
    Calculate Halpha luminosity for a galaxy using Halpha flux and luminosity distance
    -----------------------------------------
    Note that Halpha flux is assumed in this function. If using Hbeta flux instead, the value should be multiplied by 2.86
    beforehand to get a rough Halpha flux equivalent
    -----------------------------------------
    Args:
        distance: numpy array of luminosity distance values in Mpc
        Ha_flux: numpy array of Halpha flux values in ergs*s^-1*cm^-2
    -----------------------------------------
    Returns:
        numpy array of Halpha luminosity values in ergs*sec^-1
    """
    
    Ha_luminosity = Ha_flux * (4 * np.pi * ((distance * 3.08567758128 * 10**24)**2))  
    
    return(Ha_luminosity)