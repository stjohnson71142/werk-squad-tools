from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)

def Distance(z):  
    """
    Calculate luminosity distance from redshift
    --------------------------------------------------
    Uses FlatLambdaCDM and .luminosity_distance() function from astropy.cosmology to make calculation
    --------------------------------------------------
    Args:
        z: numpy array of redshift values
    --------------------------------------------------
    Returns:
        numpy array of luminosity distance values in Mpc
    """
    
    distance = cosmo.luminosity_distance(z).value
    
    return(distance)