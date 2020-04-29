from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np

def findR23(OII, OIII, NII, Halpha, HB, mass):

    '''

    A method for finding galaxy metallicity based on R23 flux ratio[1].

    Parameters
    ----------

    OII : [float of some kind]
        OII doublet emission line flux.
    OIII : [float of some kind]
        OIII [wavelength] emission line flux.
    NII : [float of some kind]
        NII [wavelength] emission line flux.
    Halpha : [float of some kind]
        Halpha 6563 emission line flux.
    HB : [float of some kind]
        Hbeta 4860 emission line flux.
    mass : [float of some kind]
        Stellar mass of galaxy, in units of log(M_sun).

    Returns
    -------

    Z_R23 : [float of some kind]
        R23 metallicity, units of 12+log(O/H). Or nan on failure.
    flag : str
        Indicates branch/method used, if any.

    Notes
    -----

    A metallicity determined this way is twofold degenerate. The bulk of this method
    is concerned with breaking that degeneracy based on comparisons of flux ratios,
    or failing that by using the galaxy's stellar mass. This method of galaxy metallicity
    determination is based on McGaugh 1991.

    [1] McGaugh S.S., 1991, ApJ 380, 140.
    '''

    # PART 1: defining things

    # line flux ratios
    OIIHB = OII/HB

    OIIIHB = OIII/HB

    N2 = NII/Halpha

    O3N2 = np.log10(OIIIHB/N2)

    N2O2 = NII/OII

    R23 = OIIHB+OIIIHB

    # actual arguments for metallicity functions
    x = np.log10(R23)
    y = np.log10(1.3*OIII/OII)

    # PART 2: getting metallicities

    # Programming Hard Mode: Write this without if statements (to be clear: I
    # have selected programming easy mode).

    # N2O2 method (if present at all):
    if np.log10(N2O2) > -1.2:
        Z_R23 = upper_branch(x,y)
        flag = 'upper'
    # N2 index:
    elif N2 > 0.1:
        Z_R23 = upper_branch(x,y)
        flag = 'upper'
    # mass:
    elif mass > 9 or mass < 8:
        if mass < 8:
            Z_R23 = lower_branch(x,y)
            flag = 'lower'
        if mass > 9:
            Z_R23 = upper_branch(x,y)
            flag = 'upper w/ mass'
    # remains degenerate:
    else:
        Z_R23 = np.NaN

    if np.isnan(Z_R23) == True:
        flag = 'none'

    return(Z_R23, flag)

# Our two R23 metallicity functions (w/ conditions determining their respective uses given above):

def lower_branch(x,y):
    '''

    Function called when degeneracy is broken by galaxy mass of < 8 log(M_sun).

    Parameters
    ----------

    x : [float of some kind]
        log(R23) from findR23.
    y : [float of some kind]
        log(1.3*`OIII`/`OII`) from findR23.

    Returns
    -------

    Z : [float of some kind]
        Metallicity with units of (12 + log(O/H)).
    '''
    Z = 12.0 - 4.944 + 0.767*x + 0.602*x**2 - y*(0.29+0.332*x-0.331*x**2)

    return(Z)
def upper_branch(x,y):
    '''

    Called when degeneracy is broken by N2O2, N2 index, or mass > 9 log(M_sun).

    Parameters
    ----------

    x : [float of some kind]
        log(R23) from findR23.
    y : [float of some kind]
        log(1.3*`OIII`/`OII`) from findR23.

    Returns
    -------

    Z : [float of some kind]
        Metallicity with units of (12 + log(O/H)).
    '''
    Z = 12.0 - 2.939 - 0.2*x - 0.237*x**2 - 0.305*x**3 - 0.0283*x**4 - y*(
    0.0047 - 0.0221*x - 0.102*x**2 - 0.0817*x**3 - 0.00717*x**4)
    return(Z)
