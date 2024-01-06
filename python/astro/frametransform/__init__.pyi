"""
Transformations between coordinate frames, and associated utility functions
"""

from __future__ import annotations
import typing
import numpy.typing as npt
import numpy as np

import astro
import datetime

def gmst(tm: astro.time|npt.ArrayLike[astro.time]) -> float|npt.ArrayLike[np.float]:
    """    
    Greenwich Mean Sidereal Time
    
    Vallado algorithm 15:
    
    GMST = 67310.5481 + (876600h + 8640184.812866) * tᵤₜ₁ * (0.983104 + tᵤₜ₁ * −6.2e−6)
    
    Input is astro.time object or list or numpy array of astro.time objects.
    
    Output is float or numpy array of floats with GMST in radians matched element-wise
    to the input times.    
    """
    
def gmst(tm: astro.time|npt.ArrayLike[astro.time]) -> float|npt.ArrayLike[np.float]:
    """    
    Greenwich Apparent Sidereal Time
    
    Input is astro.time object or list or numpy array of astro.time objects.
    
    Output is float or numpy array of floats with GAST in radians matched element-wise
    to the input times.    
    """    

def earth_rotation_angle(tm: astro.time|npt.ArrayLike[astro.time]) -> float|npt.ArrayLike[np.float]:
    """    
    Earth rotation angle    
    
    See:
    https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36_043.pdf?__blob=publicationFile&v=1

    Equation 5.15
    
    Input is astro.time object or list or numpy array of astro.time objects.
    
    Output is float or numpy array of floats with Earth Rotation Angle in radians
    matched element-wise to the input times.    
    """    


def qitrf2tirs(tm: astro.time|npt.ArrayLike[astro.time]) -> astro.quaternion|npt.ArrayLike[astro.quaternion]:
    """
    Rotation from International Terrestrial Reference Frame
    (ITRF) to the Terrestrial Intermediate Reference System (TIRS)
    represented as astro.quaterinion object

    Input is astro.time object or list or numpy array of astro.time objects.
    
    Output is astro.quaternion or numpy array of astro.quaternion representiong
    rotations from itrf to tirs matched element-wise to the input times       
    """

def qtirs2cirs(tm: astro.time|npt.ArrayLike[astro.time]) -> astro.quaternion|npt.ArrayLike[astro.quaternion]:
    """
    Rotation from Terrestrial Intermediate Reference System (TIRS)
    to the Celestial Intermediate Reference System (CIRS)


    Input is astro.time object or list or numpy array of astro.time objects.
    
    Output is astro.quaternion or numpy array of astro.quaternion representiong
    rotations from itrf to tirs matched element-wise to the input times       
    """


def qitrf2gcrf(tm: astro.time|npt.ArrayLike[astro.time]) -> astro.quaternion|npt.ArrayLike[astro.quaternion]:
    """
    Quaternion representing rotation from the
    International Terrestrial Reference Frame (ITRF)
    to the Geocentric Celestial Reference Frame (GCRF)

    Uses full IAU2006 Reduction
    See IERS Technical Note 36, Chapter 5

    but does not include solid tides, ocean tides

    Note: Very computationally expensive    
    
    Input is astro.time object or list or numpy array of astro.time objects.
    
    Output is astro.quaternion or numpy array of astro.quaternion representiong
    rotations from itrf to gcrf matched element-wise to the input times  
    """

def qteme2itrf(tm: astro.time|npt.ArrayLike[astro.time]) -> astro.quaternion|npt.ArrayLike[astro.quaternion]:
    """
    Quaternion representing rotation from the
    True Equator Mean Equinox (TEME) frame
    to the International Terrestrial Reference Frame (ITRF)
    
    This is equation 3-90 in Vallado
    
    Note: TEME is the output frame of the SGP4 propagator used to 
    compute position from two-line element sets.
    
    Input is astro.time object or list or numpy array of astro.time objects.
    
    Output is astro.quaternion or numpy array of astro.quaternion representiong
    rotations from TEME to ITRF matched element-wise to the input times  
    """