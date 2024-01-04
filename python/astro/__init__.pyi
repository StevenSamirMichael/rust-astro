"""
Toolkit containing functions and classes used in satellite dynamics
calculations.
"""

from __future__ import annotations
import typing
import numpy.typing as npt
import numpy as np

import astro

__all__ = [
    "TLE",
    "duration",
    "frametransform",
    "gravity",
    "gravmodel",
    "itrfcoord",
    "jplephem",
    "lpephem",
    "nrlmsise00",
    "quaternion",
    "satprop",
    "sgp4",
    "solarsystem",
    "time",
    "timescale",
    "univ",
    "utils"
]


__all__ = ['time', 'duration', 'timescale', 'quaternion', 'sgp4', 'gravmodel', 'gravity', 'nrlmsise00', 'univ', 'solarsystem', 'TLE', 'itrfcoord', 'frametransform', 'jplephem', 'lpephem', 'satprop', 'utils']

class itrfcoord():
    """
    Representation of a coordinate in the
    International Terrestrial Reference Frame (ITRF)

    This coordinate object can be created from and also
    output to Geodetic coordinates (latitude, longitude,
    height above ellipsoid).

    Functions are also available to provide rotation
    quaternions to the East-North-Up frame
    and North-East-Down frame at this coordinate
    """

    @typing.overload
    def __init__(self) -> None:
        """
        Represent a coordinate in the ITRF (International Terrestrial Reference Frame)
        Inputs are 3 separate floats representing ITRF Cartesian position in meters
        """ 
        
    def from_geodetic() -> itrfcoord:
        """            
        Create coordinate from input geodetic
        Optional inputs, in order:
        
        latitude, radians
        longitude, radians
        heigth above ellipsoid, meters
        
        
        Optional kwargs:
        
        latitude_deg: latitude, degrees
        longitude_deg: longitude, degrees
        latitude_rad: latitude, radians
        longitude_rad: longitude, radians
        altitude: height above ellipsoid, meters    
        """
    
    @property
    def latitude_deg(self) -> float:
        """
        Latitude in degrees
        """
    
    @property
    def longitude_deg(self) -> float:
        """
        Longitude in degrees
        """
        
    @property
    def latitude_rad(self) -> float:
        """
        Latitude in radians
        """
        
    @property
    def longitude_rad(self) -> float:
        """
        Longitude in radians
        """
        
    @property
    def altitude(self) -> float:
        """
        Altitude above ellipsoid, in meters
        """
        
    @property
    def geodetic_rad(self) -> (float, float, float):
        """
        Tuple with: (latitude_rad, longitude_rad, altitude)
        """
        
    @property
    def geodetic_deg(self) -> (float, float, float):
        """
        Tuple with (latitude_deg, longitude_deg, altitude)
        """
        
    @property
    def vec(self) -> npt.NDArray[np.float64]:
        """
        Cartesian ITRF coord as numpy array
        """
        
    def qned2itrf(self) -> astro.quaternion:
        """
        Quaternion representing rotation from North-East-Down (NED)
        to ITRF at this location
        """
        
    def qenu2itrf(self) -> astro.quaternion:
        """
        Quaternion representiong rotation from East-North-Up (ENU)
        to ITRF at this location
        """