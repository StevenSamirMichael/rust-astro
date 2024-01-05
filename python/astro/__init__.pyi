"""
Toolkit containing functions and classes used in satellite dynamics
calculations.
"""

from __future__ import annotations
import typing
import numpy.typing as npt
import numpy as np

import astro
import datetime

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

class timescale():
    """
    Specify time scale used to represent or convert between the "astro.time"
    representation of time     
    
    Most of the time, these are not needed directly, but various time scales
    are needed to compute precise rotations between various inertial and 
    Earth-fixed coordinate frames
    
    For details, see:
    https://stjarnhimlen.se/comp/time.html
    """

    @property
    def Invalid() -> int:
        """
        Invalid time scale
        """
    
    @property
    def UTC() -> int:
        """
        Universal Time Coordinate
        """
        
    def TT() -> int:
        """
        Terrestrial Time
        """
        
    def UT1() -> int:
        """
        UT1
        """
        
    def TAI() -> int:
        """
        International Atomic Time
        (nice because it is monotonically increasing)
        """
        
    def GPS() -> int:
        """
        Global Positioning System (GPS) time
        """
        
    def TDB() -> int:
        """
        Barycentric Dynamical Time
        """

class time():
    """
    Representation of an instant in time
    
    This has functionality similar to the "datetime" object, and in fact has 
    the ability to convert to an from the "datetime" object.  However, a separate
    time representation is needed as the "datetime" object does not allow for 
    conversion between various time epochs (GPS, TAI, UTC, UT1, etc...)
    """
        
    def __init__(self, *args) -> astro.time:
        """
        Create a "time" object. 
        
        
        If no arguments are passed in, the created object represents
        the current time, i.e. the time at which the function was called
        
        
        If 3 integers are passed in, they represent a UTC date specified
        by the standard Gregorian year, month (1-based), and day of month
        (1-based)
        
        if 5 integers and a float are passed in, they represent a UTC 
        date and time.  The 1st 3 numbers represent the standard
        Gregorian year, month, and day as above.  The last 3 represent the
        hour of the day [0,23], the minute of the hour [0,59], and the 
        second (including fractional component) of the minute    
        
        Example 1: 
        print(astro.time(2023, 3, 5, 11, 3,45.453))
        2023-03-05 11:03:45.453Z
        
        Example 2:
        print(astro.time(2023, 3, 5))
        2023-03-05 00:00:00.000Z
                    
        """
        
    @staticmethod
    def now() -> astro.time:
        """
        Create a "time" object representing the instant of time at the 
        calling of the function.
        """
        
    @staticmethod
    def from_date(year: int, month: int, day: int) -> astro.time:
        """
        Returns a time object representing the start of the day (midnight)
        on the provided date, specified by standard Gregorian year, month
        (1-based), and day of month (1-based)
        """

    def to_date(self) -> (int, int, int):
        """
        Return tuple representing as UTC Gegorian date of the 
        time object.  Tuple has 6 elements:
        1 : Gregorian Year
        2 : Gregorian month (1 = January, 2 = February, ...)
        3 : Day of month, beginning with 1
        
        Fractional day components are neglected
        """
        
    @staticmethod
    def from_gregorian(self, year: int, month: int,
                       day: int, hour: int,
                       min: int, sec: float) -> astro.time:
        """
        Create time object from 6 input arguments representing
        UTC Gregorian time.  
        
        Inputs are:    
        1 : Gregorian Year
        2 : Gregorian month (1 = January, 2 = February, ...)
        3 : Day of month, beginning with 1
        4 : Hour of day, in range [0,23]
        5 : Minute of hour, in range [0,59]
        6 : floating point second of minute, in range [0,60) 
        
        Example:
        print(astro.time.from_gregorian(2023, 3, 5, 11, 3,45.453))
        2023-03-05 11:03:45.453Z        
        """    
    
    
    def to_gregorian(self) -> (int, int, int, int, int, float):
        """
        Return tuple representing as UTC Gegorian date and time of the 
        time object.  Tuple has 6 elements:
        1 : Gregorian Year
        2 : Gregorian month (1 = January, 2 = February, ...)
        3 : Day of month, beginning with 1
        4 : Hour of day, in range [0,23]
        5 : Minute of hour, in range [0,59]
        6 : floating point second of minute, in range [0,60)    
        """

    @staticmethod
    def from_datetime(dt: datetime.datetime) -> astro.time:
        """
        Convert input "datetime.datetime" object to an 
        "astro.time" object represenging the same 
        instant in time
        """
        
    def datetime(self, utc: bool=True) -> datetime.datetime:
        """
        Convert object to "datetime.datetime" object representing
        same instant in time.
        
        The optional boolean input "utc" specifies wither to make the
        "daettime.datetime" object represent time in the local timezone, 
        or whether to have the "datetime.datetime" object be in "UTC" time.
        Default is true
        
        Example: (from Easterm Standard Time time zone)
        dt = astro.time(2023, 6, 3, 6, 19, 34).datetime(True)
        print(dt)
        2023-06-03 06:19:34+00:00

        dt = astro.time(2023, 6, 3, 6, 19, 34).datetime(False)
        print(dt)
        2023-06-03 02:19:34
        """
        
    def to_mjd(self, scale: astro.timescale) -> float:
        """
        Represent time instance as a Modified Julian Date 
        with the provided time scale
        """
        
    def to_jd(self, scale: astro.timescale) -> float:
        """
        Represent time instance as Julian Date with 
        the provided time scale
        """
    
    def to_unixtime(self) -> float:
        """
        Represent time as unixtime
        
        (seconds since Jan 1, 1970 UTC)
        
        Includes fractional comopnent of seconds
        """
        
    def __add__(self, 
            other: astro.duration|npt.ArrayLike[float]|
            float|list[float]|npt.ArrayLike[astro.duration]) -> astro.time|npt.ArrayLike[astro.time]:
        """
        Return an astro.time object or nunpy array of astro.time objects 
        representing the input "added to" the current object
        
        Possible inputs and corresponding outputs:
        
        1. float: return astro.time object incremented by input number of days
        2. astro.duration: return astro.time object incremented by duration
        3. list[float]: return numpy array of astro.time objects, representing
           an element-wise addition of days to the "self"
        4. list[duration]: reuturn numpy array of astro.time objects, with each
           object representing an element-wise addition of "duration" to the "self"
        5. numpy.array(float): return numpy array of astro.time objects, with each
           object representing an element-wise addition of days to the "self"        
        """
    
    def __sub__(self,
                other: astro.duration|astro.time|npt.ArrayLoke[float]|
                npt.ArrayLike[astro.duration]|
                list[float]) -> astro.time | astro.duration | npt.ArrayLike[astro.time]:
        """
        Return an astro.time object or numpy array of astro.time objects
        representing the input "subtracted from" the current object
        
        Possible inputs and corresponding outputs:
        
        1. astro.time: output is duration representing the difference
           between the "other" time and "self"
        2. astro.duration: output is astro.time object representing time minus
           the input duration
        3. list[float]: return numpy array of astro.time objects, representing
           an element-wise subtraction of days to the "self"
        4. list[duration]: return numpy array of astro.time objects, representing
           an element-wise subtraction of "duration" from the "self"
        5. numpy.array(float): return numpy array of astro.time objects, with 
           each object representing an element-wise subtraction of days from 
           the "self".
           
        """

class duration():
    """
    Representation of a time duration
    """
    
    @staticmethod
    def from_days(d: float) -> duration:
        """
        Create duration object given input number of days
        Note: a day is defined as 86,400 seconds
        """
        
    @staticmethod
    def from_seconds(d: float) -> duration:
        """
        Create duration object representing input number of seconds
        """
    
    @staticmethod
    def from_minutes(d: float) -> duration:
        """
        Create duration object representing input number of minutes
        """
        
    @staticmethod
    def from_hours(d: float) -> duration:
        """
        Create duration object representing input number of hours
        """

    def __add__(self, other: duration|astro.time) -> duration|astro.time:
        """
        Add a duration to either another duration or a time
        
        if "other" is a duration, output is a duration representing the
        sum, or concatenation, of both durations
        
        if "other" is a time, output is a time representing the input 
        time plus the duration    
        
        Example 1: 
        print(duration.from_hours(1) + duration.from_minutes(1))
        Duration: 1 hours, 1 minutes, 0.000 seconds
        
        Example 2: 
        print(duration.from_hours(1) + astro.time(2023, 6, 4, 11,30,0))
        2023-06-04 13:30:00.000Z
        
        """
        
    def __sub__(self, other: duration) -> duration:
        
        """
        Take the difference between two durations
        example:
        
        print(duration.from_hours(1) - duration.from_minutes(1))
        Duration: 59 minutes, 0.000 seconds
        
        """

    def __mul__(self, other: float) -> duration:
        """
        Multiply (or scale) duration by given value
        
        Example:
        print(duration.from_days(1) * 2.5)
        Duration: 2 days, 12 hours, 0 minutes, 0.000 seconds    
        """
        
    def days(self) -> float:
        """
        Floating point number of days represented by duration
        """
        
    def hours(self) -> float:
        """
        Floating point number of hours represented by duration
        """
        
    def minutes(self) -> float:
        """
        Floating point number of minutes represented by duration
        """
        
    def seconds(self) -> float:
        """
        Floating point number of seconds represented by duration
        """



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