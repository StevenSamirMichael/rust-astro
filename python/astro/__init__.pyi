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
from . import jplephem
from . import frametransform

__all__ = ['time', 'duration', 'timescale', 'quaternion', 'sgp4', 'gravmodel', 'gravity', 'nrlmsise00', 'univ', 'solarsystem', 'TLE', 'itrfcoord', 'frametransform', 'lpephem', 'satprop', 'jplephem', 'utils']


class TLE():
    """    
    Stucture representing a Two-Line Element Set (TLE), a satellite
    ephemeris format from the 1970s that is still somehow in use
    today and can be used to calculate satellite position and
    velcocity in the "TEME" frame (not-quite GCRF) using the
    "Simplified General Perturbations-4" (SGP-4) mathemematical
    model that is also included in this package.

    For details, see: https://en.wikipedia.org/wiki/Two-line_element_set

    The TLE format is still commonly used to represent satellite
    ephemerides, and satellite ephemerides catalogs in this format
    are publicly availalble at www.space-track.org (registration
    required)

    TLEs sometimes have a "line 0" that includes the name of the satellite

    
 
    """

    @typing.staticmethod
    def from_lines(lines: list[str]) -> list[astro.TLE]:
        """
        Return a list of TLEs from input lines, represented as a 
        list of strings
        """

    @typing.staticmethod
    def single_from_lines(lines: list[str]) -> astro.TLE:
        """
        Return a single TLE from a 2 or 3-element list of lines
        
        If additional lines are included, they are ignored.
        """
        
    @property
    def satnum(self) -> int:
        """
        Satellite number, or equivalently the NORAD ID
        """
        
    @property
    def eccen(self) -> float:
        """
        Satellite eccentricity, in range [0,1]
        """
        
    @property
    def mean_anomaly(self) -> float:
        """
        Mean anomaly in degrees
        """
        
    @property
    def mean_motion(self) -> float:
        """
        Mean motion in revs / day
        """
        
    @property
    def inclination(self) -> float:
        """
        Inclination, in degrees
        """
    @property
    def epoch(self) -> astro.time:
        """
        TLE epoch
        """
        
    @property
    def arg_of_perigee(self) -> astro.time:
        """
        Argument of Perigee, in degrees
        """
        
    @property
    def mean_motion_dot(self) -> float:
        """
        1/2 of first derivative of mean motion, in revs/day^2
        
        the "1/2" is because that is how number is stored in the TLE
        """
        
    @property
    def mean_motion_dot_dot(self) -> float:
        """
        1/6 of 2nd derivative of mean motion, in revs/day^3
        
        the "1/6" is because that is how number is stored in the TLE
        """
    
    @property
    def name(self) -> str:
        """
        The name of the satellite
        """

    @property
    def bstar(self) -> str:
        """
        "B Star" or drag of the satellite
        
        should be rho0 * Cd * A / 2 / m
        
        Units (which are strange) is multiples of 
        1 / Earth radius
        """

def sgp4(
    tle: astro.TLE, 
    tm: astro.time|list[astro.time]|npt.ArrayLike[astro.time]
    ) -> (npt.ArrayLike[np.float64], npt.ArrayLike[np.float64]):
    """
    Run Simplified General Perturbations (SGP)-4 propagator on
    Two-Line Element Set to
    output satellite position and velocity at given time
    in the "TEME" coordinate system

    A detailed description is at:
    https://celestrak.org/publications/AIAA/2008-6770/AIAA-2008-6770.pdf


    # Arguments

    tle: The TLE on which top operate.
    
    tm: astro.time object or list of objects or numpy array of 
        objects representimg time(s) at which to compute
        position and velocity
          

    # Return

    tuple with the following elements:
    
    0 : a Nx3 numpy array representing position in meters in the TEME frame at 
        each of the "N" input times (or 3-element array for a single time)
        
    1 : a Nx3 numpy array representing velocity in meters / second in the TEME
        frame at each of the "N" input times
        (or 3-element array for a single time)

    Example usage: show Geodetic position of satellite at TLE epoch
    
    lines = [
        "0 INTELSAT 902",
        "1 26900U 01039A   06106.74503247  .00000045  00000-0  10000-3 0  8290",
        "2 26900   0.0164 266.5378 0003319  86.1794 182.2590  1.00273847 16981   9300."
    ]


    tle = astro.TLE.single_from_lines(lines)

    # Compute TEME position & velocity at epoch
    pteme, vteme = astro.sgp4(tle, tle.epoch)

    # Rotate to ITRF frame
    q = astro.frametransform.qteme2itrf(tm)
    pitrf = q * pteme
    vitrf = q * vteme - np.cross(np.array([0, 0, astro.univ.omega_earth]), pitrf)

    # convert to ITRF coordinate object
    coord = astro.itrfcoord.from_vector(pitrf)
    # Print ITRF coordinate object location
    print(coord)

    Output:
    
    ITRFCoord(lat:  -0.0363 deg, lon:  -2.2438 deg, hae: 35799.51 km)

    """    
    

class gravmodel():
    """
    Earth gravity models available for use
    
    For details, see: http://icgem.gfz-potsdam.de/    
    """
    
    @property
    def jgm3() -> int:
        """
        The "JGM3" gravity model
        
        This model is used by default in the orbit propagators
        """
        
    @property
    def jgm2() -> int:
        """
        The "JGM2" gravity model
        """
        
    @property
    def egm96() -> int:
        """
        The "EGM96" gravity model
        """

    @property
    def itugrace16() -> int:
        """
        the ITU Grace 16 gravity model
        """

def gravity(pos: list[float]|astro.itrfcoord|npt.ArrayLike[np.float],
            **kwargs) -> npt.ArrayLike[np.float]:
    """
    gravity(pos)
    --

    Return acceleration due to Earth gravity at the input position. The
    acceleration does not include the centrifugal force, and is output
    in m/s^2 in the International Terrestrial Reference Frame (ITRF)

    Inputs:

        pos:   Position as ITRF coordinate (astro.itrfcoord) or numpy
                3-vector representing ITRF position in meters or 
                list 3-vector representing ITRF position in meters

    Kwargs:
        
        model:   The gravity model to use.  Options are:
                    astro.gravmodel.jgm3
                    astro.gravmodel.jgm2
                    astro.gravmodel.egm96
                    astro.gravmodel.itugrace16

                Default is astro.gravmodel.jgm3

                For details of models, see:
                http://icgem.gfz-potsdam.de/tom_longtime

        order:    The order of the gravity model to use.
                Default is 6, maximum is 16


                For details of calculation, see Chapter 3.2 of:
                "Satellite Orbits: Models, Methods, Applications",
                O. Montenbruck and B. Gill, Springer, 2012.
        
        """    

class solarsystem():
    """
    Solar system bodies for which high-precision ephemeris can be computed    
    """
    
    @property
    def Mercury() -> int:
        """
        Mercury
        """
    
    @property
    def Venus() -> int:
        """
        Venus
        """
        
    @property
    def EMB() -> int:
        """
        Earth-Moon Barycenter
        """        
        
    @property
    def Mars() -> int:
        """
        Mars
        """
        
    @property
    def Jupiter() -> int:
        """
        Jupter
        """
        
    @property
    def Saturn() -> int:
        """
        Saturn
        """
        
    @property
    def Uranus() -> int:
        """
        Uranus
        """
        
    @property
    def Neptune() -> int:
        """
        Neptune
        """
    
    @property
    def Pluto() -> int:
        """
        Pluto
        """
        
    @property
    def Moon() -> int:
        """
        Moon
        """
        
    @property
    def Sun() -> int:
        """
        Sun
        """

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

class quaternion():
    """
    Quaternion representing rotation of 3D Cartesian axes

    Quaternion is right-handed rotation of a vector,
    e.g. rotation of +xhat 90 degrees by +zhat give +yhat

    This is different than the convention used in Vallado, but
    it is the way it is commonly used in mathematics and it is
    the way it should be done.

    For the uninitiated: quaternions are a more-compact and
    computationally efficient way of representing 3D rotations.  
    They can also be multipled together and easily renormalized to
    avoid problems with floating-point precision eventually causing
    changes in the rotated vecdtor norm.

    For details, see:

    https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation    
    
    Under the hood, this is using the "UnitQuaternion" object in the 
    rust "nalgebra" crate.
    """
    
    def new() -> astro.quaternion:
        """
        Return unit quaternion (no rotation)
        """
      
    @staticmethod
    def from_axis_angle(axis: npt.ArrayLike[np.float64], angle: float) -> astro.quaternion:
        """
        Return quaternion representing right-handed rotation by 
        "angle" degrees about the given axis.  The axis does not 
        have to be normalized.
        """
      
    @staticmethod  
    def rotx(theta) -> astro.quaternion:
        """
        Return quaternion representing right-handed rotation of vector
        by "theta" degrees about the xhat unit vector
        
        Equivalent rotation matrix:
        | 1             0            0|
        | 0    cos(theta)   sin(theta)|
        | 0   -sin(theta)   cos(theta)|        
        """
        
      
    @staticmethod  
    def roty(theta) -> astro.quaternion:
        """
        Return quaternion representing right-handed rotation of vector
        by "theta" degrees about the yhat unit vector
        
        Equivalent rotation matrix:
        | cos(theta)     0   -sin(theta)|
        |          0     1             0|
        | sin(theta)     0    cos(theta)|          
        """        
    
    @staticmethod  
    def rotz(theta) -> astro.quaternion:
        """
        Return quaternion representing right-handed rotation of vector
        by "theta" degrees about the zhat unit vector
        
        Equivalent rotation matrix:
        |  cos(theta)      sin(theta)   0|
        | -sin(theta)      cos(theta)   0|
        |           0               0   1|
        """ 
    
    def to_rotation_matrix(self) -> npt.ArrayLike[np.float64]:
        """
        Return 3x3 rotation matrix representing equivalent rotation
        """
        
    def to_euler(self) -> (float, float, float):
        """
        Return equivalent rotation angle represented as rotation angles:
        ("roll", "pitch", "yaw") in radians
        roll = rotation about x axis
        pitch = rotation about y axis
        yaw = rotation about z axis
        """
        
    def angle(self) -> float:
        """
        Return the angle in radians of the rotation
        """
        
    def axis(self) -> npt.ArrayLike(np.float64):
        """
        Return the axis of rotation as a unit vector
        """
        
    def conj(self) -> astro.quaternion:
        """
        Return conjucate or inverse of the rotation
        """
        
    def conjugate(self) -> astro.quaternion:
        """
        Return conjugate or inverse of the rotation
        """

    @typing.overload
    def __mul__(self, other: astro.quaternion) -> astro.quaternion:
        """
        Multiply represents concatenation of two rotations representing 
        the quaternions.  The left value rotation is applied after 
        the right value, per the normal convention
        """
        
    @typing.overload
    def __mul__(self, other: npt.ArrayLike[np.float64]) -> npt.ArrayLike[np.float64]:
        """
        Multply by a vector to rotate the vector
        
        The vector is represented as a numpy array
        
        If the array is 1 demensional it must have 3 elements
        
        If the array is 2 dimensionsl and the dimensions are Nx3, 
        each of the "N" vectors is rotated by the quaternion and a
        Nx3 array is returned
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

    def __init__(x: float, y: float, z: float) -> astro.itrfcoord:
        """
        Represent a coordinate in the ITRF (International Terrestrial Reference Frame)
        Inputs are 3 separate floats representing ITRF Cartesian position in meters
        """ 
        
    @typing.static_method
    def from_vector(v: list[float]|npt.ArrayLike[np.float64]) -> astro.itrfcoord:
        """
        Create ITRF coordinate from 3-element list or numpy array
        representing ITRF cartesian position in meters
        """
            
    @typing.staticmethod                
    def from_geodetic() -> astro.itrfcoord:
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
