"""
High-precision satellite propagation using Runga-Kutta methods for
differential equation solving

Force models are mostly pulled from Montenbruck & Gill:
https://link.springer.com/book/10.1007/978-3-642-58351-3

The propagator can also compute the state transition matrix, meaning
position and velocity covariances can be propagated as well.

The default propagator uses a Runga-Kutta 9(8) integrator 
with coefficient computed by Verner:
https://www.sfu.ca/~jverner/

This works much better than lower-order Runga-Kutta solvers
such as Dormund-Prince, and I don't know why it isn't more
popular in numerical packages

This module includes a function to propagate a position and time directly,
and a convenience "satstate" object that represents satellite position, velocity,
and optionally covariance and can propagate itself to different times 

Forces included in the propagator:

1. Earth gravity with higher-order zonal terms
2. Sun, Moon gravity
3. Radiation pressure
4. Atmospheric drag: NRL-MISE 2000 density model, with option
to include space weather effects (can be large)

"""

from __future__ import annotations
import typing
import numpy.typing as npt
import numpy as np

import astro
import datetime

class PropStats(typing.TypedDict):
    num_eval: int
    accepted_steps: int
    rejected_steps: int

class PropResult(typing.TypedDict):
    time: npt.ArrayLike[astro.time]
    pos: npt.ArrayLike[np.float64]
    vel: npt.ArrayLike[np.float64]
    Phi: typing.NotRequired[npt.ArrayLike[np.float64]]
    stats: PropStats
    

class propsettings():
    
    @property
    def abs_error() -> float:
        """
        Maxmum absolute value of error for any element
        in propagated state following ODE integration
        
        Default: 1e-8
        """
    
    @property
    def rel_error() -> float:
        """
        Maximum relative error of any element in 
        propagated state following ODE integration
        
        Default: 1e-8
        """
        
    @property
    def gravity_order() -> int:
        """
        Earth gravity order to use in ODE integration
        
        Default: 4
        """
        
    @property
    def use_spaceweather() -> bool:
        """
        Use historical space weather data when computing
        atmospheric density for drag forces
        
        Default: true
        """
        
    @property
    def use_jplephem() -> bool:
        """
        Use high-precision but computationally expensive JPL
        ephemerides for sun and mun when computing their 
        gravitational force
        """
        
def propagate(pos: npt.ArrayLike(np.float64),
            vel: npt.ArrayLike(np.float64),
            start: astro.time,
            **kwargs
            ) -> PropResult:
    """
    High-precision orbit propagator
    
    Propagator uses advanced Runga-Kutta integrators and includes the following
    forces:
    
        1) Earth gravity, with zonal gravity up to order 16 (default is 4)
        2) Gravitational force of moon
        3) Gravitational force of sun
        4) Solar radiation pressure (with user-specified satellite model)
        5) Atmospheric drag, with correction for space wither
            (with user-specified satellite model)
    
    
    Propagate statellite ephemeris (position, velocity in gcrs & time) to new time
    and output new position and velocity
    
    Inputs and outputs are all in the Geocentric Celestial Reference Frame (GCRF)
    
    Inputs:
        
        pos:   3-element numpy array representing satellite GCRF position in meters
        vel:   3-element numpy array representing satellite GCRF velocity in m/s
            tm:   astro.time object representing instant at which satellite is at "pos" & "vel"
    
    Optional keyword arguments:
    
    
    4 ways of setting propagation end:
    (one of these must be used)
        
            stoptime: astro.time object representing instant at
                        which new position and velocity will be computed
        duration_secs: duration in seconds from "tm" for at which new
                        position and velocity will be computed.  
        duration_days: duration in days from "tm" at which new position and
                        velocity will be computed.  
            duration: An astro.duration object setting duration from "tm"
                        at which new position & velocity will be computed.
    
    
        3 ways of setting smaller interval over which to compute solution:
        (defualt is none, i.e. solution only computed at propagation end)
    
                dt_secs: Interval in seconds between "starttime" and "stoptime"
                        at which solution will also be computed
                dt_days: Interval in days between "starttime" and "stoptime" at which
                        solution will also be computed
                    dt: astro.duration representing interval over which
                        new position & velocity will be computed
    
    
        Other keywords:
    
    
            output_phi: Output 6x6 state transition matrix between "starttime" and
                        "stoptime" (and at intervals, if specified)
                        default is False
        propsettings: "propsettings" object with input settings for
                        the propagation. if left out, default will be used.
        satproperties: "SatPropertiesStatic" object with drag and
                        radiation pressure succeptibility of satellite.
                        If left out, drag and radiation pressure are neglected
                        Dynamic drag & radiation pressure models are not
                        yet implemented
    
    
    Output: Python dictionary with the following elements:
        
        "time": list of astro.time objects at which solution is computed
            "pos": GCRF position in meters at "time".  Output is a Nx3 numpy
                matrix, where N is the length of the output "time" list
            "vel": GCRF velocity in meters / second at "time".  Output is a
                Nx3 numpy matrix, where N is the length of the output
                "time" list
            "Phi": 6x6 State transition matrix corresponding to each time.
                Output is Nx6x6 numpy matrix, where N is the lenght of
                the output "time" list. Not included if output_phi
                kwarg is set to false (the default)
        "stats": Python dictionary with statistics for the propagation.  
                This includes:
                        "num_eval": Number of function evaluations of the force model
                                    required to get solution with desired accuracy
                "accepted_steps": Accepted steps in the adpative Runga-Kutta solver
                "rejected_steps": Rejected steps in the adaptive Runga-Kutta solver
    """
