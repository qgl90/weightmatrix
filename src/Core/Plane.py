"""
Plane definition: measurement + material + acceptance.

A `Plane` represents a z-location where:

- the track state can be evaluated
- optional measurement information can be applied (weight matrix `WMeas`)
- passive material contributes to multiple scattering via an x/X0 thickness

Planes are created from:

- CSV detector parameterization (`detectionlayers`)
- GDML/TGeo-derived segments inserted by the material locator

Geometry acceptance
------------------
For CSV-defined planes, an `xyshape` tag selects a pickled shapely polygon under
the `Geometries/` folder. The `PicklePlane.InAcceptance(x,y)` query is used to
decide if the track has a valid state at that plane.

Measurement model
-----------------
Two cases are supported:

- strip plane: `sigmaX` and `angle` (stereo angle in degrees)
- pixel plane: `sigmaX` and `sigmaY` (optional `angle` overrides)

The measurement information is stored as a 5x5 weight matrix `WMeas` in the
state basis `(x, y, tx, ty, q/p)`.
"""

from logging import Logger
from multiprocessing import Value
from hepunits.units import centimeter, centimeter2, gram , millimeter
from termcolor  import colored
#from typing     import Self
from .Matrix    import Matrix
from .Precision import array_d
import math 
from .TrackState import TrackState
from shapely import Polygon, Point
from shapely.geometry import Polygon, Point
from shapely import affinity
# from shapely.plotting import plot_polygon, plot_points
from shapely.ops import unary_union
from shapely import geometry
from shapely.ops import unary_union
import matplotlib.pyplot as plt
import numpy as np
import pickle 
import pandas as pd 


import os
from pathlib import Path
def _find_repo_root(start: Path) -> Path | None:
    for candidate in [start, *start.parents]:
        if (candidate / ".git").is_dir():
            return candidate
    return None


def _default_geometry_folder() -> str:
    # Prefer explicit overrides; support historical typo as well.
    env = os.environ.get("WEIGHTMATRIX_GEOMETRY_DIR") or os.environ.get("WEIGTHMATRIX_GEOMETRY_DIR")
    if env:
        return str(Path(env).expanduser())
    here = Path(__file__).resolve()
    repo_root = _find_repo_root(here)
    if repo_root is not None:
        return str(repo_root / "Geometries" )
    # Fallback when running from an unpacked tree without `.git`
    return str(here.parents[2] / "Geometries")


_GEOMETRY_FOLDER_ = _default_geometry_folder()
from typing import TypeVar
TPlane = TypeVar("TPlane", bound =  "Plane")
"""Simple Velo Geometry with an inner and outer radius configurable"""


### names with GDML are effectively 'any' since the x-y acceptance if from GDML crossing
"""
A bit of naming : 
VeloGDML    = GDML loaded Velo active regions 
FTGDML      = GDML loaded SciFi active regions 
UPGDML      = GDML loaded UP active regions 
UTGDMLInner = GDML loaded UT active regions inner having a better resolution
UTGDMLOuter = GDML loaded UT active regions outer having a worse resolution
MPGDML      = GDML loaded MightyPixel active region
GDML        = Generic Dead material from GDML or Measurement found 
"""
_AVAILABLE_GEOMETRIES_GDMLS_ = [ "VeloGDML", "FTGDML", "UPGDML","UTGDMLInner","UTGDMLOuter", "MPGDML", "GDML" ]
_AVAILABLE_GEOMETRIES_ = [ 
    "VeloGDML", "FTGDML", "UPGDML","UTGDMLInner","UTGDMLOuter", "MPGDML", "GDML",
    "VeloClose",
    "VeloLeft", 
    "VeloRight",  
    "VeloLeftOpen", 
    "VeloRightOpen", 
    "VeloOpen", 
    "SciFi_x"       , "SciFi_u"       , "SciFi_v",
    "FTDR_Fibre_x"  , "FTDR_Fibre_u"  , "FTDR_Fibre_v"  ,  "FTDR_Pixel_x",
    "Modest_Fibre_x", "Modest_Fibre_u", "Modest_Fibre_v",  "Modest_Pixel_x", 
    "Frugal_Fibre_x", "Frugal_Fibre_u", "Frugal_Fibre_v",  "Frugal_Pixel_x", 
    "Blake_Fibre_x", "Blake_Fibre_u", "Blake_Fibre_v",  "Blake_Pixel_x",     
    "any" , "UTaX", "UTaU", "UTbX", "UTbV", "UT_U2", "UT_U2_Wide", "UT_U2_BorderLess"
]
_VELO_GEOMETRIES_ = [  "VeloClose",  "VeloLeft",  "VeloRight",  "VeloLeftOpen",  "VeloRightOpen", "VeloOpen", "VeloGDML" ]
_UT_GEOMETRIES_   = [   "UTaX", "UTaU", "UTbX", "UTbV", "UT_U2", "UT_U2_Wide" , "UT_U2_BorderLess",
                        "UPGDML","UTGDMLInner", "UTGDMLOuter" ]
_DT_GEOMETRIES_   = [   "SciFi_x",      "SciFi_u"       , "SciFi_v",
                        "FTDR_Fibre_x"  ,"FTDR_Fibre_u" , "FTDR_Fibre_v", "FTDR_Pixel_x",
                        "Modest_Fibre_x", "Modest_Fibre_u", "Modest_Fibre_v",  "Modest_Pixel_x", 
                        "Frugal_Fibre_x", "Frugal_Fibre_u", "Frugal_Fibre_v",  "Frugal_Pixel_x", 
                        "Blake_Fibre_x", "Blake_Fibre_u", "Blake_Fibre_v",  "Blake_Pixel_x",                             
                        "MPGDML" , "FTGDML"]


def _GEOMETRY_PICKLE_LOAD_(name : str) : 
    fileLoc = Path(_GEOMETRY_FOLDER_) / f"{name}.pickle"
    with open(fileLoc, "rb") as fileIn:
        shape = pickle.load( fileIn)        
    return shape

class PicklePlane : 
    def __init__(self, zval, name):
        if name not in _AVAILABLE_GEOMETRIES_: 
            Logger.fatal( f"Cannot load pickle file of geometry {name}, available = {_AVAILABLE_GEOMETRIES_}")
            raise ValueError("Unavailable Geometry")
        self.z        = zval 
        self.TAG     = name 
        if name not in ["any", "VeloGDML", "FTGDML","UPGDML","UTGDMLInner","UTGDMLOuter","UPGDML", "MPGDML"]: 
            # print(f"AT z= {zval} , Loading Geometry from Pickle with name = {name}")
            self.geometry = _GEOMETRY_PICKLE_LOAD_( self.TAG )
        else : 
            # print(f"AT z= {zval} , Loading Geometry from Pickle with name = {name}")
            self.geometry = None 
    def InAcceptance( self, x : float , y: float ): 
        """Return wether a point is in the acceptance"""
        if self.geometry != None    : return self.geometry.contains( Point(x,y))        
        elif self.TAG == "any"      : return True
        elif self.TAG in _AVAILABLE_GEOMETRIES_GDMLS_ : return True        

        raise ValueError("Cannot call InAcceptance, TAG is invalid!")        
    def IsVelo( self): 
        return self.TAG in _VELO_GEOMETRIES_
    def IsUpstreamTracker(self):
        return self.TAG in _UT_GEOMETRIES_
    def IsDownstreamTracker(self):
        return self.TAG in _DT_GEOMETRIES_ 

def _CREATEWMEAS_( sigma, angle) : 
    wmeas = 1./( sigma*sigma)
    angle_rad  = np.deg2rad( angle)
    cos_stereo = np.cos( angle_rad)
    sin_stereo = np.sin( angle_rad)            
    if angle == 0. : cos_stereo, sin_stereo = 1. , 0. #sometimes the function gives no sense close to 0.
    if angle ==90. : cos_stereo, sin_stereo = 0. , 1.  #sometimes the function gives no sense close to 0.
    if angle ==180.: cos_stereo, sin_stereo =-1. , 0. #sometimes the function gives no sense close to 0.
    if angle ==270.: cos_stereo, sin_stereo = 0. ,-1.  #sometimes the function gives no sense close to 0.
    if sigma == 0. : wmeas = 0.            
    WMeas = Matrix( array_d( [                  [ wmeas*cos_stereo*cos_stereo, wmeas*cos_stereo*sin_stereo,0,0,0],
                                                [ wmeas*cos_stereo*sin_stereo, wmeas*sin_stereo*sin_stereo,0,0,0],
                                                [0,0,0,0,0], 
                                                [0,0,0,0,0], 
                                                [0,0,0,0,0 ] ] ) )
    return WMeas
    

class Plane:
    def __init__(self , zval      : float, \
                        thickness : float, \
                        dzMat     : float, \
                        angle     : float=None ,  
                        sigmaX    : float=None ,  
                        sigmaY    : float = None, #either supply sigma + angle, or  sigmaX, sigmaY, angle                         
                        geomLayer : str  ="any", 
                        hitefficiency = 1.0,
                        name = "" , 
                        dzMatAssignZero = False 
                        ):
        geomLayer = geomLayer.replace(" ","")
        """
        Plane defining the position, thickness of the materia, and shape in XY for that specific layer 
        and storing the Covariance of measurements done and what the Stateless Kalman does
        Basic inputs here...
        """
        self.name = name.replace(" ","")
        self.Z            = zval
        self.dzMat = dzMat
        # if dzMat == 0. and not(dzMatAssignZero) : self.dzMat = None 
        # else :           self.dzMat = dzMat
        self.Thickness        = thickness #in radiation length .... ( maybe density material to pass ? ) 
        self.WMeas            = Matrix()
        self.WPredictForward  = Matrix()
        self.WUpdateForward   = Matrix()
        self.WPredictBackward = Matrix()
        self.WUpdateBackward  = Matrix()
        self.IsMeasurement    = False
        self.associatedState  = None 
        """
        sigma X is always wrt the y vertical line
        sigma Y is always orthogonal to the y vertical line, it's as sigma X + 90 deg rotation
        If both are set, we have a pixel measurement, else it's a strip one with angle, sigmaX defined angles used
        """
        self.sigmaX           = sigmaX if not(pd.isnull( sigmaX)) else  None 
        self.sigmaY           = sigmaY if not(pd.isnull( sigmaY)) else  None 
        self.angle            = angle  if not(pd.isnull(  angle)) else  None 
        self.geomLayer        = geomLayer
        if self.geomLayer in _AVAILABLE_GEOMETRIES_ : 
            self.geom = PicklePlane( zval = self.Z, name = self.geomLayer)        
        else :
            print(self.geomLayer )
            Logger.fatal(f"Cannot create detector plane with name {self.name}, available are {_AVAILABLE_GEOMETRIES_+['any']}")
            raise ValueError("Invaid xyshape passed, fix me!")
        self.hiteff= hitefficiency
        
        # either a strip or a direct pixel
        # if we have sigmaX and Y is none, must have a angle passed
        # if we have sigmaX and sigmaY, angle is not needed since we can deduce it from the 2 sigmas, 
        # but if angle is passed, we will use it to compute the WMeas matrix and not the one of a 0,90 degree rotation
        self.IsMeasurement = (sigmaX != None and sigmaY == None and angle != None ) or \
                             (sigmaX != None and sigmaY!=None and angle ==   None ) or \
                             (sigmaX != None and sigmaY!=None and angle !=   None )
        ### TODO : work out the mathatic of Wmeas + Wmeas of a x-y 0,90 degree
        ### CONTINUE FROM HERE! 
        # If we have sigmaX and an angle passed (strip measurement treatment)
        if self.sigmaX != None and self.angle != None  and self.sigmaY == None : 
            #strip detector case, we have an angle and a sigmaX, we can compute the WMeas matrix with the angle and sigmaX
            self.WMeas = _CREATEWMEAS_( sigma= sigmaX, angle = angle)                     
            self.IsMeasurement = True
        elif self.sigmaX != None and self.sigmaY != None :
            #pixel detector case, we have sigmaX and sigmaY, 
            # we can compute the WMeas matrix with the angle if passed, or deduce it from the 2 sigmas if not passed (angle = arctan(sigmaY/sigmaX))
            if angle is None : 
                self.WMeas = _CREATEWMEAS_( sigma = self.sigmaX , angle = 0 ) + \
                             _CREATEWMEAS_( sigma = self.sigmaY,  angle = 90.)
            else :
                self.WMeas = _CREATEWMEAS_( sigma= self.sigmaX, angle = angle) + \
                             _CREATEWMEAS_( sigma= self.sigmaY, angle = angle + 90.)
            self.IsMeasurement = True 
        else : 
            self.IsMeasurement = False 
            self.WMeas = None 
    def __repr__(self):
        if self.IsMeasurement :             
            return colored(f"""{'Pixel' if self.sigmaX != None and self.sigmaY != None else 'Strip'} <MEASUREMENT PROVIDER: Plane z:{self.Z:.2f} mm, x/X0 :{self.Thickness*100.}% , dzMat :{self.dzMat} mm, sigmax :{self.sigmaX} mm, sigmay :{self.sigmaY} mm, hiteff :{self.hiteff *100:.2f} % , geomLayer :{self.geomLayer}, geomObjAssociated: {self.geom!=None} state {self.associatedState}>""","blue")
        else : 
            return colored(f"""<PASSIVE     PROVIDER: Plane z:{self.Z:.2f} mm, x/x0 :{self.Thickness*100.}% , dzMat :{self.dzMat} mm, sigmax :{self.sigmaX} mm, sigmay :{self.sigmaY} um, hiteff :{self.hiteff *100:.2f} % , geomLayer :{self.geomLayer}, geomObjAssociated: {self.geom!=None} >""","yellow")
    def __str__(self):
        if self.IsMeasurement :
            return colored(f"""{'Pixel' if self.sigmaX != None and self.sigmaY != None else 'Strip'} MEASUREMENT PROVIDER: <Plane z:{self.Z:.2f} mm, x/X0 :{self.Thickness*100.}% , dzMat :{self.dzMat} mm, sigmax :{self.sigmaX} mm, sigmay :{self.sigmaY} mm, hiteff :{self.hiteff *100:.2f} % , geomLayer :{self.geomLayer}, geomObjAssociated: {self.geom!=None}, state {self.associatedState}>""","blue")
        else : 
            return colored(f"PASSIVE     PROVIDER: <Plane z:{self.Z:.2f} mm, x/x0 :{self.Thickness*100.}% , dzMat :{self.dzMat} mm, sigmax :{self.sigmaX} mm, sigmay :{self.sigmaY} um, hiteff :{self.hiteff *100:.2f} % , geomLayer :{self.geomLayer}, geomObjAssociated: {self.geom!=None} >","yellow")
    def IsPixelMeasurement(self): 
        return self.IsMeasurement and self.sigmaX != None and self.sigmaY != None 
    def IsVelo(self): 
        return self.geom.IsVelo()
    def IsUT(self): 
        return self.geom.IsUpstreamTracker()
    def IsDownStreamTracker(self):
        return self.geom.IsDownstreamTracker()    
    def IsDownStreamFibreTracker(self):
        return self.geom.IsDownstreamTracker() and self.IsMeasurement and not(self.IsPixelMeasurement())
    def IsDownStreamPixelTracker(self):
        return self.geom.IsDownstreamTracker() and self.IsMeasurement and (self.IsPixelMeasurement())
    
    def PointInAcceptance(self, x : float, y: float) :
        return self.geom.InAcceptance(x,y)
    def HasStateInAcceptance(self ): 
        """Should be always true for "ANY" geometry """
        return self.geom.InAcceptance( self.associatedState.X(), self.associatedState.Y() )
    def associateState( self, trState : TrackState) :
        self.associatedState = trState    
    def GetName(self):
        return self.name    
    def __lt__(self, other : TPlane) -> bool :
        """override it so that we can sort the planes"""           
        if self.Z == other.Z :
            # if self.Thickness == other.Thickness : 
            if self.dzMat != None and other.dzMat != None : 
                return self.dzMat < other.dzMat 
            else : 
                return self.Z <= other.Z
            # return self.Thickness < other.Thickness
        return self.Z <= other.Z

    def GetPolygon( self): 
        return self.geom.geometry
