"""
Magnetic field map loader and interpolator.

The field map is provided as a CSV on a regular (x,y,z) grid with columns:

- `x`, `y`, `z` (grid coordinates)
- `Bx`, `By`, `Bz` (field components)

The loader builds three ROOT `TH3D` histograms and provides a single query API:

    `BField.GetFieldVector(x_mm, y_mm, z_mm) -> (Bx, By, Bz)`

Notes
-----
- Positions are treated as **mm**.
- Field components are treated as **Tesla**.
- `GridFactor` rescales coordinates, `VectorFactor` rescales field values.
- The map must be regular with the same spacing in x, y and z (enforced).
"""

import pandas as pd 
import ROOT   as r
import numpy  as np 
from Utils.Logger import Logger

class BField( ):
    def __init__(self, csvfilemap= "../", GridFactor=1., VectorFactor=1., startScalingAfterZ = None ): 
        """constructor passing a csv file"""
        self.gridFactor = GridFactor
        self.vectFactor = VectorFactor
        

        self.filemap  = csvfilemap
        self.bfield   = pd.read_csv( csvfilemap)
        
        self.x_min_max =  ( np.min(self.bfield["x"]*self.gridFactor) , np.max(self.bfield["x"]*self.gridFactor))
        self.y_min_max =  ( np.min(self.bfield["y"]*self.gridFactor) , np.max(self.bfield["y"]*self.gridFactor))
        self.z_min_max =  ( np.min(self.bfield["z"]*self.gridFactor) , np.max(self.bfield["z"]*self.gridFactor))       
        self.grid_vals =  ( np.unique(self.bfield["x"]*self.gridFactor) ,\
                            np.unique(self.bfield["y"]*self.gridFactor) ,\
                            np.unique(self.bfield["z"]*self.gridFactor))
        
        self.grid_step = (  (self.grid_vals[0][-1]-self.grid_vals[0][0])/float( len(self.grid_vals[0])-1)  , \
                            (self.grid_vals[1][-1]-self.grid_vals[1][0])/float( len(self.grid_vals[1])-1)  , \
                            (self.grid_vals[2][-1]-self.grid_vals[2][0])/float( len(self.grid_vals[2])-1)  )
        Logger.debug(f"BField::GridFactor {self.gridFactor}")        
        Logger.debug(f"BField::VectFactor {self.vectFactor}")                
        Logger.debug(f"BField::GridSteps  {self.grid_step}")
        Logger.debug(f"BField::XMinMax    {self.x_min_max}")
        Logger.debug(f"BField::YMinMax    {self.y_min_max}")
        Logger.debug(f"BField::ZMinMax    {self.z_min_max}")
        Logger.debug(f"BField::XGrid      {self.grid_vals[0]}")
        Logger.debug(f"BField::YGrid      {self.grid_vals[1]}")
        Logger.debug(f"BField::ZGrid      {self.grid_vals[2]}")

        if( self.grid_step[0] != self.grid_step[1] ): raise Exception("Field map must be regular in x,y")
        if( self.grid_step[1] != self.grid_step[2] ): raise Exception("Field map must be regular in y,z")
        if( self.grid_step[0] != self.grid_step[2] ): raise Exception("Field map must be regular in x,z")
        
        self.hBx = r.TH3D( "FieldMap_Bx", "FieldMap LHCb Bx", len(self.grid_vals[0]) ,
                                             self.grid_vals[0][0]  - self.grid_step[0]/2.,
                                             self.grid_vals[0][-1] + self.grid_step[0]/2.,
                                             len(self.grid_vals[1]) , 
                                             self.grid_vals[1][0]  - self.grid_step[1]/2., 
                                             self.grid_vals[1][-1] + self.grid_step[1]/2.,
                                             len(self.grid_vals[2]) , 
                                             self.grid_vals[2][0]  - self.grid_step[2]/2., 
                                             self.grid_vals[2][-1] + self.grid_step[2]/2.)                                            
        self.hBy = r.TH3D( "FieldMap_By", "FieldMap LHCb By", len(self.grid_vals[0]) ,
                                             self.grid_vals[0][0]  - self.grid_step[0]/2.,
                                             self.grid_vals[0][-1] + self.grid_step[0]/2.,
                                             len(self.grid_vals[1]) , 
                                             self.grid_vals[1][0]  - self.grid_step[1]/2., 
                                             self.grid_vals[1][-1] + self.grid_step[1]/2.,
                                             len(self.grid_vals[2]) , 
                                             self.grid_vals[2][0]  - self.grid_step[2]/2., 
                                             self.grid_vals[2][-1] + self.grid_step[2]/2.)     
        self.hBz = r.TH3D( "FieldMap_Bz", "FieldMap LHCb Bz", len(self.grid_vals[0]) ,
                                             self.grid_vals[0][0]  - self.grid_step[0]/2.,
                                             self.grid_vals[0][-1] + self.grid_step[0]/2.,
                                             len(self.grid_vals[1]) , 
                                             self.grid_vals[1][0]  - self.grid_step[1]/2., 
                                             self.grid_vals[1][-1] + self.grid_step[1]/2.,
                                             len(self.grid_vals[2]) , 
                                             self.grid_vals[2][0]  - self.grid_step[2]/2., 
                                             self.grid_vals[2][-1] + self.grid_step[2]/2.)
        
        for x,y,z, bx, by, bz in zip( self.bfield["x"] , 
                                      self.bfield["y"] , 
                                      self.bfield["z"] , 
                                      self.bfield["Bx"] , 
                                      self.bfield["By"] ,                               
                                      self.bfield["Bz"] ):
            self.hBx.Fill( x*self.gridFactor,y*self.gridFactor,z*self.gridFactor , bx*self.vectFactor)
            self.hBy.Fill( x*self.gridFactor,y*self.gridFactor,z*self.gridFactor , by*self.vectFactor)
            self.hBz.Fill( x*self.gridFactor,y*self.gridFactor,z*self.gridFactor , bz*self.vectFactor)
        return
    def GetFieldVector(self, x: float, y: float, z: float) -> tuple[float,float,float]:
        if x < self.x_min_max[0] or x > self.x_min_max[1] : return ( 0.,0.,0.)
        if y < self.y_min_max[0] or y > self.y_min_max[1] : return ( 0.,0.,0.)
        if z < self.z_min_max[0] or z > self.z_min_max[1] : return ( 0.,0.,0.)
        
        Bx_ValBin = self.hBx.GetBinContent( self.hBx.FindFixBin( x,y,z))
        By_ValBin = self.hBy.GetBinContent( self.hBy.FindFixBin( x,y,z))
        Bz_ValBin = self.hBz.GetBinContent( self.hBz.FindFixBin( x,y,z))
        
        #handle edges of 3-d grid map        
        if x < self.x_min_max[0] or    x > self.x_min_max[1] : return (0.,0.,0.)
        if x < self.grid_vals[0][0] or x > self.grid_vals[0][-1] : ( Bx_ValBin, By_ValBin, Bz_ValBin)
        
        if y < self.y_min_max[0] or y > self.y_min_max[1] : return (0.,0.,0.)
        if y < self.grid_vals[1][0] or y > self.grid_vals[1][-1] : ( Bx_ValBin, By_ValBin, Bz_ValBin)
        if z < self.z_min_max[0] or z    > self.z_min_max[1] : return (0.,0.,0.)
        if z < self.grid_vals[2][0] or y > self.grid_vals[2][-1] : ( Bx_ValBin, By_ValBin, Bz_ValBin)        
        
        # if bigger than 11.45 meters 
        if z > 11450. : 
            return (0,0,0)
        
        # Logger.info(f"Interpolating BField at x={x:.1f} mm, y={y:.1f} mm, z={z:.1f} mm --> Bx={Bx_ValBin:.3f} T, By={By_ValBin:.3f} T, Bz={Bz_ValBin:.3f} T")
        return ( self.hBx.Interpolate( x,y,z) , self.hBy.Interpolate( x,y,z), self.hBz.Interpolate( x,y,z))
