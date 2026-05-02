"""
Parabolic extrapolator (field-dependent propagation).

This implements a simple parabolic track extrapolation in a magnetic field,
roughly following the LHCb `TrackParabolicExtrapolator` style:

- propagate a 5D state `(x, y, tx, ty, q/p)` from `zOld` to `zNew`
- compute the corresponding 5x5 transport matrix

It is used by `MomentumResolutionCalculator` to:

- create intermediate states along a trajectory
- build transport matrices between successive planes

Notes
-----
- Coordinates are in **mm**.
- The code assumes a constant `c_light = 0.3 mm/ps` and uses a ROOT-interpolated
  field in **Tesla**. Momentum unit conventions are embedded in the numerical
  factors and should be audited before any physics publication usage.
"""

from logging import Logger

from hepunits.units import  meter, second , megaelectronvolt
from .BField import BField
import numpy as np
from Core.TrackState import TrackState
from Core.Matrix import Matrix 
from Core.Precision import array_d
from Utils.Logger import Logger
import math 

CURL_ANGLE = 80.
class ParabolicExtrapolator :
    def __init__(self, FieldMap : BField):
        """
        Args:
            FieldMap (BField): The Magnetic Field map to use for 
            bending tracks in transportation
            B field map must be regular in x,y,z and have the same step in all dimensions, otherwise it will raise an exception
            Also, B field intensity should be in correct units. 
            C_light here is expressed as 0.3 mm/ps 
            We use mm, MeV energy for momentum!            
        """
        self.c_light   = 0.3 #0.3 mm / ps
        self.eplus     = 1.
        self.BFieldMap = FieldMap
    def GetTransportMatrix( self, dz : float, trState : TrackState , BVec : np.array, ax : float , ay : float , noCharge: bool = False  ) -> Matrix  :     
        if dz == 0 : 
            transMat = Matrix( array_d([[ 1., 0 , 0 , 0 , 0] , 
                                        [  0, 1., 0 , 0 , 0] , 
                                        [  0 ,0 , 1., 0 , 0] , 
                                        [  0, 0 , 0 , 1., 0] , 
                                        [  0, 0 , 0 , 0 , 1.]]) )
            return transMat
        Tx    = trState.TX()
        Ty    = trState.TY()
        norm2 = 1. + Tx * Tx + Ty * Ty
        norm  = math.sqrt( norm2 )
        
        
        dAx_dTx = ( Tx * ax / norm2 ) + norm * (  Ty * BVec[0] - ( 2. * Tx * BVec[1] ) )
        dAx_dTy = ( Ty * ax / norm2 ) + norm * (  Tx * BVec[0] + BVec[2] )
        dAy_dTx = ( Tx * ay / norm2 ) + norm * ( -Ty * BVec[1] - BVec[2] )
        dAy_dTy = ( Ty * ay / norm2 ) + norm * ( -Tx * BVec[1] + ( 2. * Ty * BVec[0] ) )

        #1 MeV/c = 5.34428599 × 10-19 (Kg * mm)/s
        #      c = 0.3 mm / ps  = 0.3 * 10^12 mm /  s
        # q/p = 1/ ( X * 5.3 
        fac  = self.eplus * self.c_light * dz
        if noCharge : fac = 0.
        fact = fac * trState.QOP()
        ### p is in MeV/c = 5.36 x 10-22 kg-m/s
        """
        See https://gitlab.cern.ch/lhcb/Rec/-/blob/master/Tr/TrackExtrapolators/src/TrackParabolicExtrapolator.cpp
        """
        transMat = Matrix( array_d([   [  1.,   0   , dz + 0.5 * dAx_dTx * fact * dz, 0.5 * dAx_dTy * fact * dz,   0.5 * ax * fac * dz],
                                       [  0 ,   1.  , 0.5 * dAy_dTx * fact * dz,  dz + 0.5 * dAy_dTy * fact * dz , 0.5 * ay * fac * dz] ,
                                       [  0 ,   0   , 1.0 + dAx_dTx * fact              ,  dAx_dTy * fact                 ,  ax * fac          ] , 
                                       [  0 ,   0   ,       dAy_dTx * fact              ,  1.0 + dAy_dTy * fact           ,  ay * fac          ] ,                     
                                       [  0 ,   0   ,                    0              ,  0                              ,  1.                 ] ] ) )
                        
        return transMat
    def Propagate( self, trState : TrackState , zOld : float , zNew : float , noCharge : bool = False ) :
        if trState.Z() != zOld: 
            raise ValueError("Invalid Propagate Track state old->new must match old to State z")
        
        if zOld == zNew : 
            TrackStateNew    = TrackState( x  =  trState.X()  , \
                                           y  =  trState.Y()  , \
                                           tx =  trState.TX() , \
                                           ty =  trState.TY() , \
                                           qop=  trState.QOP() , \
                                           z  =  zNew )   
            transMatrix = self.GetTransportMatrix( dz=0, trState = TrackStateNew, BVec = np.array([0,0,0]), ax=0, ay=0 )
        else: 
            dz = zNew - zOld
            xMid = trState.X()        + ( 0.5 * trState.TX() * dz)
            yMid = trState.Y()        + ( 0.5 * trState.TY() * dz)        
            zMid = (zOld + zNew)/2.
            point = ( xMid, yMid, zMid)
            Field_At_Point = self.BFieldMap.GetFieldVector( point[0], point[1], point[2])
            tx = trState.TX()
            ty = trState.TY()
            ntx2 = 1.+ tx*tx
            nty2 = 1.+ ty*ty
            norm = math.sqrt( ntx2 + nty2 -1.)
            Bx = Field_At_Point[0]
            By = Field_At_Point[1]
            Bz = Field_At_Point[2]            
            BVec = np.array( [Bx, By, Bz])
            #calculate the A factors
            ax = norm * (  ty * (tx * Bx + Bz) - ntx2 * By)
            ay = norm * ( -tx * (ty * By + Bz) + nty2 * Bx)            
            fac  = abs(self.eplus)   * dz  * self.c_light
            if noCharge : fac = 0.
            fact = fac * trState.QOP()
            #create new state vector
            TrackStateNew    = TrackState( x  =  trState.X() + dz*( tx + 0.5 * ax * fact) , \
                                            y  =  trState.Y() + dz*( ty + 0.5 * ay * fact) , \
                                            tx =  trState.TX() + ax*fact , \
                                            ty =  trState.TY() + ay*fact , \
                                            qop=  trState.QOP() , \
                                            z  =  zNew )
            #update TranspMatrix                
            transMatrix  = self.GetTransportMatrix( dz, TrackStateNew, BVec, ax, ay )        
        return TrackStateNew, transMatrix
    def CreateStatesAlongTrajectory( self, stateOrigin: TrackState ,  zFinal : float = 18000., step_size: float = 100 , noCharge : bool = False, intermediate_z_locs = [] ):
        all_states = [ stateOrigin]
        z = all_states[-1].Z()
        
        z_locs = sorted([z + 100, 900, 1200, 1500, 1700, 1900, 2100, 2300, 2450, 2600, 2750, 3000, 3200, 3500, 3800, 4100, 4400, 4700, 5000, 5300, 5600, 5900, 6200, 6500, 6800, 7100, 7400, 7700, 8000, 8300, 8600, 8900, 9200, 9500, 12000, 15000, 18000] )
        curling = False        
        # z_locs  = [ _ for _ in np.linspace( start = z, stop = zFinal + step_size, num = int( (zFinal - z) / step_size) + 2 ) if _ > z ]
        z_locs.extend( intermediate_z_locs )
        if z_locs[0] < z : 
            raise ValueError(f"First z location to propagate is upstream of the state origin, this is not allowed. State origin z = {z} , first z location = {z_locs[0]}")
        for zNew in z_locs :
            state , trasp = self.Propagate( 
                trState = all_states[-1],
                zOld    = all_states[-1].Z(),
                zNew    = zNew,
            )
            if np.rad2deg(abs(state.TX()) ) > 80.:
                Logger.info(f"""
                            Large angle TX reached, stopping trajectory creation. TX = {state.TX()} TY = {state.TY()}
                            """)
                all_states.append(state) 
                curling = True
                break 
            if np.rad2deg(abs(state.TY()) ) > 80.:
                Logger.info(f"""
                            Large angle TY reached, stopping trajectory creation. TX = {state.TX()} TY = {state.TY()}
                            """)
                all_states.append(state) 
                curling = True
                break            
            all_states.append( state)            
        return all_states,curling 
            
        
        
    def CreateStatesAlongTrajectory_ZVALUES( self, stateOrigin: TrackState ,  Z_VALUES):
        all_states = [ stateOrigin]
        z = all_states[-1].Z()
        curling = False
        z_locs  = [ _ for _ in Z_VALUES if _ > z ]
        for zNew in z_locs :
            state , trasp = self.Propagate( 
                trState = all_states[-1],
                zOld    = all_states[-1].Z(),
                zNew    = zNew,
            )
            if np.rad2deg(abs(state.TX()) ) > 80.:
                Logger.info(f"""
                            Large angle TX reached, stopping trajectory creation. TX = {state.TX()} TY = {state.TY()}
                            """)
                all_states.append(state) 
                curling = True
                break 
            if np.rad2deg(abs(state.TY()) ) > 80.:
                Logger.info(f"""
                            Large angle TY reached, stopping trajectory creation. TX = {state.TX()} TY = {state.TY()}
                            """)
                all_states.append(state) 
                curling = True
                break
        return all_states,curling 
                    
        
    
