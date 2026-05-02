from alive_progress import alive_bar

import ROOT as r 
from ROOT import TGeoManager
r.gSystem.Load("libGeom")
from Core.TrackState import TrackState
import numpy as np 
from Utils.Logger import Logger
import math 
from termcolor import colored 
import json 
from Core.Plane import Plane 
import fnmatch 
import os 
class MaterialAndMeasurementLocator: 
    def __init__(self, gdmlFile = "geometries/LHCb_Upgrade_Simple.gdml", settings_measurements : dict = None   ) -> None:
        Logger.debug("MaterialAndMeasurementLocator Initialization....")
        if r.gGeoManager  : 
            print(r.gGeoManager)
            self.geoManager = r.gGeoManager
        else: 
            TGeoManager.LockDefaultUnits( False )
            if TGeoManager.GetDefaultUnits() != 0 : 
                TGeoManager.SetDefaultUnits(TGeoManager.kG4Units)
            Logger.debug(f"TGeoManager Units = {TGeoManager.GetDefaultUnits()}")
            if not( os.path.exists( gdmlFile)): 
                Logger.error( f"{gdmlFile} do not exists")
                raise ValueError("Cannot open gdmlFile, not existing")
            self.geoManager = TGeoManager.Import(gdmlFile)
        # Logger.debug( settings_measurements.keys())
        self.skipVelo  = True 
        self.skipUT    = True
        self.skipSciFi = True
        self.skipRICH  = True
        self.skipWorld = True
        self.VeloMeasurementFromGDML = False        
        self.ActiveElementsCreator = { 
            "TVSensorSilicon:lvDet_shape_": { 
                "sigmaX" : 0.0158 , 
                "sigmaY" : 0.0158 ,
                "type"   : "VeloGDML",
                "thickScaleFactor" : 1 ,
                "eff" : 0.99
            },
            # TODO : add wrapper to load stuff! 
            # "Silicon:lvDet_shape_"                       : { 
            #     "sigmaX" : 0.0158 , 
            #     "sigmaY" : 0.0158 ,
            #     "type"   : "VeloGDML",
            #     "thickScaleFactor" : 1 ,
            #     "eff" : 0.99
            # },
            # "FTSciFibre:lvFibreMat*_shape_"               : { 
            #     "sigmaX" : 0.130  ,  
            #     "alphaFromMother" : { "lvFTModuleHoleRightU" :-5 , 
            #                         "lvFTModuleHoleLeftU"  : 5 ,
            #                         "lvFTModuleHoleLeftX1" : 0 , 
            #                         "lvFTModuleHoleRightX1": 0},
            #     "type"   : "FTGDML",
            #     "thickScaleFactor" : 1,
            #     "eff" : 0.98
            # },
            # "MPSiliconTracker:ChipBox*_shape_"                  : { 
            #     "sigmaX" : 0.0144 ,
            #     "sigmaY" : 0.0433 ,
            #     "type"   : "MPGDML",
            #     "thickScaleFactor" : 1,
            #     "eff" : 0.99
            # },        
            # "Silicon:upSensor*_shape_" : { 
            #     "sigmaX" : 0.0144, 
            #     "sigmaY" : 0.0433, 
            #     "type"   : "UPGDML", 
            #     "thickScaleFactor" : 1,
            #     "eff" : 0.99
            # }
        }    

    def IDNODE( self,node) :
        volName  = node.GetVolume().GetName()
        shapeName=node.GetVolume().GetShape().GetName()        
        matName  = node.GetVolume().GetMedium().GetMaterial().GetName() 
        return f"{matName}:{shapeName}"
        # return f"{volName}:{shapeName}"

    def IsMeasurementNode( self, node):
        volName  = node.GetVolume().GetName()
        shapeName=node.GetVolume().GetShape().GetName()
        ID = f"{volName}:{shapeName}"
        Logger.debug(f"Element ID in geometry = {ID}")
        if ID in self.ActiveElementsCreator.keys() : return True 
        return False
    def makeArr(self,  cPoint): 
        return np.array( [cPoint[0], cPoint[1], cPoint[2]])
    def MakePlanesListFromTo( self, startp : np.array, endp : np.array, track_qop : float ):
        if endp[2] >  startp[2]:
            forward = True
        else:
            forward = False
        dirvec = (endp - startp)
        dirvec /= math.sqrt( dirvec[0]**2 + dirvec[1]**2  + dirvec[2]**2 )
        try : 
            self.geoManager.InitTrack( startp[0], startp[1], startp[2] , dirvec[0], dirvec[1],dirvec[2])
        except :
            Logger.error(f"Cannot initTrack startp = {startp} , endp = {endp}")
            raise ValueError("Cannot initTrack Geomanager")        
        entry_exit_points  = [     ]
        radLen_midPoints   = [     ]        
        volIDStep          = [     ]        
        MeasurementsPlanes = [     ]
        PassivePlanes      = [     ]
        node = self.geoManager.GetCurrentNode()
        while(node):
            entryP        = self.makeArr( self.geoManager.GetCurrentPoint() )        
            # LInt          = node.GetVolume().GetMedium().GetMaterial().GetIntLen()
            # MaterialName. = node.GetVolume().GetMedium().GetMaterial().GetName()
            IDNODE = self.IDNODE(node)   
            Logger.debug(f"Stepping into node with ID = {IDNODE} at entryP = {entryP} with track dirvec = {dirvec}")
            """
            ActiveElementsCreator is a dict with keys being the ID of the node to be considered as measurement, and values being a dict with the following keys
                - sigmaX : the sigma of the measurement in x direction in mm
                - sigmaY : the sigma of the measurement in y direction in mm
                - type   : the type of the measurement, used for categorization and plotting (e.g. VeloGDML, UPGDML, MPGDML, FTGDML)
                - thickScaleFactor : the factor to apply to the thickness of the material to be considered as measurement, used to tune the contribution of the measurement to the scattering (e.g. for a silicon sensor, we can consider only 10% of the thickness for scattering)
                - eff : the hit efficiency of the measurement, used to tune the contribution of the measurement to the scattering (e.g. for a silicon sensor, we can consider 99% of the efficiency)            
            """
            #########################################
            #baseline for any passive elements.
            ########################################
            X0Use         = node.GetVolume().GetMedium().GetMaterial().GetRadLen()        
            sigmaX = None 
            sigmaY = None
            angle  = None
            hitefficiency = 1.0
            # name          = IDNODE
            geomLayer = "any"
            ############################################################
            # Figure out if it's a measurement layer or not 
            ############################################################
            def check_regexp( node_id , keys_elements):
                matched = False 
                keymatch = None 
                for k in keys_elements: 
                    if fnmatch.fnmatch(  node_id, k):
                        matched = True 
                        keymatch = k 
                        break     
                if matched : 
                    Logger.debug( f"CHECKING : node_id = {node_id} for measurement among {keys_elements} Matched!")                
                return matched , keymatch
            if IDNODE in self.ActiveElementsCreator.keys() :
                Logger.debug(f"Exact match for measurement node found with ID = {IDNODE} Must Be true")
                IsMeasurement , KMATCH= check_regexp( IDNODE , self.ActiveElementsCreator.keys())
            else: 
                IsMeasurement = False
                                
            if IsMeasurement :
                Logger.debug(f"Measurement node found with ID = {IDNODE}")
                name          = IDNODE                
                X0Use         = X0Use/self.ActiveElementsCreator[ KMATCH ]["thickScaleFactor"]                
                sigmaX        = self.ActiveElementsCreator[    KMATCH ]["sigmaX"]
                sigmaY        = 0.
                angle         = 45.0
                hitefficiency = self.ActiveElementsCreator[    KMATCH ]["eff"]
                geomLayer     = self.ActiveElementsCreator[   KMATCH ]["type"]
                if   geomLayer in ["VeloGDML", "UPGDML", "MPGDML"]: 
                    sigmaY , angle = self.ActiveElementsCreator[KMATCH]["sigmaY"] , 0.
                elif geomLayer in ["UTGDMLInner" , "UTGDMLOuter"] : 
                    ###strip stuff
                    sigmaY = None
                    angle  = None
                    for m in [ f"{self.geoManager.GetMother(i).GetName()}" for i in range(1,9,1)]:
                        if    ( "UT" in m and "X" in m) : angle = 0.
                        elif  ( "UT" in m and "V" in m) : angle =-5.
                        elif  ( "UT" in m and "U" in m) : angle = 5.
                        if angle != None : break
                    if angle == None : raise ValueError("Invalid UT angle discovery")
                elif geomLayer == "FTGDML":  
                    ###fibre stuff
                    sigmaY = None 
                    angle  = None
                    for m in [ f"{self.geoManager.GetMother(i).GetName()}" for i in range(1,9,1)]: 
                        if    ( "lvLayer" in m  and ("X1_" in m or "X2_" in m)) : angle = 0.
                        elif  ( "lvLayer" in m  and ("U_" in m )) : angle =  5.
                        elif  ( "lvLayer" in m  and ("V_" in m )) : angle = -5.
                        if angle != None : break                        
                        # Logger.info(f"FT discovering U/V --> {m}")                        
                    # Logger.info(f"Angle assigned = {angle}")
                    if angle == None : raise ValueError("impossible, angle has to be found")
            else: 
                Logger.debug(f"Passive node found with ID = {IDNODE}")
                
            node = self.geoManager.FindNextBoundaryAndStep()
            breakIt = False 
            if node : 
                exitP = self.makeArr( self.geoManager.GetCurrentPoint() )              
                if exitP[2] > endp[2]  and forward : 
                    breakIt = True
                    exitP = endp #set exitPoint to endPoint and break afterwards
                elif exitP[2] < endp[2] and not forward     : 
                    breakIt = True
                    exitP = endp #set exitPoint to endPoint and break afterwards
                else: 
                    breakIt = False 
            else : 
                exitP = endp
                breakIt = True
            ##### we want to store the Plane to be at z = entryZ 
            ##### Exit is at                          z = z+dzMat 
            ##### z Used for scattering               z = z+dzMat/2. [?]
            if IsMeasurement : 
                MeasurementsPlanes.append( 
                    Plane( zval = entryP[2],
                           dzMat= (exitP-entryP)[2],                     
                           dzMatAssignZero=True,                            
                           thickness= (exitP - entryP)[2]/X0Use,
                           sigmaX = sigmaX ,
                           sigmaY = sigmaY ,
                           angle  = float(angle),
                           geomLayer= geomLayer , 
                           hitefficiency =hitefficiency,
                           name = IDNODE, 
                    )
                )
                # State where ? entry or exit point? In the middle?
                # MeasurementsPlanes[-1].associateState( TrackState(x =  (exitP+entryP)[0]/2. ,
                #                                                   y =  (exitP+entryP)[1]/2. ,
                #                                                   z=   (exitP+entryP)[2]/2. ,
                MeasurementsPlanes[-1].associateState( TrackState(x =   entryP[0] ,
                                                                  y =   entryP[1] ,
                                                                  z=    entryP[2] ,                                                                  
                                                                  tx=   dirvec[0]/dirvec[2] ,
                                                                  ty=   dirvec[1]/dirvec[2] ,
                                                                  qop=  track_qop))
            else: 
                # print( geomLayer)
                PassivePlanes.append(
                    Plane( zval = entryP[2],
                           dzMat=  (exitP - entryP)[2],
                           dzMatAssignZero= True,     
                           thickness= (exitP - entryP)[2]/X0Use, #not slope corrected
                           #thickness     =np.sqrt( np.sum( exitP - entryP))/X0Use,
                           sigmaX        =sigmaX,
                           sigmaY        =sigmaY,
                           angle         =float(angle) if angle else angle,
                           geomLayer     =geomLayer,
                           hitefficiency =hitefficiency,
                           name = IDNODE, 
                    )
                )
                PassivePlanes[-1].associateState(      TrackState(x =  entryP[0],
                                                                  y =  entryP[1],
                                                                  z =  entryP[2],
                                                                  tx=  dirvec[0]/dirvec[2],
                                                                  ty=  dirvec[1]/dirvec[2],
                                                                  qop= track_qop ))
            if breakIt : break                    
        return MeasurementsPlanes, PassivePlanes
    def GetMeasurementsAndPassiveFromStateToState( self,  stateVecIn :TrackState,  stateVecOut :TrackState , KeepOnly = None  ):
        endp   = np.array( [ stateVecOut.X(), stateVecOut.Y(), stateVecOut.Z()] )
        startp = np.array( [ stateVecIn.X(),  stateVecIn.Y() , stateVecIn.Z()]  )
        meas , passive = self.MakePlanesListFromTo( startp=startp, endp=endp, track_qop=stateVecIn.QOP() )
        return meas, passive
    
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_nTV_vs_eta(parquet_file: str = "test_planes_list.parquet", 
                    save_fig: bool = False, 
                    fig_name: str = "nTV_vs_eta.png"):
    """
    Read parquet file and create publication-ready plots of nTV vs eta.
    
    Parameters:
        parquet_file: path to the saved parquet file
        save_fig: whether to save the figure to disk
        fig_name: filename when saving
    """
    
    # Read data
    df = pd.read_parquet(parquet_file)
    print(f"Loaded {len(df):,} tracks from {parquet_file}")
    print(f"Unique nTV values: {sorted(df['nTV'].unique())}")
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Number of Measurements (nTV) vs Pseudorapidity (η)', fontsize=16)
    
    # 1. Scatter plot
    axes[0, 0].scatter(df['eta'], df['nTV'], alpha=0.25, s=6, color='steelblue')
    axes[0, 0].set_xlabel('η (Pseudorapidity)')
    axes[0, 0].set_ylabel('nTV')
    axes[0, 0].set_title('Scatter Plot')
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. 2D Histogram (Hexbin)
    hb = axes[0, 1].hexbin(df['eta'], df['nTV'], gridsize=60, cmap='viridis', mincnt=1)
    plt.colorbar(hb, ax=axes[0, 1], label='Number of tracks')
    axes[0, 1].set_xlabel('η (Pseudorapidity)')
    axes[0, 1].set_ylabel('nTV')
    axes[0, 1].set_title('Density (Hexbin)')
    
    # 3. Mean nTV vs eta (Profile)
    eta_bins = np.linspace(-6, 6, 80)
    df['eta_bin'] = pd.cut(df['eta'], bins=eta_bins)
    profile = df.groupby('eta_bin')['nTV'].agg(['mean', 'std', 'count'])
    
    mid_eta = profile.index.map(lambda x: x.mid)
    axes[1, 0].errorbar(mid_eta, profile['mean'], yerr=profile['std'],
                        fmt='o-', capsize=3, linewidth=2, color='darkorange')
    axes[1, 0].axhline(y=3, color='red', linestyle='--', alpha=0.5)
    axes[1, 0].set_xlabel('η (Pseudorapidity)')
    axes[1, 0].set_ylabel('Mean nTV ± Std')
    axes[1, 0].set_title('Mean nTV vs η')
    axes[1, 0].grid(True, alpha=0.3)
    
    # 4. Boxplot by eta regions
    bins = [-6, -4, -2, -1, 0, 1, 2, 4, 6]
    labels = ['-6..-4', '-4..-2', '-2..-1', '-1..0', '0..1', '1..2', '2..4', '4..6']
    df['eta_region'] = pd.cut(df['eta'], bins=bins, labels=labels)
    
    sns.boxplot(x='eta_region', y='nTV', data=df, ax=axes[1, 1], palette='Set2')
    axes[1, 1].set_xlabel('η Region')
    axes[1, 1].set_ylabel('nTV')
    axes[1, 1].set_title('Distribution by η Region')
    axes[1, 1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    
    if save_fig:
        plt.savefig(fig_name, dpi=300, bbox_inches='tight')
        print(f"Figure saved as {fig_name}")
        
    return df  # return dataframe in case you want to use it further
    
if __name__ == "__main__":
    
    # if False: 
    df= {
        "nTV" : [] , 
        "eta" : [], 
        "phi" : []
    }
    Mat = MaterialAndMeasurementLocator( gdmlFile = "gdmls/LHCb_23Jan2026_run5_proto-v00.03.gdml")    
    NPHI = 50 
    NETA = 50    
    with alive_bar(NPHI * NETA, title="Processing bin for nTV vs ETA", bar="smooth", spinner="dots", enrich_print=False) as bar:
            
        for phi in np.linspace(-math.pi,math.pi,NPHI):
            for eta in np.linspace(-5,6,NETA):
                x , y ,z = 0, 0, 0 
                v = r.TVector3()
                v.SetPtEtaPhi(500.,eta , phi)
                P = float(10 * 1000. )
                v.SetMag(P)                                
                tx = v.X()/v.Z()
                ty = v.Y()/v.Z()                
                if eta > 0 : 
                    xend = x + tx*800.
                    yend = y + ty*800.
                    zend = 800.
                else: 
                    xend = x + tx*(-800.)
                    yend = y + ty*(-800.)
                    zend = -800.
                stateVecIn =  TrackState( x=x,y=y,z=z, tx=v.X()/v.Z(), ty=v.Y()/v.Z(), qop=1/P)
                stateVecOut = TrackState(x=xend,y=yend,z=zend, tx=v.X()/v.Z(), ty=v.Y()/v.Z(), qop=1/P)
                            
                Meas, Passive = Mat.GetMeasurementsAndPassiveFromStateToState(stateVecIn= stateVecIn, stateVecOut=    stateVecOut)
                nTV = len(Meas)                    
                df["nTV"].append( nTV )
                df["eta"].append( eta )
                df["phi"].append( phi )
                bar()   # update progress bar                
    df = pd.DataFrame(df)
    df.to_parquet("test_planes_list.parquet")       
    df = plot_nTV_vs_eta(
        parquet_file="test_planes_list.parquet",
        save_fig=True,
        fig_name="nTV_vs_eta_analysis.png"
    )
