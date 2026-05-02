"""
Material locator using ROOT TGeo navigation.

This module wraps `TGeoManager` stepping between two 3D points and returns a list
of crossed volumes (segments) with:

    [entry_z, exit_z, radLen, dz, node_id]

where:
- `radLen` is the material radiation length returned by ROOT for the current volume
- `dz` is the traveled z-distance in that volume segment (in mm)

The caller typically converts `dz/radLen` into an x/X0 contribution and inserts
corresponding passive planes into the tracking sequence.

Important
---------
`TGeoManager.Import(...)` supports both GDML and ROOT geometry files depending
on how the geometry was produced. In this project settings often use a `.root`
geometry file.
"""

import ROOT as r 
from ROOT import TGeoManager
r.gSystem.Load("libGeom")
from Core.TrackState import TrackState
from Utils.Logger import Logger
import numpy as np 
import math 
r.gInterpreter.Declare( 
    """
    const int MeanExcEnergy_NELEMENTS = 93; // 0 = vacuum, 1 = hydrogen, 92 = uranium
    const float MeanExcEnergy_vals[] = {1.e15, 19.2, 41.8, 40.0, 63.7, 76.0, 78., 82.0, 95.0, 115.0, 137.0, 149.0, 156.0, 166.0, 173.0, 173.0, 180.0, 174.0, 188.0, 190.0, 191.0, 216.0, 233.0, 245.0, 257.0, 272.0, 286.0, 297.0, 311.0, 322.0, 330.0, 334.0, 350.0, 347.0, 348.0, 343.0, 352.0, 363.0, 366.0, 379.0, 393.0, 417.0, 424.0, 428.0, 441.0, 449.0, 470.0, 470.0, 469.0, 488.0, 488.0, 487.0, 485.0, 491.0, 482.0, 488.0, 491.0, 501.0, 523.0, 535.0, 546.0, 560.0, 574.0, 580.0, 591.0, 614.0, 628.0, 650.0, 658.0, 674.0, 684.0, 694.0, 705.0, 718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0, 800.0, 810.0, 823.0, 823.0, 830.0, 825.0, 794.0, 827.0, 826.0, 841.0, 847.0, 878.0, 890.0};

    float MeanExcEnergy_get(int Z)
    {
    assert(Z >= 0 && Z < MeanExcEnergy_NELEMENTS);
    return MeanExcEnergy_vals[Z];
    }

    float MeanExcEnergy_get(TGeoMaterial* mat)
    {
    if (mat->IsMixture()) {
        double logMEE = 0.;
        double denom  = 0.;
        TGeoMixture* mix = (TGeoMixture*)mat;
        for (int i = 0; i < mix->GetNelements(); ++i) {
        int index = int(floor((mix->GetZmixt())[i]));
        logMEE += 1. / (mix->GetAmixt())[i] * (mix->GetWmixt())[i] * (mix->GetZmixt())[i] * log(MeanExcEnergy_get(index));
        denom  += (mix->GetWmixt())[i] * (mix->GetZmixt())[i] * 1. / (mix->GetAmixt())[i];
        }
        logMEE /= denom;
        return exp(logMEE);
    } else { // not a mixture
        int index = int(floor(mat->GetZ()));
        return MeanExcEnergy_get(index);
    }
    }    
    """
)
import os 
class MaterialLocator: 
    def __init__(self, gdmlFile = "geometries/LHCb_Upgrade_Simple.gdml") -> None:
        print("MaterialLocator Initialization....")
        if r.gGeoManager  : 
            print(r.gGeoManager)
            self.geoManager = r.gGeoManager
        else: 
            TGeoManager.LockDefaultUnits( False )
            if TGeoManager.GetDefaultUnits() != 0 : 
                TGeoManager.SetDefaultUnits(TGeoManager.kG4Units)
            Logger.info(f"TGeoManager Units = {TGeoManager.GetDefaultUnits()}")
            if not( os.path.exists( gdmlFile)): 
                Logger.error( f"{gdmlFile} do not exists")
                raise ValueError("Cannot open gdmFile, not existing")
            Logger.info(f"Importing GDML file {gdmlFile} ...")
            self.geoManager = TGeoManager.Import(gdmlFile)
            self.geoManager.SetVisLevel( 10000 )                 
    def makeArr(self,  cPoint): 
        return np.array( [cPoint[0], cPoint[1], cPoint[2]])
    def node_id(self, node) -> str :
        volName    = node.GetVolume().GetName()
        shapeName  = node.GetVolume().GetShape().GetName()
        # mediumName = node.GetVolume().GetMedium().GetName()
        matName    = node.GetVolume().GetMedium().GetMaterial().GetName() 
        return f"{volName}:{shapeName}:{matName}"
    # def FineMaterialFromPositionToPositionXYZ( self, startp :np.array , endp :np.array , KeepOnly = None  ): 
    #     if startp[2] == endp[2] : 
    #         Logger.fatal("StartP and EndP have same z location, no material returned")
    #         return []
    #     if startp[2] > endp[2] : 
    #         Logger.fatal("StartP Z is greater than EndP Z, no material returned")
    #         return []
        
    #     startP = r.TVector3(startp[0], startp[1], startp[2]) 
    #     endP   = r.TVector3(endp[0], endp[1], endp[2])
        
    #     dirvec = endP - startP
    #     dirvec.Unit()
    #     try : 
    #         self.geoManager.InitTrack( startP.X(), startP.Y(), startP.Z() , dirvec.X(), dirvec.Y(),dirvec.Z())
    #     except :
    #         Logger.error(f"""
    #                      Cannot initTrack startp = ({startP.X()}, {startP.Y()}, {startP.Z()}) , 
    #                                       endp =   ({endP.X()}, {endP.Y()}, {endP.Z()})""")
    #         raise ValueError("Cannot initTrack Geomanager")
    #     node = self.geoManager.GetCurrentNode()
    #     entry_point         = self.makeArr(self.geoManager.GetCurrentPoint())
    #     entry_x_over_x0     = node.GetVolume().GetMedium().GetMaterial().GetRadLen()
    #     entry_IDNODE        = self.node_id( node)        
    #     nextNode            = self.geoManager.FindNextBoundary((endP-startP).Mag())
    #     exit_point          = self.makeArr(self.geoManager.GetCurrentPoint())
    #     break_condition = False
    #     if exit_point[2] > endP.Z() : 
    #         break_condition = True
    #         exit_point = self.makeArr( [endP.X(), endP.Y(), endP.Z()])
    #     Layers_To_Add_In_Step = [] 
    #     Layers_To_Add_In_Step.append(  [entry_point[2] , exit_point[2], entry_x_over_x0, exit_point[2]-entry_point[2], entry_IDNODE])        
    #     if break_condition :
    #         return Layers_To_Add_In_Step        
    #     while nextNode is not None: 
            
    #         node = nextNode 
    #         entry_point         = self.makeArr( self.geoManager.GetCurrentPoint())
    #         entry_x_over_x0     =  node.GetVolume().GetMedium().GetMaterial().GetRadLen()
    #         entry_IDNODE        = self.node_id( node)
    #         nextNode            = self.geoManager.FindNextBoundary()
    #         exit_point          = self.makeArr( self.geoManager.GetCurrentPoint())            
    #         if exit_point[2] > endP.Z() : 
    #             exit_point = self.makeArr([endP.X(), endP.Y(), endP.Z()])
    #             break_condition = True
    #         # entry, exit, x/x0 in the layer, dz in the layer, nodeID
    #         Layers_To_Add_In_Step.append([entry_point[2] , exit_point[2], entry_x_over_x0, exit_point[2]-entry_point[2], entry_IDNODE])
    #         if break_condition :
    #             break
        
    #     return Layers_To_Add_In_Step

    def FineMaterialFromPositionToPositionXYZ(self, startp: np.array, endp: np.array, KeepOnly=None):
        if abs(startp[2] - endp[2]) < 1e-8:
            Logger.fatal(f"""
                         StartP and EndP have same z location, no material returned
                         {startp[2]} --> {endp[2]}
                         """)
            return []

        if startp[2] > endp[2]:
            Logger.fatal("StartP Z is greater than EndP Z, no material returned")
            return []

        # Direction and total length
        startP = r.TVector3(startp[0], startp[1], startp[2])
        endP   = r.TVector3(endp[0],   endp[1],   endp[2])
        dirvec = endP - startP
        total_length = dirvec.Mag()

        if total_length < 1e-8:
            return []

        dirvec = dirvec.Unit()

        # Initialize track
        try:
            self.geoManager.InitTrack(startP.X(), startP.Y(), startP.Z(),
                                    dirvec.X(), dirvec.Y(), dirvec.Z())
        except Exception as e:
            Logger.error(f"Cannot InitTrack: start=({startP.X():.3f}, {startP.Y():.3f}, {startP.Z():.3f}) "
                            f"end=({endP.X():.3f}, {endP.Y():.3f}, {endP.Z():.3f})")
            raise

        Layers_To_Add_In_Step = []
        traveled = 0.0
        target_z = endP.Z()

        while True:
            current_point = self.makeArr(self.geoManager.GetCurrentPoint())
            current_node  = self.geoManager.GetCurrentNode()
            
            if not current_node or not current_node.GetVolume():
                break

            vol   = current_node.GetVolume()
            med   = vol.GetMedium()
            mat   = med.GetMaterial() if med else None
            radLen = mat.GetRadLen() if mat else 0.0

            nodeID = self.node_id(current_node)

            # Find next boundary, but do not go beyond remaining distance
            remaining = total_length - traveled
            next_node = self.geoManager.FindNextBoundaryAndStep(remaining)
            exit_point = self.makeArr(self.geoManager.GetCurrentPoint())
            step       = self.geoManager.GetStep()

            traveled += step

            # Cap at target Z if we overshoot
            break_condition = False
            if exit_point[2] > target_z + 1e-6:
                exit_point = np.array([endP.X(), endP.Y(), endP.Z()])
                break_condition = True

            dz = exit_point[2] - current_point[2]

            # Add the layer
            Layers_To_Add_In_Step.append([
                current_point[2],      # entry z
                exit_point[2],         # exit z
                radLen,                # x/x0 = RadLen
                dz,                    # thickness in this volume
                nodeID                 # volume:shape:material
            ])

            if break_condition or self.geoManager.IsOutside() or step < 1e-8 or traveled >= total_length - 1e-6:
                break

            # Continue into next volume (daughters are automatically handled by the navigator)

        return Layers_To_Add_In_Step    
    def FindMaterialFromPositionToPosition( self,  stateVecIn :TrackState,  stateVecOut :TrackState , KeepOnly = None  ):
        endp   = np.array( [ stateVecOut.X(), stateVecOut.Y(), stateVecOut.Z()] )
        startp = np.array( [ stateVecIn.X(),  stateVecIn.Y() , stateVecIn.Z()]  )
        Logger.debug(f"FindMaterialFromPositionToPosition: {startp} → {endp}")
        return self.FineMaterialFromPositionToPositionXYZ( startp=startp, endp=endp, KeepOnly = KeepOnly)
    def TabulateMaterials( self, Entry_Exit_Infos: list, printIt=False): 
        """
        Extry_Exit_Infos = FineMaterialFromPositionToPositionXYZ output return
        """
        Logger.info("Tabulating.....")
        from prettytable import PrettyTable
        x = PrettyTable()
        HEADERS = [ "VolumeName" , "EntryZ [mm]", "ExitZ [mm]", "x/X0 [integrated]", "x/IntLength [integrated]", "DZ in material [mm]"]
        x.field_names = HEADERS
        x.align['VolumeName'] = 'l'
        total_thickness = 0
        total_intlen=0
        for i in range(len(Entry_Exit_Infos)):
            entryZ       = Entry_Exit_Infos[i][0]
            exitZ        = Entry_Exit_Infos[i][1]
            radLen       = Entry_Exit_Infos[i][2]
            dzInMat      = Entry_Exit_Infos[i][3]
            NodeID       = Entry_Exit_Infos[i][4]
            thickness    = 100*(exitZ - entryZ)/radLen
            print(Entry_Exit_Infos)
            intLen       = Entry_Exit_Infos[i][5]
            IntLength = 100 * (exitZ - entryZ)/intLen
            total_thickness += thickness
            total_intlen+=IntLength
            ROW = [f"{NodeID}" , 
                   f"{entryZ:.2f} mm", 
                   f"{exitZ:.2f}", 
                   f"{thickness:.2f} [{total_thickness:.2f}]%" , 
                   f"{IntLength:.2f} [{total_intlen:.2f}]%" ,f"{dzInMat:.2f}"]
            x.add_row( ROW)
        if printIt:
            print(x)
        return x
