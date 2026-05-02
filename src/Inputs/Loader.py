
from   Utils.Logger import Logger
import pandas as pd
import ROOT as r
import numpy as np
def LoadInfoTracks(kinematic_file : str , start_entry : int , stop_entry : int ) -> dict :
    Logger.info(f"Loading entries from {kinematic_file}")
    inputFileLoad = kinematic_file
    startEntry    = start_entry
    endEntry      = stop_entry    
    if ".root" in inputFileLoad: 
        dfROOT = r.RDataFrame( "tree", f"{inputFileLoad}")
        
        df = pd.DataFrame(dfROOT.AsNumpy())
        df.to_parquet(f'{inputFileLoad.replace(".root", ".parquet")}')
        del df 
    elif ".parquet" in inputFileLoad:
        inputFileLoad = inputFileLoad.replace(".root", ".parquet")
    df = pd.read_parquet(inputFileLoad)
    inputs = { 
        "pvPos"    : [],
        "originPos": [],
        "PxPyPzPE" : [],
        "charge"   : [] 
    }
    entry = 0 
    for orig_x, orig_y, orig_z , px, py, pz, pe , charge, pvx,pvy,pvz  in zip( df["K_TRUEORIGINVERTEX_X"], 
                                                                               df["K_TRUEORIGINVERTEX_Y"], 
                                                                               df["K_TRUEORIGINVERTEX_Z"], 
                                                                               df["K_TRUEP_X"], 
                                                                               df["K_TRUEP_Y"], 
                                                                               df["K_TRUEP_Z"], 
                                                                               df["K_TRUEP_E"], 
                                                                               df["K_Charge"], 
                                                                               df["PV_X"], 
                                                                               df["PV_Y"],
                                                                               df["PV_Z"]):
        
        if entry < startEntry : 
            entry += 1
            continue
        if entry > endEntry : break
        entry += 1
        #shift observed on PV pos , to apply on both pvPos and originPos to be at (0,0,0) as mean value
        inputs["pvPos"].append( np.array(     [pvx   ,pvy   ,pvz   ]))
        inputs["originPos"].append( np.array( [orig_x,orig_y,orig_z]))
        inputs["PxPyPzPE"].append(  np.array( [px    ,py    ,pz, pe]))
        inputs["charge"].append( charge)    
    return inputs