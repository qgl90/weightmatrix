from __future__ import annotations
from unicodedata import name

def preliminary() : 
    import sys
    from pathlib import Path
    # def add_src_to_path() -> None:
    """
    Allow running scripts directly (e.g. `python3 scripts/run_it.py`) without
    requiring an editable install by prepending `<repo>/src` to `sys.path`.
    """
    repo_root = Path(__file__).resolve().parents[2]    
    src_dir = repo_root / "src"
    if src_dir.is_dir():
        print(f"Append to sys.path {src_dir}")
        sys.path.insert(0, str(src_dir))
preliminary()


from alive_progress import alive_bar
import ROOT as r
from ROOT import TLorentzVector
r.DisableImplicitMT() 
from Utils.Logger import Logger
from Settings.Settings import load_settings
from Core.MomentumResolutionCalculator import MomentumResolutionCalculator
from weightmatrix.config import Settings

from argparse import ArgumentParser
if __name__ == "__main__":
    Logger.info("Running script to compute momentum resolution on sampled tracks")    
    Logger.warning("This code is an adaption of Pierre's Billoir work")
    arg_ = ArgumentParser()
    arg_.add_argument("--kinematic_file", default="samples/kaons_B2Kmumu.parquet", help="Path to the kinematic file")
    arg_.add_argument("--settings_file", default="layouts/u2_baseline.yaml", help="Path to the settings file")
    arg_.add_argument("--output_file", default="test.root", help="Path to the output file")    
    arg_.add_argument("--start_entry",type = int, default = 0, help="Number of events to run")        
    arg_.add_argument("--stop_entry",type = int, default = 1000, help="Number of events to run")        
    args = arg_.parse_args()
    
    Logger.info("Loading settings and kinematic file")
    Logger.info(f'kinematic_file           = {args.kinematic_file}')
    Logger.info(f'settings_file            = {args.settings_file}') 
    Logger.info(f'output_file              = {args.output_file}')
    Logger.info(f'start_entry              = {args.start_entry}')
    Logger.info(f'stop_entry               = {args.stop_entry}')
    
    
    kinematic_file  = args.kinematic_file
    # Use refactored settings resolver (YAML-relative paths, env expansion)
    settings_confs = Settings.from_yaml(args.settings_file).as_dict()
    Logger.info(f'pmin                     = {settings_confs["pmin"]}')
    Logger.info(f'pmax                     = {settings_confs["pmax"]}')
    Logger.info(f'steps                    = {settings_confs["steps"]}')
    Logger.info(f'path_geom_pickle         = {settings_confs["path_geom_pickle"]}')
    Logger.info(f'detectionlayers          = {settings_confs["detectionlayers"]}')
    Logger.info(f'fieldMap                 = {settings_confs["fieldMap"]}')
    Logger.info(f'storing_planes           = {settings_confs["storing_planes"]}')
    Logger.info(f'useGDML                  = {settings_confs["useGDML"]}')
    Logger.info(f'skipVelo                 = {settings_confs["skipVelo"]}')
    Logger.info(f'skipRICH1                = {settings_confs["skipRICH1"]}')
    Logger.info(f'skipUT                   = {settings_confs["skipUT"]}')
    Logger.info(f'skipFT                   = {settings_confs["skipFT"]}')
    Logger.info(f'skipWorld                = {settings_confs["skipWorld"]}')
    Logger.info(f'fieldShrinkFactor        = {settings_confs["fieldShrinkFactor"]}')
    Logger.info(f'fieldMagFactor           = {settings_confs["fieldMagFactor"]}')
    Logger.info(f'correctSciFiAngles       = {settings_confs["correctSciFiAngles"]}')
    Logger.info(f'correctUTInnerResolution = {settings_confs["correctUTInnerResolution"]}')
    Logger.info(f'VeloMeasurementsFromGDML = {settings_confs["VeloMeasurementsFromGDML"]}')
    Logger.info(f'UPMeasurementsFromGDML   = {settings_confs["UPMeasurementsFromGDML"]}')
    Logger.info(f'DownMeasurementsFromGDML = {settings_confs["DownMeasurementsFromGDML"]}')
    Logger.info(f'useHelium = {settings_confs["useHelium"]}')
    start_entry = args.start_entry
    stop_entry  = args.stop_entry
    CalculatorResolution = MomentumResolutionCalculator(
        settings_confs,
        outFileName= args.output_file.replace(".root", f"_{start_entry}_{stop_entry}.root"), 
        flyupto = None
    )    
    import pandas as pd  
    if not kinematic_file.endswith(".parquet"):
        raise ValueError(
            f"--kinematic_file must be a parquet file compatible with this runner, got: {kinematic_file}"
        )
    df = pd.read_parquet(kinematic_file)
    df["M1_CHARGE"] = df["M1_ID"].apply(lambda x : -1 if x > 0 else +1 )
    df["M2_CHARGE"] = df["M2_ID"].apply(lambda x : -1 if x > 0 else +1 )    
    
    tot_loops = stop_entry - start_entry
    import numpy as np 
    entries = np.arange( start_entry, stop_entry )
    with alive_bar( tot_loops) as bar :                
        for i in entries :
            entry = df.iloc[i]
            for fstate in ["M1", "M2"] :
                mu_pv_x   = float(entry[f'{fstate}_TRUEPRIMARYVERTEX_X'])
                mu_pv_y   = float(entry[f'{fstate}_TRUEPRIMARYVERTEX_Y'])
                mu_pv_z   = float(entry[f'{fstate}_TRUEPRIMARYVERTEX_Z'])            
                mu_ovtx_x = float(entry[f'{fstate}_TRUEORIGINVERTEX_X'])
                mu_ovtx_y = float(entry[f'{fstate}_TRUEORIGINVERTEX_Y'])
                mu_ovtx_z = float(entry[f'{fstate}_TRUEORIGINVERTEX_Z'])
                mu_px = float(entry[f'{fstate}_PX'])
                mu_py = float(entry[f'{fstate}_PY'])
                mu_pz = float(entry[f'{fstate}_PZ'])
                mu_e  = float(entry[f'{fstate}_ENERGY'])
                mu_charge = int(entry[f'{fstate}_CHARGE'])                
                
                pv = (mu_pv_x, mu_pv_y, mu_pv_z)
                origin_mu1 = (mu_ovtx_x, mu_ovtx_y, mu_ovtx_z)            
                mu_vect = TLorentzVector()
                mu_vect.SetPxPyPzE(mu_px, mu_py, mu_pz, mu_e)            
                originX = origin_mu1[0]
                originY = origin_mu1[1]
                originZ = origin_mu1[2]
                pvPos_X  = pv[0]
                pvPos_Y  = pv[1]
                pvPos_Z  = pv[2]
                # if mu_vect.PseudoRapidity() >3 and mu_vect.PseudoRapidity() < 5.0 :
                CalculatorResolution.ComputeResolutions(
                    fourvector= mu_vect,
                    charge = mu_charge,
                    originPosition=(originX,originY,originZ),
                    PVPosition    =(pvPos_X, pvPos_Y, pvPos_Z),
                )
            bar()
        CalculatorResolution.SaveTupleToTFile()
    Logger.warning("It's over")
