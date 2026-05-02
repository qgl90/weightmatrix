from __future__ import annotations
from unicodedata import name
from tabulate import tabulate

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
from ROOT import TVector3
r.DisableImplicitMT() 
from Utils.Logger import Logger
from Settings.Settings import load_settings
from Core.ParabolicExtrapolator import ParabolicExtrapolator    
from Core.MaterialLocator import MaterialLocator 
from Core.BField import BField
from Core.TrackState import TrackState
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep 
hep.style.use("LHCb2")


if __name__ == "__main__":     
    settings_confs  = load_settings("layouts/u2_baseline.yaml")
    Extrap = ParabolicExtrapolator( FieldMap = settings_confs["fieldMap"])        
    Logger.info(f"""
                Field Map = {settings_confs['fieldMap']},
                ShrinkFact= {settings_confs['fieldShrinkFactor']},
                MagnifyFac= {settings_confs['fieldMagFactor']}
                """)
    
    TrajectoryCreator = ParabolicExtrapolator(
         BField( settings_confs["fieldMap"],
                settings_confs["fieldShrinkFactor"] ,
                settings_confs["fieldMagFactor"])    
    )
    
    MatLocator = MaterialLocator( gdmlFile=  settings_confs["gdmlInput"] )
    
    # =================== Setup color scales ====================
    from matplotlib.gridspec import GridSpec
    import matplotlib.cm as cm

    momenta    = [1000, 3000, 5000, 10000, 15000, 20000, 30000, 50000, 100000]  # MeV/c
    # momenta    = [20000]  # MeV/c    
    n          = len(momenta)+3
    colors_pos = [cm.Blues(0.35 + 0.65 * i / (n - 1)) for i in range(n)]
    colors_neg = [cm.Reds (0.35 + 0.65 * i / (n - 1)) for i in range(n)]

    # =================== 3D Plot ====================
    fig = plt.figure(figsize=(40, 16))
    gs  = GridSpec(3, 2, figure=fig, height_ratios=[6, 6, 1], hspace=0.35, wspace=0.3)
    ax     = fig.add_subplot(gs[0:2, 0], projection='3d')
    ax2    = fig.add_subplot(gs[0, 1])   # YZ plane
    ax3    = fig.add_subplot(gs[1, 1])   # XZ plane
    ax_leg = fig.add_subplot(gs[2, :])   # legend panel
    ax_leg.axis('off')
    legend_handles = []
    for i, p in enumerate(momenta):    
        Logger.info(f"Creating trajectory for p = {p} MeV/c")        
        MomentumOrigin = TVector3(0, 0, p)
        UnitDir = MomentumOrigin.Unit()
        tx = UnitDir.X() / UnitDir.Z()
        ty = UnitDir.Y() / UnitDir.Z()
        # 7 meters, must go in MP region
        OriginState_PosCharge = TrackState( x  =  0.0 , \
                                            y  =  0.0 , \
                                            tx =  0.040 ,
                                            ty =  0.050 , \
                                            qop=  1./MomentumOrigin.Mag(), \
                                            z  =  0. )        
        trajectory_pos, curling_pos = TrajectoryCreator.CreateStatesAlongTrajectory( stateOrigin=OriginState_PosCharge, zFinal=10000. )
        Material_Crossed = []
        for idx in range(len(trajectory_pos)-1):
            # current to next state material found 
            start = trajectory_pos[idx]
            stop  = trajectory_pos[idx+1]
            Logger.info(f"""Finding material crossed between start = ({start.X():.1f}, {start.Y():.1f}, {start.Z():.1f}) mm and stop = ({stop.X():.1f}, {stop.Y():.1f}, {stop.Z():.1f}) mm""")
            layers = MatLocator.FindMaterialFromPositionToPosition( stateVecIn =start, stateVecOut =stop )
            Material_Crossed.extend(layers)
        total = 0.0
        table_data = []
        for mat in Material_Crossed:
            entry_z = mat[0]
            exit_z  = mat[1]
            dz      = mat[3]
            x0      = mat[2]
            node_id = mat[4]            
            measurement = False 
            # if "lvMP" in node_id or "lvT:" in node_id or "lvFT" in node_id and exit_z > 6000.:
            #     x0 = 305288.4 #air x0 ! 
            #     node_id = "AIR-REPLACEMENT-MightyTracker"
            if ("MP:SiliconTracker" in node_id) : 
                # MP interceptor for hits
                Logger.critical(f"Found a silicon sensor crossed in MP → Measurement to create")
                node_id = f"MP-MEASUREMENT_{node_id}"
                measurement = True
            if ("lvUP" in node_id or "up" in node_id) and ":Silicon" in node_id:
                Logger.critical(f"Found a silicon sensor crossed in lvUP → Measurement to create")
                node_id = f"UP-MEASUREMENT_{node_id}"                
                measurement = True                
            if ("lvDet" in node_id and "TVSensorSilicon" in node_id) or ("TVSubstrateSilicon" in node_id) or ("TV:SubstrateSilicon" in node_id) or ("TV:SensorSilicon" in node_id) : 
                Logger.critical(f"Found a silicon sensor crossed in TV → Measurement to create")
                node_id = f"TV-MEASUREMENT_{node_id}"
                measurement = True
            if "lvFTFibre" in node_id and "FT:SciFibre" in node_id:
                Logger.critical(f"Found a scintillating fibre crossed in FT → Measurement   to create")
                node_id = f"FT-MEASUREMENT_{node_id}"                                
                measurement = True
            x_over_x0 = dz / x0
            total   += x_over_x0                        
            def red(text, measurement=measurement):
                return f"\033[91m{text}\033[0m" if measurement else text
            table_data.append(
                [ red(_, measurement) for _ in 
                    [
                    f"{entry_z:.1f}",
                    f"{exit_z:.1f}",
                    f"{dz:.2f}",
                    f"{x_over_x0*100:.2f}",
                    f"{total*100:.2f}",
                    f"{x0/10.:.2f}",
                    node_id
                    ]       
            ])
        # Print beautiful table
        # headers = ["Entry Z (mm)", "Exit Z (mm)", "dz (mm)", "x/X0 (%)", "Total x/X0 (%)", "X0 (cm)", "Node ID"]
        # Logger.info("\n" + "="*140)
        # Logger.info("Material Crossed Summary")
        # Logger.info(tabulate(table_data, headers=headers, tablefmt="pretty", floatfmt=".2f"))
        # Logger.info("="*140)
        # Logger.info(f"Total radiation length traversed: {total*100:.2f}% of X0")
        # Logger.info("="*140)
        OriginState_NegCharge = TrackState( x  =  0.0 , \
                                            y  =  0.0 , \
                                            tx =  -0.040 ,
                                            ty =  -0.050 , \
                                            qop=  -1./MomentumOrigin.Mag(), \
                                            z  =  0. )
        trajectory_neg, curling_neg = TrajectoryCreator.CreateStatesAlongTrajectory( stateOrigin=OriginState_NegCharge, zFinal=10000. )

        xyz_pos = np.array([(st.X(), st.Y(), st.Z()) for st in trajectory_pos])
        x_pos, y_pos, z_pos = xyz_pos[:, 0]/1000., xyz_pos[:, 1]/1000., xyz_pos[:, 2]/1000.

        xyz_neg = np.array([(st.X(), st.Y(), st.Z()) for st in trajectory_neg])
        x_neg, y_neg, z_neg = xyz_neg[:, 0]/1000., xyz_neg[:, 1]/1000., xyz_neg[:, 2]/1000.

        c_pos   = colors_pos[i]
        c_neg   = colors_neg[i]
        p_label = f'p = {p/1000.:.1f} GeV/c'

        # 3D
        ax.plot(z_pos, x_pos, y_pos, '-',  color=c_pos, linewidth=2)
        ax.plot(z_neg, x_neg, y_neg, '--', color=c_neg, linewidth=2)

        # YZ
        ax2.plot(z_pos, y_pos, '-',  color=c_pos, linewidth=2)
        ax2.plot(z_neg, y_neg, '--', color=c_neg, linewidth=2)

        # XZ – collect one pair of handles per momentum for the shared legend
        h_pos, = ax3.plot(z_pos, x_pos, '-',  color=c_pos, linewidth=2, label=f'q=+1  {p_label}')
        h_neg, = ax3.plot(z_neg, x_neg, '--', color=c_neg, linewidth=2, label=f'q=\u22121  {p_label}')
        legend_handles.extend([h_pos, h_neg])

    # ==================== Axis formatting (outside loop) ====================
    ax.set_xlabel('Z [m]', labelpad=20.0)
    ax.set_ylabel('X [m]', labelpad=20.0)
    ax.set_zlabel('Y [m]', labelpad=20.0)
    ax.set_xlim(-0.2, 18)
    ax.set_ylim(-3, 3)
    ax.set_zlim(-3, 3)
    ax.view_init(elev=20, azim=-120)
    ax.set_title('3D Trajectory')
    ax.grid(True)

    ax2.set_xlabel('Z [m]')
    ax2.set_ylabel('Y [m]')
    ax2.set_xlim(-0.2, 18)
    ax2.set_ylim(-3, 3)
    ax2.set_title('YZ Projection')
    ax2.grid(True)

    ax3.set_xlabel('Z [m]')
    ax3.set_ylabel('X [m]')
    ax3.set_xlim(-0.2, 18)
    ax3.set_ylim(-3, 3)
    ax3.set_title('XZ Projection')
    ax3.grid(True)

    # ==================== Shared legend pad ====================
    ax_leg.legend(handles=legend_handles, loc='center', ncol=n,
                  frameon=False, fontsize=13, handlelength=2.5)

    plt.savefig("plots/trajectory_checks/test.pdf", bbox_inches='tight')