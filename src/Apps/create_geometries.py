from __future__ import annotations
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

import os
from Utils.Logger import Logger
from argparse import ArgumentParser
from shapely import Polygon, Point
from shapely.geometry import Polygon, Point
from shapely import affinity
from shapely import geometry
from shapely.plotting import plot_polygon, plot_points
from shapely.ops import unary_union
import matplotlib.pyplot as plt
import pickle
import mplhep as hep
from math import sqrt
import numpy as np 
import pickle 
import math 
from adjustText import adjust_text

plt.style.use( hep.style.LHCb2)

def alabel(poly, name, decimals=2):
    """Quick area label in m²"""
    area_m2 = poly.area / 1_000_000
    return f"{name} (area = {area_m2:.{decimals}f} m²)"
def add_corners_text(thepolygon, ax, OFFSET=20, interiorsOnly=False):
    texts = []          # collect Text artists
    points = []         # optional: the (x,y) points for adjust_text
    
    if interiorsOnly:
        corners = []
        for interior in thepolygon.interiors:
            corners.extend(interior.coords)
    else:
        corners = list(thepolygon.exterior.coords)
    
    for corner in corners:
        x, y = corner[0], corner[1]
        
        # Only label in the quadrant you want, or remove this if you want all
        if not (x > 0 and y > 0):
            continue
            
        txt = ax.annotate(
            f'({x:.2f}, {y:.2f})',
            xy=(x, y),
            xytext=(x + OFFSET, y + OFFSET),   # initial guess
            arrowprops=dict(arrowstyle='->', color='blue', lw=0.8),
            ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7)  # optional nice box
        )
        texts.append(txt)
        points.append((x, y))
    
    # === AUTOMATIC ADJUSTMENT ===
    adjust_text(
        texts,
        x=[p[0] for p in points],
        y=[p[1] for p in points],
        ax=ax,
        arrowprops=dict(arrowstyle='->', color='blue', lw=0.8),
        expand_points=(1.5, 1.5),      # how much space to give around points
        expand_text=(1.2, 1.2),        # space around texts
        force_text=(0.1, 0.2),
        force_points=(0.1, 0.2),
        only_move={'points': 'xy', 'text': 'xy'},
        autoalign=True,
        # time_lim=2,                  # safety timeout if needed
    )
    
# def add_corners_text( thepolygon , ax, OFFSET, interiorsOnly= False):
#     corners = thepolygon.exterior.coords[:]
#     if interiorsOnly : 
#         corners = [] 
#         for interior in thepolygon.interiors : 
#             xy = interior 
#             print(xy)
#             corners+= list(xy.coords)
#     for corner in corners:
#         if corner[0] < 0. : OFFSETX = OFFSET -50
#         if corner[0] > 0. : OFFSETX = OFFSET +50
#         if corner[1] < 0. : OFFSETY = OFFSET -100
#         if corner[1] > 0. : OFFSETY = OFFSET +100
#         # if not(interiorsOnly) : 
#         #     OFFSETX = OFFSETY = OFFSET
#         if corner[0] >0 and corner[1] >0 : 
#             ax.annotate(f'({corner[0]:.2f}, {corner[1]:.2f})', xy=corner, xytext=(corner[0] + OFFSETX, corner[1] + OFFSETY),
#                         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5'), ha='center', va='bottom') #, auto=True)
        
def VeloGeometries(out_dir): 
    ## Velo plane creator 
    def VeloClosed(): 
        print("1. VeloClosed")
        coordsIn  = ( ( -4.7 , -4.7 ),   (-4.7, 4.7)   , (4.7, 4.7) , (4.7 , -4.7) )
        coordsOut = ( ( -35.4 , -35.4 ), (-35.4, 35.4) , (35.4, 35.4) , (35.4 , -35.4) )
        inner = Polygon(coordsIn)
        outer = Polygon(coordsOut)
        inner_rot = affinity.rotate(inner, 45, (0, 0))
        outer_rot = affinity.rotate(outer, 45, (0, 0))
        geometry = outer_rot.difference(inner_rot)
        print("1. done")
        return geometry 
    def VeloLeftSide(shift = 0.): 
        print(f"2. VeloLeftSide shift = {shift}")
        coordsIn  = ( ( -4.7 , -4.7 ),   (-4.7, 4.7)   , (0, 4.7) ,  (0 , -4.7) )
        coordsOut = ( ( -35.4 , -35.4 ), (-35.4, 35.4) , (0, 35.4) , (0 , -35.4) )
        inner = Polygon(coordsIn)
        outer = Polygon(coordsOut)
        inner_rot = affinity.rotate(inner, 45, (0, 0))
        outer_rot = affinity.rotate(outer, 45, (0, 0))
        geometry = outer_rot.difference(inner_rot)
        print("2. done")

        return geometry 
    def VeloRightSide(shift = 0.): 
        print(f"3. VeloRightSide shift = {shift}")
        coordsIn  = ( (  0. , -4.7 ),   ( 0., 4.7)   , (4.7, 4.7) , (4.7 , -4.7) )
        coordsOut = ( (  0. , -35.4 ),  ( 0., 35.4) , (35.4, 35.4) , (35.4 , -35.4) )
        inner = Polygon(coordsIn)
        outer = Polygon(coordsOut)
        inner_rot = affinity.rotate(inner, 45, (0, 0))
        outer_rot = affinity.rotate(outer, 45, (0, 0))
        geometry = outer_rot.difference(inner_rot)
        print(f"3. Done")
        return geometry 
    fig, axs = plt.subplots( 1,2, figsize = (20,10), sharex = True, sharey = True)
    vClose     = VeloClosed()
    vLeft      = VeloLeftSide()
    vRight     = VeloRightSide()
    vLeftOpen  = affinity.translate(vLeft,  xoff=-4.9/sqrt(2),yoff=-4.9/sqrt(2), zoff=0.0)
    vRightOpen = affinity.translate(vRight, xoff=+4.9/sqrt(2),yoff=+4.9/sqrt(2), zoff=0.0)
    print("Velo unary union")
    vOpen = unary_union( [vLeftOpen, vRightOpen])
    plot_polygon(vClose    , label = alabel(vClose, 'VeloClose'),       add_points = False, facecolor = "blue", edgecolor = "black",alpha = 0.1, ax = axs[0], hatch = 'xxx')
    plot_polygon(vLeft     , label = alabel(vLeft, 'VeloLeft'),        add_points = False,  facecolor = "none", edgecolor = "blue",  alpha = 0.8, ax = axs[0]) #, hatch = 'x')
    plot_polygon(vRight    , label = alabel(vRight, 'VeloRight'),       add_points = False,  facecolor = "none", edgecolor = "red",   alpha = 0.8, ax = axs[0]) #, hatch = 'x')
    plot_polygon(vOpen,      label = alabel(vOpen, 'VeloOpen'),        add_points = False, facecolor = "blue", edgecolor = "black",alpha = 0.1, ax = axs[1], hatch = 'xxx')
    plot_polygon(vLeftOpen , label = alabel(vLeftOpen, 'VeloOpenLeft'),    add_points = False, facecolor = "none", edgecolor = "blue",alpha = 0.8, ax = axs[1]) #, hatch = 'x')
    plot_polygon(vRightOpen ,label = alabel(vRightOpen, 'VeloOpenRight'),   add_points = False, facecolor = "none", edgecolor = "red",alpha = 0.8, ax = axs[1]) #, hatch = 'x')
    for i in range(2): axs[i].legend( loc = 'upper right')
    for i in range(2): axs[i].set_xlabel( "x [mm]", loc = 'center')
    for i in range(2): axs[i].set_ylabel( "y [mm]", loc = 'center')
    for i in range(2): axs[i].set_xlim(  -100,100) 
    for i in range(2): axs[i].set_ylim(  -100,100) 
    geometries = { 
        "VeloClose" : VeloClosed(), 
        "VeloLeft"  : VeloLeftSide(), 
        "VeloRight" : VeloRightSide(), 
        "VeloOpen"  : vOpen, 
        "VeloLeftOpen" : vLeftOpen, 
        "VeloRightOpen" : vRightOpen
    }
    fig.savefig( f"{out_dir}/Velo_Geometries.pdf")
    return geometries
def DownstreamTrackerGeometries(out_dir): 
    print("Donwstramtrackers")
    def getalpha( layer): 
        if layer not in ["x","u","v"] : raise ValueError( "Invalid flag for layer type, xuv accepted only")    
        alpha =  0. if layer == "x" else None 
        alpha = -5. if layer == "u" else alpha 
        alpha = +5. if layer == "v" else alpha    
        if alpha == None : raise ValueError( "Invalid alpha, fix me")
        return np.deg2rad(alpha)

    def MightyTracker( flag = "FullSciFi", layer = "x" , fibrePart = True):
        print( f"MightyTracker({flag}) layer = {layer} , fibrePart = {fibrePart}")
        if flag not in ["Frugal", "Modest", "FTDR", "FullSciFi", "Blake"] : raise ValueError("Invalid flag for SciFi region")
        MODULE_SIZE  =  528.0 
        yMax         = +2500.0
        hole_coords  = [(-130, -130)  , (130, -130)  , ( 130, 130) , (-130, 130)  ]   
        alpha = getalpha( layer = layer)    
        edges         = np.array( [  MODULE_SIZE * i  for i in [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]])
        edges_centres = [ (edges[i] + edges[i+1] )/2. for i in range(len(edges)-1)]
        # TODO : Frugal add 
        if flag == "Blake" :
            coords_q1 = (   (0.0 , 0.0),
                            (1525.60, 0.0), 
                            (1525.60, 104.3),
                            (1049.20, 104.3),
                            (1049.20, 312.9),
                            (127.45, 312.9),
                            (127.45, 338.10),
                            (0.0   , 338.10),
                            (0.0   , 300.0) )
            hole_coords = [ 
                (-127.45 , -127.45)  , (127.45, -127.45)  , (127.45, 127.45) , (-127.45, 127.45)
            ]
        elif flag == "Modest": 
            coords_q1 = (   (0.0  , 0.0),
                         (MODULE_SIZE*3, 0.0),
                         (MODULE_SIZE*3, 200), 
                         (MODULE_SIZE*2, 200), 
                         (MODULE_SIZE*2, 300), 
                         (MODULE_SIZE  , 300),
                         (MODULE_SIZE  , 500),
                         (0    , 500))
        elif flag == "Frugal":
            coords_q1 = ( (0.0 , 0.0),
                        (MODULE_SIZE*3, 0.0), 
                        (MODULE_SIZE*3, 100.0),
                        (MODULE_SIZE*2, 100.0),
                        (MODULE_SIZE*2, 200.0),
                        (MODULE_SIZE  , 200.0),
                        (MODULE_SIZE  , 300.0),
                        (0.0          , 300.0) )           
        elif flag == "FTDR" : 
            coords_q1 = (   (0.0  , 0.0),
                        (MODULE_SIZE*4, 0.0),
                        (MODULE_SIZE*4, 200),
                        (MODULE_SIZE*3, 200),
                        (MODULE_SIZE*3, 300),
                        (MODULE_SIZE*2, 300),
                        (MODULE_SIZE*2, 400),
                        (MODULE_SIZE  , 400),
                        (MODULE_SIZE  , 500),
                        (0    , 500))
        elif flag == "FullSciFi": 
            coords_q1 = ( ( -135.0, -120.0,), 
                          (  135.0, -120.0), 
                          (  135.0, 120.0),
                          ( -135.0, 120.0))      
         
        BEAMHOLE  = Polygon(hole_coords)
        poly_q1 = Polygon(coords_q1)
        poly_q2 = affinity.scale( poly_q1, -1, 1, 1 , Point(0,0))
        poly_q3 = affinity.rotate(poly_q1, 180, (0, 0))
        poly_q4 = affinity.scale( poly_q3, -1, 1, 1 , Point(0,0))   
        PIXEL     = unary_union( [poly_q1, poly_q2, poly_q3, poly_q4])
        BEAMHOLE  = Polygon(hole_coords)
        MODULES_UP_Fibre  = [] 
        MODULES_DOWN_Fibre= [] 

        MODULES_UP_Pixels  = [] 
        MODULES_DOWN_Pixels= [] 

        for i in range( len( edges) -1): 
            xMin = edges[i]
            xMax = edges[i+1]
            MODULE_UP   = Polygon( [ ( xMin , 0) , (xMax,0), ( xMax,  yMax) , (xMin,  yMax)])
            MODULE_DOWN = Polygon( [ ( xMin , 0) , (xMax,0), ( xMax, -yMax) , (xMin, -yMax)])        
            """Cut in non-rotated frame the Pixel area out"""
            if flag == "FullSciFi" : 
                MODULE_UP_Fibre     = MODULE_UP   #.difference(   BEAMHOLE)  
                MODULE_DOWN_Fibre   = MODULE_DOWN #.difference( BEAMHOLE)            
            else : 
                MODULE_UP_Fibre     = MODULE_UP.difference(   PIXEL)   if MODULE_UP.overlaps( PIXEL)    else MODULE_UP
                MODULE_DOWN_Fibre   = MODULE_DOWN.difference( PIXEL) if MODULE_DOWN.overlaps( PIXEL)  else MODULE_DOWN

            MODULE_UP_Pixel     = MODULE_UP.difference(   MODULE_UP_Fibre).difference( BEAMHOLE)
            MODULE_DOWN_Pixel   = MODULE_DOWN.difference( MODULE_DOWN_Fibre).difference( BEAMHOLE)

            MODULES_UP_Pixels.append( affinity.rotate( MODULE_UP_Pixel  , alpha , (edges_centres[i],0.), use_radians = True))
            MODULES_DOWN_Pixels.append( affinity.rotate( MODULE_DOWN_Pixel  , alpha , (edges_centres[i],0.), use_radians = True))


            MODULES_UP_Fibre.append(    affinity.rotate( MODULE_UP_Fibre,    alpha,   (edges_centres[i],0.), use_radians=True ))
            MODULES_DOWN_Fibre.append(  affinity.rotate( MODULE_DOWN_Fibre,  alpha,   (edges_centres[i],0.), use_radians=True ))

        SciFi_Up   =  unary_union( MODULES_UP_Fibre)
        SciFi_Down =  unary_union( MODULES_DOWN_Fibre)        
        Pixel_Up   =  unary_union( MODULES_UP_Pixels)
        Pixel_Down =  unary_union( MODULES_DOWN_Pixels)
        SciFi_polygon   =  unary_union( [ SciFi_Up, SciFi_Down])    
        if flag == "FullSciFi" :
            SciFi_polygon = SciFi_polygon.difference( affinity.rotate( BEAMHOLE, alpha, ( 0, 0), use_radians=True))
        Pixel_ploygon   =  unary_union( [ Pixel_Up, Pixel_Down])    
        if fibrePart : return SciFi_polygon
        else : return Pixel_ploygon
        return None 

    fig, axs = plt.subplots( ncols=4,nrows=5, figsize = (40,40), sharex = True, sharey = True)
    shapes = { 
         "Modest_Fibre_x" : MightyTracker(flag = "Modest" , layer = "x", fibrePart = True), 
         "Modest_Pixel_x" : MightyTracker(flag = "Modest" , layer = "x", fibrePart = False),         
         "Modest_Fibre_u" : MightyTracker(flag = "Modest" , layer = "u", fibrePart = True), 
         "Modest_Fibre_v" : MightyTracker(flag = "Modest" , layer = "v", fibrePart = True), 
         "FTDR_Fibre_x" : MightyTracker(flag = "FTDR"   , layer = "x", fibrePart = True), 
         "FTDR_Fibre_u" : MightyTracker(flag = "FTDR"   , layer = "u", fibrePart = True),
         "FTDR_Fibre_v" : MightyTracker(flag = "FTDR"   , layer = "v", fibrePart = True), 
         "FTDR_Pixel_x" : MightyTracker(flag = "FTDR"   , layer = "x", fibrePart = False),
         "Frugal_Fibre_x" : MightyTracker(flag = "Frugal" , layer = "x", fibrePart = True), 
         "Frugal_Pixel_x" : MightyTracker(flag = "Frugal" , layer = "x", fibrePart = False),         
         "Frugal_Fibre_u" : MightyTracker(flag = "Frugal" , layer = "u", fibrePart = True), 
         "Frugal_Fibre_v" : MightyTracker(flag = "Frugal" , layer = "v", fibrePart = True),
         "SciFi_x" : MightyTracker(flag = "FullSciFi"  , layer = "x", fibrePart = True),
         "SciFi_u" : MightyTracker(flag = "FullSciFi"  , layer = "u", fibrePart = True),
         "SciFi_v" : MightyTracker(flag = "FullSciFi"  , layer = "v", fibrePart = True),
         "Blake_Fibre_x" : MightyTracker(flag = "Blake" , layer = "x", fibrePart = True), 
         "Blake_Pixel_x" : MightyTracker(flag = "Blake" , layer = "x", fibrePart = False),         
         "Blake_Fibre_u" : MightyTracker(flag = "Blake" , layer = "u", fibrePart = True), 
         "Blake_Fibre_v" : MightyTracker(flag = "Blake" , layer = "v", fibrePart = True),
         
    }
    plot_polygon( shapes["FTDR_Fibre_x"], add_points= False,label = alabel(shapes["FTDR_Fibre_x"], "FTDR_Fibre_x"), ax = axs[0][0] ,      edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")
    plot_polygon( shapes["FTDR_Pixel_x"], add_points= False,label = alabel(shapes["FTDR_Pixel_x"], "FTDR_Pixel_x"), ax = axs[0][1] ,      edgecolor = 'teal', facecolor = 'blue', alpha = 0.5, hatch = "+++")
    plot_polygon( shapes["FTDR_Fibre_u"], add_points= False,label = alabel(shapes["FTDR_Fibre_u"], "FTDR_Fibre_u"), ax = axs[0][2] ,      edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")
    plot_polygon( shapes["FTDR_Fibre_v"], add_points= False,label = alabel(shapes["FTDR_Fibre_v"], "FTDR_Fibre_v"), ax = axs[0][3] ,      edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")

    plot_polygon( shapes["Modest_Fibre_x"], label = alabel(shapes["Modest_Fibre_x"], "Modest_Fibre_x"), add_points= False, ax = axs[1][0] , edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")
    plot_polygon( shapes["Modest_Pixel_x"], label = alabel(shapes["Modest_Pixel_x"], "Modest_Pixel_x"), add_points= False, ax = axs[1][1] , edgecolor = 'teal', facecolor = 'blue', alpha = 0.5, hatch = "+++")
    plot_polygon( shapes["Modest_Fibre_u"], label = alabel(shapes["Modest_Fibre_u"], "Modest_Fibre_u"), add_points= False, ax = axs[1][2] , edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")
    plot_polygon( shapes["Modest_Fibre_v"], label = alabel(shapes["Modest_Fibre_v"], "Modest_Fibre_v"), add_points= False, ax = axs[1][3] , edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")
  
    plot_polygon( shapes["Frugal_Fibre_x"], add_points= False,label = alabel(shapes["Frugal_Fibre_x"], "Frugal_Fibre_x"), ax = axs[2][0] ,      edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")
    plot_polygon( shapes["Frugal_Pixel_x"], add_points= False,label = alabel(shapes["Frugal_Pixel_x"], "Frugal_Pixel_x"), ax = axs[2][1] ,      edgecolor = 'teal', facecolor = 'blue', alpha = 0.5, hatch = "+++")
    plot_polygon( shapes["Frugal_Fibre_u"], add_points= False,label = alabel(shapes["Frugal_Fibre_u"], "Frugal_Fibre_u"), ax = axs[2][2] ,      edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")
    plot_polygon( shapes["Frugal_Fibre_v"], add_points= False,label = alabel(shapes["Frugal_Fibre_v"], "Frugal_Fibre_v"), ax = axs[2][3] ,      edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")

  
    plot_polygon( shapes["Blake_Fibre_x"], add_points= False,label = alabel(shapes["Blake_Fibre_x"], "Blake_Fibre_x"), ax = axs[3][0] ,      edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")
    plot_polygon( shapes["Blake_Pixel_x"], add_points= False,label = alabel(shapes["Blake_Pixel_x"], "Blake_Pixel_x"), ax = axs[3][1] ,      edgecolor = 'teal', facecolor = 'blue', alpha = 0.5, hatch = "+++")
    plot_polygon( shapes["Blake_Fibre_u"], add_points= False,label = alabel(shapes["Blake_Fibre_u"], "Blake_Fibre_u"), ax = axs[3][2] ,      edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")
    plot_polygon( shapes["Blake_Fibre_v"], add_points= False,label = alabel(shapes["Blake_Fibre_v"], "Blake_Fibre_v"), ax = axs[3][3] ,      edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "||")

    # Blake shape pixel
    fig2, axs2 = plt.subplots(1,4, figsize=(50,10) )
    ax2 = axs2.flatten() 
    plot_polygon( shapes["Blake_Pixel_x"], add_points= False,label = alabel(shapes["Blake_Pixel_x"], "Blake_Pixel_x"), ax =  ax2[0] , edgecolor = 'teal', facecolor = 'blue', alpha = 0.5, hatch = "+++")    
    add_corners_text( thepolygon= shapes["Blake_Pixel_x"], ax = ax2[0], OFFSET= 0 , interiorsOnly=False )
    ax2[0].set_xlim(-2000,2000)
    ax2[0].set_ylim(-600,600)

    plot_polygon( shapes["Blake_Fibre_x"], add_points= False,label = alabel(shapes["Blake_Fibre_x"], "Blake_Fibre_x"), ax =  ax2[1] , edgecolor = 'teal', facecolor = 'blue', alpha = 0.5, hatch = "+++")    
    add_corners_text( thepolygon= shapes["Blake_Fibre_x"], ax = ax2[1], OFFSET= 0 , interiorsOnly=False )
    ax2[1].set_xlim(-4000,4000)
    ax2[1].set_ylim(-4000,4000)    

    plot_polygon( shapes["Blake_Fibre_u"], add_points= False,label = alabel(shapes["Blake_Fibre_u"], "Blake_Fibre_u"), ax =  ax2[2] , edgecolor = 'teal', facecolor = 'blue', alpha = 0.5, hatch = "+++")    
    add_corners_text( thepolygon= shapes["Blake_Fibre_u"], ax = ax2[2], OFFSET= 0 , interiorsOnly=False )
    ax2[2].set_xlim(-4000,4000)
    ax2[2].set_ylim(-4000,4000)    
      
    plot_polygon( shapes["Blake_Fibre_u"], add_points= False,label = alabel(shapes["Blake_Fibre_u"], "Blake_Fibre_u"), ax =  ax2[3] , edgecolor = 'teal', facecolor = 'blue', alpha = 0.5, hatch = "+++")    
    add_corners_text( thepolygon= shapes["Blake_Fibre_u"], ax = ax2[3], OFFSET= 0 , interiorsOnly=False )
    ax2[3].set_xlim(-4000,4000)
    ax2[3].set_ylim(-4000,4000)    
    fig2.tight_layout()  
    fig2.savefig(f"{out_dir}/Blake_Geometry.pdf")
    # # Blake shape pixel
    # fig, ax = plt.subplots(1,1, figsize=(50,20) )
    # plot_polygon( shapes["Blake_Pixel_x"], add_points= False,label = "Blake_Pixel_x", ax =  ax ,      edgecolor = 'teal', facecolor = 'blue', alpha = 0.5, hatch = "+++")    
    # add_corners_text( thepolygon= shapes["Blake_Pixel_x"], ax = ax, OFFSET= 30 , interiorsOnly=True )
    # add_corners_text( thepolygon= shapes["Blake_Pixel_x"], ax = ax, OFFSET= 30 , interiorsOnly=False )
    # ax.set_xlim(-2000,2000)
    # ax.set_ylim(-600,600)
    # fig.savefig(f"{out_dir}/Blake_Pixel_x_ZOOM.pdf")
    


    # fig, ax = plt.subplots(1,1, figsize=(50,20) )
    # plot_polygon( shapes["Frugal_Pixel_x"], add_points= False,label = "Frugal_Pixel_x", ax =  ax ,      edgecolor = 'teal', facecolor = 'blue', alpha = 0.5, hatch = "+++")    
    # add_corners_text( thepolygon= shapes["Frugal_Pixel_x"], ax = ax, OFFSET= 30 , interiorsOnly=True )
    # add_corners_text( thepolygon= shapes["Frugal_Pixel_x"], ax = ax, OFFSET= 30 , interiorsOnly=False )
    # ax.set_xlim(-2000,2000)
    # ax.set_ylim(-600,600)
    # fig.savefig(f"{out_dir}/Frugal_Pixel_x_ZOOM.pdf")

    plot_polygon( shapes["SciFi_x"], add_points= False,label = alabel(shapes["SciFi_x"], "SciFi_x"), ax = axs[4][0] , edgecolor = 'blue', facecolor = 'teal',                alpha = 0.5, hatch = "||")
    plot_polygon( shapes["SciFi_u"], add_points= False,label = alabel(shapes["SciFi_u"], "SciFi_u"), ax = axs[4][1] , edgecolor = 'blue', facecolor = 'teal',                alpha = 0.5, hatch = "||")
    plot_polygon( shapes["SciFi_v"], add_points= False,label = alabel(shapes["SciFi_v"], "SciFi_v"), ax = axs[4][2] , edgecolor = 'blue', facecolor = 'teal',                alpha = 0.5, hatch = "||")


    for i in range(5):
        for j in range(4): 
            axs[i][j].set_xlim( -4200,4200)
            axs[i][j].set_ylim( -4200,4200)
            axs[i][j].legend( loc = 'best')
            axs[i][0].set_xlabel( "x [mm]" , loc = 'center')
            axs[i][0].set_ylabel( "y [mm]" , loc = 'center')        
    fig.savefig(f"{out_dir}/Downstream_Trackers.pdf")
    return shapes 
def UpstreamTrackerGeometries(out_dir):
    print("UpstreamTrackerGeometries")
    def UTGeometries( name ,layer = "x",  _hole_radius_ = 66.8/2. , _nstaves_ = 16, _stave_size_ = 95.5, yMax = 1338./2): 
        if "UTU2" in name :
            FULL_MODULE_X = 1672.0
            FULL_MODULE_Y = 1355.0
            if name == "UTU2_BorderLess" : 
                FULL_MODULE_X = FULL_MODULE_X * 10./12.
                FULL_MODULE_Y = FULL_MODULE_Y * 32./36. 
            HOLE_SIZE_X = 39. 
            HOLE_SIZE_Y = 37.  
            hole_coords   = [(-HOLE_SIZE_X , -HOLE_SIZE_Y)  ,
                            ( HOLE_SIZE_X , -HOLE_SIZE_Y)   , 
                            ( HOLE_SIZE_X ,  HOLE_SIZE_Y)   , 
                            (-HOLE_SIZE_X ,  HOLE_SIZE_Y)    ]
            BEAMHOLE   = Polygon(hole_coords)
            coords = [   (-FULL_MODULE_X/2., -FULL_MODULE_Y/2.),
                ( FULL_MODULE_X/2., -FULL_MODULE_Y/2.),
                ( FULL_MODULE_X/2.,  FULL_MODULE_Y/2.),
                (-FULL_MODULE_X/2.,  FULL_MODULE_Y/2.)]
            External      = Polygon(coords)
            UP  = External.difference(BEAMHOLE)
            return UP
            
        _stave_size_ = 95.5 
        alpha =  0.  if layer == "x" else None 
        alpha =  5. if layer == "u" else alpha
        alpha = -5  if layer == "v" else alpha
        if alpha is None : raise ValueError("Invalid angle")
        edges            = np.array(  [  _stave_size_ * i  for i in np.arange( -_nstaves_/2, _nstaves_/2+1, 1) ])
        edges_centres = [ (edges[i] + edges[i+1] )/2. for i in range(len(edges)-1)]
        BEAM_HOLE = Point(0, 0).buffer(_hole_radius_)
        STAVES = [] 
        for i in range( len( edges) -1): 
            xMin = edges[i]
            xMax = edges[i+1]
            STAVE_FULL  = Polygon( [ ( xMin ,-yMax) , (xMax,-yMax), ( xMax,  yMax) , (xMin,  yMax)])     
            STAVES.append( affinity.rotate( STAVE_FULL  , np.deg2rad(alpha) , (edges_centres[i],0.), use_radians = True))
        UT = unary_union( STAVES)
        UT = UT.difference( BEAM_HOLE)
        return UT

    UTDet = {      
        "UTaX"        : UTGeometries( "UTaX", layer ="x"   , _hole_radius_ = 66.8/2.     , _nstaves_ = 16 , _stave_size_ = 95.5),
        "UTaU"        : UTGeometries( "UTaU", layer ="u"   , _hole_radius_ = 66.8/2.     , _nstaves_ = 16 , _stave_size_ = 95.5),
        "UTbV"        : UTGeometries( "UTbV", layer ="u"   , _hole_radius_ = 66.8/2.     , _nstaves_ = 19 , _stave_size_ = 95.5),
        "UTbX"        : UTGeometries( "UTbX", layer ="x"   , _hole_radius_ = 66.8/2.     , _nstaves_ = 19 , _stave_size_ = 95.5),    
        "UT_U2"       : UTGeometries( "UTU2", layer ="x"   , _hole_radius_ = 66.8/2.     , _nstaves_ = 20 , _stave_size_ = 95.5),
        "UT_U2_BorderLess": UTGeometries( "UTU2_BorderLess", layer ="x", _hole_radius_ = 66.8/2.     , _nstaves_ = 20 , _stave_size_ = 95.5),
        # "UT_U2_Wide" : UTGeometries( "UTU2Wide", layer ="x", _hole_radius_ = 66.8/2. , _nstaves_ = 30 , _stave_size_ = 95.5),    
    }
    fig, axs = plt.subplots( 2,3, figsize = (30,20), sharex = True, sharey = True)

    # ax.text(corner[0], corner[1], f'({corner[0]:.2f}, {corner[1]:.2f})', ha='center', va='bottom')        
    plot_polygon( UTDet["UTaX"],       label = "UTaX", add_points= False, ax = axs[0][0] , edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "+")
    plot_polygon( UTDet["UTaU"],       label = "UTaU", add_points= False, ax = axs[0][1] , edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "+")
    plot_polygon( UTDet["UTbV"],       label = "UTbV", add_points= False, ax = axs[1][1] , edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "+")
    plot_polygon( UTDet["UTbX"],       label = "UTbX", add_points= False, ax = axs[1][0] , edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "+")

    plot_polygon( UTDet["UT_U2"],      label = "UT_U2",      add_points= False, ax = axs[0][2] , edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "+")
    add_corners_text( thepolygon= UTDet["UT_U2"], ax = axs[0][2], OFFSET= 200 , interiorsOnly=False )
   
    plot_polygon( UTDet["UT_U2"],      label = "UT_U2 (base)", add_points= False, ax = axs[1][2] , edgecolor = 'blue', facecolor = 'teal', alpha = 0.1, hatch = "+")
    plot_polygon( UTDet["UT_U2_BorderLess"], label = "UT_U2_BorderLess", add_points= False, ax = axs[1][2] , edgecolor = 'blue', facecolor = 'teal', alpha = 0.5, hatch = "+")
    add_corners_text( thepolygon= UTDet["UT_U2"], ax = axs[1][2], OFFSET = +200, interiorsOnly = True )
    add_corners_text( thepolygon= UTDet["UT_U2"], ax = axs[1][2], OFFSET = -200, interiorsOnly = False )
    add_corners_text( thepolygon= UTDet["UT_U2_BorderLess"], ax = axs[1][2], OFFSET = +200)

   
    for i in range(2):
        for j in range(3): 
            axs[i][j].set_xlim( -1500,1500)
            axs[i][j].set_ylim( -1500,1500)
            axs[i][j].legend( loc = 'best')
            axs[i][0].set_xlabel( "x [mm]" , loc = 'center')
            axs[i][0].set_ylabel( "y [mm]" , loc = 'center')
            axs[i][j].grid(True)
    plt.tight_layout() 
    fig.savefig(f"{out_dir}/Upstream_Tracker_Geometries.pdf")
    return UTDet


if __name__ == '__main__':
    parse = ArgumentParser( description = "Create geometries for the VELO and the trackers, and save them as pickle files")
    parse.add_argument("--output-dir", "-o", default = "Geometries", help = "Output directory to save the geometries")
    args = parse.parse_args()
    os.makedirs( args.output_dir, exist_ok = True)
    os.makedirs( f"{args.output_dir}/plots", exist_ok = True)
    
    Logger.info("Starting geometry creation... VELO ")
    VeloTracker    = VeloGeometries( out_dir = f"{args.output_dir}/plots")
    DownTracker    = DownstreamTrackerGeometries( out_dir = f"{args.output_dir}/plots")
    UpstreamTracker= UpstreamTrackerGeometries( out_dir = f"{args.output_dir}/plots") 


    for label in VeloTracker:
        foutname = f"{args.output_dir}/{label}.pickle"
        with open(foutname,"wb") as fout:
            pickle.dump( VeloTracker[label], fout)
    for label in DownTracker:
        foutname = f"{args.output_dir}/{label}.pickle"
        with open(foutname,"wb") as fout:
            pickle.dump( DownTracker[label], fout)
    for label in UpstreamTracker:
        foutname = f"{args.output_dir}/{label}.pickle"
        with open(foutname,"wb") as fout:
            pickle.dump( UpstreamTracker[label], fout)            
    Logger.info("Creation of geometries and saving done")
