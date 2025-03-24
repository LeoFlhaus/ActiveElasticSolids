# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 13:25:44 2025

@author: flehi

Perform a Voronoi tessellation on a flat surface
- Input are two .txt files from the initialization
- Output are two files with a) the positions and b) the neighbors
"""

import numpy as np
import scipy as sp
import pandas as pd
from scipy.spatial import SphericalVoronoi



def VoronoiFlat(SimulationFile, InformationFile, NeighborFile, PositionFile):
    # ----------------------------------------------------------------------------------------
    info = pd.read_csv(InformationFile, sep = '\t') # Load files
    data = pd.read_csv(SimulationFile, sep = '\t')
    # ----------------------------------------------------------------------------------------
    data_final = data[data.t == data.t.max()] # Load positions and radii
    x = np.array(data_final.x)
    y = np.array(data_final.y)
    N = info.N[0] # Number of cells
    R = info.R[0] # Disk radius
    # -----------------------------------------------------------------------------------------
    points = [] # Perform Voronoi tessellation
    for i in range(N):
        points.append([x[i],y[i]])
        sv = sp.spatial.Voronoi(points)
        sv.regions.pop(0)
    polygons = [] # Calculate polygons for the Voronoi tessellation
    for i in range(len(sv.regions)):
        polygons.append([])
        for j in range(len(sv.regions[i])):
            if (sv.regions[i][j] != -1):
                if ((sv.vertices[sv.regions[i][j]][0]**2 + sv.vertices[sv.regions[i][j]][1]**2) < R**2):
                    polygons[i].append(np.array([sv.vertices[sv.regions[i][j]][0],sv.vertices[sv.regions[i][j]][1]]))
                if ((sv.vertices[sv.regions[i][j]][0]**2 + sv.vertices[sv.regions[i][j]][1]**2) > R**2):
                    phi = np.arctan2(sv.vertices[sv.regions[i][j]][1],sv.vertices[sv.regions[i][j]][0])
                    polygons[i].append(np.array([R*np.cos(phi),R*np.sin(phi)]))
    neighbors = [] # Estimate the closest neighbors
    for i in range(len(points)):
        if ((points[i][0]**2 + points[i][1]**2) < R**2):
            neighbors.append([])
            for j in range(len(sv.ridge_points)):
                if ((points[sv.ridge_points[j][1]][0]**2 + points[sv.ridge_points[j][1]][1]**2) < R**2):
                    if ((sv.ridge_points[j][0] == i) and (sv.ridge_points[j][1] not in neighbors[i])):
                        neighbors[i].append(sv.ridge_points[j][1])
                if ((points[sv.ridge_points[j][0]][0]**2 + points[sv.ridge_points[j][0]][1]**2) < R**2):
                    if ((sv.ridge_points[j][1] == i) and (sv.ridge_points[j][0] not in neighbors[i])):
                        neighbors[i].append(sv.ridge_points[j][0])  
                                 
    df1 = pd.DataFrame({x[0] : x[1:], y[0]: y[1:]}) # Store positions
    df1.to_csv(PositionFile, sep='\t', index=False)
    NeighbourNumber = [len(neighbors[i]) for i in range(len(neighbors))] # Store neighbors, number of cells, disk radius
    Neighbours = [[] for i in range(10)]
    for i in range(10):
        for j in range(N):
            if (NeighbourNumber[j] > i):
                Neighbours[i].append(neighbors[j][i]+1)
            else:
                Neighbours[i].append(0)
    df2 = pd.DataFrame({R : NeighbourNumber[0:],
                        N : Neighbours[0][0:],
                        '2' : Neighbours[1][0:],
                        '3' : Neighbours[2][0:],
                        '4' : Neighbours[3][0:],
                        '5' : Neighbours[4][0:],
                        '6' : Neighbours[5][0:],
                        '7' : Neighbours[6][0:],
                        '8' : Neighbours[7][0:],
                        '9' : Neighbours[8][0:],
                        '10' : Neighbours[9][0:]})
    df2.to_csv(NeighborFile, sep='\t', index=False)

def VoronoiSphere(SimulationFile, InformationFile, NeighborFile, PositionFile):
    # Load data
    data = pd.read_csv(SimulationFile, sep = '\t')
    info = pd.read_csv(InformationFile, sep = '\t')
    
    data.t.max()
    data_final = data[data.t == data.t.max()]
    x = np.array(data_final.x)
    y = np.array(data_final.y)
    z = np.array(data_final.z)
    R = info.R[0]
    N = info.N[0]

    # Calculate Voronoi tessellation
    points = []
    for i in range(N):
        points.append([x[i],y[i],z[i]])
    radius = info.R[0]
    center = np.array([0, 0, 0])
    sv = SphericalVoronoi(points, radius, center,threshold=1e-03)
    sv.sort_vertices_of_regions()
    # Calculate polygons
    polygons = []
    for i in range(len(sv.regions)):
        polygons.append([])
        for j in range(len(sv.regions[i])):
            polygons[i].append(sv.vertices[sv.regions[i][j]])
    neighbors = []
    for i in range(N):
        neighbors.append([])
        for j in sv.regions[i]:
            for k in range(N):
                if (j in sv.regions[k]):
                    if (k not in neighbors[i]):
                        if (k != i):
                            neighbors[i].append(k)

    df3 = pd.DataFrame({x[0] : x[1:], y[0]: y[1:],z[0]:z[1:]})
    df3.to_csv(PositionFile, sep='\t', index=False)
    NeighbourNumber = [len(neighbors[i]) for i in range(len(neighbors))]
    Neighbours = [[] for i in range(10)]
    for i in range(10):
        for j in range(N):
            if (NeighbourNumber[j] > i):
                Neighbours[i].append(neighbors[j][i]+1)
            else:
                Neighbours[i].append(0)
    df4 = pd.DataFrame({R : NeighbourNumber[0:],
                        N : Neighbours[0][0:],
                        '2' : Neighbours[1][0:],
                        '3' : Neighbours[2][0:],
                        '4' : Neighbours[3][0:],
                        '5' : Neighbours[4][0:],
                        '6' : Neighbours[5][0:],
                        '7' : Neighbours[6][0:],
                        '8' : Neighbours[7][0:],
                        '9' : Neighbours[8][0:],
                        '10' : Neighbours[9][0:]})
    df4.to_csv(NeighborFile, sep='\t', index=False)
