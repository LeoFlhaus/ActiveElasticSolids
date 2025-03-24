# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 16:16:45 2024

@author: flehi
"""

import numpy as np
import pandas as pd



###############################################################################
# Functions to analyse the data
# 1) Calculate the displacement field
# 2) Calculate the velocity field
# 3) Calculate the velocity field in a grid
# 4) Calculate the angular velocity of the polarity field
# 5) Calculate the angular velocity of the polarity field in a grid
# 6) Calculate the angle between two vectors
# 7) Calculate topological defects in the polarity field
# 8) Calculate topological defects in the velocity field
# 9) Calculate the velocity field in 3d
# 10) Calculate the vector product
# 11) Calculate the angular velocity on a sphere
###############################


# 1)
def CalcDisplacementField(simulation, information, tmin):
    data = pd.read_csv(simulation, sep = '\t') #11.11 16:03, 12.11 16:52
    info = pd.read_csv(information, sep = '\t')
    N = info.N[0]
    samples = info.samples[0]-1
    dt = info.dt[0]
    steps = info.stepsinsample[0]
    K = info.K[0]
    x = np.array([data['x'][i*N:(i*N+N)] for i in range(samples)])
    y = np.array([data['y'][i*N:(i*N+N)] for i in range(samples)])
    t = np.linspace(0,K * dt*samples*steps,samples)
    xi = []
    yi = []
    uxi = []
    uyi = []
    ri = []
    thetai = []
    uri = []
    uthetai = []
    ti = t[tmin:]
    for i in range(N):
        x1 = [x[j+tmin][i] for j in range(samples-tmin)]
        y1 = [y[j+tmin][i] for j in range(samples-tmin)] 
        xi.append(x1[0])
        yi.append(y1[0])
        ri.append(np.sqrt(x1[0]**2+y1[0]**2))
        thetai.append(np.arctan2(y1[0] , x1[0]))
        uxi.append([])
        uyi.append([])
        uri.append([])
        uthetai.append([])
        for j in range(samples-tmin):
            uxi[i].append(x1[j]-x1[0])
            uyi[i].append(y1[j]-y1[0])
            uri[i].append(uxi[i][j] * np.cos(thetai[i]) + uyi[i][j] * np.sin(thetai[i]))
            uthetai[i].append(-uxi[i][j] * np.sin(thetai[i]) + uyi[i][j] * np.cos(thetai[i]))    
    return [ri,thetai,uri,uthetai, ti]


# 2)
def CalcVel(simulation, information):
    data = pd.read_csv(simulation, sep = '\t') #Load data
    info = pd.read_csv(information, sep = '\t')
    N = info.N[0]
    samples = info.samples[0]-1
    dt = info.dt[0]
    K = info.K[0]
    steps = info.stepsinsample[0]
    x = np.array([data['x'][i*N:(i*N+N)] for i in range(samples)])
    y = np.array([data['y'][i*N:(i*N+N)] for i in range(samples)])
    vx = []
    vy = []
    t = np.linspace(0,K * dt * steps * samples, samples)
    for j in range(N):
        vx.append(np.gradient([x[i][j] for i in range(samples)],t))
        vy.append(np.gradient([y[i][j] for i in range(samples)],t))
    return [vx,vy]
# 3)
def CalcVelGrid(x,y,vx, vy, GridDist,l,R):
    N = len(vx)
    Grid = [[]for i in range(int(R / GridDist + 1))]#[[]for i in range(int(R))]
    for i in range(int(R / GridDist + 1)):
        for j in range(int(R / GridDist + 1)):
            Grid[i].append([-R + 2 * GridDist * i, -R + 2 * GridDist * j])
    vxgrid = [[[] for j in range(int(R / GridDist + 1))] for i in range(int(R / GridDist + 1))]
    vygrid = [[[] for j in range(int(R / GridDist + 1))] for i in range(int(R / GridDist + 1))]
    for n in range(N):
        for i in range(int(R / GridDist) + 1):
            for j in range(int(R / GridDist) + 1):
                if ((x[l][n] > (-R + 2 * GridDist * i)) and 
                    (x[l][n] < (-R + 2 * GridDist * (i + 1))) and 
                    (y[l][n] > (-R + 2 * GridDist * j)) and 
                    (y[l][n] < (-R + 2 * GridDist * (j + 1)))):
                    vxgrid[i][j].append(vx[n][l])
                    vygrid[i][j].append(vy[n][l])
    vxgridav = [[[] for j in range(int(R / GridDist + 1))] for i in range(int(R / GridDist + 1))]
    vygridav = [[[] for j in range(int(R / GridDist + 1))] for i in range(int(R / GridDist + 1))]
    for i in range(int(R / GridDist + 1)):
        for j in range(int(R / GridDist + 1)):
            if (vxgrid[i][j] != []):
                vxgridav[i][j].append(np.mean(vxgrid[i][j]))
                vygridav[i][j].append(np.mean(vygrid[i][j]))
            else:
                vxgridav[i][j].append(0)
                vygridav[i][j].append(0)
    return [Grid, vxgridav, vygridav]

# 4)
def CalcDqDt(q,t): # Calculate velocity of the polarity phase
    N = len(q[0])
    samples = len(q)
    omega = []
    dtq = []
    for j in range(N):
        qj = [q[i][j] for i in range(samples)]
        dtq.append(np.gradient(qj,t))
        omega.append(np.mean(dtq[j][200:]))
    return omega

# 5)
def CalcPolGrid(x,y,px, py, GridDist,l,R):
    N = len(px[0])
    Grid = [[]for i in range(int(2 * R / GridDist + 1))]#[[]for i in range(int(R))]
    for i in range(int(2 * R / GridDist + 1)):
        for j in range(int(2 * R / GridDist + 1)):
            Grid[i].append([-R + GridDist * i, -R + GridDist * j])
    pxgrid = [[[] for j in range(int(2 * R / GridDist + 1))] for i in range(int(2 * R / GridDist + 1))]
    pygrid = [[[] for j in range(int(2 * R / GridDist + 1))] for i in range(int(2 * R / GridDist + 1))]
    for n in range(N):
        for i in range(int(2 * R / GridDist + 1)):
            for j in range(int(2 * R / GridDist + 1)):
                if ((x[l][n] > (-R + GridDist * i - GridDist / 2)) and 
                    (x[l][n] < (-R + GridDist * i + GridDist / 2)) and 
                    (y[l][n] > (-R + GridDist * j - GridDist / 2)) and 
                    (y[l][n] < (-R + GridDist * j + GridDist / 2))):
                    pxgrid[i][j].append(px[l][n])
                    pygrid[i][j].append(py[l][n])
    pxgridav = [[[] for j in range(int(2 * R / GridDist + 1))] for i in range(int(2 * R / GridDist + 1))]
    pygridav = [[[] for j in range(int(2 * R / GridDist + 1))] for i in range(int(2 * R / GridDist + 1))]
    for i in range(int(2 * R / GridDist + 1)):
        for j in range(int(2 * R / GridDist + 1)):
            if (pxgrid[i][j] != []):
                pxgridav[i][j].append(np.mean(pxgrid[i][j]))
                pygridav[i][j].append(np.mean(pygrid[i][j]))
            else:
                pxgridav[i][j].append(0)
                pygridav[i][j].append(0)
    return [Grid, pxgridav, pygridav]

# 6)
def CalcAngle(a1, a2, b1,b2):
    dot = a1 * b1 + a2 * b2
    det = a1 * b2 - b1 * a2
    angle = np.arctan2(det, dot)
    return angle
# 7)
def CalcTopDefects(x,y,px,py,GridDist,l,R):
    Grid, pxgrid, pygrid = CalcPolGrid(x,y,px,py,GridDist,l,R)
    defects = [[],[],[]]
    for i in range(len(Grid)-18):
        i += 9
        for j in range(len(Grid[0])-18):
            j += 9
            sig = 0
            sig += CalcAngle(pxgrid[i-1][j][0],  pygrid[i-1][j][0],  pxgrid[i-1][j+1][0],pygrid[i-1][j+1][0])
            sig += CalcAngle(pxgrid[i-1][j+1][0],pygrid[i-1][j+1][0],pxgrid[i][j+1][0],  pygrid[i][j+1][0])
            sig += CalcAngle(pxgrid[i][j+1][0],  pygrid[i][j+1][0],  pxgrid[i+1][j+1][0],pygrid[i+1][j+1][0])
            sig += CalcAngle(pxgrid[i+1][j+1][0],pygrid[i+1][j+1][0],pxgrid[i+1][j][0],  pygrid[i+1][j][0])
            sig += CalcAngle(pxgrid[i+1][j][0],  pygrid[i+1][j][0],  pxgrid[i+1][j-1][0],pygrid[i+1][j-1][0])
            sig += CalcAngle(pxgrid[i+1][j-1][0],pygrid[i+1][j-1][0],pxgrid[i][j-1][0],  pygrid[i][j-1][0])
            sig += CalcAngle(pxgrid[i][j-1][0],  pygrid[i][j-1][0],  pxgrid[i-1][j-1][0],pygrid[i-1][j-1][0])
            sig += CalcAngle(pxgrid[i-1][j-1][0],pygrid[i-1][j-1][0],pxgrid[i-1][j][0],  pygrid[i-1][j][0])
            if (np.abs(sig) > 4):
                defects[0].append(i)
                defects[1].append(j)
                defects[2].append(sig)
    defectpos = [-R + np.mean(defects[0]) * GridDist, 
                 - R + GridDist * np.mean(defects[1])]
    defects.append(defectpos)
    return defects

# 8)
def CalcTopDefectsVel(x,y,px,py,GridDist,l,R):
    Grid, pxgrid, pygrid = CalcVelGrid(x,y,px,py,GridDist,l,R)
    defects = [[],[],[]]
    for i in range(len(Grid)-10):
        i += 5
        for j in range(len(Grid[0])-10):
            j += 5
            sig = 0
            sig += CalcAngle(pxgrid[i-1][j][0],  pygrid[i-1][j][0],  pxgrid[i-1][j+1][0],pygrid[i-1][j+1][0])
            sig += CalcAngle(pxgrid[i-1][j+1][0],pygrid[i-1][j+1][0],pxgrid[i][j+1][0],  pygrid[i][j+1][0])
            sig += CalcAngle(pxgrid[i][j+1][0],  pygrid[i][j+1][0],  pxgrid[i+1][j+1][0],pygrid[i+1][j+1][0])
            sig += CalcAngle(pxgrid[i+1][j+1][0],pygrid[i+1][j+1][0],pxgrid[i+1][j][0],  pygrid[i+1][j][0])
            sig += CalcAngle(pxgrid[i+1][j][0],  pygrid[i+1][j][0],  pxgrid[i+1][j-1][0],pygrid[i+1][j-1][0])
            sig += CalcAngle(pxgrid[i+1][j-1][0],pygrid[i+1][j-1][0],pxgrid[i][j-1][0],  pygrid[i][j-1][0])
            sig += CalcAngle(pxgrid[i][j-1][0],  pygrid[i][j-1][0],  pxgrid[i-1][j-1][0],pygrid[i-1][j-1][0])
            sig += CalcAngle(pxgrid[i-1][j-1][0],pygrid[i-1][j-1][0],pxgrid[i-1][j][0],  pygrid[i-1][j][0])
            if (np.abs(sig) > 4):
                defects[0].append(i)
                defects[1].append(j)
                defects[2].append(sig)
    if (defects[0] != []):
        defectpos = [-R + np.mean(defects[0]) * GridDist, 
                     - R + GridDist * np.mean(defects[1])]
        defects.append(defectpos)
    else:
        defects.append([0,0])
    return defects

# 9)
def CalcVelSphere(simulation, information):
    data = pd.read_csv(simulation, sep = '\t') #Load data
    info = pd.read_csv(information, sep = '\t')
    N = info.N[0]
    samples = info.samples[0]-1
    dt = info.dt[0]
    K = info.K[0]
    steps = info.stepsinsample[0]
    x = np.array([data['x'][i*N:(i*N+N)] for i in range(samples)])
    y = np.array([data['y'][i*N:(i*N+N)] for i in range(samples)])
    z = np.array([data['z'][i*N:(i*N+N)] for i in range(samples)])
    vx = []
    vy = []
    vz = []
    t = np.linspace(0,K * dt * steps * samples, samples)
    for j in range(N):
        vx.append(np.gradient([x[i][j] for i in range(samples)],t))
        vy.append(np.gradient([y[i][j] for i in range(samples)],t))
        vz.append(np.gradient([z[i][j] for i in range(samples)],t))
    return [vx,vy,vz]
# 10)
def vector_product(x1,x2,x3,y1,y2,y3):
    z1 = x2*y3-x3*y2
    z2 = x3*y1-x1*y3
    z3 = x1*y2-x2*y1
    return [z1,z2,z3]
# 11)
def CalcAngVelSphere(simulation, information):
    ###########################################################################
    # Load simulation and settings
    ###########################################################################
    data = pd.read_csv(simulation, sep = '\t')
    info = pd.read_csv(information, sep = '\t')
    samples = info.samples[0]
    N       = info.N[0]
    ###########################################################################
    # Save arrays with positions, polarizations and forces
    ###########################################################################
    x = np.array([data['x'][i*N:(i*N+N)] for i in range(samples)])
    y = np.array([data['y'][i*N:(i*N+N)] for i in range(samples)])
    z = np.array([data['z'][i*N:(i*N+N)] for i in range(samples)])
    vx,vy,vz = CalcVelSphere(simulation, information)
    ###########################################################################
    # Calculate the tensor of intertia and its mean values
    ###########################################################################
    M = []
    for i in range(samples):
        M.append([])
        for j in range(N):
            M[i].append([])
            rij = np.sqrt(x[i][j]**2+y[i][j]**2+z[i][j]**2)
            M[i][j] = [[rij**2-x[i][j]**2,-x[i][j]*y[i][j],-x[i][j]*z[i][j]],
                       [-x[i][j]*y[i][j],rij**2-y[i][j]**2,-y[i][j]*z[i][j]],
                       [-x[i][j]*z[i][j],-y[i][j]*z[i][j],rij**2-z[i][j]**2]]
    Mmean = []
    for i in range(samples):
        Mmean.append([[0,0,0],[0,0,0],[0,0,0]])
        for j in range(N):
            for k in range(3):
                for l in range(3):
                    Mmean[i][k][l] += M[i][j][k][l]
    for i in range(samples):
        for k in range(3):
            for l in range(3):
                Mmean[i][k][l] = Mmean[i][k][l]/N   
    ###########################################################################
    # Calculate the angular momentum and its mean values
    ###########################################################################
    L = []
    for i in range(samples):
        L.append([])
        for j in range(N):
            L[i].append(vector_product(x[i][j],y[i][j],z[i][j],vx[i][j],vy[i][j],vz[i][j]))
    Lmean = []
    for i in range(samples):
        Lmean.append([0,0,0])
        for j in range(N):
            for k in range(3):
                Lmean[i][k] += L[i][j][k]
    for i in range(samples):
        for k in range(3):
            Lmean[i][k] = Lmean[i][k]/N
    ###########################################################################
    # Invert the tensor of inertia and calculate the angular velocity
    ###########################################################################
    M_inv = []
    for i in range(samples):
        M_inv.append([np.linalg.inv(Mmean[i])])
    Omega = []
    for i in range(samples):
        Omega.append(np.matmul(M_inv[i],Lmean[i]))
    return Omega