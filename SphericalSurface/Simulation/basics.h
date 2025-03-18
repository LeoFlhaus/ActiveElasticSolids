#ifndef _BASICS_H_
#define _BASICS_H_

// -------------------------------------------------------------------------------------------------------------------------
// Define here constants and basic settings for the simulation
// -------------------------------------------------------------------------------------------------------------------------


#define NP 2000 // Number of cells
#define PositionFile "./ConfigFiles2/xyzN2000.txt" // File with initial positions
#define NeighborFile "./ConfigFiles2/connectN2000.txt" // File with all neighbors from Voronoi tessellation


#define PI 3.14159265
#define DT 0.0001 // Integration time step

#define Samples 500 // samples written to file
#define StepsInSamples 5000//number of steps between each sample

#define gamma1 1.0 // turning rate with force
#define gamma2 1.0 // turning rate  on spherical surface
#define K1 80.0 // Linear term in spring force
#define K3 1.0 // normal force to keep particles on sphere
#define KWall 800.0 // Boundary strength
#define FP 1.0 // Active force


#endif
