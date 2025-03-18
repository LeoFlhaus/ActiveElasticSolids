#ifndef _BASICS_H_
#define _BASICS_H_

// -------------------------------------------------------------------------------------------------------------------------
// Define here simulation settings
// -------------------------------------------------------------------------------------------------------------------------
#define NP 5000             // Number of cells
#define PI 3.14159265
#define v0 3.0              // Initial active force
#define DT 0.0001           // Integration time step

#define samples 50          // Samples written to file
#define stepsInSample 50000 // Number of steps to run between each sample
#define gamma1 1.0          // Turning rate
#define gamma2 1.0          // Second torque to ensure P remain in the tangent plane
#define EPS 1000.0          // Interaction strength
#define sig 666.7           // Wall strength
#define K3 1.0                // Normal force to keep cells on the sphere
#endif
