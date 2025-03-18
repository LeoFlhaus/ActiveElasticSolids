#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "initfnc.h"
#include "basics.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


// -------------------------------------------------------------------------------------------------------------------------
// * Program to initialize the simulation
// * Read the positions from the initialization 
// * Calculate (random) initial polarizations
// -------------------------------------------------------------------------------------------------------------------------

void initconfig(double* x, double* y, double* z, double* px, double* py, double* pz) { // Read positions
	FILE* config;
	config = fopen(PositionFile, "r");
	int i;
	double alpha;
	for (i = 0; i < NP; i++) {
		fscanf(config, "%lf  %lf  %lf", &x[i], &y[i], &z[i]); // Scan file
		// Calculate polarity vectors
		//alpha = 0;
		alpha = 2 * PI * (float)rand() / (float)(RAND_MAX);
		px[i] = sin(alpha) * cos(phiinit(x[i], y[i])) * cos(thetainit(x[i], y[i], z[i])) - cos(alpha) * sin(phiinit(x[i], y[i]));
		py[i] = sin(alpha) * sin(phiinit(x[i], y[i])) * cos(thetainit(x[i], y[i], z[i])) + cos(alpha) * cos(phiinit(x[i], y[i]));
		pz[i] = -sin(alpha) * sin(thetainit(x[i], y[i], z[i]));
	}
	printf("x = %f \t y = %f\t z = %f\n", x[0], y[0], z[0]);
	fclose(config);
}

// Functions to calculate spherical coordinates from Cartesian 
double thetainit(double x, double y, double z) {
	return acos(z / sqrt(x * x + y * y + z * z));
}
double phiinit(double x, double y) {
	return atan2(y, x);
}

double* nrvector(int NN)
{
	double* v;
	v = (double*)malloc(NN * sizeof(double));
	return v;
}
