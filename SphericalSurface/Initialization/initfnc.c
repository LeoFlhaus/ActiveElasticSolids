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
// Create randomly distributed cells on a sphere
// Sphere radius: R
// Initial cell radii: RPev
// Number of cells: NP
// -------------------------------------------------------------------------------------------------------------------------

void initxyz(double* x, double* y, double* z, double* px, double* py, double* pz, double* RPev, 
	double R, double RP) {
	int i;
	double phi, theta,alpha;

	for (i = 0; i < NP; i++) {
		phi = ((float)(rand()) / (float)(RAND_MAX)) * 2 * PI; 
		theta = ((float)(rand()) / (float)(RAND_MAX)) * PI;
		x[i] = R * sin(theta) * cos(phi);
		y[i] = R * sin(theta) * sin(phi);
		z[i] = R * cos(theta);
		alpha = 2 * PI * (float)rand() / (float)(RAND_MAX); // Initial random polarization 
		px[i] = sin(alpha) * cos(phiinit(x[i], y[i])) * cos(thetainit(x[i], y[i], z[i])) - cos(alpha) * sin(phiinit(x[i], y[i]));
		py[i] = sin(alpha) * sin(phiinit(x[i], y[i])) * cos(thetainit(x[i], y[i], z[i])) + cos(alpha) * cos(phiinit(x[i], y[i]));
		pz[i] = -sin(alpha) * sin(thetainit(x[i], y[i], z[i]));
		RPev[i] = RP + 0.1 * (float)rand() / (float)(RAND_MAX);
	}
}

// -------------------------------------------------------------------------------------------------------------------------
// Calculate spherical coordinates
// -------------------------------------------------------------------------------------------------------------------------
double thetainit(double x, double y, double z) {
	return acos(z / sqrt(x * x + y * y + z * z));
}
double phiinit(double x, double y) {
	return atan2(y, x);
}