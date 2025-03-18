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

void initxyz(double* x, double* y, double* px, double* py, double* RPev,  double R, double RP) {
	int i;
	double phi, r;

	for (i = 0; i < NP; i++) {
		phi = ((float)(rand()) / (float)(RAND_MAX)) * 2 * PI;
		r = ((float)(rand()) / (float)(RAND_MAX)) * (R - RP);
		x[i] = r * sin(phi);
		y[i] = r * cos(phi);
		phi = ((float)(rand()) / (float)(RAND_MAX)) * 2 * PI;
		px[i] = cos(phi);
		py[i] = sin(phi);
		RPev[i] = RP + RP * (float)rand() / (float)(RAND_MAX);
	}
}
