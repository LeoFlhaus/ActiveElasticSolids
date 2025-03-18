#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "update.h"
#include "basics.h"
#include "initfnc.h"

// -------------------------------------------------------------------------------------------------------------------------
// Program to create an initial configuration on a spherical surface
// - Initialize cell in initfnc.c
// - Calculate forces and torques and integrate equations of motion in update.c
// - Basic settings in basics.h
// -------------------------------------------------------------------------------------------------------------------------


int main() {

	double R, RP; // Radius of sphere and cells
	int i, j, k;
	double* x, * y, * px, * py, * fx, * fy, * tpx, * tpy, * RPev, T, vp;
	double rsum;
	// ---------------------------------------------------------------------------------------------------------------------
	// Create files
	// ---------------------------------------------------------------------------------------------------------------------
	char fileName_fptr[100];
	char fileName_information[100];
	time_t now = time(NULL);
	struct tm* t = localtime(&now);
	strftime(fileName_fptr, sizeof(fileName_fptr) - 1, "./simulation_%d_%m_%Y_%H_%M.txt", t);
	strftime(fileName_information, sizeof(fileName_information) - 1, "./information_%d_%m_%Y_%H_%M.txt", t);
	FILE* fptr;
	FILE* information;
	fptr = fopen(fileName_fptr, "w");
	information = fopen(fileName_information, "w");

	// ---------------------------------------------------------------------------------------------------------------------
	// Calculate disk and initial cell radii
	// ---------------------------------------------------------------------------------------------------------------------
	R  = 0.49805873147651974 * sqrt(NP) / 2;  // Same spacing per cell
	RP = 0.4 * R / sqrt(NP);
	// ---------------------------------------------------------------------------------------------------------------------
	// Initialize arrays
	// ---------------------------------------------------------------------------------------------------------------------
	x = malloc(NP * sizeof(double));
	y = malloc(NP * sizeof(double));
	px = malloc(NP * sizeof(double));
	py = malloc(NP * sizeof(double));
	fx = malloc(NP * sizeof(double));
	fy = malloc(NP * sizeof(double));
	tpx = malloc(NP * sizeof(double));
	tpy = malloc(NP * sizeof(double));
	RPev = malloc(NP * sizeof(double));
	// ---------------------------------------------------------------------------------------------------------------------
	T = 0;   // Time
	vp = v0; // Active force
	// ---------------------------------------------------------------------------------------------------------------------
	// Print informations 
	// ---------------------------------------------------------------------------------------------------------------------
	fprintf(information, "N\tv0\tdt\tK\tR\tRP\tsamples\tstepsinsample\tgamma\tepsilon\n");
	fprintf(information, " %d \t %f \t %f \t%f \t%f \t %f \t %d      \t %d              \t  %f    \t %f        ",
		NP, v0, DT, K, R, RP, samples, stepsInSample, gamma, EPS);
	fprintf(fptr, "t\tx\ty\tpx\tpy\tfx\tfy\ttpx\ttpy\trp\n");
	initxyz(x, y, px, py, RPev, NP, R, RP);
	for (j = 0; j < NP; j++) {
		fprintf(fptr, "%f \t %f \t %f\t %f \t %f \t %f \t %f \t %f \t  %f \t  %f \n",
			T, x[j], y[j], px[j], py[j], 0.0, 0.0, 0.0, 0.0, RPev[j]);

	}
	// ---------------------------------------------------------------------------------------------------------------------
	// Initialize cells
	// ---------------------------------------------------------------------------------------------------------------------
	initxyz(x, y, z, px, py, pz, RPev, R, RP);
	for (j = 0; j < NP; j++) {
		fprintf(fptr, "%f \t %f \t %f\t  %f\t %f \t %f \t %f \t %f \t %f \t %f \t %f \t  %f \t  %f \t  %f \n",
			T, x[j], y[j], z[j], px[j], py[j], pz[j], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, RPev[j]);
	}
	// ---------------------------------------------------------------------------------------------------------------------
	// Integrate equations of motion
	// ---------------------------------------------------------------------------------------------------------------------
	for (i = 0; i < samples - 1; i++) {
		for (j = 0; j < stepsInSample; j++) {
			updateeuler(x, y, px, py, fx, fy, tpx, tpy, RPev, vp, T, R);
			T += DT;
			rsum = 0;
			for (k = 0; k < NP; k++) {
				rsum += RPev[k] * RPev[k];
			}
			if ((rsum / (R * R)) < 0.88) {
				for (k = 0; k < NP; k++) {
					RPev[k] += 0.000001;
				}
				vp -= 0.000003;
			}
			else {
				vp = 0;
			}
		}
		for (j = 0; j < NP; j++) {
			fprintf(fptr, "%f \t %f \t %f\t %f \t %f \t %f \t %f \t  %f \t  %f \t  %f \n",
				T, x[j], y[j], px[j], py[j], fx[j], fy[j], tpx[j], tpy[j], RPev[j]);

		}
		fprintf(fptr, "\n");
		printf("Packing %f \n ", rsum / (R * R));
		printf("Status: Sample %d / %d \n", i, samples);
	}

	fclose(fptr); // Close files
	fclose(information); 
	return 0;
}
