#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "update.h"
#include "basics.h"
#include "initfnc.h"


// -------------------------------------------------------------------------------------------------------------------------
// Simulate active elastic solids on a spherical surface
// - Initialize the cells in initfnc.c
// - Integrate the equations of motion in update.c
// - Basic settings in basics.h
// -------------------------------------------------------------------------------------------------------------------------

int main() {
	// ---------------------------------------------------------------------------------------------------------------------
	//  Declaration
	// ---------------------------------------------------------------------------------------------------------------------
	int i, j, k, ii, ind; // Inidices
	double T, R, r, xi; 
	double* x, * y, * z, * px, * py, * pz;
	double * fx, * fy, * fz, * tpx, * tpy, * tpz; // Vectors with positions, polarizations, forces and torques
	int list[NP][10],numb[NP]; // Define neighbors for force calculation
	float dist[NP][10];

	x = nrvector(NP);
	y = nrvector(NP);
	z = nrvector(NP);
	fx = nrvector(NP);
	fy = nrvector(NP);
	fz = nrvector(NP);
	px = nrvector(NP);
	py = nrvector(NP);
	pz = nrvector(NP);
	tpx = nrvector(NP);
	tpy = nrvector(NP);
	tpz = nrvector(NP);
	// -------------------------------------
	initconfig(x, y, z, px, py, pz); // Initialize positions and polarizations in initfnc.c
	// -------------------------------------
	for (i = 0; i < NP; ++i) {
		numb[i] = 0; // Number of neighbours to particle i
		for (j = 0; j < 10; ++j) {
			list[i][j] = 0;
			dist[i][j] = 0;
		}
	}
	// -------------------------------------
	FILE* ptr2 = fopen(NeighborFile, "r");
	fscanf(ptr2, "%lf %d %d %d %d %d %d %d %d %d %d", &R, &j, &j, &j, &j, &j, &j, &j, &j, &j, &j);
	for (i = 0; i < NP; ++i) {
		fscanf(ptr2, "%d %d %d %d %d %d %d %d %d %d %d", &numb[i], &list[i][0], &list[i][1], &list[i][2], &list[i][3], &list[i][4], &list[i][5], &list[i][6], &list[i][7], &list[i][8], &list[i][9]);
	}
	fclose(ptr2);
	printf("%lf \n", R);
	// ------------------------------------
	for (i = 0; i < NP; ++i) {
		for (j = 0; j < 10; ++j) {
			list[i][j] += -1; // particle number for 0 to N-1
		}
	}
	// ------------------------------------
	for (i = 0; i < NP; ++i) {
		for (j = 0; j < numb[i]; ++j) {
			k = list[i][j];
			r = (x[k] - x[i]) * (x[k] - x[i]) + (y[k] - y[i]) * (y[k] - y[i]) + (z[k] - z[i]) * (z[k] - z[i]);
			r = sqrt(r);
			dist[i][j] = r;
		}
	}
	// -------------------------------------
	// Create files with a time stamp
	char fileName_fp[100];
	char fileName_information[100];
	time_t now = time(NULL);
	struct tm* tt = localtime(&now);
	strftime(fileName_fp, sizeof(fileName_fp) - 1,  "./R/simulation_%d_%m_%Y_%H_%M.txt", tt);
	strftime(fileName_information, sizeof(fileName_information) - 1, "./R/information_%d_%m_%Y_%H_%M.txt", tt);
	
	FILE* fp;
	FILE* information;
	fp = fopen(fileName_fp, "w"); 	// open simulation and information file
	information = fopen(fileName_information, "w");
	fprintf(information, "N\tv0\tdt\tK\tR\tsamples\tstepsinsample\tgamma\tgamma2\n");
	fprintf(information, " %d\t%f\t%f\t%f\t%f  \t%d    \t%d          \t%f\t%f",
		NP, FP, DT, K1, R, Samples, StepsInSamples, gamma1, gamma2);
	fclose(information);
	T = 0.0;
	fprintf(fp, "t\tx\ty\tz\tpx\tpy\tpz\n");
	for (i = 0; i < NP; ++i) {
		fprintf(fp, "%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",
			T, x[i], y[i], z[i], px[i], py[i], pz[i]);
	}
	// --------------------------------------------------------
	// Integrate the equations of motion
	// --------------------------------------------------------
	for (k = 1; k < Samples; ++k) {
		ii = 250;
		for (i = 1; i <= StepsInSamples; ++i) {
			updateeuler(x, y, z, px, py, pz, fx, fy, fz, tpx, tpy, tpz, numb,list, dist, R);
			T += DT;
		}
		for (i = 0; i < NP; ++i) {
			fprintf(fp, "%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", T, x[i],y[i],z[i],px[i],py[i],pz[i]);
		}
		printf("%lf %lf %lf \n", x[ii], px[ii], fx[ii]);
	}
	fclose(fp);
	return 0;
}
