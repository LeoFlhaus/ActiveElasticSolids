#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "update.h"
#include "basics.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


// -------------------------------------------------------------------------------------------------------------------------
// Calculate forces and torques
// Integrate equations of motion
// -------------------------------------------------------------------------------------------------------------------------

void updateeuler(double* x, double* y, double* z, double* px, double* py, double* pz,
	double* fx, double* fy, double* fz, double* tpx, double* tpy, double* tpz, double* RPev, double vp, double T, double R) {
	// ---------------------------------------------------------------------------------------------------------------------
	// Euler integration of position and polarity vectors
	// ---------------------------------------------------------------------------------------------------------------------
	int i;
	double fn,norm;
	forceUpdate(x, y, z, px, py, pz, fx, fy, fz, tpx, tpy, tpz, RPev, vp, R); // Calculate forces and torques
	for (i = 0; i < NP; i++) {
		fn = K3 * (x[i] * x[i] + y[i] * y[i] + z[i] * z[i] - R * R); // Normal force
		x[i] += DT * (vp * px[i] + fx[i] - x[i] * fn / R);
		y[i] += DT * (vp * py[i] + fy[i] - y[i] * fn / R);
		z[i] += DT * (vp * pz[i] + fz[i] - z[i] * fn / R);
		px[i] += DT * tpx[i];
		py[i] += DT * tpy[i];
		pz[i] += DT * tpz[i];
		norm = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]); // Normalize P
		px[i] = px[i] / norm;
		py[i] = py[i] / norm;
		pz[i] = pz[i] / norm;
	}
}


void forceUpdate(double* x, double* y, double* z, double* px, double* py, double* pz,
	double* fx, double* fy, double* fz, double* tpx, double* tpy, double* tpz, double* RPev,
	double vp, double R) {
	// ---------------------------------------------------------------------------------------------------------------------
	// Calculate forces and torques for all cells
	// ---------------------------------------------------------------------------------------------------------------------
	int i, j;
	double felix, feliy, feliz; // Elastic force
	double ftimesp, rtimesp;
	for (i = 0; i < NP; i++) {  // Calculate the force for all N cells
		felix = 0; // First set the force to 0
		feliy = 0;
		feliz = 0;
		for (j = 0; j < NP; j++) { //Interaction between particles, Sum elastic forces force for all particles i != j
			if (i != j) {
				felix += Felx(x[i], y[i], z[i], x[j], y[j], z[j], RPev[i], RPev[j]);
				feliy += Fely(x[i], y[i], z[i], x[j], y[j], z[j], RPev[i], RPev[j]);
				feliz += Felz(x[i], y[i], z[i], x[j], y[j], z[j], RPev[i], RPev[j]);
			}
		}
		fx[i] = felix;
		fy[i] = feliy;
		fz[i] = feliz;

		ftimesp = px[i] * fx[i] + py[i] * fy[i] + pz[i] * fz[i];
		rtimesp = px[i] * x[i] + py[i] * y[i] + pz[i] * z[i];
		tpx[i] = gamma1 * (fx[i] - ftimesp * px[i]) + gamma2 * rtimesp / R * (rtimesp - x[i]);
		tpy[i] = gamma1 * (fy[i] - ftimesp * py[i]) + gamma2 * rtimesp / R * (rtimesp - y[i]);
		tpz[i] = gamma1 * (fz[i] - ftimesp * pz[i]) + gamma2 * rtimesp / R * (rtimesp - z[i]);
	}
}


// -------------------------------------------------------------------------------------------------------------------------
// Calculate forces based on Hertzian contact mechanics
// -------------------------------------------------------------------------------------------------------------------------

double Felx(double x1, double y1, double z1, double x2, double y2, double z2, double RPev1, double RPev2) {
	double rij = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
	if (rij < (RPev1 + RPev2)) {
		double c = sig * sqrt(RPev1 * RPev2 / (RPev1 + RPev2));
		return c * pow(RPev1 + RPev2 - rij, 1.5) * (x1 - x2) / rij;
	}
	else {
		return 0;
	}
}
double Fely(double x1, double y1, double z1, double x2, double y2, double z2, double RPev1, double RPev2) {
	double rij = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
	if (rij < (RPev1 + RPev2)) {
		double c = sig * sqrt(RPev1 * RPev2 / (RPev1 + RPev2));
		return c * pow(RPev1 + RPev2 - rij, 1.5) * (y1 - y2) / rij;
	}
	else {
		return 0;
	}
}
double Felz(double x1, double y1, double z1, double x2, double y2, double z2, double RPev1, double RPev2) {
	double rij = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
	if (rij < (RPev1 + RPev2)) {
		double c = sig * sqrt(RPev1 * RPev2 / (RPev1 + RPev2));
		return c * pow(RPev1 + RPev2 - rij, 1.5) * (z1 - z2) / rij;
	}
	else {
		return 0;
	}
}
