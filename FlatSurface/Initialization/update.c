#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "update.h"
#include "basics.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void updateeuler(double* x, double* y, double* px, double* py,
	double* fx, double* fy, double* tpx, double* tpy, double* RPev, double vp, double T, double R) {

	int i;
	double temp; 

	forceUpdate(x, y, px, py, fx, fy, tpx, tpy, RPev, vp, R);

	for (i = 0; i < NP; i++) {

		x[i] += DT * fx[i];
		y[i] += DT * fy[i];
		px[i] += DT * tpx[i];
		py[i] += DT * tpy[i];

	}

	for (i = 0; i < NP; i++) {
		temp = sqrt(px[i] * px[i] + py[i] * py[i]);
		px[i] = px[i] / temp;
		py[i] = py[i] / temp;
	}
}

void forceUpdate(double* x, double* y, double* px, double* py,
	double* fx, double* fy, double* tpx, double* tpy, double* RPev, 
	double vp, double R) {
	int i, j;
	double felix, feliy; // Elastic force
	for (i = 0; i < NP; i++) {  // Calculate the force for all N particles
		felix = 0; // First set the force to 0
		feliy = 0;
		for (j = 0; j < NP; j++) { //Interaction between particles, Sum elastic forces force for all particles i != j
			if (i != j) {
				felix += Felx(x[i], y[i], x[j], y[j], RPev[i], RPev[j]);
				feliy += Fely(x[i], y[i], x[j], y[j], RPev[i], RPev[j]);			}
		}
		felix += FWallX(x[i], y[i], RPev[i], R);
		feliy += FWallY(x[i], y[i], RPev[i], R);

		// Calculate the force on each particle as the sum of the active force, the elastic force and the normal force
		fx[i] = vp * px[i] + felix;
		fy[i] = vp * py[i] + feliy;

		// Calculate the torques changing the polarity
		tpx[i] = gamma1 * (fx[i] - (px[i] * fx[i] + py[i] * fy[i]) * px[i]); 
		tpy[i] = gamma1 * (fy[i] - (px[i] * fx[i] + py[i] * fy[i]) * py[i]);
	}
}


double Felx(double x1, double y1, double x2, double y2, double RPev1, double RPev2) {
	double rij = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	double c = 4.0 / (6.0 * sig) * sqrt(RPev1 * RPev2 / (RPev1 + RPev2));
	if (rij < (RPev1 + RPev2)) {
		return c * pow(RPev1 + RPev2 - rij, 1.5) * (x1 - x2) / rij;
	}
	else {
		return 0;
	}
}
double Fely(double x1, double y1, double x2, double y2, double RPev1, double RPev2) {
	double rij = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	double c = 4.0 / (6.0 * sig) * sqrt(RPev1 * RPev2 / (RPev1 + RPev2));
	if (rij < (RPev1 + RPev2)) {
		return c * pow(RPev1 + RPev2 - rij, 1.5) * (y1 - y2) / rij;
	}
	else {
		return 0;
	}
}
double FWallX(double x, double y, double RPev, double R) {
	double r = sqrt(x * x + y * y) + RPev;
	if (r > R ) {
		return -EPS * (r - R) * x;
	}
	else {
		return 0;
	}
}
double FWallY(double x, double y, double RPev, double R) {
	double r  = sqrt(x * x + y * y) + RPev;
	if (r > (R)) {
		return -EPS * (r - R) * y;
	}
	else {
		return 0;
	}
}
