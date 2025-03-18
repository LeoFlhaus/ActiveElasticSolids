#ifndef _UPDATE_H
#define _UPDATE_H_

void updateeuler(double* x, double* y, double* z, double* px, double* py, double* pz,
	double* fx, double* fy, double* fz, double* tpx, double* tpy, double* tpz, double* RPev, 
	double vp, double T, double R);
void forceUpdate(double* x, double* y, double* z, double* px, double* py, double* pz,
	double* fx, double* fy, double* fz, double* tpx, double* tpy, double* tpz, double* RPev,
	double vp, double R);
double CalcTheta(double z, double R);
double CalcPhi(double x, double y);
double Felx(double x1, double y1, double z1, double x2, double y2, double z2, double RPev1, double RPev2);
double Fely(double x1, double y1, double z1, double x2, double y2, double z2, double RPev1, double RPev2);
double Felz(double x1, double y1, double z1, double x2, double y2, double z2, double RPev1, double RPev2);
#endif

