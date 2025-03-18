#ifndef _UPDATE_H
#define _UPDATE_H_
void updateeuler(double* x, double* y, double* px, double* py,
	double* fx, double* fy, double* tpx, double* tpy,  double* RPev, double vp, double T, double R);

void forceUpdate(double* x, double* y, double* px, double* py,
	double* fx, double* fy, double* tpx, double* tpy, double* RPev, double vp, double R);
double theta(double x, double y, double z);
double phi(double x, double y);
double Felx(double x1, double y1, double x2, double y2, double RPev1, double RPev2);
double Fely(double x1, double y1, double x2, double y2, double RPev1, double RPev2);
double FWallX(double x, double y, double RPev, double R);
double FWallY(double x, double y, double RPev, double R);
#endif

