#ifndef _UPDATE_H
#define _UPDATE_H_
#include "basics.h"
void updateeuler(double* x, double* y, double* z, double* px, double* py, double* pz,
	double* fx, double* fy, double* fz, double* tpx, double* tpy, double* tpz, 
	int numb[NP], int list[NP][10], float dist[NP][10], double R);

void forceUpdate(double* x, double* y, double* z, double* px, double* py, double* pz,
	double* fx, double* fy, double* fz, double* tpx, double* tpy, double* tpz,
	int numb[NP], int list[NP][10], float dist[NP][10] double R);

double CalcTheta(double x, double y, double z);
double CalcPhi(double x, double y);
#endif

