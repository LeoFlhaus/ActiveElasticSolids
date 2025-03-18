#ifndef _INITFNC_H
#define _INITFNC_H


void initconfig(double* x, double* y, double* z, double* px, double* py, double* pz);
//void readneighbors(int* Nparticle, int* Nneighbor, double* distance);
double thetainit(double x, double y, double z);
double phiinit(double x, double y);
double* nrvector(int NN);
#endif 
