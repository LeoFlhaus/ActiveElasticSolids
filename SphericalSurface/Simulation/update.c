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
// Calculate elastic force
// Calculate and integrate equations of motion
// -------------------------------------------------------------------------------------------------------------------------

void updateeuler(double* x, double* y, double* z, double* px, double* py, double* pz,
	double* fx, double* fy, double* fz, double* tpx, double* tpy, double* tpz, 
	int numb[NP], int list[NP][10], float dist[NP][10], double R) {

	int i;
	double fn, pnorm, vx, vy, vz, nx, ny, nz, nnorm, Dpx, Dpy, Dpz, xold, yold, zold, CosAlpha, SinAlpha, DotProd;
	// ---------------------------------------------------------------------------------------------------------------------
	// Calculate forces and torques
	forceUpdate(x, y, z, px, py, pz, fx, fy, fz, tpx, tpy, tpz, numb, list, dist, xi, R);
	// ---------------------------------------------------------------------------------------------------------------------
	// Integrate  vectors
	for (i = 0; i < NP; i++) {
		xold = x[i]; // Old positions
		yold = y[i];
		zold = z[i];
		fn = K3 * (x[i] * x[i] + y[i] * y[i] + z[i] * z[i] - R * R); // Normal force
		vx = FP * px[i] + fx[i] - x[i] * fn / R; // Calculate velocities
		vy = FP * py[i] + fy[i] - y[i] * fn / R;
		vz = FP * pz[i] + fz[i] - z[i] * fn / R;
		}
		nx = y[i] * vz - z[i] * vy; // Calculate rotation axis and rotate polarity vector, integrate vectors
		ny = z[i] * vx - x[i] * vz;
		nz = x[i] * vy - y[i] * vx;
		nnorm = sqrt(nx * nx + ny * ny + nz * nz);
		nx = nx / nnorm;
		ny = ny / nnorm;
		nz = nz / nnorm;
		x[i] += DT * vx;
		y[i] += DT * vy;
		z[i] += DT * vz;
		DotProd = ((x[i] * xold + y[i] * yold + z[i] * zold ) / (sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]) * sqrt(xold * xold + yold * yold + zold * zold)));
		if (DotProd > 1.0) {
			CosAlpha = 1.0;
		}
		else {
			CosAlpha = DotProd;
		}
		SinAlpha = sqrt(1 - CosAlpha * CosAlpha);
		Dpx = px[i] + DT * tpx[i];
		Dpy = py[i] + DT * tpy[i];
		Dpz = pz[i] + DT * tpz[i];
		px[i] += tpx[i] * DT - SinAlpha * nz * py[i] + SinAlpha * ny * pz[i];
		py[i] += tpy[i] * DT - SinAlpha * nx * pz[i] + SinAlpha * nz * px[i];
		pz[i] += tpz[i] * DT - SinAlpha * ny * px[i] + SinAlpha * nx * py[i];

		pnorm = sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
		if (pnorm != 0) {
			px[i] = px[i] / pnorm;
			py[i] = py[i] / pnorm;
			pz[i] = pz[i] / pnorm;
		}
	}
}



void forceUpdate(double* x, double* y, double* z, double* px, double* py, double* pz,
	double* fx, double* fy, double* fz, double* tpx, double* tpy, double* tpz, 
	int numb[NP], int list[NP][10], float dist[NP][10], double R) {

	int i, j, jj;
	double felix, feliy, feliz, rr, dV, U;

	for (i = 0; i < NP; i++) {
		felix = 0; // First set the force to 0
		feliy = 0;
		feliz = 0;
		dV = 0;
		
		// ---------------------------------------------------------------------------------------------------------------
		// Calculate elastic forces
		for (j = 0; j < numb[i]; ++j) {
			jj = list[i][j];
			rr = (x[jj] - x[i]) * (x[jj] - x[i]) + (y[jj] - y[i]) * (y[jj] - y[i]) + (z[jj] - z[i]) * (z[jj] - z[i]);
			rr = sqrt(rr);
			U = (rr - dist[i][j]);
	
			if (rr != 0) {
				dV = K1 * (U + 20 * U * U * U) / rr;
			}
			felix += -dV * (x[i] - x[jj]);
			feliy += -dV * (y[i] - y[jj]);
			feliz += -dV * (z[i] - z[jj]);
		}
		fx[i] = felix;
		fy[i] = feliy;
		fz[i] = feliz;
		// -----------------------------------------------------------------------------------------------------------------
		// Calculate the torques changing the polarity
		double ftimesp = fx[i] * px[i] + fy[i] * py[i] + fz[i] * pz[i];
		double tx = fx[i] - ftimesp * px[i];
		double ty = fy[i] - ftimesp * py[i];
		double tz = fz[i] - ftimesp * pz[i];
		double tr = tx * x[i] + ty * y[i] + tz * z[i];
		double rtimesp = px[i] * x[i] + py[i] * y[i] + pz[i] * z[i];

		tpx[i] = gamma1 * (tx - tr * x[i] / R) + gamma2 * rtimesp / R * (rtimesp * px[i] - x[i]);
		tpy[i] = gamma1 * (ty - tr * y[i] / R) + gamma2 * rtimesp / R * (rtimesp * py[i] - y[i]);
		tpz[i] = gamma1 * (tz - tr * z[i] / R) + gamma2 * rtimesp / R * (rtimesp * pz[i] - z[i]);
	}
}


double CalcTheta(double x, double y, double z) {
	return acos(z/ sqrt(x * x + y * y + z * z));
}
double CalcPhi(double x, double y) {
	return atan2(y, x);
}

