#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>



#define NP 5000                                       // Number of cells, 5000 or 20000
#define PositionFile "./ConfigFiles/xyN5000.txt"      // File with initial positions
#define NeighborFile "./ConfigFiles/connectN5000.txt" // File with Delauney triangulation
#define DT 0.0002                                     // integration time step
#define KK 80.0                                       // Spring constant
#define FP 1.0                                        // Active force
#define EPS 1.0                                       // Turning rate
#define KWall 800.0                                   // Boundary strength
#define PI 3.14159265

#define Samples 500
#define StepsInSamples 1500
#define Boundary 1                                    // Boundary condition: 0 for rotation, 1 for oscillation

int list[NP][10], numb[NP];                           // List of particles and neighbors
float dist[NP][10];                                   // List of equilibrium distances

// -------------------------------------------------------------------------------------------------------------------------
// Define functions
// -------------------------------------------------------------------------------------------------------------------------
double* nrvector(int NN);
void   velocities(double* x, double* y, double* q, double* vx, double* vy, double* omega, int* BoundaryParticle, double *BoundaryX, double* BoundaryY, double RR);
// -------------------------------------------------------------------------------------------------------------------------
// Main program
// // ----------------------------------------------------------------------------------------------------------------------
int main() {
	int i, j, k, ii, *BoundaryParticle, ind;
	double RR, rr, * x, * y, * q, * vx, * vy, * omega, * xm, * ym, * qm, t, *BoundaryX, *BoundaryY, xi;
	
	// ---------------------------------------------------------------------------------------------------------------------
	// Initialize arrays
	// ---------------------------------------------------------------------------------------------------------------------
	BoundaryParticle = malloc(NP * sizeof(int)); // Outermost cells
	BoundaryX = nrvector(NP);
	BoundaryY = nrvector(NP);
	x = nrvector(NP);    // Coordinates
	y = nrvector(NP);
	vx = nrvector(NP);   // Velocities
	vy = nrvector(NP);
	omega = nrvector(NP);// Torques
	q = nrvector(NP);    // Polarities
	xm = nrvector(NP);   // For RK2 (midpoints)
	ym = nrvector(NP);
	qm = nrvector(NP);
	// ---------------------------------------------------------------------------------------------------------------------
	// Read initial positions and create random angles
	// ---------------------------------------------------------------------------------------------------------------------
	FILE* ptr = fopen(PositionFile, "r");
	for (i = 0; i < NP; ++i) {
		fscanf(ptr, "%lf %lf", &x[i], &y[i]);
	}
	fclose(ptr);
	for (j = 0; j < NP; ++j) q[j] = 2 * PI * (float)rand() / (float)(RAND_MAX);
	// ---------------------------------------------------------------------------------------------------------------------
	// Establish neighbor list, all entries 0
	// ---------------------------------------------------------------------------------------------------------------------
	for (i = 0; i < NP; ++i){
		numb[i] = 0; // number of neighbours to particle i
		for (j = 0; j < 10; ++j){
			list[i][j] = 0;
			dist[i][j] = 0;
		}
	}
	// ---------------------------------------------------------------------------------------------------------------------
	// Read neighbors from file
	// ---------------------------------------------------------------------------------------------------------------------
	FILE* ptr2 = fopen(NeighborFile, "r");
	fscanf(ptr2, "%lf %d %d %d %d %d %d %d %d %d %d", &RR, &j,&j,&j,&j,&j,&j,&j,&j,&j,&j); 
	for (i = 0; i < NP; ++i) {
		fscanf(ptr2, "%d %d %d %d %d %d %d %d %d %d %d", &numb[i], &list[i][0], &list[i][1], &list[i][2], &list[i][3], &list[i][4], &list[i][5], &list[i][6], &list[i][7], &list[i][8], &list[i][9]);
	}
	fclose(ptr2);
	printf("%lf \n", RR);
	// ---------------------------------------------------------------------------------------------------------------------
	// Establish boundary particle lists
	// ---------------------------------------------------------------------------------------------------------------------
	for (i = 0; i < NP; ++i) {
		BoundaryParticle[i] = 0;
		if ((x[i] * x[i] + y[i] * y[i]) >= (RR - 0.6) * (RR - 0.6)) {
			BoundaryParticle[i] = 1;
		}
		BoundaryX[i] = 0.0;
		BoundaryY[i] = 0.0;
		if (BoundaryParticle[i] == 1) {
			BoundaryX[i] = x[i];
			BoundaryY[i] = y[i];
		}
	}
	for (i = 0; i < NP; ++i) {
		for (j = 0; j < 10; ++j) {
			list[i][j] += -1;      //particle number 0 til NP-1
		}
	}
	// ---------------------------------------------------------------------------------------------------------------------
	// Calculate equlibrium distances
	// ---------------------------------------------------------------------------------------------------------------------
	for (i = 0; i < NP; ++i) {
		for (j = 0; j < numb[i]; ++j) {
			k = list[i][j];
			rr = (x[k] - x[i]) * (x[k] - x[i]) + (y[k] - y[i]) * (y[k] - y[i]);
			rr = sqrt(rr);
			dist[i][j] = rr;
		}
	}
	// ---------------------------------------------------------------------------------------------------------------------
	// Create files with a time stamp
	// - One file for positions and polarities
	// - One for simulation parameters
	// ---------------------------------------------------------------------------------------------------------------------
	char fileName_fp[100];
	char fileName_information[100];
	time_t now = time(NULL);
	struct tm* tt = localtime(&now);
	strftime(fileName_fp, sizeof(fileName_fp) - 1, "./N5000/simulation_%d_%m_%Y_%H_%M.txt", tt);
	strftime(fileName_information, sizeof(fileName_information) - 1, "./N5000/information_%d_%m_%Y_%H_%M.txt", tt);
	FILE* fp;
	FILE* information;
	fp = fopen(fileName_fp, "w"); 	// open simulation and information file
	information = fopen(fileName_information, "w");

	fprintf(information, "N\tv0\tdt\tK\tR\tsamples\tstepsinsample\tgamma\tBoundary\tKWall\n");
	fprintf(information, " %d\t%f\t%f\t%f\t%f  \t%d    \t%d          \t%f   \t%d\t%f",
			NP, FP, DT, KK, RR, Samples, StepsInSamples, EPS, Boundary, KWall);
	fclose(information);
	fprintf(fp, "t\tx\ty\tq\n");
	// ---------------------------------------------------------------------------------------------------------------------
	// time integration
	// - RK2 method for polarity and position vectors
	// ---------------------------------------------------------------------------------------------------------------------
	t = 0.0; // Time
	for (k = 1; k < Samples; ++k) {
		ii = 250;
		for (i = 1; i <= StepsInSamples; ++i) {
			velocities(x, y, q, vx, vy, omega, BoundaryParticle, BoundaryX, BoundaryY, RR);  // calculate velocties
			for (j = 0; j < NP; ++j) {
				xm[j] = x[j] + 0.5 * vx[j] * DT;
				ym[j] = y[j] + 0.5 * vy[j] * DT;
				qm[j] = q[j] + 0.5 * omega[j] * DT;
			}
			velocities(xm, ym, qm, vx, vy, omega, BoundaryParticle, BoundaryX, BoundaryY, RR);
			for (j = 0; j < NP; ++j) {
				x[j] = x[j] + vx[j] * DT;
				y[j] = y[j] + vy[j] * DT;
				q[j] = q[j] + omega[j] * DT;
			}
			t += DT;
		} 
		for (i = 0; i < NP; ++i) fprintf(fp, "%f\t%.2lf\t%.2lf\t%.2lf\n", t, x[i], y[i], q[i]); // Print to file
		printf("%f %lf %lf %lf %lf %lf\n", t, vx[ii], vy[ii], q[ii], -sin(q[ii]) * vx[ii] + cos(q[ii]) * vy[ii], cos(q[ii]) * vx[ii] + sin(q[ii]) * vy[ii]);
	}
	fclose(fp); // Close file
	return 0;   
}

// -------------------------------------------------------------------------------------------------------------------------
// Calculate forces and torques
// -------------------------------------------------------------------------------------------------------------------------
void  velocities(double* x, double* y, double* q, double* vx, double* vy, double* omega, int* BoundaryParticle, double* BoundaryX, double *BoundaryY, double RR) {
	int i, j, jj;
	double rr, dV, U, U3;
	for (i = 0; i < NP; ++i) {
		vx[i] = 0.0;
		vy[i] = 0.0;
		omega[i] = 0.0;
		if (Boundary == 0) { // For a slip-boundary condition
			rr = x[i] * x[i] + y[i] * y[i];   //force from wall
			rr = sqrt(rr);
			if (rr > RR - 0.6) {
				dV = KWall * (rr - RR) * (rr - RR + 0.6);  //zero wall force // at one, keep particles in a potential
				vx[i] += -dV * x[i] / rr;
				vy[i] += -dV * y[i] / rr;
			}
		}
		if (Boundary == 1) { // For a stick-boundary condition
			if (BoundaryParticle[i] == 1) {
				rr = sqrt((BoundaryX[i] - x[i]) * (BoundaryX[i] - x[i]) + (BoundaryY[i] - y[i]) * (BoundaryY[i] - y[i]));
				U = rr;
				U3 = U * U * U;
				dV = KWall * (U);// + 20 * U3);
				vx[i] += -dV * (x[i] - BoundaryX[i]);
				vy[i] += -dV * (y[i] - BoundaryY[i]);
			}
		}

		for (j = 0; j < numb[i]; ++j) {  // force between particles
			jj = list[i][j];
			rr = (x[jj] - x[i]) * (x[jj] - x[i]) + (y[jj] - y[i]) * (y[jj] - y[i]);
			rr = sqrt(rr);
			U = (rr - dist[i][j]);
			U3 = U * U * U;
			dV = KK * (U + 20 * U3) / rr; 

			vx[i] += -dV * (x[i] - x[jj]);
			vy[i] += -dV * (y[i] - y[jj]);
		}
		omega[i] = EPS * (-sin(q[i]) * vx[i] + cos(q[i]) * vy[i]);  //turning torque in the direction of forces or here velocities
		vx[i] = vx[i] + FP * cos(q[i]); 
		vy[i] = vy[i] + FP * sin(q[i]);
	}
}

double* nrvector(int NN)
{
	double* v;
	v = (double*)malloc(NN * sizeof(double));
	return v;
}
