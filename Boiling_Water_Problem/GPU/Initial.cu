#include <stdlib.h>
#include <math.h>
#include "Compute_Cv.h"

__global__ void GPU_Initial(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6,
			    float *d_mass, float *d_momentum_X, float *d_momentum_Y, float *d_energy, float *d_mass_fraction,
			    float R_dry, float R_v, int NX, int NY, int N_CELLS){
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+2);
	int j = (int)INDEX - i * (NY+2);
	float T_i, Y_i, Cv_mix, R_mix;
	if (INDEX < N_CELLS){
		if (i >= 1 && i < NX + 1 && j >= 1 && j < NY + 1){
			d_p1[INDEX] = 5.0; // u = 5 m/s
			d_p2[INDEX] = 0.0; // v = 0 m/s
			d_p3[INDEX] = 300.0; // T = 300 K
			d_p4[INDEX] = 101325.0; // P = 1 atm
			d_p5[INDEX] = 0.0 ; // phi = 0
			d_p6[INDEX] = 0.0; // RH = 0
			R_mix = R_dry * (1 - d_p5[INDEX]) + R_v * d_p5[INDEX];
			T_i = d_p3[INDEX];
			Y_i = d_p5[INDEX];
			Cv_mix = Compute_Cv(T_i, Y_i);
			d_p0[INDEX] = d_p4[INDEX] / (R_mix * d_p3[INDEX]);

			d_mass[INDEX] = d_p0[INDEX];
			d_momentum_X[INDEX] = d_p0[INDEX] * d_p1[INDEX];
			d_momentum_Y[INDEX] = d_p0[INDEX] * d_p2[INDEX];
			d_energy[INDEX] = 0.5 * d_p0[INDEX] * (d_p1[INDEX] * d_p1[INDEX] + d_p2[INDEX] * d_p2[INDEX]) + d_p0[INDEX] * Cv_mix * d_p3[INDEX];
			d_mass_fraction[INDEX] = d_p0[INDEX] * d_p5[INDEX];
		}
	}
}

void Initial(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6,
	     float *d_mass, float *d_momentum_X, float *d_momentum_Y, float *d_energy, float *d_mass_fraction,
	     float R_dry, float R_v, int NX, int NY, int N_CELLS){
	     int TPB = 128;
	     int GPB = (TPB + N_CELLS - 1) / TPB;
	     GPU_Initial<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6,
				       d_mass, d_momentum_X, d_momentum_Y, d_energy, d_mass_fraction,
				       R_dry,  R_v,  NX,  NY,  N_CELLS);
}
