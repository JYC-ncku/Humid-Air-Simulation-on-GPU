#include <stdlib.h>
#include <math.h>
#include "Compute_Cv.h"

void Initial(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6,
	     float *mass, float *momentum_X, float *momentum_Y, float *energy, float *mass_fraction,
	     float *MAX_Freq, float *MAX_CFL, float R_dry, float R_v, int NX, int NY, int N_CELLS){
/*
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+4);
	int j = (int)INDEX - i * (NY+4);
*/
	float R_mix;
//	if (INDEX < N_CELLS){
//		if (i >= 2 && i < NX + 2 && j >= 2 && j < NY + 2){
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY+4) + j;
			p1[INDEX] = 5.0; // u = 5 m/s
			p2[INDEX] = 0.0; // v = 0 m/s
			p3[INDEX] = 300.0; // T = 300 K
			p4[INDEX] = 101325.0; // P = 1 atm
			p5[INDEX] = 0.0 ; // phi = 0
			p6[INDEX] = 0.0; // RH = 0
			R_mix = R_dry * (1 - p5[INDEX]) + R_v * p5[INDEX];
			p0[INDEX] = p4[INDEX] / (R_mix * p3[INDEX]);
			MAX_Freq[INDEX] = 0.0;
			MAX_CFL[INDEX] = 0.0;

			mass[INDEX] = p0[INDEX];
			momentum_X[INDEX] = p0[INDEX] * p1[INDEX];
			momentum_Y[INDEX] = p0[INDEX] * p2[INDEX];

			float Cv_mix_init = Compute_Cv(p3[INDEX], p5[INDEX]);
			float e_init = p3[INDEX] * Cv_mix_init;
			energy[INDEX] = mass[INDEX] * (e_init + 0.5 * (p1[INDEX] * p1[INDEX] + p2[INDEX] * p2[INDEX]));

			mass_fraction[INDEX] = mass[INDEX] * p5[INDEX];
		}
	}
}
/*
void Initial(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6,
	     float *mass, float *momentum_X, float *momentum_Y, float *energy, float *mass_fraction,
	     float *MAX_Freq, float *MAX_CFL, float R_dry, float R_v, int NX, int NY, int N_CELLS){
	     int TPB = 128;
	     int GPB = (TPB + N_CELLS - 1) / TPB;
	     GPU_Initial<<<GPB, TPB>>>(p0, p1, p2, p3, p4, p5, p6,
				       mass, momentum_X, momentum_Y, energy, mass_fraction,
				       MAX_Freq, MAX_CFL, R_dry,  R_v,  NX,  NY,  N_CELLS);
}
*/
