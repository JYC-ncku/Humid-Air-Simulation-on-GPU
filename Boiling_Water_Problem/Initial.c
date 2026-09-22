#include <stdlib.h>
#include <math.h>
#include "Compute_Cv.h"

void Initial(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6,
	     float *mass, float *momentum_X, float *momentum_Y, float *energy, float *mass_fraction,
	     float R_dry, float R_v, int NX, int NY){
	float P_sat, P_v, T_i, Y_i, Cv_mix, R_mix, phi_max;
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY+2) + j;
			p1[INDEX] = 5.0; // u = 5 m/s
			p2[INDEX] = 0.0; // v = 0 m/s
			p3[INDEX] = 300.0; // T = 300 K
			p4[INDEX] = 101325.0; // P = 1 atm
			p5[INDEX] = 0.0 ; // phi = 0
			R_mix = R_dry * (1 - p5[INDEX]) + R_v * p5[INDEX];
			Cv_mix = Compute_Cv(T_i, Y_i);
			p0[INDEX] = p4[INDEX] / (R_mix * p3[INDEX]);
			T_i = p3[INDEX];
			Y_i = p5[INDEX];
			mass[INDEX] = p0[INDEX];
			momentum_X[INDEX] = p0[INDEX] * p1[INDEX];
			momentum_Y[INDEX] = p0[INDEX] * p2[INDEX];
			energy[INDEX] = 0.5 * p0[INDEX] * (p1[INDEX] * p1[INDEX] + p2[INDEX] * p2[INDEX]) + p0[INDEX] * Cv_mix * p2[INDEX];
			mass_fraction[INDEX] = p0[INDEX] * p5[INDEX];
			p6[INDEX] = 0.0;
		}
	}
}
