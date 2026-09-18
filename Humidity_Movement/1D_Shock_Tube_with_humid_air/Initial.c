#include <stdlib.h>
#include "Compute_Cv.h"

void Initial(float *x, float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *P_sat, float *P_v, float *mass, float *momentum, float *energy, float *mass_fraction,
	     float R_dry, float R_v, float dx, int N_CELLS){
	for (int i = 1; i < N_CELLS + 1; i++){
		if ( i < N_CELLS/2){
			p1[i] = 0.0;
			p2[i] = 300.0; //unit: K
			p3[i] = 5 * 101325.0; //unit: Pa
			p4[i] = 0.0045; //phi_max
			float R_mix_L = R_dry * (1 - p4[i]) + R_v * p4[i];
			p0[i] = p3[i] / (R_mix_L * p2[i]);
		} else{
			p1[i] = 0.0;
			p2[i] = 300.0;
			p3[i] = 0.5 * 101325.0;
			p4[i] = 0.0;
			float R_mix_R = R_dry * (1 - p4[i]) + R_v * p4[i];
			p0[i] = p3[i] / (R_mix_R * p2[i]);
		}
		float T_i = p2[i];
		float Y_i = p4[i];
		float Cv_mix = Compute_Cv(T_i, Y_i);
		mass[i] = p0[i];
		momentum[i] = p0[i] * p1[i];
		energy[i] = 0.5 * p0[i] * p1[i] * p1[i] + p0[i] * Cv_mix * p2[i];
		mass_fraction[i] = p0[i] * p4[i];
		p5[i] = 0.0;
		P_sat[i] = 0.0;
		P_v[i] = 0.0;

	}
}
