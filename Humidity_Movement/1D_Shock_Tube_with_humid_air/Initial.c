#include <stdlib.h>

void Initial(float *x, float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *P_sat, float *P_v, float *mass, float *momentum, float *energy, float *mass_fraction, float dx, float GAMMA, int N_CELLS){
	for (int i = 1; i < N_CELLS + 1; i++){
		x[i] = (i+0.5) * dx;
		if ( i < N_CELLS/2){
			p0[i] = 10.0;
			p1[i] = 0.0;
			p2[i] = 1.0;
			p3[i] = 10.0;
			p4[i] = 0.0;
		} else{
			p0[i] = 1.0;
			p1[i] = 0.0;
			p2[i] = 1.0;
			p3[i] = 1.0;
			p4[i] = 0.0;
		}
		mass[i] = p0[i];
		momentum[i] = p0[i] * p1[i];
		energy[i] = 0.5 * p0[i] * p1[i] * p1[i] + (p3[i] / (GAMMA - 1.0));
		mass_fraction[i] = p0[i] * p4[i];
		p5[i] = 0.0;
		P_sat[i] = 0.0;
		P_v[i] = 0.0;
	}
}
