#include <stdlib.h>

void Initial(float *x, float *p0, float *p1, float *p2, float *p3, float *mass, float *momentum, float *energy, float GAMMA, float dx, int N_CELLS){
	for (int i = 0; i < N_CELLS; i++){
		x[i] = (i + 0.5) * dx;
		if(i < N_CELLS/2){
			p0[i] = 10.0;
			p1[i] = 0.0;
			p2[i] = 1.0;
			p3[i] = 10.0;
		} else{
			p0[i] = 1.0;
			p1[i] = 0.0;
			p2[i] = 1.0;
			p3[i] = 1.0;
		}
		mass[i] = p0[i];
		momentum[i] = p0[i] * p1[i];
		energy[i] = 0.5 * p0[i] * p1[i] * p1[i] + (p3[i] / (GAMMA - 1));
	}
}
