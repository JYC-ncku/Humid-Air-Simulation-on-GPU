#include <stdlib.h>

void Calc_variable(float *p0, float *p1, float *p2, float *p3, float *p4, float *mass, float *momentum_X, float *momentum_Y, float *energy, float GAMMA, float R, int NX, int NY){
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY + 2) + j;
			p0[INDEX] = mass[INDEX];
			p1[INDEX] = momentum_X[INDEX] / mass[INDEX];
			p2[INDEX] = momentum_Y[INDEX] / mass[INDEX];
			float internal_e = (energy[INDEX] / p0[INDEX]) - 0.5 * (p1[INDEX] * p1[INDEX] + p2[INDEX] * p2[INDEX]);
			float CV = R / (GAMMA - 1.0);
			p3[INDEX] = internal_e / CV;
			p4[INDEX] = p0[INDEX] * R *p3[INDEX];
		}
	}
}
