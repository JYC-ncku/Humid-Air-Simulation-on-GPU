#include <stdlib.h>

__global__ void GPU_Calc_variable(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_mass, float *d_momentum, float *d_energy, float GAMMA, float R, int N_CELLS){
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < N_CELLS){
		if ( i >= 0 && i < N_CELLS){
			d_p0[i] = d_mass[i];
			d_p1[i] = d_momentum[i] / d_mass[i];
			d_p3[i] = (GAMMA - 1) * (d_energy[i] - 0.5 * d_p0[i] * d_p1[i] * d_p1[i]);
			d_p2[i] = d_p3[i] / (d_p0[i] * R);
		}
	}
}

void Calc_variable(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_mass, float *d_momentum, float *d_energy, float GAMMA, float R, int N_CELLS){
	int TPB = 128;
	int GPB = (N_CELLS + TPB - 1) / TPB;
	GPU_Calc_variable<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_mass, d_momentum, d_energy, GAMMA, R, N_CELLS);
}
