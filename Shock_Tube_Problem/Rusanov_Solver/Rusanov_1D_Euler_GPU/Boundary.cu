#include <stdlib.h>

__global__ void GPU_Boundary(float *d_mass_flux, float *d_momentum_flux, float *d_energy_flux, int N_CELLS){
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i <= N_CELLS){
		if (i == 0){
			d_mass_flux[i] = d_mass_flux[i+1];
			d_momentum_flux[i] = d_momentum_flux[i+1];
			d_energy_flux[i] = d_energy_flux[i+1];
		} else if (i == N_CELLS){
			d_mass_flux[i] = d_mass_flux[i-1];
			d_momentum_flux[i] = d_momentum_flux[i-1];
			d_energy_flux[i] = d_energy_flux[i-1];
		}
	}
}

void Boundary(float *d_mass_flux, float *d_momentum_flux, float *d_energy_flux, int N_CELLS){
	int TPB = 128;
	int GPB = (N_CELLS + TPB - 1) / TPB;
	GPU_Boundary<<<GPB, TPB>>>(d_mass_flux, d_momentum_flux, d_energy_flux, N_CELLS);
}
