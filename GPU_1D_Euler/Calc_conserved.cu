#include <stdlib.h>

__global__ void GPU_Calc_conserved(float *d_mass, float *d_momentum, float *d_energy, float *d_mass_flux, float *d_momentum_flux, float *d_energy_flux,
				   float dt, float dx, int N_CELLS){
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i < N_CELLS){
		// Use FVM to get new conservation values
		if ( i >= 0 && i < N_CELLS){
			d_mass[i] = d_mass[i] - (dt / dx)*(d_mass_flux[i+1] - d_mass_flux[i]);
			d_momentum[i] = d_momentum[i] - (dt / dx)*(d_momentum_flux[i+1] - d_momentum_flux[i]);
			d_energy[i] = d_energy[i] - (dt / dx)*(d_energy_flux[i+1] - d_energy_flux[i]);
		}
	}
}
void Calc_conserved(float *d_mass, float *d_momentum, float *d_energy, float *d_mass_flux, float *d_momentum_flux, float *d_energy_flux, float dt, float dx, int N_CELLS){
	int TPB = 128;
	int GPB = (N_CELLS + TPB - 1) / TPB;
	GPU_Calc_conserved<<<TPB, GPB>>>(d_mass, d_momentum, d_energy, d_mass_flux, d_momentum_flux, d_energy_flux, dt, dx, N_CELLS);
}
