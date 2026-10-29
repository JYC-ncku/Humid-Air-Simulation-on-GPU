#include <stdlib.h>

void Calc_conserved(float *mass, float *momentum_X, float *momentum_Y, float *energy,
		    float *mass_flux_X, float *momentum_X_flux_X, float *momentum_Y_flux_X, float *energy_flux_X,
		    float *mass_flux_Y, float *momentum_X_flux_Y, float *momentum_Y_flux_Y, float *energy_flux_Y,
		    float dx, float dy, float dt, int NX, int NY){
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY + 2) + j;
			int INDEX_L = (i - 1) * (NY + 2) + j;
			int INDEX_B = i * (NY + 2) + (j - 1);
			mass[INDEX] = mass[INDEX] - (dt / dx) * (mass_flux_X[INDEX] - mass_flux_X[INDEX_L]) - (dt / dy) * (mass_flux_Y[INDEX] - mass_flux_Y[INDEX_B]);
			momentum_X[INDEX] = momentum_X[INDEX] - (dt / dx) * (momentum_X_flux_X[INDEX] - momentum_X_flux_X[INDEX_L])
							      - (dt / dy) * (momentum_X_flux_Y[INDEX] - momentum_X_flux_Y[INDEX_B]);
			momentum_Y[INDEX] = momentum_Y[INDEX] - (dt / dx) * (momentum_Y_flux_X[INDEX] - momentum_Y_flux_X[INDEX_L])
							      - (dt / dy) * (momentum_Y_flux_Y[INDEX] - momentum_Y_flux_Y[INDEX_B]);
			energy[INDEX] = energy[INDEX] - (dt / dx) * (energy_flux_X[INDEX] - energy_flux_X[INDEX_L])
						      - (dt / dy) * (energy_flux_Y[INDEX] - energy_flux_Y[INDEX_B]);
		}
	}
}
