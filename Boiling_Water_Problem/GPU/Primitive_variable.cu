#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "compute_T.h"

__global__ void GPU_Calc_primitive_variable(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6,
					    float *d_mass, float *d_momentum_X, float *d_momentum_Y, float *d_energy, float *d_mass_fraction,
					    float *d_mass_flux_X, float *d_momentum_X_flux_X, float *d_momentum_Y_flux_X, float *d_energy_flux_X, float *d_mass_fraction_flux_X,
					    float *d_mass_flux_Y, float *d_momentum_X_flux_Y, float *d_momentum_Y_flux_Y, float *d_energy_flux_Y, float *d_mass_fraction_flux_Y,
					    float R_dry, float R_v, float dx, float dy, float dt, int NX, int NY, int N_CELLS){
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+2);
	int j = (int)INDEX - i * (NY+2);
	float e_target, T_old, phi_new, phi_max,  R_mix, P_sat, P_v;
	if (INDEX < N_CELLS){
		if (i >= 1 && i < NX + 1 && j >= 1 && j < NY + 1){
			int INDEX = i * (NY+2) + j;
			int INDEX_R = (i+1) * (NY+2) + j;
			int INDEX_T = i * (NY+2) + (j+1);
		        // Use FVM to get new conservation values
			d_mass[INDEX] = d_mass[INDEX] - (dt / dx) * (d_mass_flux_X[INDEX_R] - d_mass_flux_X[INDEX])
						  - (dt / dy) * (d_mass_flux_Y[INDEX_T] - d_mass_flux_Y[INDEX]);
			d_momentum_X[INDEX] = d_momentum_X[INDEX] - (dt / dx) * (d_momentum_X_flux_X[INDEX_R] - d_momentum_X_flux_X[INDEX])
								  - (dt / dy) * (d_momentum_X_flux_Y[INDEX_T] - d_momentum_X_flux_Y[INDEX]);
			d_momentum_Y[INDEX] = d_momentum_Y[INDEX] - (dt / dx) * (d_momentum_Y_flux_X[INDEX_R] - d_momentum_Y_flux_X[INDEX])
								  - (dt / dy) * (d_momentum_Y_flux_Y[INDEX_T] - d_momentum_Y_flux_Y[INDEX]);
			d_energy[INDEX] = d_energy[INDEX] - (dt / dx) * (d_energy_flux_X[INDEX_R] - d_energy_flux_X[INDEX])
							  - (dt / dy) * (d_energy_flux_Y[INDEX_T] - d_energy_flux_Y[INDEX]);
			d_mass_fraction[INDEX] = d_mass_fraction[INDEX] - (dt / dx) * (d_mass_fraction_flux_X[INDEX_R] - d_mass_fraction_flux_X[INDEX])
									- (dt / dy) * (d_mass_fraction_flux_Y[INDEX_T] - d_mass_fraction_flux_Y[INDEX]);

			//Get new variable
			d_p0[INDEX] = d_mass[INDEX];
			d_p1[INDEX] = d_momentum_X[INDEX] / d_mass[INDEX];
			d_p2[INDEX] = d_momentum_Y[INDEX] / d_mass[INDEX];
			d_p5[INDEX] = d_mass_fraction[INDEX] / d_p0[INDEX];
			e_target = (d_energy[INDEX] / d_p0[INDEX]) - 0.5 * (d_p1[INDEX] * d_p1[INDEX] + d_p2[INDEX] * d_p2[INDEX]);
			T_old = d_p3[INDEX];
			phi_new = d_p5[INDEX];
			d_p3[INDEX] = compute_T(T_old, phi_new, e_target);
			R_mix = R_dry * (1.0 - d_p5[INDEX]) + R_v * d_p5[INDEX];
			d_p4[INDEX] = d_p0[INDEX] * R_mix * d_p3[INDEX];

			P_sat = 611.0 * exp((17.27 * (d_p3[INDEX] - 273.15)) / ((d_p3[INDEX] - 273.15) + 237.3)); // Tetens equation
			P_v = (d_p0[INDEX] * d_p5[INDEX]) * R_v * d_p3[INDEX];
			d_p6[INDEX] = P_v / P_sat;
			if (d_p4[INDEX] <= P_sat) {
				phi_max = 1.0; // 壓力低於飽和蒸氣壓，允許 100% 水氣
			} else {
				phi_max = (P_sat / R_v) / (((d_p4[INDEX] - P_sat) / R_dry) + (P_sat / R_v));
			}

			if (d_p5[INDEX] > phi_max || d_p6[INDEX] > 1.0){
				d_p5[INDEX] = phi_max;
				d_p6[INDEX] = 1.0; // 100%!
				d_mass_fraction[INDEX] = d_p5[INDEX] * d_p0[INDEX];
			}
			if (d_p5[INDEX] < 0.0 || d_p6[INDEX] < 0.0){
				d_p5[INDEX] = 0.0;
				d_p6[INDEX] = 0.0;
				d_mass_fraction[INDEX] = 0.0;
			}

		}
	}
}

void Calc_primitive_variable(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6,
			     float *d_mass, float *d_momentum_X, float *d_momentum_Y, float *d_energy, float *d_mass_fraction,
			     float *d_mass_flux_X, float *d_momentum_X_flux_X, float *d_momentum_Y_flux_X, float *d_energy_flux_X, float *d_mass_fraction_flux_X,
			     float *d_mass_flux_Y, float *d_momentum_X_flux_Y, float *d_momentum_Y_flux_Y, float *d_energy_flux_Y, float *d_mass_fraction_flux_Y,
			     float R_dry, float R_v, float dx, float dy, float dt, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Calc_primitive_variable<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6,
						  d_mass, d_momentum_X, d_momentum_Y, d_energy, d_mass_fraction,
						  d_mass_flux_X, d_momentum_X_flux_X, d_momentum_Y_flux_X, d_energy_flux_X, d_mass_fraction_flux_X,
						  d_mass_flux_Y, d_momentum_X_flux_Y, d_momentum_Y_flux_Y, d_energy_flux_Y, d_mass_fraction_flux_Y,
						  R_dry, R_v, dx, dy, dt, NX, NY, N_CELLS);
}
