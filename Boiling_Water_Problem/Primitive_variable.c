#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "compute_T.h"

void Calc_primitive_variable(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6,
			     float *mass, float *momentum_X, float *momentum_Y, float *energy, float *mass_fraction,
			     float *mass_flux_X, float *momentum_X_flux_X, float *momentum_Y_flux_X, float *energy_flux_X, float *mass_fraction_flux_X,
			     float *mass_flux_Y, float *momentum_X_flux_Y, float *momentum_Y_flux_Y, float *energy_flux_Y, float *mass_fraction_flux_Y,
			     float R_dry, float R_v, float dx, float dy, float dt, int NX, int NY){
	float e_target, T_old, phi_new, phi_max,  R_mix, P_sat, P_v;
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY+2) + j;
			int INDEX_R = (i+1) * (NY+2) + j;
			int INDEX_T = i * (NY+2) + (j+1);
		        // Use FVM to get new conservation values
			mass[INDEX] = mass[INDEX] - (dt / dx) * (mass_flux_X[INDEX_R] - mass_flux_X[INDEX])
						  - (dt / dy) * (mass_flux_Y[INDEX_T] - mass_flux_Y[INDEX]);
			momentum_X[INDEX] = momentum_X[INDEX] - (dt / dx) * (momentum_X_flux_X[INDEX_R] - momentum_X_flux_X[INDEX])
							      - (dt / dy) * (momentum_X_flux_Y[INDEX_T] - momentum_X_flux_Y[INDEX]);
			momentum_Y[INDEX] = momentum_Y[INDEX] - (dt / dx) * (momentum_Y_flux_X[INDEX_R] - momentum_Y_flux_X[INDEX])
							      - (dt / dy) * (momentum_Y_flux_Y[INDEX_T] - momentum_Y_flux_Y[INDEX]);
			energy[INDEX] = energy[INDEX] - (dt / dx) * (energy_flux_X[INDEX_R] - energy_flux_X[INDEX])
						      - (dt / dy) * (energy_flux_Y[INDEX_T] - energy_flux_Y[INDEX]);
			mass_fraction[INDEX] = mass_fraction[INDEX] - (dt / dx) * (mass_fraction_flux_X[INDEX_R] - mass_fraction_flux_X[INDEX])
								    - (dt / dy) * (mass_fraction_flux_Y[INDEX_T] - mass_fraction_flux_Y[INDEX]);

			//Get new variable
			p0[INDEX] = mass[INDEX];
			p1[INDEX] = momentum_X[INDEX] / mass[INDEX];
			p2[INDEX] = momentum_Y[INDEX] / mass[INDEX];
			p5[INDEX] = mass_fraction[INDEX] / p0[INDEX];
			e_target = (energy[INDEX] / p0[INDEX]) - 0.5 * (p1[INDEX] * p1[INDEX] + p2[INDEX] * p2[INDEX]);
			T_old = p3[INDEX];
			phi_new = p5[INDEX];
			p3[INDEX] = compute_T(T_old, phi_new, e_target);
			R_mix = R_dry * (1.0 - p5[INDEX]) + R_v * p5[INDEX];
			p4[INDEX] = p0[INDEX] * R_mix * p3[INDEX];

			P_sat = 611.0 * exp((17.27 * (p3[INDEX] - 273.15)) / ((p3[INDEX] - 273.15) + 237.3)); // Tetens equation
			P_v = (p0[INDEX] * p5[INDEX]) * R_v * p2[INDEX];
			p6[INDEX] = P_v / P_sat;

			phi_max = (P_sat / R_v) / (((p4[INDEX] - P_sat) / R_dry) + (P_sat / R_v));

			if (p5[INDEX] > phi_max || p6[INDEX] > 1.0){
				p5[INDEX] = phi_max;
				p6[INDEX] = 1.0; // 100%!
			}
		}
	}
}
