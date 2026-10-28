#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "compute_T.h"

void Calc_primitive_variable(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *mass, float *momentum, float *energy, float *mass_fraction,
			     float *mass_flux, float *momentum_flux, float *energy_flux, float *mass_fraction_flux, float R_dry, float R_v, float dx, float dt,
			     int N_CELLS){
	float phi_max;
	for (int i = 2; i < N_CELLS + 2; i++){
	        // Use FVM to get new conservation values
		mass[i] = mass[i] - (dt / dx) * (mass_flux[i] - mass_flux[i-1]);
		momentum[i] = momentum[i] - (dt / dx) * (momentum_flux[i] - momentum_flux[i-1]);
		energy[i] = energy[i] - (dt / dx) * (energy_flux[i] - energy_flux[i-1]);
		mass_fraction[i] = mass_fraction[i] - (dt / dx) * (mass_fraction_flux[i] - mass_fraction_flux[i-1]);
		//Get new variable
		p0[i] = mass[i];
		p1[i] = momentum[i] / mass[i];
		p4[i] = mass_fraction[i] / p0[i];
		float e_target = (energy[i] / p0[i]) - 0.5 * p1[i] * p1[i];
		float T_old = p2[i];
		float phi_old = p4[i];
		p2[i] = compute_T(T_old, phi_old, e_target);
		float R_mix = R_dry * (1.0 - p4[i]) + R_v * p4[i];
		p3[i] = p0[i] * R_mix * p2[i];

		float P_sat = 611.0 * exp((17.27 * (p2[i] - 273.15)) / ((p2[i] - 273.15) + 237.3)); // Tetens equation
		float P_v = (p0[i] * p4[i]) * R_v * p2[i];
		p5[i] = P_v / P_sat;

		if (p3[i] <= P_sat) {
			phi_max = 1.0; // 壓力低於飽和蒸氣壓，允許 100% 水氣
		} else {
			phi_max = (P_sat / R_v) / (((p3[i] - P_sat) / R_dry) + (P_sat / R_v));
		}

		if (p4[i] > phi_max || p5[i] > 1.0){
			p4[i] = phi_max;
			p5[i] = 1.0; // 100%!
			mass_fraction[i] = p4[i] * p0[i];
		}

		if (p4[i] < 0.0 || p5[i] < 0.0){
			p4[i] = 0.0;
			p5[i] = 0.0;
			mass_fraction[i] = 0.0;

		}
	}
}
