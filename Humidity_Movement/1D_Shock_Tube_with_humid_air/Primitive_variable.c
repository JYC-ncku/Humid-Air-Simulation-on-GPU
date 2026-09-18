#include <stdlib.h>
#include <math.h>
#include "compute_T.h"

void Calc_primitive_variable(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *P_sat, float *P_v, float *mass, float *momentum, float *energy, float *mass_fraction,
			     float *mass_flux, float *momentum_flux, float *energy_flux, float *mass_fraction_flux, float *R_mix, float *Cv_mix, float R_v, float dx, float dt,
			     int N_CELLS){
	for (int i = 1; i < N_CELLS + 1; i++){
	        // Use FVM to get new conservation values
		mass[i] = mass[i] - (dt / dx) * (mass_flux[i+1] - mass_flux[i]);
		momentum[i] = momentum[i] - (dt / dx) * (momentum_flux[i+1] - momentum_flux[i]);
		energy[i] = energy[i] - (dt / dx) * (energy_flux[i+1] - energy_flux[i]);
		mass_fraction[i] = mass_fraction[i] - (dt / dx) * (mass_fraction_flux[i+1] - mass_fraction_flux[i]);
		//Get new variable
		p0[i] = mass[i];
		p1[i] = momentum[i] / mass[i];
		float e_target = (energy[i] / p0[i]) - 0.5 * p1[i] * p1[i];
		float T_old = p2[i];
		float phi_old = p4[i];
		p2[i] = compute_T(T_old, phi_old, e_target);
		p3[i] = p0[i] * R_mix[i] * p2[i];
		p4[i] = mass_fraction[i] / p0[i];

		P_sat[i] = 0.611 * exp((17.27 * (p2[i] - 273.15)) / ((p2[i] - 273.15) + 237.3)); // Tetens equation
		P_v[i] = (p0[i] * p4[i]) * R_v * p2[i];
		p5[i] = P_v[i] / P_sat[i];
	}
}
