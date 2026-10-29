#include <stdlib.h>
#include <math.h>
void Calc_HLL_flux(float rho_L, float rho_R, float u_L, float u_R, float T_L, float T_R, float P_L, float P_R, float e_L, float e_R, float a_L, float a_R,
		   float *mass_flux, float *momentum_flux, float *energy_flux, int i){
	float W_L = fmin(u_L - a_L, u_R - a_R);
	float W_R = fmax(u_L + a_L, u_R + a_R);
	float mass_L = rho_L;
	float momentum_L = rho_L * u_L;
	float energy_L = e_L;
	float mass_R = rho_R;
	float momentum_R = rho_R * u_R;
	float energy_R = e_R;

	float mass_flux_L = rho_L * u_L;
	float mass_flux_R = rho_R * u_R;
	float momentum_flux_L = rho_L * u_L * u_L + P_L;
	float momentum_flux_R = rho_R * u_R * u_R + P_R;
	float energy_flux_L = (e_L + P_L) * u_L;
	float energy_flux_R = (e_R + P_R) * u_R;
	if (W_L >= 0.0){
		mass_flux[i] = mass_flux_L;
		momentum_flux[i] = momentum_flux_L;
		energy_flux[i] = energy_flux_L;
	} else if ( W_R <= 0.0){
		mass_flux[i] = mass_flux_R;
		momentum_flux[i] = momentum_flux_R;
		energy_flux[i] = energy_flux_R;
	} else {
		mass_flux[i] = (W_R * mass_flux_L - W_L * mass_flux_R + W_L * W_R * (mass_R - mass_L)) / (W_R - W_L);
		momentum_flux[i] = (W_R * momentum_flux_L - W_L * momentum_flux_R + W_L * W_R * (momentum_R - momentum_L)) / (W_R - W_L);
		energy_flux[i] = (W_R * energy_flux_L - W_L * energy_flux_R + W_L * W_R * (energy_R - energy_L)) / (W_R - W_L);
	}
}
