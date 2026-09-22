#include <stdlib.h>
#include <math.h>
void Calc_HLL_flux(float rho_L, float rho_R, float u_L, float u_R, float T_L, float T_R, float P_L, float P_R, float Y_L, float Y_R, float E_L, float E_R, float a_L, float a_R,
		   float *mass_flux, float *momentum_X_flux, float *momentum_Y_flux, float *energy_flux, float *mass_fraction_flux, int INDEX){
	float W_L = fmin(u_L - a_L, u_R - a_R);
	float W_R = fmax(u_L + a_L, u_R + a_R);
	float mass_L = rho_L;
	float momentum_X_L = rho_L * u_L;
	float energy_L = rho_L * E_L;
	float mass_R = rho_R;
	float momentum_X_R = rho_R * u_R;
	float energy_R = rho_R * E_R;
	float mass_fraction_L = mass_L * Y_L;
	float mass_fraction_R = mass_R * Y_R;

	float mass_flux_L = rho_L * u_L;
	float mass_flux_R = rho_R * u_R;
	float momentum_X_flux_L = rho_L * u_L * u_L + P_L;
	float momentum_X_flux_R = rho_R * u_R * u_R + P_R;
	float energy_flux_L = (energy_L + P_L) * u_L;
	float energy_flux_R = (energy_R + P_R) * u_R;
	float mass_fraction_flux_L = mass_flux_L * Y_L;
	float mass_fraction_flux_R = mass_flux_R * Y_R;
	if (W_L >= 0.0){
		mass_flux[INDEX] = mass_flux_L;
		momentum_X_flux[INDEX] = momentum_X_flux_L;
		energy_flux[INDEX] = energy_flux_L;
		mass_fraction_flux[INDEX] = mass_fraction_flux_L;
	} else if ( W_R <= 0.0){
		mass_flux[INDEX] = mass_flux_R;
		momentum_X_flux[INDEX] = momentum_X_flux_R;
		energy_flux[INDEX] = energy_flux_R;
		mass_fraction_flux[INDEX] = mass_fraction_flux_R;
	} else {
		mass_flux[INDEX] = (W_R * mass_flux_L - W_L * mass_flux_R + W_L * W_R * (mass_R - mass_L)) / (W_R - W_L);
		momentum_X_flux[INDEX] = (W_R * momentum_X_flux_L - W_L * momentum_X_flux_R + W_L * W_R * (momentum_X_R - momentum_X_L)) / (W_R - W_L);
		energy_flux[INDEX] = (W_R * energy_flux_L - W_L * energy_flux_R + W_L * W_R * (energy_R - energy_L)) / (W_R - W_L);
		mass_fraction_flux[INDEX] = (W_R * mass_fraction_flux_L - W_L * mass_fraction_flux_R + W_L * W_R * (mass_fraction_R - mass_fraction_L)) / (W_R - W_L);
	}
}
