#include <stdlib.h>
#include <math.h>
void Calc_HLL_X_flux(float rho_L, float rho_R, float u_L, float u_R, float v_L, float v_R, float T_L, float T_R, float P_L, float P_R, float Y_L, float Y_R,
		     float E_L, float E_R, float a_L, float a_R,
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
	float momentum_Y_flux_L = rho_L * u_L * v_L;
	float momentum_Y_flux_R = rho_R * u_R * v_R;
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

void Calc_HLL_Y_flux(float rho_B, float rho_T, float u_B, float u_T, float v_B, float v_T, float T_B, float T_T, float P_B, float P_T, float Y_B, float Y_T,
		     float E_B, float E_T, float a_B, float a_T,
		     float *mass_flux, float *momentum_X_flux, float *momentum_Y_flux, float *energy_flux, float *mass_fraction_flux, int INDEX){
	float W_B = fmin(v_B - a_B, v_T - a_T);
	float W_T = fmax(v_B + a_B, v_T + a_T);
	float mass_B = rho_B;
	float momentum_X_B = rho_B * u_B;
	float energy_B = rho_B * E_B;
	float mass_T = rho_T;
	float momentum_X_T = rho_T * u_T;
	float energy_T = rho_T * E_T;
	float mass_fraction_B = mass_B * Y_B;
	float mass_fraction_T = mass_T * Y_T;

	float mass_flux_B = rho_B * u_B;
	float mass_flux_T = rho_T * u_T;
	float momentum_X_flux_B = rho_B * u_B * u_B + P_B;
	float momentum_X_flux_T = rho_T * u_T * u_T + P_T;
	float momentum_Y_flux_B = rho_B * u_B * v_B;
	float momentum_Y_flux_T = rho_T * u_T * v_T;
	float energy_flux_B = (energy_B + P_B) * u_B;
	float energy_flux_T = (energy_T + P_T) * u_T;
	float mass_fraction_flux_B = mass_flux_B * Y_B;
	float mass_fraction_flux_T = mass_flux_T * Y_T;
	if (W_B >= 0.0){
		mass_flux[INDEX] = mass_flux_B;
		momentum_X_flux[INDEX] = momentum_X_flux_B;
		energy_flux[INDEX] = energy_flux_B;
		mass_fraction_flux[INDEX] = mass_fraction_flux_B;
	} else if ( W_T <= 0.0){
		mass_flux[INDEX] = mass_flux_T;
		momentum_X_flux[INDEX] = momentum_X_flux_T;
		energy_flux[INDEX] = energy_flux_T;
		mass_fraction_flux[INDEX] = mass_fraction_flux_T;
	} else {
		mass_flux[INDEX] = (W_T * mass_flux_B - W_B * mass_flux_T + W_B * W_T * (mass_T - mass_B)) / (W_T - W_B);
		momentum_X_flux[INDEX] = (W_T * momentum_X_flux_B - W_B * momentum_X_flux_T + W_B * W_T * (momentum_X_T - momentum_X_B)) / (W_T - W_B);
		energy_flux[INDEX] = (W_T * energy_flux_B - W_B * energy_flux_T + W_B * W_T * (energy_T - energy_B)) / (W_T - W_B);
		mass_fraction_flux[INDEX] = (W_T * mass_fraction_flux_B - W_B * mass_fraction_flux_T + W_B * W_T * (mass_fraction_T - mass_fraction_B)) / (W_T - W_B);
	}
}
