#include <stdlib.h>
#include <math.h>

float MAX_WAVE_SPEED(float u_L, float u_R, float a_L, float a_R){
	float W_L = fabs(u_L) + a_L;
	float W_R = fabs(u_R) + a_R;
	float W_LOCAL_MAX;
	if (W_L > W_R){
		W_LOCAL_MAX = W_L;
	}else {
		W_LOCAL_MAX = W_R;
	}
	return W_LOCAL_MAX;
}

void Calc_flux(float rho_L, float rho_R, float u_L, float u_R, float v_L, float v_R, float T_L, float T_R, float P_L, float P_R, float e_L, float e_R, float W_LOCAL_MAX,
	       float *mass_flux, float *momentum_X_flux, float *momentum_Y_flux, float *energy_flux, int INDEX){
	float mass_L = rho_L;
	float mass_R = rho_R;
	float momentum_X_L = rho_L * u_L;
	float momentum_X_R = rho_R * u_R;
	float momentum_Y_L = rho_L * v_L;
	float momentum_Y_R = rho_R * v_R;
	float energy_L = e_L;
	float energy_R = e_R;

	float mass_flux_L = (rho_L * u_L) + (rho_L * v_L);
	float mass_flux_R = (rho_R * u_R) + (rho_R * v_R);
	float momentum_X_flux_L = (rho_L * u_L * u_L + P_L) + (rho_L * u_L * v_L);
	float momentum_X_flux_R = (rho_R * u_R * u_R + P_R) + (rho_R * u_R * v_R);
	float momentum_Y_flux_L = (rho_L * u_L * v_L) + (rho_L * v_L * v_L + P_R);
	float momentum_Y_flux_R = (rho_L * u_L * v_L) + (rho_R * v_R * v_R + P_R);
	float energy_flux_L = u_L * (e_L + P_L) + v_L * (e_L + P_L);
	float energy_flux_R = u_R * (e_R + P_R) + v_R * (e_R + P_R);

	mass_flux[INDEX] = 0.5 * (mass_flux_L + mass_flux_R) - 0.5 * W_LOCAL_MAX * (mass_R - mass_L);
	momentum_X_flux[INDEX] = 0.5 * (momentum_X_flux_L + momentum_X_flux_R) - 0.5 * W_LOCAL_MAX * (momentum_X_R - momentum_X_L);
	momentum_Y_flux[INDEX] = 0.5 * (momentum_Y_flux_L + momentum_Y_flux_R) - 0.5 * W_LOCAL_MAX * (momentum_Y_R - momentum_Y_L);
	energy_flux[INDEX] = 0.5 * (energy_flux_L + energy_flux_R) - 0.5 * W_LOCAL_MAX * (energy_R - energy_L);
}
