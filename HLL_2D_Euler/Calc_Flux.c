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

void Calc_flux_X(float rho_L, float rho_R, float u_L, float u_R, float v_L, float v_R, float T_L, float T_R, float P_L, float P_R, float e_L, float e_R, float W_LOCAL_MAX,
		 float *mass_flux_X, float *momentum_X_flux_X, float *momentum_Y_flux_X, float *energy_flux_X, int INDEX){
	float mass_L = rho_L;
	float mass_R = rho_R;
	float momentum_X_L = rho_L * u_L;
	float momentum_X_R = rho_R * u_R;
	float momentum_Y_L = rho_L * v_L;
	float momentum_Y_R = rho_R * v_R;
	float energy_L = e_L;
	float energy_R = e_R;

	float mass_flux_L = (rho_L * u_L);
	float mass_flux_R = (rho_R * u_R);
	float momentum_X_flux_L = (rho_L * u_L * u_L + P_L);
	float momentum_X_flux_R = (rho_R * u_R * u_R + P_R);
	float momentum_Y_flux_L = (rho_L * u_L * v_L);
	float momentum_Y_flux_R = (rho_R * u_R * v_R);
	float energy_flux_L = u_L * (e_L + P_L);
	float energy_flux_R = u_R * (e_R + P_R);

	mass_flux_X[INDEX] = 0.5 * (mass_flux_L + mass_flux_R) - 0.5 * W_LOCAL_MAX * (mass_R - mass_L);
	momentum_X_flux_X[INDEX] = 0.5 * (momentum_X_flux_L + momentum_X_flux_R) - 0.5 * W_LOCAL_MAX * (momentum_X_R - momentum_X_L);
	momentum_Y_flux_X[INDEX] = 0.5 * (momentum_Y_flux_L + momentum_Y_flux_R) - 0.5 * W_LOCAL_MAX * (momentum_Y_R - momentum_Y_L);
	energy_flux_X[INDEX] = 0.5 * (energy_flux_L + energy_flux_R) - 0.5 * W_LOCAL_MAX * (energy_R - energy_L);
}

void Calc_flux_Y(float rho_B, float rho_T, float u_B, float u_T, float v_B, float v_T, float T_B, float T_T, float P_B, float P_T, float e_B, float e_T, float W_LOCAL_MAX,
		 float *mass_flux_Y, float *momentum_X_flux_Y, float *momentum_Y_flux_Y, float *energy_flux_Y, int INDEX){
	float mass_B = rho_B;
	float mass_T = rho_T;
	float momentum_X_B = rho_B * u_B;
	float momentum_X_T = rho_T * u_T;
	float momentum_Y_B = rho_B * v_B;
	float momentum_Y_T = rho_T * v_T;
	float energy_B = e_B;
	float energy_T = e_T;

	float mass_flux_B = (rho_B * v_B);
	float mass_flux_T = (rho_T * v_T);
	float momentum_X_flux_B = rho_B * u_B * v_B;
	float momentum_X_flux_T = rho_T * u_T * v_T;
	float momentum_Y_flux_B = (rho_B * v_B * v_B) + P_B;
	float momentum_Y_flux_T = (rho_T * v_T * v_T) + P_T;
	float energy_flux_B = v_B * (e_B + P_B);
	float energy_flux_T = v_T * (e_T + P_T);

	mass_flux_Y[INDEX] = 0.5 * (mass_flux_B + mass_flux_T) - 0.5 * W_LOCAL_MAX * (mass_T - mass_B);
	momentum_X_flux_Y[INDEX] = 0.5 * (momentum_X_flux_B + momentum_X_flux_T) - 0.5 * W_LOCAL_MAX * (momentum_X_T - momentum_X_B);
	momentum_Y_flux_Y[INDEX] = 0.5 * (momentum_Y_flux_B + momentum_Y_flux_T) - 0.5 * W_LOCAL_MAX * (momentum_Y_T - momentum_Y_B);
	energy_flux_Y[INDEX] = 0.5 * (energy_flux_B + energy_flux_T) - 0.5 * W_LOCAL_MAX * (energy_T - energy_B);
}
