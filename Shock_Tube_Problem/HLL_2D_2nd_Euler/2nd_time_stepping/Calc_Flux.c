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

void Calc_flux_X(float rho_L, float rho_R, float u_L, float u_R, float v_L, float v_R, float P_L, float P_R, float e_L, float e_R, float a_L, float a_R,
		 float *mass_flux_X, float *momentum_X_flux_X, float *momentum_Y_flux_X, float *energy_flux_X, int INDEX){
	float W_L = fmin(u_L - a_L, u_R - a_R);
	float W_R = fmax(u_L + a_L, u_R + a_R);
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
	if (W_L >= 0.0){
		mass_flux_X[INDEX] = mass_flux_L;
		momentum_X_flux_X[INDEX] = momentum_X_flux_L;
		momentum_Y_flux_X[INDEX] = momentum_Y_flux_L;
		energy_flux_X[INDEX] = energy_flux_L;
	} else if ( W_R <= 0.0){
		mass_flux_X[INDEX] = mass_flux_R;
		momentum_X_flux_X[INDEX] = momentum_X_flux_R;
		momentum_Y_flux_X[INDEX] = momentum_Y_flux_R;
		energy_flux_X[INDEX] = energy_flux_R;
	} else {
		mass_flux_X[INDEX] = (W_R * mass_flux_L - W_L * mass_flux_R + W_L * W_R * (mass_R - mass_L)) / (W_R - W_L);
		momentum_X_flux_X[INDEX] = (W_R * momentum_X_flux_L - W_L * momentum_X_flux_R + W_L * W_R * (momentum_X_R - momentum_X_L)) / (W_R - W_L);
		momentum_Y_flux_X[INDEX] = (W_R * momentum_Y_flux_L - W_L * momentum_Y_flux_R + W_L * W_R * (momentum_Y_R - momentum_Y_L)) / (W_R - W_L);
		energy_flux_X[INDEX] = (W_R * energy_flux_L - W_L * energy_flux_R + W_L * W_R * (energy_R - energy_L)) / (W_R - W_L);
	}
}

void Calc_flux_Y(float rho_B, float rho_T, float u_B, float u_T, float v_B, float v_T, float P_B, float P_T, float e_B, float e_T, float a_B, float a_T,
		 float *mass_flux_Y, float *momentum_X_flux_Y, float *momentum_Y_flux_Y, float *energy_flux_Y, int INDEX){
	float W_B = fmin(v_B - a_B, v_T - a_T);
	float W_T = fmax(v_B + a_B, v_T + a_T);
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
	if (W_B >= 0.0){
		mass_flux_Y[INDEX] = mass_flux_B;
		momentum_X_flux_Y[INDEX] = momentum_X_flux_B;
		momentum_Y_flux_Y[INDEX] = momentum_Y_flux_B;
		energy_flux_Y[INDEX] = energy_flux_B;
	} else if ( W_T <= 0.0){
		mass_flux_Y[INDEX] = mass_flux_T;
		momentum_X_flux_Y[INDEX] = momentum_X_flux_T;
		momentum_Y_flux_Y[INDEX] = momentum_Y_flux_T;
		energy_flux_Y[INDEX] = energy_flux_T;
	} else {
		mass_flux_Y[INDEX] = (W_T * mass_flux_B - W_B * mass_flux_T + W_B * W_T * (mass_T - mass_B)) / (W_T - W_B);
		momentum_X_flux_Y[INDEX] = (W_T * momentum_X_flux_B - W_B * momentum_X_flux_T + W_B * W_T * (momentum_X_T - momentum_X_B)) / (W_T - W_B);
		momentum_Y_flux_Y[INDEX] = (W_T * momentum_Y_flux_B - W_B * momentum_Y_flux_T + W_B * W_T * (momentum_Y_T - momentum_Y_B)) / (W_T - W_B);
		energy_flux_Y[INDEX] = (W_T * energy_flux_B - W_B * energy_flux_T + W_B * W_T * (energy_T - energy_B)) / (W_T - W_B);
	}
}
