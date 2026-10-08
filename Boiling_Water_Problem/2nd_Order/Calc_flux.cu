#include <stdlib.h>
#include <math.h>
#include "Compute_Cv.h"

__device__ float MINMOD(float U_L, float U_C, float U_R, float dx){
	float dU_dx;
	float Forward = (U_R - U_C) / dx;
	float Backward = (U_C - U_L) / dx;
	if (Backward * Forward < 0){
		dU_dx = 0;
	} else if ( fabs(Forward) < fabs(Backward) ){
			dU_dx = Forward;
		} else {
			dU_dx = Backward;
	}
	return dU_dx;
}

__device__ void Calc_HLL_X_flux(float rho_L, float rho_R, float u_L, float u_R, float v_L, float v_R, float T_L, float T_R, float P_L, float P_R, float Y_L, float Y_R,
				float E_L, float E_R, float a_L, float a_R,
				float *d_mass_flux_X, float *d_momentum_X_flux_X, float *d_momentum_Y_flux_X, float *d_energy_flux_X, float *d_mass_fraction_flux_X,
				float D, float dx, int INDEX){
	float W_L = fmin(u_L - a_L, u_R - a_R);
	float W_R = fmax(u_L + a_L, u_R + a_R);
	float mass_L = rho_L;
	float mass_R = rho_R;
	float momentum_X_L = rho_L * u_L;
	float momentum_X_R = rho_R * u_R;
	float momentum_Y_L = rho_L * v_L;
	float momentum_Y_R = rho_R * v_R;
	float energy_L = rho_L * E_L;
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
		d_mass_flux_X[INDEX] = mass_flux_L;
		d_momentum_X_flux_X[INDEX] = momentum_X_flux_L;
		d_momentum_Y_flux_X[INDEX] = momentum_Y_flux_L;
		d_energy_flux_X[INDEX] = energy_flux_L;
		d_mass_fraction_flux_X[INDEX] = mass_fraction_flux_L;
	} else if ( W_R <= 0.0){
		d_mass_flux_X[INDEX] = mass_flux_R;
		d_momentum_X_flux_X[INDEX] = momentum_X_flux_R;
		d_momentum_Y_flux_X[INDEX] = momentum_Y_flux_R;
		d_energy_flux_X[INDEX] = energy_flux_R;
		d_mass_fraction_flux_X[INDEX] = mass_fraction_flux_R;
	} else {
		d_mass_flux_X[INDEX] = (W_R * mass_flux_L - W_L * mass_flux_R + W_L * W_R * (mass_R - mass_L)) / (W_R - W_L);
		d_momentum_X_flux_X[INDEX] = (W_R * momentum_X_flux_L - W_L * momentum_X_flux_R + W_L * W_R * (momentum_X_R - momentum_X_L)) / (W_R - W_L);
		d_momentum_Y_flux_X[INDEX] = (W_R * momentum_Y_flux_L - W_L * momentum_Y_flux_R + W_L * W_R * (momentum_Y_R - momentum_Y_L)) / (W_R - W_L);
		d_energy_flux_X[INDEX] = (W_R * energy_flux_L - W_L * energy_flux_R + W_L * W_R * (energy_R - energy_L)) / (W_R - W_L);
		d_mass_fraction_flux_X[INDEX] = (W_R * mass_fraction_flux_L - W_L * mass_fraction_flux_R + W_L * W_R * (mass_fraction_R - mass_fraction_L)) / (W_R - W_L);
	}

	float rho_face_X = 0.5 * (rho_L + rho_R);
	d_mass_fraction_flux_X[INDEX] -= rho_face_X * D * ((Y_R - Y_L) / dx); //Central difference
}

__device__ void Calc_HLL_Y_flux(float rho_B, float rho_T, float u_B, float u_T, float v_B, float v_T, float T_B, float T_T, float P_B, float P_T, float Y_B, float Y_T,
				float E_B, float E_T, float a_B, float a_T,
				float *d_mass_flux_Y, float *d_momentum_X_flux_Y, float *d_momentum_Y_flux_Y, float *d_energy_flux_Y, float *d_mass_fraction_flux_Y,
				float D, float dy, int INDEX){
	float W_B = fmin(v_B - a_B, v_T - a_T);
	float W_T = fmax(v_B + a_B, v_T + a_T);
	float mass_B = rho_B;
	float mass_T = rho_T;
	float momentum_X_B = rho_B * u_B;
	float momentum_X_T = rho_T * u_T;
	float momentum_Y_B = rho_B * v_B;
	float momentum_Y_T = rho_T * v_T;
	float energy_B = rho_B * E_B;
	float energy_T = rho_T * E_T;
	float mass_fraction_B = mass_B * Y_B;
	float mass_fraction_T = mass_T * Y_T;

	float mass_flux_B = rho_B * v_B;
	float mass_flux_T = rho_T * v_T;
	float momentum_X_flux_B = rho_B * u_B * v_B;
	float momentum_X_flux_T = rho_T * u_T * v_T;
	float momentum_Y_flux_B = rho_B * v_B * v_B + P_B;
	float momentum_Y_flux_T = rho_T * v_T * v_T + P_T;
	float energy_flux_B = (energy_B + P_B) * v_B;
	float energy_flux_T = (energy_T + P_T) * v_T;
	float mass_fraction_flux_B = mass_flux_B * Y_B;
	float mass_fraction_flux_T = mass_flux_T * Y_T;

	if (W_B >= 0.0){
		d_mass_flux_Y[INDEX] = mass_flux_B;
		d_momentum_X_flux_Y[INDEX] = momentum_X_flux_B;
		d_momentum_Y_flux_Y[INDEX] = momentum_Y_flux_B;
		d_energy_flux_Y[INDEX] = energy_flux_B;
		d_mass_fraction_flux_Y[INDEX] = mass_fraction_flux_B;
	} else if ( W_T <= 0.0){
		d_mass_flux_Y[INDEX] = mass_flux_T;
		d_momentum_X_flux_Y[INDEX] = momentum_X_flux_T;
		d_momentum_Y_flux_Y[INDEX] = momentum_Y_flux_T;
		d_energy_flux_Y[INDEX] = energy_flux_T;
		d_mass_fraction_flux_Y[INDEX] = mass_fraction_flux_T;
	} else {
		d_mass_flux_Y[INDEX] = (W_T * mass_flux_B - W_B * mass_flux_T + W_B * W_T * (mass_T - mass_B)) / (W_T - W_B);
		d_momentum_X_flux_Y[INDEX] = (W_T * momentum_X_flux_B - W_B * momentum_X_flux_T + W_B * W_T * (momentum_X_T - momentum_X_B)) / (W_T - W_B);
		d_momentum_Y_flux_Y[INDEX] = (W_T * momentum_Y_flux_B - W_B * momentum_Y_flux_T + W_B * W_T * (momentum_Y_T - momentum_Y_B)) / (W_T - W_B);
		d_energy_flux_Y[INDEX] = (W_T * energy_flux_B - W_B * energy_flux_T + W_B * W_T * (energy_T - energy_B)) / (W_T - W_B);
		d_mass_fraction_flux_Y[INDEX] = (W_T * mass_fraction_flux_B - W_B * mass_fraction_flux_T + W_B * W_T * (mass_fraction_T - mass_fraction_B)) / (W_T - W_B);
	}

	float rho_face_Y = 0.5 * (rho_B + rho_T);
	d_mass_fraction_flux_Y[INDEX] -= rho_face_Y * D * ((Y_T - Y_B) / dy);
}


__global__ void GPU_Calc_Tot_Flux(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5,
				  float *d_mass_flux_X, float *d_momentum_X_flux_X, float *d_momentum_Y_flux_X, float *d_energy_flux_X, float *d_mass_fraction_flux_X,
				  float *d_mass_flux_Y, float *d_momentum_X_flux_Y, float *d_momentum_Y_flux_Y, float *d_energy_flux_Y, float *d_mass_fraction_flux_Y,
				  float R_dry, float R_v, float D, float dx, float dy, int NX, int NY, int N_CELLS){
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int) INDEX / (NY+4);
	int j = (int) INDEX - i * (NY+4);
	int INDEX_R = (i + 1) * (NY + 4) + j;
	int INDEX_RR = (i + 2) * (NY + 4) + j;
	int INDEX_L = (i - 1) * (NY + 4) + j;
	int INDEX_B = i * (NY + 4) + (j - 1);
	int INDEX_T = i * (NY + 4) + (j + 1);
	int INDEX_TT = i * (NY + 4) + (j + 2);
	if (INDEX < N_CELLS){
		//X-dir
		if (i >= 1 && i < NX + 2 && j >= 2 && j < NY + 2){
			float rho_L = d_p0[INDEX_L];
			float rho_C = d_p0[INDEX];
			float rho_R = d_p0[INDEX_R];
			float rho_RR = d_p0[INDEX_RR];
			float drho_dx_L = MINMOD(rho_L, rho_C, rho_R, dx);
			float drho_dx_R = MINMOD(rho_C, rho_R, rho_RR, dx);
			float rho_L_star = rho_C + 0.5 * dx * drho_dx_L;
			float rho_R_star = rho_R - 0.5 * dx * drho_dx_R;

			float u_L = d_p1[INDEX_L];
			float u_C = d_p1[INDEX];
			float u_R = d_p1[INDEX_R];
			float u_RR = d_p1[INDEX_RR];
			float du_dx_L = MINMOD(u_L, u_C, u_R, dx);
			float du_dx_R = MINMOD(u_C, u_R, u_RR, dx);
			float u_L_star = u_C + 0.5 * dx * du_dx_L;
			float u_R_star = u_R - 0.5 * dx * du_dx_R;

			float v_L = d_p2[INDEX_L];
			float v_C = d_p2[INDEX];
			float v_R = d_p2[INDEX_R];
			float v_RR = d_p2[INDEX_RR];
			float dv_dx_L = MINMOD(v_L, v_C, v_R, dx);
			float dv_dx_R = MINMOD(v_C, v_R, v_RR, dx);
			float v_L_star = v_C + 0.5 * dx * dv_dx_L;
			float v_R_star = v_R - 0.5 * dx * dv_dx_R;

			float T_L = d_p3[INDEX_L];
			float T_C = d_p3[INDEX];
			float T_R = d_p3[INDEX_R];
			float T_RR = d_p3[INDEX_RR];
			float dT_dx_L = MINMOD(T_L, T_C, T_R, dx);
			float dT_dx_R = MINMOD(T_C, T_R, T_RR, dx);
			float T_L_star = T_C + 0.5 * dx * dT_dx_L;
			float T_R_star = T_R - 0.5 * dx * dT_dx_R;

			float Y_L = d_p5[INDEX_L];
			float Y_C = d_p5[INDEX];
			float Y_R = d_p5[INDEX_R];
			float Y_RR = d_p5[INDEX_RR];
			float dY_dx_L = MINMOD(Y_L, Y_C, Y_R, dx);
			float dY_dx_R = MINMOD(Y_C, Y_R, Y_RR, dx);
			float Y_L_star = Y_C + 0.5 * dx * dY_dx_L;
			float Y_R_star = Y_R - 0.5 * dx * dY_dx_R;

			float R_mix_L = R_dry * (1 - Y_L) + R_v * Y_L;
			float R_mix_R = R_dry * (1 - Y_R) + R_v * Y_R;
			float Cv_mix_L = Compute_Cv(T_L, Y_L);
			float Cv_mix_R = Compute_Cv(T_R, Y_R);
			float Gamma_L = 1 + R_mix_L / Cv_mix_L;
			float Gamma_R = 1 + R_mix_R / Cv_mix_R;

			float P_L_star = rho_L_star * R_mix_L * T_L_star;
			float P_R_star = rho_R_star * R_mix_R * T_R_star;

			float E_L_star = 0.5 * (u_L_star * u_L_star + v_L_star * v_L_star) + Cv_mix_L * T_L_star;
			float E_R_star = 0.5 * (u_R_star * u_R_star + v_R_star * v_R_star) + Cv_mix_R * T_R_star;
			float a_L_star = sqrt(Gamma_L * R_mix_L * T_L_star); // Sound speed a = (R*T)^0.5
			float a_R_star = sqrt(Gamma_R * R_mix_R * T_R_star);
			Calc_HLL_X_flux(rho_L_star, rho_R_star, u_L_star, u_R_star, v_L_star, v_R_star, T_L_star, T_R_star, P_L_star, P_R_star, Y_L_star, Y_R_star, E_L_star, E_R_star, a_L_star, a_R_star,
				        d_mass_flux_X, d_momentum_X_flux_X, d_momentum_Y_flux_X, d_energy_flux_X, d_mass_fraction_flux_X,
				        D, dx, INDEX);
		}
		//Y-dir
		if (i >= 2 && i < NX + 2 && j >= 1 && j < NY + 2){
			float rho_B = d_p0[INDEX_B];
			float rho_C = d_p0[INDEX];
			float rho_T = d_p0[INDEX_T];
			float rho_TT = d_p0[INDEX_TT];
			float drho_dx_B = MINMOD(rho_B, rho_C, rho_T, dx);
			float drho_dx_T = MINMOD(rho_C, rho_T, rho_TT, dx);
			float rho_B_star = rho_C + 0.5 * dy * drho_dx_B;
			float rho_T_star = rho_T - 0.5 * dy * drho_dx_T;

			float u_B = d_p1[INDEX_B];
			float u_C = d_p1[INDEX];
			float u_T = d_p1[INDEX_T];
			float u_TT = d_p1[INDEX_TT];
			float du_dx_B = MINMOD(u_B, u_C, u_T, dx);
			float du_dx_T = MINMOD(u_C, u_T, u_TT, dx);
			float u_B_star = u_C + 0.5 * dy * du_dx_B;
			float u_T_star = u_T - 0.5 * dy * du_dx_T;

			float v_B = d_p2[INDEX_B];
			float v_C = d_p2[INDEX];
			float v_T = d_p2[INDEX_T];
			float v_TT = d_p2[INDEX_TT];
			float dv_dx_B = MINMOD(v_B, v_C, v_T, dx);
			float dv_dx_T = MINMOD(v_C, v_T, v_TT, dx);
			float v_B_star = v_C + 0.5 * dy * dv_dx_B;
			float v_T_star = v_T - 0.5 * dy * dv_dx_T;

			float T_B = d_p3[INDEX_B];
			float T_C = d_p3[INDEX];
			float T_T = d_p3[INDEX_T];
			float T_TT = d_p3[INDEX_TT];
			float dT_dx_B = MINMOD(T_B, T_C, T_T, dx);
			float dT_dx_T = MINMOD(T_C, T_T, T_TT, dx);
			float T_B_star = T_C + 0.5 * dy * dT_dx_B;
			float T_T_star = T_T - 0.5 * dy * dT_dx_T;

			float Y_B = d_p5[INDEX_B];
			float Y_C = d_p5[INDEX];
			float Y_T = d_p3[INDEX_T];
			float Y_TT = d_p3[INDEX_TT];
			float dY_dx_B = MINMOD(Y_B, Y_C, Y_T, dx);
			float dY_dx_T = MINMOD(Y_C, Y_T, Y_TT, dx);
			float Y_B_star = Y_C + 0.5 * dy * dY_dx_B;
			float Y_T_star = Y_T - 0.5 * dy * dY_dx_T;

			float R_mix_B = R_dry * (1 - Y_B) + R_v * Y_B;
			float R_mix_T = R_dry * (1 - Y_T) + R_v * Y_T;
			float Cv_mix_B = Compute_Cv(T_B, Y_B);
			float Cv_mix_T = Compute_Cv(T_T, Y_T);
			float Gamma_B = 1 + R_mix_B / Cv_mix_B;
			float Gamma_T = 1 + R_mix_T / Cv_mix_T;

			float P_B_star = rho_B_star * R_mix_B * T_B_star;
			float P_T_star = rho_T_star * R_mix_T * T_T_star;

			float E_B_star = 0.5 * (u_B_star * u_B_star + v_B_star * v_B_star) + Cv_mix_B * T_B_star;
			float E_T_star = 0.5 * (u_T_star * u_T_star + v_T_star * v_T_star) + Cv_mix_T * T_T_star;
			float a_B_star = sqrt(Gamma_B * R_mix_B * T_B_star); // Sound speed a = (R*T)^0.5
			float a_T_star = sqrt(Gamma_T * R_mix_T * T_T_star);
			Calc_HLL_Y_flux(rho_B_star, rho_T_star, u_B_star, u_T_star, v_B_star, v_T_star, T_B_star, T_T_star, P_B_star, P_T_star, Y_B_star, Y_T_star, E_B_star, E_T_star, a_B_star, a_T_star,
				        d_mass_flux_Y, d_momentum_X_flux_Y, d_momentum_Y_flux_Y, d_energy_flux_Y, d_mass_fraction_flux_Y,
				        D, dy, INDEX);
		}
	}
}

void Calc_Tot_Flux(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5,
		   float *d_mass_flux_X, float *d_momentum_X_flux_X, float *d_momentum_Y_flux_X, float *d_energy_flux_X, float *d_mass_fraction_flux_X,
		   float *d_mass_flux_Y, float *d_momentum_X_flux_Y, float *d_momentum_Y_flux_Y, float *d_energy_flux_Y, float *d_mass_fraction_flux_Y,
		   float R_dry, float R_v, float D, float dx, float dy, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Calc_Tot_Flux<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5,
					d_mass_flux_X, d_momentum_X_flux_X, d_momentum_Y_flux_X, d_energy_flux_X, d_mass_fraction_flux_X,
					d_mass_flux_Y, d_momentum_X_flux_Y, d_momentum_Y_flux_Y, d_energy_flux_Y, d_mass_fraction_flux_Y,
					R_dry, R_v, D, dx, dy, NX, NY, N_CELLS);
}
