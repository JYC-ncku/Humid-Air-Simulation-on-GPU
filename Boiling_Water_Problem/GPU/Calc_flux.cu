#include <stdlib.h>
#include <math.h>
#include "Compute_Cv.h"

__device__ float MAX_Wave_Speed(float u_L, float u_R, float a_L, float a_R){
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
	d_mass_fraction_flux_X[INDEX] -= rho_face_X * D * ((Y_R - Y_L) / dx);
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
				  float *W_GLOBAL_MAX, float R_dry, float R_v, float D, float dx, float dy, int NX, int NY, int N_CELLS){
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int) INDEX / (NY+2);
	int j = (int) INDEX - i * (NY+2);
	int INDEX_L = (i-1) * (NY+2) + j;
	int INDEX_B = i * (NY+2) + (j-1);
	if (INDEX < N_CELLS){
		//X-dir
		if (i >= 1 && i < NX + 2 && j >= 1 && j < NY + 1){
			float rho_L = d_p0[INDEX_L];
			float rho_R = d_p0[INDEX];
			float u_L = d_p1[INDEX_L];
			float u_R = d_p1[INDEX];
			float v_L = d_p2[INDEX_L];
			float v_R = d_p2[INDEX];
			float T_L = d_p3[INDEX_L];
			float T_R = d_p3[INDEX];
			float P_L = d_p4[INDEX_L];
			float P_R = d_p4[INDEX];
			float Y_L = d_p5[INDEX_L];
			float Y_R = d_p5[INDEX];
			float R_mix_L = R_dry * (1 - Y_L) + R_v * Y_L;
			float R_mix_R = R_dry * (1 - Y_R) + R_v * Y_R;
			float Cv_mix_L = Compute_Cv(T_L, Y_L);
			float Cv_mix_R = Compute_Cv(T_R, Y_R);
			float Gamma_L = 1 + R_mix_L / Cv_mix_L;
			float Gamma_R = 1 + R_mix_R / Cv_mix_R;
			float E_L = 0.5 * (u_L * u_L + v_L * v_L) + Cv_mix_L * T_L;
			float E_R = 0.5 * (u_R * u_R + v_R * v_R) + Cv_mix_R * T_R;
			float a_L = sqrt(Gamma_L * R_mix_L * T_L); // Sound speed a = (R*T)^0.5
			float a_R = sqrt(Gamma_R * R_mix_R * T_R);
			float W_LOCAL_MAX = MAX_Wave_Speed(u_L, u_R, a_L, a_R);
			Calc_HLL_X_flux(rho_L, rho_R, u_L, u_R, v_L, v_R, T_L, T_R, P_L, P_R, Y_L, Y_R, E_L, E_R, a_L, a_R,
				        d_mass_flux_X, d_momentum_X_flux_X, d_momentum_Y_flux_X, d_energy_flux_X, d_mass_fraction_flux_X,
				        D, dx, INDEX);
			if (W_LOCAL_MAX > W_GLOBAL_MAX[INDEX]){
				W_GLOBAL_MAX[INDEX] = W_LOCAL_MAX;
			}
		}
		//Y-dir
		if (i >= 1 && i < NX + 1 && j >= 1 && j < NY + 2){
			float rho_B = d_p0[INDEX_B];
			float rho_T = d_p0[INDEX];
			float u_B = d_p1[INDEX_B];
			float u_T = d_p1[INDEX];
			float v_B = d_p2[INDEX_B];
			float v_T = d_p2[INDEX];
			float T_B = d_p3[INDEX_B];
			float T_T = d_p3[INDEX];
			float P_B = d_p4[INDEX_B];
			float P_T = d_p4[INDEX];
			float Y_B = d_p5[INDEX_B];
			float Y_T = d_p5[INDEX];
			float R_mix_B = R_dry * (1 - Y_B) + R_v * Y_B;
			float R_mix_T = R_dry * (1 - Y_T) + R_v * Y_T;
			float Cv_mix_B = Compute_Cv(T_B, Y_B);
			float Cv_mix_T = Compute_Cv(T_T, Y_T);
			float Gamma_B = 1 + R_mix_B / Cv_mix_B;
			float Gamma_T = 1 + R_mix_T / Cv_mix_T;
			float E_B = 0.5 * (u_B * u_B + v_B * v_B) + Cv_mix_B * T_B;
			float E_T = 0.5 * (u_T * u_T + v_T * v_T) + Cv_mix_T * T_T;
			float a_B = sqrt(Gamma_B * R_mix_B * T_B); // Sound speed a = (R*T)^0.5
			float a_T = sqrt(Gamma_T * R_mix_T * T_T);
			float W_LOCAL_MAX = MAX_Wave_Speed(v_B, v_T, a_B, a_T);
			Calc_HLL_Y_flux(rho_B, rho_T, u_B, u_T, v_B, v_T, T_B, T_T, P_B, P_T, Y_B, Y_T, E_B, E_T, a_B, a_T,
				        d_mass_flux_Y, d_momentum_X_flux_Y, d_momentum_Y_flux_Y, d_energy_flux_Y, d_mass_fraction_flux_Y,
				        D, dy, INDEX);
			if (W_LOCAL_MAX > W_GLOBAL_MAX[INDEX]){
				W_GLOBAL_MAX[INDEX] = W_LOCAL_MAX;
			}
		}
	}
}

void Calc_Tot_Flux(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5,
		   float *d_mass_flux_X, float *d_momentum_X_flux_X, float *d_momentum_Y_flux_X, float *d_energy_flux_X, float *d_mass_fraction_flux_X,
		   float *d_mass_flux_Y, float *d_momentum_X_flux_Y, float *d_momentum_Y_flux_Y, float *d_energy_flux_Y, float *d_mass_fraction_flux_Y,
		   float *W_GLOBAL_MAX, float R_dry, float R_v, float D, float dx, float dy, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Calc_Tot_Flux<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5,
					d_mass_flux_X, d_momentum_X_flux_X, d_momentum_Y_flux_X, d_energy_flux_X, d_mass_fraction_flux_X,
					d_mass_flux_Y, d_momentum_X_flux_Y, d_momentum_Y_flux_Y, d_energy_flux_Y, d_mass_fraction_flux_Y,
					W_GLOBAL_MAX, R_dry, R_v, D, dx, dy, NX, NY, N_CELLS);
}
