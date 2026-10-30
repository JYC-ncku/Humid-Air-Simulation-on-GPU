#include <stdlib.h>
#include <math.h>
#include <stdio.h>
//#include <thrust/extrema.h>
//#include <thrust/execution_policy.h>
#include "compute_T.h"
#include "Compute_Cv.h"

// Return the maximum CFL number across all cells
void Compute_MAX_CFL(float *MAX_Freq, float *p0, float *p1, float *p2, float *p3, float *p5, float dx, float dy, float R_dry, float R_v, int NX, int NY, int N_CELLS){
/*
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+4);
	int j = (int)INDEX - i * (NY+4);
*/
	float R_mix, Cv_mix, Gamma_mix;
//	if (INDEX < N_CELLS){
	// Only care inner cells.
//		if (i >= 2 && i < NX + 2 && j >= 2 && j < NY + 2){
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY+2) + j;
			float u = p1[INDEX];
			float v = p2[INDEX];
			float T = p3[INDEX];
			float Y = p5[INDEX];
			R_mix = R_dry * (1.0 - p5[INDEX]) + R_v * p5[INDEX];
			Cv_mix = Compute_Cv(T, Y);
			Gamma_mix = 1.0 + R_mix / Cv_mix;
			if (T<0){
				printf("Error: Negative temperature in cell %d: T = %f\n. Aborting.", INDEX, T);
			}
			float a = sqrt(Gamma_mix * R_mix * T); // GAMMA = 1.4, R = 1.0
			float MAX_WAVE_SPEED_X = (fabs(u) + a) / dx;
			float MAX_WAVE_SPEED_Y = (fabs(v) + a) / dy;
			MAX_Freq[INDEX] = MAX_WAVE_SPEED_X + MAX_WAVE_SPEED_Y;
		}
	}
}


void Calc_primitive_variable(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6,
			     float *mass, float *momentum_X, float *momentum_Y, float *energy, float *mass_fraction,
			     float *flux_X, float *flux_Y,
			     float R_dry, float R_v, float dx, float dy, float dt, int NX, int NY, int N_CELLS){
/*
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+4);
	int j = (int)INDEX - i * (NY+4);
*/
	float e_target, T_old, phi_new, phi_max,  R_mix, P_sat, P_v;
//	if (INDEX < N_CELLS){
//		if (i >= 2 && i < NX + 2 && j >= 2 && j < NY + 2){
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){

			int INDEX = i * (NY+2) + j;

			// 我們是 i*5，所以左界面是 (i-1)*5，右界面是 i*5
			int L_interface = ((i - 1) * (NY + 2) + j) * 6;
			int R_interface = INDEX * 6; // T_intewrface = R_interface
			int B_interface = (i * (NY + 2) + (j - 1)) * 6;
			int T_interface = INDEX * 6; // T_intewrface = R_interface

		        // Use FVM to get new conservation values
			mass[INDEX] = mass[INDEX] - (dt / dx) * (flux_X[R_interface + 0] - flux_X[L_interface + 0])
						      - (dt / dy) * (flux_Y[T_interface + 0] - flux_Y[B_interface + 0]);
			momentum_X[INDEX] = momentum_X[INDEX] - (dt / dx) * (flux_X[R_interface + 1] - flux_X[L_interface + 1])
								  - (dt / dy) * (flux_Y[T_interface + 1] - flux_Y[B_interface + 1]);
			momentum_Y[INDEX] = momentum_Y[INDEX] - (dt / dx) * (flux_X[R_interface + 2] - flux_X[L_interface + 2])
								  - (dt / dy) * (flux_Y[T_interface + 2] - flux_Y[B_interface + 2]);
			energy[INDEX] = energy[INDEX] - (dt / dx) * (flux_X[R_interface + 4] - flux_X[L_interface + 4])
							  - (dt / dy) * (flux_Y[T_interface + 4] - flux_Y[B_interface + 4]);
			mass_fraction[INDEX] = mass_fraction[INDEX] - (dt / dx) * (flux_X[R_interface + 5] - flux_X[L_interface + 5])
									- (dt / dy) * (flux_Y[T_interface + 5] - flux_Y[B_interface + 5]);

			//Get new variable
			p0[INDEX] = mass[INDEX];
			p1[INDEX] = momentum_X[INDEX] / mass[INDEX];
			p2[INDEX] = momentum_Y[INDEX] / mass[INDEX];
			p5[INDEX] = mass_fraction[INDEX] / p0[INDEX];
			e_target = (energy[INDEX] / p0[INDEX]) - 0.5 * (p1[INDEX] * p1[INDEX] + p2[INDEX] * p2[INDEX]);
			T_old = p3[INDEX];
			phi_new = p5[INDEX];
			p3[INDEX] = compute_T(T_old, phi_new, e_target);
			R_mix = R_dry * (1.0 - p5[INDEX]) + R_v * p5[INDEX];
			p4[INDEX] = p0[INDEX] * R_mix * p3[INDEX];

			P_sat = 611.0 * exp((17.27 * (p3[INDEX] - 273.15)) / ((p3[INDEX] - 273.15) + 237.3)); // Tetens equation
			P_v = (p0[INDEX] * p5[INDEX]) * R_v * p3[INDEX];
			p6[INDEX] = P_v / P_sat;
			if (p4[INDEX] <= P_sat) {
				phi_max = 1.0; // 壓力低於飽和蒸氣壓，允許 100% 水氣
			} else {
				phi_max = (P_sat / R_v) / (((p4[INDEX] - P_sat) / R_dry) + (P_sat / R_v));
			}

			if (p5[INDEX] > phi_max || p6[INDEX] > 1.0){
				p5[INDEX] = phi_max;
				p6[INDEX] = 1.0; // 100%!
				mass_fraction[INDEX] = p5[INDEX] * p0[INDEX];
			}
			if (p5[INDEX] < 0.0 || p6[INDEX] < 0.0){
				p5[INDEX] = 0.0;
				p6[INDEX] = 0.0;
				mass_fraction[INDEX] = 0.0;
			}

		}
	}
}
/*
void Calc_primitive_variable(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6,
			     float *mass, float *momentum_X, float *momentum_Y, float *energy, float *mass_fraction,
			     float *flux_X, float *flux_Y,
			     float R_dry, float R_v, float dx, float dy, float dt, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Calc_primitive_variable<<<GPB, TPB>>>(p0, p1, p2, p3, p4, p5, p6,
						  mass, momentum_X, momentum_Y, energy, mass_fraction,
						  flux_X, flux_Y,
						  R_dry, R_v, dx, dy, dt, NX, NY, N_CELLS);
}
*/
float Compute_dt(float *MAX_Freq, float *p0, float *p1, float *p2, float *p3, float *p5, float dx, float dy, float R_dry, float R_v, float CFL, int NX, int NY, int N_CELLS){
/*
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Compute_MAX_CFL<<<GPB, TPB>>>(MAX_Freq, p0, p1, p2, p3, p5, dx, dy, R_dry, R_v, NX, NY, N_CELLS);
	MAX_CFL = thrust::max_element(thrust::device, MAX_Freq, MAX_Freq + N_CELLS);
	float h_MAX_CFL;
	cudaMemcpy(&h_MAX_CFL, MAX_CFL, sizeof(float), cudaMemcpyDeviceToHost);
*/
	float MAX_CFL;
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY+2) + j;
			MAX_CFL = MAX_Freq[INDEX];
		}
	}
	float dt = CFL / MAX_CFL;
	return dt;
}
