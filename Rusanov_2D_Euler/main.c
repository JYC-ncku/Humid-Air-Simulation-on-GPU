#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "Calc_Flux.h"
#include "Calc_conserved.h"
#include "Calc_variable.h"
#include "Boundary.h"
int main(){
	int NX = 1000;
	int NY = 5;
//	int NX = 5;
//	int NY = 1000;
	int N_CELLS = (NX+2) * (NY+2);
	float L = 1.0;
	float H = 0.005;
//	float L = 0.005;
//	float H = 1.0;
	float dx = L/NX;
	float dy = H/NY;
	float *x, *y, *p0, *p1, *p2, *p3, *p4, *mass, *momentum_X, *momentum_Y, *energy,
	      *mass_flux_X, *momentum_X_flux_X, *momentum_Y_flux_X, *energy_flux_X,
	      *mass_flux_Y, *momentum_X_flux_Y, *momentum_Y_flux_Y, *energy_flux_Y;
	float t = 0.0;
	float t_FINAL = 0.2;
	float R = 1.0;
	float GAMMA = 1.4;
	float CFL = 0.5;
	Allocate_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &mass, &momentum_X, &momentum_Y, &energy,
			&mass_flux_X, &momentum_X_flux_X, &momentum_Y_flux_X, &energy_flux_X,
			&mass_flux_Y, &momentum_X_flux_Y, &momentum_Y_flux_Y, &energy_flux_Y,
			N_CELLS);
	//Initial condition (p0 is density, p1 is X-direction veloctiy, p2 is Y-direction veloctiy, p3 is temperature, p4 is pressure)
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY+2) + j;
			if (i < NX/2){
//			if (j < NY/2){
				p0[INDEX] = 10.0;
				p1[INDEX] = 0.0;
				p2[INDEX] = 0.0;
				p3[INDEX] = 1.0;
			} else {
				p0[INDEX] = 1.0;
				p1[INDEX] = 0.0;
				p2[INDEX] = 0.0;
				p3[INDEX] = 1.0;
			}
			p4[INDEX] = p0[INDEX] * R * p3[INDEX];
			mass[INDEX] = p0[INDEX];
			momentum_X[INDEX] = p0[INDEX] * p1[INDEX];
			momentum_Y[INDEX] = p0[INDEX] * p2[INDEX];
			energy[INDEX] = 0.5 * p0[INDEX] * (p1[INDEX] * p1[INDEX] + p2[INDEX] * p2[INDEX]) +  (p4[INDEX] / (GAMMA - 1.0));
		}
	}

	while(t < t_FINAL){
		float W_GLOBAL_MAX = 1e-10;
		Boundary(p0, p1, p2, p3, p4, NX, NY);
		//X-direction flux
		for (int i = 0; i < NX + 1; i++){
			for (int j = 1; j < NY + 1; j++){
				int INDEX = i * (NY+2) + j;
				int INDEX_R = (i+1) * (NY+2) + j;
				float rho_L = p0[INDEX];	//LEFT = BOTTOM
				float rho_R = p0[INDEX_R];
				float u_L = p1[INDEX];
				float u_R = p1[INDEX_R];
				float v_L = p2[INDEX];
				float v_R = p2[INDEX_R];
				float T_L = p3[INDEX];
				float T_R = p3[INDEX_R];
				float P_L = p4[INDEX];
				float P_R = p4[INDEX_R];
				float e_L = 0.5 * rho_L * (u_L * u_L + v_L * v_L) + P_L / (GAMMA - 1);
				float e_R = 0.5 * rho_R * (u_R * u_R + v_R * v_R) + P_R / (GAMMA - 1);
				float a_L = sqrt(GAMMA * R * T_L);
				float a_R = sqrt(GAMMA * R * T_R);
				float W_LOCAL_MAX_X = MAX_WAVE_SPEED(u_L, u_R, a_L, a_R);
				if (W_LOCAL_MAX_X > W_GLOBAL_MAX){
					W_GLOBAL_MAX = W_LOCAL_MAX_X;
				}
				Calc_flux_X(rho_L, rho_R, u_L, u_R, v_L, v_R, T_L, T_R, P_L, P_R, e_L, e_R, W_LOCAL_MAX_X,
					    mass_flux_X, momentum_X_flux_X, momentum_Y_flux_X, energy_flux_X, INDEX);
			}
		}
		//Y-direction flux
		for (int i = 1; i < NX + 1; i++){
			for (int j = 0; j < NY + 1; j++){
				int INDEX = i * (NY+2) + j;
				int INDEX_T = i * (NY+2) + (j+1);
				float rho_B = p0[INDEX];	//LEFT = BOTTOM
				float rho_T = p0[INDEX_T];
				float u_B = p1[INDEX];
				float u_T = p1[INDEX_T];
				float v_B = p2[INDEX];
				float v_T = p2[INDEX_T];
				float T_B = p3[INDEX];
				float T_T = p3[INDEX_T];
				float P_B = p4[INDEX];
				float P_T = p4[INDEX_T];
				float e_B = 0.5 * rho_B * (u_B * u_B + v_B * v_B) + P_B / (GAMMA - 1);
				float e_T = 0.5 * rho_T * (u_T * u_T + v_T * v_T) + P_T / (GAMMA - 1);
				float a_B = sqrt(GAMMA * R * T_B);
				float a_T = sqrt(GAMMA * R * T_T);
				float W_LOCAL_MAX_Y = MAX_WAVE_SPEED(v_B, v_T, a_B, a_T);
				if (W_LOCAL_MAX_Y > W_GLOBAL_MAX){
					W_GLOBAL_MAX = W_LOCAL_MAX_Y;
				}
				Calc_flux_Y(rho_B, rho_T, u_B, u_T, v_B, v_T, T_B, T_T, P_B, P_T, e_B, e_T, W_LOCAL_MAX_Y,
					    mass_flux_Y, momentum_X_flux_Y, momentum_Y_flux_Y, energy_flux_Y, INDEX);
			}
		}

		float dt = CFL * dx / W_GLOBAL_MAX;	//這裡用的是正方形網格(dx = dy)，所以可以把delta提出來。

		Calc_conserved(mass, momentum_X, momentum_Y, energy,
			       mass_flux_X, momentum_X_flux_X, momentum_Y_flux_X, energy_flux_X,
			       mass_flux_Y, momentum_X_flux_Y, momentum_Y_flux_Y, energy_flux_Y,
			       dx, dy, dt, NX, NY);
		Calc_variable(p0, p1, p2, p3, p4, mass, momentum_X, momentum_Y, energy, GAMMA, R, NX, NY);

		t += dt;
	}

	FILE *pFile = fopen("Results_of_5000_cells_X_direction", "w");
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
//	for (int j = 1; j < NY + 1; j++){
//		for (int i = 1; i < NX + 1; i++){
			int INDEX = i *(NY+2) + j;
			float X = (i - 0.5) * dx;
			float Y = (j - 0.5) * dy;
			fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", X, Y, p0[INDEX], p1[INDEX], p2[INDEX], p3[INDEX], p4[INDEX]);
		}
	}
	fclose(pFile);

	Free_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &mass, &momentum_X, &momentum_Y, &energy,
		    &mass_flux_X, &momentum_X_flux_X, &momentum_Y_flux_X, &energy_flux_X,
		    &mass_flux_Y, &momentum_X_flux_Y, &momentum_Y_flux_Y, &energy_flux_Y);

	return 0;
}
