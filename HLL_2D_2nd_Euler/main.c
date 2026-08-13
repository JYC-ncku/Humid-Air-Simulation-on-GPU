#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "Calc_Flux.h"
#include "Calc_conserved.h"
#include "Calc_variable.h"
#include "Boundary.h"

float MINMOD(float QL, float QC, float QR, float dx){
	float dU_dx;
	float Forward = (QR - QC) / dx;
	float Backward = (QC - QL) / dx;
	if (Backward * Forward < 0){
		dU_dx = 0;
	} else if ( fabs(Forward) < fabs(Backward) ){
			dU_dx = Forward;
		} else {
			dU_dx = Backward;
	}
	return dU_dx;
}

int main(){
	int NX = 1000;
	int NY = 5;
//	int NX = 5;
//	int NY = 1000;
	int N_CELLS = (NX+4) * (NY+4);
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
			int INDEX = i * (NY+4) + j;
			if (i < (NX/2 + 1)){
//			if (j < (NY/2 + 1)){
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
				int INDEX_L = (i-1) * (NY+4) + j;
				int INDEX = i * (NY+4) + j;
				int INDEX_R = (i+1) * (NY+4) + j;
				int INDEX_RR = (i+2) * (NY+4) + j;
				float rho_L = p0[INDEX_L];
				float rho_C = p0[INDEX];	//LEFT = BOTTOM
				float rho_R = p0[INDEX_R];
				float rho_RR = p0[INDEX_RR];
				float rho_L_star = MINMOD(rho_L, rho_C, rho_R, dx);
				float rho_R_star = MINMOD(rho_C, rho_R, rho_RR, dx);

				float u_L = p1[INDEX_L];
				float u_C = p1[INDEX];
				float u_R = p1[INDEX_R];
				float u_RR = p1[INDEX_RR];
				float u_L_star = MINMOD(u_L, u_C, u_R, dx);
				float u_R_star = MINMOD(u_C, u_R, u_RR, dx);

				float v_L = p2[INDEX_L];
				float v_C = p2[INDEX];
				float v_R = p2[INDEX_R];
				float v_RR = p2[INDEX_RR];
				float v_L_star = MINMOD(v_L, v_C, v_R, dx);
				float v_R_star = MINMOD(v_C, v_R, v_RR, dx);

				float T_L = p3[INDEX_L];
				float T_C = p3[INDEX];
				float T_R = p3[INDEX_R];
				float T_RR = p3[INDEX_RR];

				float P_L = p4[INDEX_L];
				float P_C = p4[INDEX];
				float P_R = p4[INDEX_R];
				float P_RR = p4[INDEX_RR];
				float P_L_star = MINMOD(P_L, P_C, P_R, dx);
				float P_R_star = MINMOD(P_C, P_R, P_RR, dx);

				float e_L = 0.5 * rho_L * (u_L * u_L + v_L * v_L) + P_L / (GAMMA - 1);
				float e_C = 0.5 * rho_C * (u_C * u_C + v_C * v_C) + P_C / (GAMMA - 1);
				float e_R = 0.5 * rho_R * (u_R * u_R + v_R * v_R) + P_R / (GAMMA - 1);
				float e_RR = 0.5 * rho_RR * (u_RR * u_RR + v_RR * v_RR) + P_RR / (GAMMA - 1);
				float e_L_star = MINMOD(e_L, e_C, e_R, dx);
				float e_R_star = MINMOD(e_C, e_R, e_RR, dx);

				float a_L = sqrt(GAMMA * R * T_L);
				float a_C = sqrt(GAMMA * R * T_C);
				float a_R = sqrt(GAMMA * R * T_R);
				float a_RR = sqrt(GAMMA * R * T_RR);
				float a_L_star = MINMOD(a_L, a_C, a_R, dx);
				float a_R_star = MINMOD(a_C, a_R, a_RR, dx);

				float W_LOCAL_MAX_X = MAX_WAVE_SPEED(u_L, u_C, a_C, a_C);
				if (W_LOCAL_MAX_X > W_GLOBAL_MAX){
					W_GLOBAL_MAX = W_LOCAL_MAX_X;
				}
				Calc_flux_X(rho_L_star, rho_R_star, u_L_star, u_R_star, v_L_star, v_R_star, P_L_star, P_R_star, e_L_star, e_R_star, a_L_star, a_R_star,
					    mass_flux_X, momentum_X_flux_X, momentum_Y_flux_X, energy_flux_X, INDEX);
			}
		}
		//Y-direction flux
		for (int i = 1; i < NX + 1; i++){
			for (int j = 0; j < NY + 1; j++){
				int INDEX_B = i * (NY+4) + (j-1);
				int INDEX = i * (NY+4) + j;
				int INDEX_T = i * (NY+4) + (j+1);
				int INDEX_TT = i * (NY+4) + (j+2);
				float rho_B = p0[INDEX_B];	//LEFT = BOTTOM
				float rho_C = p0[INDEX];
				float rho_T = p0[INDEX_T];
				float rho_TT = p0[INDEX_TT];
				float rho_B_star = MINMOD(rho_B, rho_C, rho_T, dy);
				float rho_T_star = MINMOD(rho_C, rho_T, rho_TT, dy);

				float u_B = p1[INDEX_B];
				float u_C = p1[INDEX];
				float u_T = p1[INDEX_T];
				float u_TT = p1[INDEX_TT];
				float u_B_star = MINMOD(u_B, u_C, u_T, dy);
				float u_T_star = MINMOD(u_C, u_T, u_TT, dy);

				float v_B = p2[INDEX_B];
				float v_C = p2[INDEX];
				float v_T = p2[INDEX_T];
				float v_TT = p2[INDEX_TT];
				float v_B_star = MINMOD(v_B, v_C, v_T, dy);
				float v_T_star = MINMOD(v_C, v_T, v_TT, dy);

				float T_B = p3[INDEX_B];
				float T_C = p3[INDEX];
				float T_T = p3[INDEX_T];
				float T_TT = p3[INDEX_TT];

				float P_B = p4[INDEX_B];
				float P_C = p4[INDEX];
				float P_T = p4[INDEX_T];
				float P_TT = p4[INDEX_TT];
				float P_B_star = MINMOD(P_B, P_C, P_T, dy);
				float P_T_star = MINMOD(P_C, P_T, P_TT, dy);

				float e_B = 0.5 * rho_B * (u_B * u_B + v_B * v_B) + P_B / (GAMMA - 1);
				float e_C = 0.5 * rho_C * (u_C * u_C + v_C * v_C) + P_C / (GAMMA - 1);
				float e_T = 0.5 * rho_T * (u_T * u_T + v_T * v_T) + P_T / (GAMMA - 1);
				float e_TT = 0.5 * rho_TT * (u_TT * u_TT + v_TT * v_TT) + P_TT / (GAMMA - 1);
				float e_B_star = MINMOD(e_B, e_C, e_T, dy);
				float e_T_star = MINMOD(e_C, e_T, e_TT, dy);

				float a_B = sqrt(GAMMA * R * T_B);
				float a_C = sqrt(GAMMA * R * T_C);
				float a_T = sqrt(GAMMA * R * T_T);
				float a_TT = sqrt(GAMMA * R * T_TT);
				float a_B_star = MINMOD(a_B, a_C, a_T, dy);
				float a_T_star = MINMOD(a_C, a_T, a_TT, dy);

				float W_LOCAL_MAX_Y = MAX_WAVE_SPEED(v_B, v_C, a_B, a_C);
				if (W_LOCAL_MAX_Y > W_GLOBAL_MAX){
					W_GLOBAL_MAX = W_LOCAL_MAX_Y;
				}
				Calc_flux_Y(rho_B_star, rho_T_star, u_B_star, u_T_star, v_B_star, v_T_star, P_B_star, P_T_star, e_B_star, e_T_star, a_B_star, a_T_star,
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
			int INDEX = i *(NY+4) + j;
			float X = (i - 1.5) * dx;
			float Y = (j - 1.5) * dy;
			fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", X, Y, p0[INDEX], p1[INDEX], p2[INDEX], p3[INDEX], p4[INDEX]);
		}
	}
	fclose(pFile);

	Free_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &mass, &momentum_X, &momentum_Y, &energy,
		    &mass_flux_X, &momentum_X_flux_X, &momentum_Y_flux_X, &energy_flux_X,
		    &mass_flux_Y, &momentum_X_flux_Y, &momentum_Y_flux_Y, &energy_flux_Y);

	return 0;
}
