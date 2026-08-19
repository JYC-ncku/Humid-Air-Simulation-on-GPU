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

	float *mass_old, *momentum_X_old, *momentum_Y_old, *energy_old,
	      *mass_now, *momentum_X_now, *momentum_Y_now, *energy_now;

	float t = 0.0;
	float t_FINAL = 0.2;
	float R = 1.0;
	float GAMMA = 1.4;
	float CFL = 0.5;
	Allocate_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &mass, &momentum_X, &momentum_Y, &energy,
			&mass_flux_X, &momentum_X_flux_X, &momentum_Y_flux_X, &energy_flux_X,
			&mass_flux_Y, &momentum_X_flux_Y, &momentum_Y_flux_Y, &energy_flux_Y,
			&mass_old, &momentum_X_old, &momentum_Y_old, &energy_old,
			&mass_now, &momentum_X_now, &momentum_Y_now, &energy_now,
			N_CELLS);
	//Initial condition (p0 is density, p1 is X-direction veloctiy, p2 is Y-direction veloctiy, p3 is temperature, p4 is pressure)
	for (int i = 2; i < NX + 2; i++){
		for (int j = 2; j < NY + 2; j++){
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

		//The conserved quantity is temporarily stored.
		for (int i = 2; i < NX + 2; i++){
			for (int j =2; j < NY + 2; j++){
				int INDEX = i * (NY+4) + j;
				mass_old[INDEX] = mass[INDEX];
				momentum_X_old[INDEX] = momentum_X[INDEX];
				momentum_Y_old[INDEX] = momentum_Y[INDEX];
				energy_old[INDEX] = energy[INDEX];
			}
		}

		//X-direction flux
		for (int i = 1; i < NX + 2; i++){
			for (int j = 2; j < NY + 2; j++){
				int INDEX_L = (i-1) * (NY+4) + j;
				int INDEX = i * (NY+4) + j;
				int INDEX_R = (i+1) * (NY+4) + j;
				int INDEX_RR = (i+2) * (NY+4) + j;
				float rho_L = p0[INDEX_L];
				float rho_C = p0[INDEX];	//LEFT = BOTTOM
				float rho_R = p0[INDEX_R];
				float rho_RR = p0[INDEX_RR];
				float drho_dx_L = MINMOD(rho_L, rho_C, rho_R, dx);
				float drho_dx_R = MINMOD(rho_C, rho_R, rho_RR, dx);
				float rho_L_star = rho_C + 0.5 * dx * drho_dx_L;
				float rho_R_star = rho_R - 0.5 * dx * drho_dx_R;

				float u_L = p1[INDEX_L];
				float u_C = p1[INDEX];
				float u_R = p1[INDEX_R];
				float u_RR = p1[INDEX_RR];
				float du_dx_L = MINMOD(u_L, u_C, u_R, dx);
				float du_dx_R = MINMOD(u_C, u_R, u_RR, dx);
				float u_L_star = u_C + 0.5 * dx * du_dx_L;
				float u_R_star = u_R - 0.5 * dx * du_dx_R;

				float v_L = p2[INDEX_L];
				float v_C = p2[INDEX];
				float v_R = p2[INDEX_R];
				float v_RR = p2[INDEX_RR];
				float dv_dx_L = MINMOD(v_L, v_C, v_R, dx);
				float dv_dx_R = MINMOD(v_C, v_R, v_RR, dx);
				float v_L_star = v_C + 0.5 * dx * dv_dx_L;
				float v_R_star = v_R - 0.5 * dx * dv_dx_R;

				float T_L = p3[INDEX_L];
				float T_C = p3[INDEX];
				float T_R = p3[INDEX_R];
				float T_RR = p3[INDEX_RR];
				float dT_dx_L = MINMOD(T_L, T_C, T_R, dx);
				float dT_dx_R = MINMOD(T_C, T_R, T_RR, dx);
				float T_L_star = T_C + 0.5 * dx * dT_dx_L;
				float T_R_star = T_R - 0.5 * dx * dT_dx_R;

				float P_L = p4[INDEX_L];
				float P_C = p4[INDEX];
				float P_R = p4[INDEX_R];
				float P_RR = p4[INDEX_RR];
				float dP_dx_L = MINMOD(P_L, P_C, P_R, dx);
				float dP_dx_R = MINMOD(P_C, P_R, P_RR, dx);
				float P_L_star = P_C + 0.5 * dx * dP_dx_L;
				float P_R_star = P_R - 0.5 * dx * dP_dx_R;

				float e_L_star = 0.5 * rho_L_star * (u_L_star * u_L_star + v_L_star * v_L_star) + P_L_star / (GAMMA - 1);
				float e_R_star = 0.5 * rho_R_star * (u_R_star * u_R_star + v_R_star * v_R_star) + P_R_star / (GAMMA - 1);
				float a_L_star = sqrt(GAMMA * R * T_L_star);
				float a_R_star = sqrt(GAMMA * R * T_R_star);

				float W_LOCAL_MAX_X = MAX_WAVE_SPEED(u_L_star, u_R_star, a_L_star, a_R_star);
				if (W_LOCAL_MAX_X > W_GLOBAL_MAX){
					W_GLOBAL_MAX = W_LOCAL_MAX_X;
				}
				Calc_flux_X(rho_L_star, rho_R_star, u_L_star, u_R_star, v_L_star, v_R_star, P_L_star, P_R_star, e_L_star, e_R_star, a_L_star, a_R_star,
					    mass_flux_X, momentum_X_flux_X, momentum_Y_flux_X, energy_flux_X, INDEX);
			}
		}
		//Y-direction flux
		for (int i = 2; i < NX + 2; i++){
			for (int j = 1; j < NY + 2; j++){
				int INDEX_B = i * (NY+4) + (j-1);
				int INDEX = i * (NY+4) + j;
				int INDEX_T = i * (NY+4) + (j+1);
				int INDEX_TT = i * (NY+4) + (j+2);
				float rho_B = p0[INDEX_B];	//LEFT = BOTTOM
				float rho_C = p0[INDEX];
				float rho_T = p0[INDEX_T];
				float rho_TT = p0[INDEX_TT];
				float drho_dy_B = MINMOD(rho_B, rho_C, rho_T, dy);
				float drho_dy_T = MINMOD(rho_C, rho_T, rho_TT, dy);
				float rho_B_star = rho_C + 0.5 * dy * drho_dy_B;
				float rho_T_star = rho_T - 0.5 * dy * drho_dy_T;

				float u_B = p1[INDEX_B];
				float u_C = p1[INDEX];
				float u_T = p1[INDEX_T];
				float u_TT = p1[INDEX_TT];
				float du_dy_B = MINMOD(u_B, u_C, u_T, dy);
				float du_dy_T = MINMOD(u_C, u_T, u_TT, dy);
				float u_B_star = u_C + 0.5 * dy * du_dy_B;
				float u_T_star = u_T - 0.5 * dy * du_dy_T;

				float v_B = p2[INDEX_B];
				float v_C = p2[INDEX];
				float v_T = p2[INDEX_T];
				float v_TT = p2[INDEX_TT];
				float dv_dy_B = MINMOD(v_B, v_C, v_T, dy);
				float dv_dy_T = MINMOD(v_C, v_T, v_TT, dy);
				float v_B_star = v_C + 0.5 * dy * dv_dy_B;
				float v_T_star = v_T - 0.5 * dy * dv_dy_T;

				float T_B = p3[INDEX_B];
				float T_C = p3[INDEX];
				float T_T = p3[INDEX_T];
				float T_TT = p3[INDEX_TT];
				float dT_dy_B = MINMOD(T_B, T_C, T_T, dy);
				float dT_dy_T = MINMOD(T_C, T_T, T_TT, dy);
				float T_B_star = T_C + 0.5 * dy * dT_dy_B;
				float T_T_star = T_T - 0.5 * dy * dT_dy_T;

				float P_B = p4[INDEX_B];
				float P_C = p4[INDEX];
				float P_T = p4[INDEX_T];
				float P_TT = p4[INDEX_TT];
				float dP_dy_B = MINMOD(P_B, P_C, P_T, dy);
				float dP_dy_T = MINMOD(P_C, P_T, P_TT, dy);
				float P_B_star = P_C + 0.5 * dy * dP_dy_B;
				float P_T_star = P_T - 0.5 * dy * dP_dy_T;

				float e_B_star = 0.5 * rho_B_star * (u_B_star * u_B_star + v_B_star * v_B_star) + P_B_star / (GAMMA - 1);
				float e_T_star = 0.5 * rho_T_star * (u_T_star * u_T_star + v_T_star * v_T_star) + P_T_star / (GAMMA - 1);
				float a_B_star = sqrt(GAMMA * R * T_B_star);
				float a_T_star = sqrt(GAMMA * R * T_T_star);

				float W_LOCAL_MAX_Y = MAX_WAVE_SPEED(v_B_star, v_T_star, a_B_star, a_T_star);
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
//======================================================================================================================================================================================================
		//X-direction flux
		for (int i = 1; i < NX + 2; i++){
			for (int j = 2; j < NY + 2; j++){
				int INDEX_L = (i-1) * (NY+4) + j;
				int INDEX = i * (NY+4) + j;
				int INDEX_R = (i+1) * (NY+4) + j;
				int INDEX_RR = (i+2) * (NY+4) + j;
				float rho_L = p0[INDEX_L];
				float rho_C = p0[INDEX];	//LEFT = BOTTOM
				float rho_R = p0[INDEX_R];
				float rho_RR = p0[INDEX_RR];
				float drho_dx_L = MINMOD(rho_L, rho_C, rho_R, dx);
				float drho_dx_R = MINMOD(rho_C, rho_R, rho_RR, dx);
				float rho_L_star = rho_C + 0.5 * dx * drho_dx_L;
				float rho_R_star = rho_R - 0.5 * dx * drho_dx_R;

				float u_L = p1[INDEX_L];
				float u_C = p1[INDEX];
				float u_R = p1[INDEX_R];
				float u_RR = p1[INDEX_RR];
				float du_dx_L = MINMOD(u_L, u_C, u_R, dx);
				float du_dx_R = MINMOD(u_C, u_R, u_RR, dx);
				float u_L_star = u_C + 0.5 * dx * du_dx_L;
				float u_R_star = u_R - 0.5 * dx * du_dx_R;

				float v_L = p2[INDEX_L];
				float v_C = p2[INDEX];
				float v_R = p2[INDEX_R];
				float v_RR = p2[INDEX_RR];
				float dv_dx_L = MINMOD(v_L, v_C, v_R, dx);
				float dv_dx_R = MINMOD(v_C, v_R, v_RR, dx);
				float v_L_star = v_C + 0.5 * dx * dv_dx_L;
				float v_R_star = v_R - 0.5 * dx * dv_dx_R;

				float T_L = p3[INDEX_L];
				float T_C = p3[INDEX];
				float T_R = p3[INDEX_R];
				float T_RR = p3[INDEX_RR];
				float dT_dx_L = MINMOD(T_L, T_C, T_R, dx);
				float dT_dx_R = MINMOD(T_C, T_R, T_RR, dx);
				float T_L_star = T_C + 0.5 * dx * dT_dx_L;
				float T_R_star = T_R - 0.5 * dx * dT_dx_R;

				float P_L = p4[INDEX_L];
				float P_C = p4[INDEX];
				float P_R = p4[INDEX_R];
				float P_RR = p4[INDEX_RR];
				float dP_dx_L = MINMOD(P_L, P_C, P_R, dx);
				float dP_dx_R = MINMOD(P_C, P_R, P_RR, dx);
				float P_L_star = P_C + 0.5 * dx * dP_dx_L;
				float P_R_star = P_R - 0.5 * dx * dP_dx_R;

				float e_L_star = 0.5 * rho_L_star * (u_L_star * u_L_star + v_L_star * v_L_star) + P_L_star / (GAMMA - 1);
				float e_R_star = 0.5 * rho_R_star * (u_R_star * u_R_star + v_R_star * v_R_star) + P_R_star / (GAMMA - 1);
				float a_L_star = sqrt(GAMMA * R * T_L_star);
				float a_R_star = sqrt(GAMMA * R * T_R_star);

				Calc_flux_X(rho_L_star, rho_R_star, u_L_star, u_R_star, v_L_star, v_R_star, P_L_star, P_R_star, e_L_star, e_R_star, a_L_star, a_R_star,
					    mass_flux_X, momentum_X_flux_X, momentum_Y_flux_X, energy_flux_X, INDEX);
			}
		}
		//Y-direction flux
		for (int i = 2; i < NX + 2; i++){
			for (int j = 1; j < NY + 2; j++){
				int INDEX_B = i * (NY+4) + (j-1);
				int INDEX = i * (NY+4) + j;
				int INDEX_T = i * (NY+4) + (j+1);
				int INDEX_TT = i * (NY+4) + (j+2);
				float rho_B = p0[INDEX_B];	//LEFT = BOTTOM
				float rho_C = p0[INDEX];
				float rho_T = p0[INDEX_T];
				float rho_TT = p0[INDEX_TT];
				float drho_dy_B = MINMOD(rho_B, rho_C, rho_T, dy);
				float drho_dy_T = MINMOD(rho_C, rho_T, rho_TT, dy);
				float rho_B_star = rho_C + 0.5 * dy * drho_dy_B;
				float rho_T_star = rho_T - 0.5 * dy * drho_dy_T;

				float u_B = p1[INDEX_B];
				float u_C = p1[INDEX];
				float u_T = p1[INDEX_T];
				float u_TT = p1[INDEX_TT];
				float du_dy_B = MINMOD(u_B, u_C, u_T, dy);
				float du_dy_T = MINMOD(u_C, u_T, u_TT, dy);
				float u_B_star = u_C + 0.5 * dy * du_dy_B;
				float u_T_star = u_T - 0.5 * dy * du_dy_T;

				float v_B = p2[INDEX_B];
				float v_C = p2[INDEX];
				float v_T = p2[INDEX_T];
				float v_TT = p2[INDEX_TT];
				float dv_dy_B = MINMOD(v_B, v_C, v_T, dy);
				float dv_dy_T = MINMOD(v_C, v_T, v_TT, dy);
				float v_B_star = v_C + 0.5 * dy * dv_dy_B;
				float v_T_star = v_T - 0.5 * dy * dv_dy_T;

				float T_B = p3[INDEX_B];
				float T_C = p3[INDEX];
				float T_T = p3[INDEX_T];
				float T_TT = p3[INDEX_TT];
				float dT_dy_B = MINMOD(T_B, T_C, T_T, dy);
				float dT_dy_T = MINMOD(T_C, T_T, T_TT, dy);
				float T_B_star = T_C + 0.5 * dy * dT_dy_B;
				float T_T_star = T_T - 0.5 * dy * dT_dy_T;

				float P_B = p4[INDEX_B];
				float P_C = p4[INDEX];
				float P_T = p4[INDEX_T];
				float P_TT = p4[INDEX_TT];
				float dP_dy_B = MINMOD(P_B, P_C, P_T, dy);
				float dP_dy_T = MINMOD(P_C, P_T, P_TT, dy);
				float P_B_star = P_C + 0.5 * dy * dP_dy_B;
				float P_T_star = P_T - 0.5 * dy * dP_dy_T;

				float e_B_star = 0.5 * rho_B_star * (u_B_star * u_B_star + v_B_star * v_B_star) + P_B_star / (GAMMA - 1);
				float e_T_star = 0.5 * rho_T_star * (u_T_star * u_T_star + v_T_star * v_T_star) + P_T_star / (GAMMA - 1);
				float a_B_star = sqrt(GAMMA * R * T_B_star);
				float a_T_star = sqrt(GAMMA * R * T_T_star);

				Calc_flux_Y(rho_B_star, rho_T_star, u_B_star, u_T_star, v_B_star, v_T_star, P_B_star, P_T_star, e_B_star, e_T_star, a_B_star, a_T_star,
					    mass_flux_Y, momentum_X_flux_Y, momentum_Y_flux_Y, energy_flux_Y, INDEX);
			}
		}
//=======================================================================================================================================================================================================

		Calc_conserved(mass, momentum_X, momentum_Y, energy,
			       mass_flux_X, momentum_X_flux_X, momentum_Y_flux_X, energy_flux_X,
			       mass_flux_Y, momentum_X_flux_Y, momentum_Y_flux_Y, energy_flux_Y,
			       dx, dy, dt, NX, NY);
		//Store the conserved quantities calculated in the second step.
		for (int i = 2; i < NX + 2; i++){
			for (int j =2; j < NY + 2; j++){
				int INDEX = i * (NY+4) + j;
				mass_now[INDEX] = mass[INDEX];
				momentum_X_now[INDEX] = momentum_X[INDEX];
				momentum_Y_now[INDEX] = momentum_Y[INDEX];
				energy_now[INDEX] = energy[INDEX];
			}
		}
		//Calculate the new conserved quantity.
		for (int i = 2; i < NX + 2; i++){
			for (int j =2; j < NY + 2; j++){
				int INDEX = i * (NY+4) + j;
				mass[INDEX] = 0.5 * (mass_old[INDEX] + mass_now[INDEX]);
				momentum_X[INDEX] = 0.5 * (momentum_X_old[INDEX] + momentum_X_now[INDEX]);
				momentum_Y[INDEX] = 0.5 * (momentum_Y_old[INDEX] + momentum_Y_now[INDEX]);
				energy[INDEX] = 0.5 * (energy_old[INDEX] + energy_now[INDEX]);
			}
		}

		Calc_variable(p0, p1, p2, p3, p4, mass, momentum_X, momentum_Y, energy, GAMMA, R, NX, NY);

		t += dt;
	}

	FILE *pFile = fopen("Results_of_5000_cells_X_direction", "w");
	for (int i = 2; i < NX + 2; i++){
		for (int j = 2; j < NY + 2; j++){
//	for (int j = 2; j < NY + 2; j++){
//		for (int i = 2; i < NX + 2; i++){
			int INDEX = i *(NY+4) + j;
			float X = (i - 1.5) * dx;
			float Y = (j - 1.5) * dy;
			fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", X, Y, p0[INDEX], p1[INDEX], p2[INDEX], p3[INDEX], p4[INDEX]);
		}
	}
	fclose(pFile);

	Free_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &mass, &momentum_X, &momentum_Y, &energy,
		    &mass_flux_X, &momentum_X_flux_X, &momentum_Y_flux_X, &energy_flux_X,
		    &mass_flux_Y, &momentum_X_flux_Y, &momentum_Y_flux_Y, &energy_flux_Y,
		    &mass_old, &momentum_X_old, &momentum_Y_old, &energy_old,
		    &mass_now, &momentum_X_now, &momentum_Y_now, &energy_now);

	return 0;
}
