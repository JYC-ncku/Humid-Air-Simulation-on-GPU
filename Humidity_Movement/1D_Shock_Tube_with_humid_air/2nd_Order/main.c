#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "Initial.h"
#include "Boundary.h"
#include "Calc_flux.h"
#include "Primitive_variable.h"
#include "Compute_Cv.h"

float MINMOD(float U_L, float U_C, float U_R, float dx){
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

float Compute_dt(float *p1, float *p2, float *p4, float dx, float CFL, float R_dry, float R_v, int N_CELLS){
	float MAX_Freq = -1.0;
	float R_mix, Cv_mix, Gamma_mix, dt;
	// Only care inner cells.
	for (int i = 2; i < N_CELLS + 2; i++){
		float u = p1[i];
		float T = p2[i];
		float Y = p4[i];
		if (T<0){
			printf("Error: Negative temperature in cell %d: T = %f\n. Aborting.", i, T);
		exit(1);
		}
		R_mix = R_dry * (1.0 - p4[i]) + R_v * p4[i];
		Cv_mix = Compute_Cv(T, Y);
		Gamma_mix = 1.0 + R_mix / Cv_mix;
		float a = sqrt(Gamma_mix * R_mix * T);
		float Freq = (fabs(u) + a) / dx;
		if (Freq>MAX_Freq){
			MAX_Freq = Freq;
		}
	}
	dt = CFL / MAX_Freq;
	return dt;
}

int main(){
	int N_CELLS = 1000;
	float *x, *p0, *p1, *p2, *p3, *p4, *p5,
	      *mass, *momentum, *energy, *mass_fraction, *mass_flux, *momentum_flux, *energy_flux, *mass_fraction_flux;
	float L = 0.01; // unit: m
	float t = 0;
	float t_FINAL = 7e-6;
	float CFL = 0.05;
	float dx = L/N_CELLS;

	float D = 1.837e-5; //Diffusivity of water vapor. unit:(m^2/s)

	float R_bar = 8.3145; // unti:J/(mol*K)
	float MW_H2O = 0.01802; // unit:kg/mol
	float MW_air = 0.02897; // unit:kg/mol
	float R_v = R_bar / MW_H2O; // unit:J/(kg*k) R = R_bar / Molecular weight
	float R_dry = R_bar / MW_air;

	Allocate_memory(&x, &p0, &p1, &p2, &p3, &p4, &p5,
			&mass, &momentum, &energy, &mass_fraction, &mass_flux, &momentum_flux, &energy_flux, &mass_fraction_flux, N_CELLS);
	Initial(x, p0, p1, p2, p3, p4, p5, mass, momentum, energy, mass_fraction, R_dry, R_v, dx, N_CELLS);
	int step = 0;
	while(t < t_FINAL){
		float dt = Compute_dt(p1, p2, p4, dx, CFL, R_dry, R_v, N_CELLS);
		Boundary(p0, p1, p2, p3, p4, N_CELLS);
		for (int i = 1; i < N_CELLS + 2; i++){
			float rho_L = p0[i-1];
			float rho_C = p0[i];
			float rho_R = p0[i+1];
			float rho_RR = p0[i+2];
			float drho_dx_L = MINMOD(rho_L, rho_C, rho_R, dx);
			float drho_dx_R = MINMOD(rho_C, rho_R, rho_RR, dx);
			float rho_L_star = rho_C + 0.5 * dx * drho_dx_L;
			float rho_R_star = rho_R - 0.5 * dx * drho_dx_R;

			float u_L = p1[i-1];
			float u_C = p1[i];
			float u_R = p1[i+1];
			float u_RR = p1[i+2];
			float du_dx_L = MINMOD(u_L, u_C, u_R, dx);
			float du_dx_R = MINMOD(u_C, u_R, u_RR, dx);
			float u_L_star = u_C + 0.5 * dx * du_dx_L;
			float u_R_star = u_R - 0.5 * dx * du_dx_R;

			float T_L = p2[i-1];
			float T_C = p2[i];
			float T_R = p2[i+1];
			float T_RR = p2[i+2];
			float dT_dx_L = MINMOD(T_L, T_C, T_R, dx);
			float dT_dx_R = MINMOD(T_C, T_R, T_RR, dx);
			float T_L_star = T_C + 0.5 * dx * dT_dx_L;
			float T_R_star = T_R - 0.5 * dx * dT_dx_R;

			float Y_L = p4[i-1];
			float Y_C = p4[i];
			float Y_R = p4[i+1];
			float Y_RR = p4[i+2];
			float dY_dx_L = MINMOD(Y_L, Y_C, Y_R, dx);
			float dY_dx_R = MINMOD(Y_C, Y_R, Y_RR, dx);
			float Y_L_star = Y_C + 0.5 * dx * dY_dx_L;
			float Y_R_star = Y_R - 0.5 * dx * dY_dx_R;

			float R_mix_L = R_dry * (1 - Y_L_star) + R_v * Y_L_star;
			float R_mix_R = R_dry * (1 - Y_R_star) + R_v * Y_R_star;
			float Cv_mix_L = Compute_Cv(T_L_star, Y_L_star);
			float Cv_mix_R = Compute_Cv(T_R_star, Y_R_star);
			float Gamma_L = 1 + R_mix_L / Cv_mix_L;
			float Gamma_R = 1 + R_mix_R / Cv_mix_R;

			float P_L_star = rho_L_star * R_mix_L * T_L_star;
			float P_R_star = rho_R_star * R_mix_R * T_R_star;

			float E_L_star = 0.5 * u_L_star * u_L_star + Cv_mix_L * T_L_star;
			float E_R_star = 0.5 * u_R_star * u_R_star + Cv_mix_R * T_R_star;
			float a_L_star = sqrt(Gamma_L * R_mix_L * T_L_star); // Sound speed a = (GAMMA*R*T)^0.5
			float a_R_star = sqrt(Gamma_R * R_mix_R * T_R_star);
			Calc_HLL_flux(rho_L_star, rho_R_star, u_L_star, u_R_star, T_L_star, T_R_star, P_L_star, P_R_star, Y_L_star, Y_R_star, E_L_star, E_R_star, a_L_star, a_R_star,
				      mass_flux, momentum_flux, energy_flux, mass_fraction_flux, i);
			float rho_face_X = 0.5 * (rho_L + rho_R);
			mass_fraction_flux[i] -= rho_face_X * D * ((Y_R - Y_C) / dx);
		}

		Calc_primitive_variable(p0, p1, p2, p3, p4, p5, mass, momentum, energy, mass_fraction,
					mass_flux, momentum_flux, energy_flux, mass_fraction_flux, R_dry, R_v, dx, dt, N_CELLS);
		t += dt;
		step++;
		if (step % 100 == 0) {
			printf("Current time = %.6f / %.2f\n", t, t_FINAL);
		}
	}

	FILE *pFile = fopen("Results_of_1000_cells.txt", "w");
	for (int i = 2; i < N_CELLS + 2; i++){
		float X = (i - 1.5) * dx;
		fprintf(pFile, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", X, p0[i], p1[i], p2[i], p3[i], p4[i], p5[i]);
	}
	fclose(pFile);

	Free_memory(&x, &p0, &p1, &p2, &p3, &p4, &p5, &mass, &momentum, &energy, &mass_fraction, &mass_flux, &momentum_flux, &energy_flux, &mass_fraction_flux);
return 0;
}

