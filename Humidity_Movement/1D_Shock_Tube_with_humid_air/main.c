#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "Initial.h"
#include "Boundary.h"
#include "Calc_flux.h"
#include "Primitive_variable.h"

float MAX_Wave_Speed(double u_L, double u_R, double a_L, double a_R){
    double W_L = fabs(u_L) + a_L;
    double W_R = fabs(u_R) + a_R;
    double W_LOCAL_MAX;
    if (W_L > W_R){
        W_LOCAL_MAX = W_L;
    }else {
        W_LOCAL_MAX = W_R;
    }
    return W_LOCAL_MAX;
}

int main(){
	int N_CELLS = 800;
	float *x, *p0, *p1, *p2, *p3, *p4, *p5, *P_sat, *P_v,
	      *mass, *momentum, *energy, *mass_fraction, *mass_flux, *momentum_flux, *energy_flux, *mass_fraction_flux;
	float L = 0.01; // unit: m
	float t = 0;
	float t_FINAL = 0.2;
//	float R = 1.0;
//	float GAMMA = 1.4;
	float CFL = 0.5;
	float dx = L/N_CELLS;

	float D = 1.837e-5; //Diffusivity of water vapor. unit:(m^2/s)

	float R_bar = 8.3145; // unti:J/(mol*K)
	float MW_H2O = 0.01802; // unit:kg/mol
	float MW_air = 0.02897; // unit:kg/mol
	float R_v = R_bar / MW_H2O; // unit:J/(kg*k) R = R_bar / Molecular weight
	float R_dry = R_bar / MW_air;

	Allocate_memory(&x, &p0, &p1, &p2, &p3, &p4, &p5, &P_sat, &P_v,
			&mass, &momentum, &energy, &mass_fraction, &mass_flux, &momentum_flux, &energy_flux, &mass_fraction_flux, N_CELLS);
	Initial(x, p0, p1, p2, p3, p4, p5, P_sat, P_v, mass, momentum, energy, mass_fraction, dx, R_dry, R_v, N_CELLS);
	int step = 0;
	while(t < t_FINAL){
		float W_GLOBAL_MAX = 1e-10;
		Boundary(p0, p1, p2, p3, p4, N_CELLS);
		for (int i = 1; i < N_CELLS + 2; i++){
			float rho_L = p0[i-1];
			float rho_R = p0[i];
			float u_L = p1[i-1];
			float u_R = p1[i];
			float T_L = p2[i-1];
			float T_R = p2[i];
			float P_L = p3[i-1];
			float P_R = p3[i];
			float Y_L = p4[i-1];
			float Y_R = p4[i];
			float R_mix_L = R_dry * (1 - Y_L) + R_v * Y_L;
			float R_mix_R = R_dry * (1 - Y_R) + R_v * Y_R;
			float Cv_mix_L = Compute_Cv(T_L, Y_L);
			float Cv_mix_R = Compute_Cv(T_R, Y_R);
			float Gamma_L = 1 + R_mix_L / Cv_mix_L;
			float Gamma_R = 1 + R_mix_R / Cv_mix_R;
			float E_L = 0.5 * u_L * u_L + Cv_mix_L * T_L;
			float E_R = 0.5 * u_R * u_R + Cv_mix_R * T_R;
			float a_L = sqrt(Gamma_L * R_mix_L * T_L); // Sound speed a = (R*T)^0.5
			float a_R = sqrt(Gamma_R * R_mix_R * T_R);
			float W_LOCAL_MAX = MAX_Wave_Speed(u_L, u_R, a_L, a_R);
			Calc_HLL_flux(rho_L, rho_R, u_L, u_R, T_L, T_R, P_L, P_R, Y_L, Y_R, E_L, E_R, a_L, a_R,
				      mass_flux, momentum_flux, energy_flux, mass_fraction_flux, i);
			mass_fraction_flux[i] -= D * ((Y_R - Y_L) / dx);
			if (W_LOCAL_MAX > W_GLOBAL_MAX){
				W_GLOBAL_MAX = W_LOCAL_MAX;
			}
		}

		float dt = CFL * (dx / W_GLOBAL_MAX);

		Calc_primitive_variable(p0, p1, p2, p3, p4, p5, P_sat, P_v, mass, momentum, energy, mass_fraction,
					mass_flux, momentum_flux, energy_flux, mass_fraction_flux, R_dry, R_v, dx, dt, N_CELLS);
		t += dt;
		step++;
		if (step % 100 == 0) {
			printf("Current time = %.6f / %.2f\n", t, t_FINAL);
		}
	}

	FILE *pFile = fopen("Results_of_800_cells.txt", "w");
	for (int i = 1; i < N_CELLS + 1; i++){
		float X = (i - 0.5) * dx;
		fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", X, p0[i], p1[i], p2[i], p3[i], p4[i], p5[i]);
	}
	fclose(pFile);

	Free_memory(&x, &p0, &p1, &p2, &p3, &p4, &p5, &P_sat, &P_v, &mass, &momentum, &energy, &mass_fraction, &mass_flux, &momentum_flux, &energy_flux, &mass_fraction_flux);
return 0;
}

