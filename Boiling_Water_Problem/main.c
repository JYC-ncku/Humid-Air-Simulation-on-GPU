#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "Initial.h"
#include "Boundary.h"
#include "Calc_flux.h"
#include "Primitive_variable.h"

float MAX_Wave_Speed(float u_L, float u_R, float a_L, float a_R){
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

int main(){
	int NX = 400;
	int NY = 200;
	int N_CELLS = (NX+2) * (NY+2); // 2 Ghost cells
	 // p0: Density (rho), p1: X-velocity (u), p2: Y-velocity (v), p3: Temperature (T), p4: Pressure (p), p5: Mass fraction (Y_v), p6: Relative humidity (RH)
	float *p0, *p1, *p2, *p3, *p4, *p5, *p6,
	      *mass, *momentum_X, *momentum_Y, *energy, *mass_fraction,
	      *mass_flux_X, *momentum_X_flux_X, *momentum_Y_flux_X, *energy_flux_X, *mass_fraction_flux_X,
	      *mass_flux_Y, *momentum_X_flux_Y, *momentum_Y_flux_Y, *energy_flux_Y, *mass_fraction_flux_Y;
	float L = 1.0; // unit: m
	float H = 0.5; // unit: m
	float t = 0;
	float t_FINAL = 5.0; // unit: s
//	float R = 1.0;
//	float GAMMA = 1.4;
	float CFL = 0.5;
	float dx = L/NX;
	float dy = H/NY;
	float D = 1.837e-5; //Diffusivity of water vapor. unit:(m^2/s)

	float R_bar = 8.3145; // unti:J/(mol*K)
	float MW_H2O = 0.01802; // unit:kg/mol
	float MW_air = 0.02897; // unit:kg/mol
	float R_v = R_bar / MW_H2O; // unit:J/(kg*k) R = R_bar / Molecular weight
	float R_dry = R_bar / MW_air;

	Allocate_memory(&p0, &p1, &p2, &p3, &p4, &p5, &p6,
			&mass, &momentum_X, &momentum_Y, &energy, &mass_fraction,
			&mass_flux_X, &momentum_X_flux_X, &momentum_Y_flux_X, &energy_flux_X, &mass_fraction_flux_X,
			&mass_flux_Y, &momentum_X_flux_Y, &momentum_Y_flux_Y, &energy_flux_Y, &mass_fraction_flux_Y,
			N_CELLS);
	Initial(p0, p1, p2, p3, p4, p5, p6, mass, momentum_X, momentum_Y, energy, mass_fraction, R_dry, R_v, NX, NY);
	int step = 0;
	while(t < t_FINAL){
		float W_GLOBAL_MAX = 1e-10;
		Boundary(p0, p1, p2, p3, p4, p5, p6, R_dry, R_v, NX, NY);
		for (int i = 1; i < NX + 2; i++){
			for (int j = 1; j < NY + 2; j ++){
				int INDEX_L = (i-1) * (NY+2) + j;
				int INDEX = i * (NY+2) + j;
				float rho_L = p0[INDEX_L];
				float rho_R = p0[INDEX];
				float u_L = p1[INDEX_L];
				float u_R = p1[INDEX];
				float v_L = p2[INDEX_L];
				float v_R = p2[INDEX];
				float T_L = p3[INDEX_L];
				float T_R = p3[INDEX];
				float P_L = p4[INDEX_L];
				float P_R = p4[INDEX];
				float Y_L = p5[INDEX_L];
				float Y_R = p5[INDEX];
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
				Calc_HLL_X_flux(rho_L, rho_R, u_L, u_R, v_L, v_R, T_L, T_R, P_L, P_R, Y_L, Y_R, E_L, E_R, a_L, a_R,
					        mass_flux_X, momentum_X_flux_X, momentum_Y_flux_X, energy_flux_X, mass_fraction_flux_X, INDEX);
				mass_fraction_flux_X[INDEX] -= D * ((Y_R - Y_L) / dx);
				if (W_LOCAL_MAX > W_GLOBAL_MAX){
					W_GLOBAL_MAX = W_LOCAL_MAX;
				}
			}
		}
		for (int i = 1; i < NX + 2; i++){
			for (int j = 1; j < NY + 2; j ++){
				int INDEX_B = i * (NY+2) + (j-1);
				int INDEX = i * (NY+2) + j;
				float rho_B = p0[INDEX_B];
				float rho_T = p0[INDEX];
				float u_B = p1[INDEX_B];
				float u_T = p1[INDEX];
				float v_B = p2[INDEX_B];
				float v_T = p2[INDEX];
				float T_B = p3[INDEX_B];
				float T_T = p3[INDEX];
				float P_B = p4[INDEX_B];
				float P_T = p4[INDEX];
				float Y_B = p5[INDEX_B];
				float Y_T = p5[INDEX];
				float R_mix_B = R_dry * (1 - Y_B) + R_v * Y_B;
				float R_mix_T = R_dry * (1 - Y_T) + R_v * Y_T;
				float Cv_mix_B = Compute_Cv(T_B, Y_B);
				float Cv_mix_T = Compute_Cv(T_T, Y_T);
				float Gamma_B = 1 + R_mix_B / Cv_mix_B;
				float Gamma_T = 1 + R_mix_T / Cv_mix_T;
				float E_B = 0.5 * u_B * u_B + Cv_mix_B * T_B;
				float E_T = 0.5 * u_T * u_T + Cv_mix_T * T_T;
				float a_B = sqrt(Gamma_B * R_mix_B * T_B); // Sound speed a = (R*T)^0.5
				float a_T = sqrt(Gamma_T * R_mix_T * T_T);
				float W_LOCAL_MAX = MAX_Wave_Speed(v_B, v_T, a_B, a_T);
				Calc_HLL_Y_flux(rho_B, rho_T, u_B, u_T, v_B, v_T, T_B, T_T, P_B, P_T, Y_B, Y_T, E_B, E_T, a_B, a_T,
					        mass_flux_Y, momentum_X_flux_Y, momentum_Y_flux_Y, energy_flux_Y, mass_fraction_flux_Y, INDEX);
				mass_fraction_flux_Y[INDEX] -= D * ((Y_T - Y_B) / dy);
				if (W_LOCAL_MAX > W_GLOBAL_MAX){
					W_GLOBAL_MAX = W_LOCAL_MAX;
				}
			}
		}
		float dt = CFL * (dx / W_GLOBAL_MAX); // dx = dy

		Calc_primitive_variable(p0, p1, p2, p3, p4, p5, p6,
					mass, momentum_X, momentum_Y, energy, mass_fraction,
					mass_flux_X, momentum_X_flux_X, momentum_Y_flux_X, energy_flux_X, mass_fraction_flux_X,
					mass_flux_Y, momentum_Y_flux_Y, momentum_Y_flux_Y, energy_flux_Y, mass_fraction_flux_Y,
					R_dry, R_v, dx, dy, dt, NX, NY);
		t += dt;
		step++;
		if (step % 100 == 0) {
			printf("Current time = %.6f / %.2f\n", t, t_FINAL);
		}
	}

	FILE *pFile = fopen("Results_of_400x200_cells.txt", "w");
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY+2) + j;
			float X = (i - 0.5) * dx;
			float Y = (j - 0.5) * dy;
		fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", X, Y, p0[INDEX], p1[INDEX], p2[INDEX], p3[INDEX], p4[INDEX], p5[INDEX]);
		}
	}
	fclose(pFile);

	Free_memory(&p0, &p1, &p2, &p3, &p4, &p5, &p6, &mass, &momentum_X, &momentum_Y, &energy, &mass_fraction,
		    &mass_flux_X, &momentum_X_flux_X, &momentum_Y_flux_X, &energy_flux_X, &mass_fraction_flux_X,
		    &mass_flux_Y, &momentum_X_flux_Y, &momentum_Y_flux_Y, &energy_flux_Y, &mass_fraction_flux_Y);
return 0;
}

