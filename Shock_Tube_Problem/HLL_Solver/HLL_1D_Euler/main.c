#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "Initial.h"
#include "Calc_flux.h"

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
	float *x, *p0, *p1, *p2, *p3, *mass, *momentum, *energy, *mass_flux, *momentum_flux, *energy_flux;
	float L = 1.0;
	float t = 0;
	float t_FINAL = 0.2;
	float R = 1.0;
	float GAMMA = 1.4;
	float CFL = 0.5;
	float dx = L/N_CELLS;
	float W_GLOBAL_MAX = 1e-10;
	Allocate_memory(&x, &p0, &p1, &p2, &p3, &mass, &momentum, &energy, &mass_flux, &momentum_flux, &energy_flux, N_CELLS);
	Initial(x, p0, p1, p2, p3, mass, momentum, energy, dx, GAMMA, N_CELLS);
	while(t < t_FINAL){
		for (int i = 1; i < N_CELLS; i++){
			float rho_L = p0[i-1];
			float rho_R = p0[i];
			float u_L = p1[i-1];
			float u_R = p1[i];
			float T_L = p2[i-1];
			float T_R = p2[i];
			float P_L = p3[i-1];
			float P_R = p3[i];
			float e_L = 0.5 * rho_L * u_L * u_L + (P_L / (GAMMA - 1));
			float e_R = 0.5 * rho_R * u_R * u_R + (P_R / (GAMMA - 1));
			float a_L = sqrt(GAMMA * R * T_L); // Sound speed a = (R*T)^0.5
			float a_R = sqrt(GAMMA * R * T_R);
			float W_LOCAL_MAX = MAX_Wave_Speed(u_L, u_R, a_L, a_R);
			Calc_HLL_flux(rho_L, rho_R, u_L, u_R, T_L, T_R, P_L, P_R, e_L, e_R, a_L, a_R,
				      mass_flux, momentum_flux, energy_flux, i);
			if (W_LOCAL_MAX > W_GLOBAL_MAX){
				W_GLOBAL_MAX = W_LOCAL_MAX;
			}
		}

		//Boundary condition
		mass_flux[0] = mass_flux[1];
		momentum_flux[0] = momentum_flux[1];
		energy_flux[0] = energy_flux[1];
		mass_flux[N_CELLS] = mass_flux[N_CELLS-1];
		momentum_flux[N_CELLS] = momentum_flux[N_CELLS-1];
		energy_flux[N_CELLS] = energy_flux[N_CELLS-1];

		float dt = CFL * (dx / W_GLOBAL_MAX);

		for (int i = 0; i < N_CELLS; i++){
		        // Use FVM to get new conservation values
			mass[i] = mass[i] - (dt / dx)*(mass_flux[i+1] - mass_flux[i]);
			momentum[i] = momentum[i] - (dt / dx)*(momentum_flux[i+1] - momentum_flux[i]);
			energy[i] = energy[i] - (dt / dx)*(energy_flux[i+1] - energy_flux[i]);
			//Get new variable
			p0[i] = mass[i];
			p1[i] = momentum[i] / mass[i];
			p3[i] = (GAMMA - 1) * (energy[i] - 0.5 * p0[i] * p1[i] * p1[i]);
			p2[i] = p3[i] / (p0[i] * R);
		}

		t += dt;
	}

	FILE *pFile = fopen("Results_of_800_cells.txt", "w");
	for (int i = 0; i < N_CELLS; i++){
		fprintf(pFile, "%g\t%g\t%g\t%g\t%g\n", x[i], p0[i], p1[i], p2[i], p3[i]);
	}
	fclose(pFile);

	Free_memory(&x, &p0, &p1, &p2, &p3, &mass, &momentum, &energy, &mass_flux, &momentum_flux, &energy_flux);
return 0;
}

