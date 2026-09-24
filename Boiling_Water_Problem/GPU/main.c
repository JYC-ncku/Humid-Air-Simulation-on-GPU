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
	int NX = 40;
	int NY = 20;
	int N_CELLS = (NX+2) * (NY+2); // 2 Ghost cells
	 // p0: Density (rho), p1: X-velocity (u), p2: Y-velocity (v), p3: Temperature (T), p4: Pressure (p), p5: Mass fraction (Y_v), p6: Relative humidity (RH)
	float *h_p0, *h_p1, *h_p2, *h_p3, *h_p4, *h_p5, *h_p6,
	      *h_mass, *h_momentum_X, *h_momentum_Y, *h_energy, *h_mass_fraction,
	      *h_mass_flux_X, *h_momentum_X_flux_X, *h_momentum_Y_flux_X, *h_energy_flux_X, *h_mass_fraction_flux_X,
	      *h_mass_flux_Y, *h_momentum_X_flux_Y, *h_momentum_Y_flux_Y, *h_energy_flux_Y, *h_mass_fraction_flux_Y,
	      *d_p0, *d_p1, *d_p2, *d_p3, *d_p4, *d_p5, *d_p6,
	      *d_mass, *d_momentum_X, *d_momentum_Y, *d_energy, *d_mass_fraction,
	      *d_mass_flux_X, *d_momentum_X_flux_X, *d_momentum_Y_flux_X, *d_energy_flux_X, *d_mass_fraction_flux_X,
	      *d_mass_flux_Y, *d_momentum_X_flux_Y, *d_momentum_Y_flux_Y, *d_energy_flux_Y, *d_mass_fraction_flux_Y,
	      *W_GLOBAL_MAX;
	float L = 1.0; // unit: m
	float H = 0.5; // unit: m
	float t = 0;
	float t_FINAL = 5.0; // unit: s
//	float R = 1.0;
//	float GAMMA = 1.4;
	float CFL = 0.25;
	float dx = L/NX;
	float dy = H/NY;
	float D = 1.837e-5; //Diffusivity of water vapor. unit:(m^2/s)

	float R_bar = 8.3145; // unti:J/(mol*K)
	float MW_H2O = 0.01802; // unit:kg/mol
	float MW_air = 0.02897; // unit:kg/mol
	float R_v = R_bar / MW_H2O; // unit:J/(kg*k) R = R_bar / Molecular weight
	float R_dry = R_bar / MW_air;

	Allocate_memory(&h_p0, &h_p1, &h_p2, &h_p3, &h_p4, &h_p5, &h_p6,
			&h_mass, &h_momentum_X, &h_momentum_Y, &h_energy, &h_mass_fraction,
			&h_mass_flux_X, &h_momentum_X_flux_X, &h_momentum_Y_flux_X, &h_energy_flux_X, &h_mass_fraction_flux_X,
			&h_mass_flux_Y, &h_momentum_X_flux_Y, &h_momentum_Y_flux_Y, &h_energy_flux_Y, &h_mass_fraction_flux_Y,
			&d_p0, &d_p1, &d_p2, &d_p3, &d_p4, &d_p5, &d_p6,
			&d_mass, &d_momentum_X, &d_momentum_Y, &d_energy, &d_mass_fraction,
			&d_mass_flux_X, &d_momentum_X_flux_X, &d_momentum_Y_flux_X, &d_energy_flux_X, &d_mass_fraction_flux_X,
			&d_mass_flux_Y, &d_momentum_X_flux_Y, &d_momentum_Y_flux_Y, &d_energy_flux_Y, &d_mass_fraction_flux_Y,
			&W_GLOBAL_MAX);

	Initial(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6, d_mass, d_momentum_X, d_momentum_Y, d_energy, d_mass_fraction, R_dry, R_v, NX, NY, N_CELLS);
	int step = 0;
	while(t < t_FINAL){
		float W_GLOBAL_MAX = 1e-10;
		Boundary(p0, p1, p2, p3, p4, p5, p6, R_dry, R_v, NX, NY);

		Calc_Tot_Flux(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5,
			      d_mass_flux_X, d_momentum_X_flux_X, d_momentum_Y_flux_X, d_energy_flux_X, d_mass_fraction_flux_X,
			      d_mass_flux_Y, d_momentum_X_flux_Y, d_momentum_Y_flux_Y, d_energy_flux_Y, d_mass_fraction_flux_Y,
			      W_GLOBAL_MAX, R_dry, R_v, D, NX, NY, N_CELLS);

		float dt = CFL * (dx / (2.0 * W_GLOBAL_MAX)); // dx = dy

		Calc_primitive_variable(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6,
					d_mass, d_momentum_X, d_momentum_Y, d_energy, d_mass_fraction,
					d_mass_flux_X, d_momentum_X_flux_X, d_momentum_Y_flux_X, d_energy_flux_X, d_mass_fraction_flux_X,
					d_mass_flux_Y, d_momentum_X_flux_Y, d_momentum_Y_flux_Y, d_energy_flux_Y, d_mass_fraction_flux_Y,
					R_dry, R_v, dx, dy, dt, NX, NY, N_CELLS);
		t += dt;
		step++;
		if (step % 100 == 0) {
			printf("Current time = %.6f / %.2f\n", t, t_FINAL);
		}
	}

	FILE *pFile = fopen("Results_of_40x20_cells.txt", "w");
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY+2) + j;
			float X = (i - 0.5) * dx;
			float Y = (j - 0.5) * dy;
		fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", X, Y, p0[INDEX], p1[INDEX], p2[INDEX], p3[INDEX], p4[INDEX], p5[INDEX], p6[INDEX]);
		}
	}
	fclose(pFile);

	Free_memory(&h_p0, &h_p1, &h_p2, &h_p3, &h_p4, &h_p5, &h_p6,
		    &h_mass, &h_momentum_X, &h_momentum_Y, &h_energy, &h_mass_fraction,
		    &h_mass_flux_X, &h_momentum_X_flux_X, &h_momentum_Y_flux_X, &h_energy_flux_X, &h_mass_fraction_flux_X,
		    &h_mass_flux_Y, &h_momentum_X_flux_Y, &h_momentum_Y_flux_Y, &h_energy_flux_Y, &h_mass_fraction_flux_Y,
		    &d_p0, &d_p1, &d_p2, &d_p3, &d_p4, &d_p5, &d_p6,
		    &d_mass, &d_momentum_X, &d_momentum_Y, &d_energy, &d_mass_fraction,
		    &d_mass_flux_X, &d_momentum_X_flux_X, &d_momentum_Y_flux_X, &d_energy_flux_X, &d_mass_fraction_flux_X,
		    &d_mass_flux_Y, &d_momentum_X_flux_Y, &d_momentum_Y_flux_Y, &d_energy_flux_Y, &d_mass_fraction_flux_Y,
		    &W_GLOBAL_MAX);
return 0;
}

