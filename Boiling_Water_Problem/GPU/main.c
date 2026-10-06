#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "Initial.h"
#include "Boundary.h"
#include "Calc_flux.h"
#include "Primitive_variable.h"

int main(){
	int NX = 400;
	int NY = 200;
	int N_CELLS = (NX+2) * (NY+2); // 2 Ghost cells
	 // p0: Density (rho), p1: X-velocity (u), p2: Y-velocity (v), p3: Temperature (T), p4: Pressure (p), p5: Mass fraction (Y_v), p6: Relative humidity (RH)
	float *h_p0, *h_p1, *h_p2, *h_p3, *h_p4, *h_p5, *h_p6,
	      *d_p0, *d_p1, *d_p2, *d_p3, *d_p4, *d_p5, *d_p6,
	      *d_mass, *d_momentum_X, *d_momentum_Y, *d_energy, *d_mass_fraction,
	      *d_mass_flux_X, *d_momentum_X_flux_X, *d_momentum_Y_flux_X, *d_energy_flux_X, *d_mass_fraction_flux_X,
	      *d_mass_flux_Y, *d_momentum_X_flux_Y, *d_momentum_Y_flux_Y, *d_energy_flux_Y, *d_mass_fraction_flux_Y,
	      *MAX_Freq, *d_MAX_CFL;
	float L = 1.0; // unit: m
	float H = 0.5; // unit: m
	float t = 0;
	float t_FINAL = 1.0; // unit: s
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
			&d_p0, &d_p1, &d_p2, &d_p3, &d_p4, &d_p5, &d_p6,
			&d_mass, &d_momentum_X, &d_momentum_Y, &d_energy, &d_mass_fraction,
			&d_mass_flux_X, &d_momentum_X_flux_X, &d_momentum_Y_flux_X, &d_energy_flux_X, &d_mass_fraction_flux_X,
			&d_mass_flux_Y, &d_momentum_X_flux_Y, &d_momentum_Y_flux_Y, &d_energy_flux_Y, &d_mass_fraction_flux_Y,
			&MAX_Freq, &d_MAX_CFL, N_CELLS);

	Initial(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6, d_mass, d_momentum_X, d_momentum_Y, d_energy, d_mass_fraction, MAX_Freq, R_dry, R_v, NX, NY, N_CELLS);
	int step = 0;
	while(t < t_FINAL){
		float dt = Compute_dt(MAX_Freq, d_MAX_CFL, d_p0, d_p1, d_p2, d_p3, d_p5, dx, dy, R_dry, R_v, NX, NY, N_CELLS);

		Boundary(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6, R_dry, R_v, dx, dy, NX, NY, N_CELLS);

		Calc_Tot_Flux(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5,
			      d_mass_flux_X, d_momentum_X_flux_X, d_momentum_Y_flux_X, d_energy_flux_X, d_mass_fraction_flux_X,
			      d_mass_flux_Y, d_momentum_X_flux_Y, d_momentum_Y_flux_Y, d_energy_flux_Y, d_mass_fraction_flux_Y,
			      R_dry, R_v, D, dx, dy, NX, NY, N_CELLS);

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
	//Get the data from device
	Get_From_Device(&h_p0, &d_p0, N_CELLS);
	Get_From_Device(&h_p1, &d_p1, N_CELLS);
	Get_From_Device(&h_p2, &d_p2, N_CELLS);
	Get_From_Device(&h_p3, &d_p3, N_CELLS);
	Get_From_Device(&h_p4, &d_p4, N_CELLS);
	Get_From_Device(&h_p5, &d_p5, N_CELLS);
	Get_From_Device(&h_p6, &d_p6, N_CELLS);

	FILE *pFile = fopen("Results_of_400x200_cells_1s.txt", "w");
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY+2) + j;
			float X = (i - 0.5) * dx;
			float Y = (j - 0.5) * dy;
		fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
				X, Y, h_p0[INDEX], h_p1[INDEX], h_p2[INDEX], h_p3[INDEX], h_p4[INDEX], h_p5[INDEX], h_p6[INDEX]);
		}
	}
	fclose(pFile);

	Free_memory(&h_p0, &h_p1, &h_p2, &h_p3, &h_p4, &h_p5, &h_p6,
		    &d_p0, &d_p1, &d_p2, &d_p3, &d_p4, &d_p5, &d_p6,
		    &d_mass, &d_momentum_X, &d_momentum_Y, &d_energy, &d_mass_fraction,
		    &d_mass_flux_X, &d_momentum_X_flux_X, &d_momentum_Y_flux_X, &d_energy_flux_X, &d_mass_fraction_flux_X,
		    &d_mass_flux_Y, &d_momentum_X_flux_Y, &d_momentum_Y_flux_Y, &d_energy_flux_Y, &d_mass_fraction_flux_Y,
		    &MAX_Freq, &d_MAX_CFL);
return 0;
}

