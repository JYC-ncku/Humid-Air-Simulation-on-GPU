#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "Initial.h"
#include "Boundary.h"
#include "Calc_flux.h"
#include "Primitive_variable.h"

int main(){
	int NX = 384;
	int NY = 620;
	int N_CELLS = (NX+4) * (NY+4); // 2 Ghost cells
	 // p0: Density (rho), p1: X-velocity (u), p2: Y-velocity (v), p3: Temperature (T), p4: Pressure (p), p5: Mass fraction (Y_v), p6: Relative humidity (RH)
	float *p0, *p1, *p2, *p3, *p4, *p5, *p6,
	      *p0, *p1, *p2, *p3, *p4, *p5, *p6,
	      *mass, *momentum_X, *momentum_Y, *energy, *mass_fraction,
	      *interface_p, *flux_X, *flux_Y,
	      *MAX_Freq, *MAX_CFL;
	float L = 0.959; // unit: m
	float H = 1.55; // unit: m
	float t = 0;
	float t_FINAL = 1.0; // unit: s
//	float R = 1.0;
//	float GAMMA = 1.4;
	float CFL = 0.25;
	float dx = L/NX; // 0.0025
	float dy = H/NY; // 0.0025
	float D = 1.837e-5; // Diffusivity of water vapor. unit:(m^2/s)
	float mu = 1.83e-5; // Dynamic Viscosity unit:(Pa * s)
	float K = 0.026; // Thermal Conductivity (W/(m * K))
	float R_bar = 8.3145; // unti:J/(mol*K)
	float MW_H2O = 0.01802; // unit:kg/mol
	float MW_air = 0.02897; // unit:kg/mol
	float R_v = R_bar / MW_H2O; // unit:J/(kg*k) R = R_bar / Molecular weight
	float R_dry = R_bar / MW_air;

	Allocate_memory(&p0, &p1, &p2, &p3, &p4, &p5, &p6,
			&p0, &p1, &p2, &p3, &p4, &p5, &p6,
			&mass, &momentum_X, &momentum_Y, &energy, &mass_fraction,
			&interface_p, &flux_X, &flux_Y,
			&MAX_Freq, &MAX_CFL, N_CELLS);

	Initial(p0, p1, p2, p3, p4, p5, p6, mass, momentum_X, momentum_Y, energy, mass_fraction, MAX_Freq, MAX_CFL, R_dry, R_v, NX, NY, N_CELLS);
	int step = 0;

	float output_interval = 0.05; // 每 0.05 秒輸出一次
	float next_output_time = 0.0; // 下一次輸出的時間點
	int frame_count = 0;          // 影格編號

	while(t < t_FINAL){
		float dt = Compute_dt(MAX_Freq, MAX_CFL, p0, p1, p2, p3, p5, dx, dy, R_dry, R_v, CFL, NX, NY, N_CELLS);

		Boundary(p0, p1, p2, p3, p4, p5, p6, R_dry, R_v, dx, dy, NX, NY, N_CELLS);

		Calc_Tot_Flux(p0, p1, p2, p3, p4, p5,
			      interface_p, flux_X, flux_Y,
			      R_dry, R_v, D, mu, K, dx, dy, NX, NY, N_CELLS);

		Calc_primitive_variable(p0, p1, p2, p3, p4, p5, p6,
					mass, momentum_X, momentum_Y, energy, mass_fraction,
					flux_X, flux_Y,
					R_dry, R_v, dx, dy, dt, NX, NY, N_CELLS);

	if (t >= next_output_time) {
		printf("Saving Frame %03d at t = %.3f...\n", frame_count, t);

		char filename[100];
		sprintf(filename, "Results_frame_%03d.txt", frame_count);

		FILE *pFile = fopen(filename, "w");
		for (int i = 2; i < NX + 2; i++){
			for (int j = 2; j < NY + 2; j++){
				int INDEX = i * (NY+4) + j;
				float X = (i - 1.5) * dx;
				float Y = (j - 1.5) * dy;
				fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\n", X, Y, p1[INDEX], p3[INDEX], p5[INDEX], p6[INDEX]);
			}
		}
		fclose(pFile);

		frame_count++;
		next_output_time += output_interval; // 設定下一次拍照時間
	}

		t += dt;
		step++;
		if (step % 100 == 0) {
			printf("Current time = %.6f / %.2f\n", t, t_FINAL);
		}
	}
/*
	FILE *pFile = fopen("Results_of_400x200_cells_1s.txt", "w");
	for (int i = 2; i < NX + 2; i++){
		for (int j = 2; j < NY + 2; j++){
			int INDEX = i * (NY+4) + j;
			float X = (i - 1.5) * dx;
			float Y = (j - 1.5) * dy;
		fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
				X, Y, p0[INDEX], p1[INDEX], p2[INDEX], p3[INDEX], p4[INDEX], p5[INDEX], p6[INDEX]);
		}
	}
	fclose(pFile);
*/
	Free_memory(&p0, &p1, &p2, &p3, &p4, &p5, &p6,
		    &p0, &p1, &p2, &p3, &p4, &p5, &p6,
		    &mass, &momentum_X, &momentum_Y, &energy, &mass_fraction,
		    &interface_p, &flux_X, &flux_Y,
		    &MAX_Freq, &MAX_CFL);
return 0;
}
