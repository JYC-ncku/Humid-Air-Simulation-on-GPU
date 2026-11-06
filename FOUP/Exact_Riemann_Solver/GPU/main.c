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
	int N_CELLS = (NX+2) * (NY+2); // 1 Ghost cells
	 // p0: Density (rho), p1: X-velocity (u), p2: Y-velocity (v), p3: Temperature (T), p4: Pressure (p), p5: Mass fraction (Y_v), p6: Relative humidity (RH)
	float *h_p0, *h_p1, *h_p2, *h_p3, *h_p4, *h_p5, *h_p6, *h_B,
	      *d_p0, *d_p1, *d_p2, *d_p3, *d_p4, *d_p5, *d_p6, *d_B,
	      *d_mass, *d_momentum_X, *d_momentum_Y, *d_energy, *d_mass_fraction,
	      *d_interface_p, *d_flux_X, *d_flux_Y,
	      *MAX_Freq, *d_MAX_CFL;
	float L = 0.959; // unit: m
	float H = 1.55; // unit: m
	float t = 0;
	float t_FINAL = 10.0; // unit: s
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

	Allocate_memory(&h_p0, &h_p1, &h_p2, &h_p3, &h_p4, &h_p5, &h_p6, &h_B,
			&d_p0, &d_p1, &d_p2, &d_p3, &d_p4, &d_p5, &d_p6, &d_B,
			&d_mass, &d_momentum_X, &d_momentum_Y, &d_energy, &d_mass_fraction,
			&d_interface_p, &d_flux_X, &d_flux_Y,
			&MAX_Freq, &d_MAX_CFL, N_CELLS);

	Initial(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6, d_B, d_mass, d_momentum_X, d_momentum_Y, d_energy, d_mass_fraction, MAX_Freq, d_MAX_CFL, R_dry, R_v, dx, dy, NX, NY, N_CELLS);

	int step = 0;

	float output_interval = 1.0; // 每 1.0 秒輸出一次
	float next_output_time = 0.0; // 下一次輸出的時間點
	int frame_count = 0;          // 影格編號

	while(t < t_FINAL){
		float dt = Compute_dt(MAX_Freq, d_MAX_CFL, d_p0, d_p1, d_p2, d_p3, d_p5, dx, dy, R_dry, R_v, CFL, NX, NY, N_CELLS);

		Boundary(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6, R_dry, R_v, dx, dy, NX, NY, N_CELLS);

		Calc_Tot_Flux(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_B,
			      d_interface_p, d_flux_X, d_flux_Y,
			      R_dry, R_v, D, mu, K, dx, dy, NX, NY, N_CELLS);

		Calc_primitive_variable(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6,
					d_mass, d_momentum_X, d_momentum_Y, d_energy, d_mass_fraction,
					d_flux_X, d_flux_Y, d_B,
					R_dry, R_v, dx, dy, dt, NX, NY, N_CELLS);

	// For make animations
	if (t >= next_output_time) {
		printf("Saving Frame %03d at t = %.3f...\n", frame_count, t);

		//Get the data from device
		Get_From_Device(&h_p0, &d_p0, N_CELLS);
		Get_From_Device(&h_p1, &d_p1, N_CELLS);
		Get_From_Device(&h_p2, &d_p2, N_CELLS);
		Get_From_Device(&h_p3, &d_p3, N_CELLS);
		Get_From_Device(&h_p4, &d_p4, N_CELLS);
		Get_From_Device(&h_p5, &d_p5, N_CELLS);
		Get_From_Device(&h_p6, &d_p6, N_CELLS);
		Get_From_Device(&h_B, &d_B, N_CELLS);

		char filename[100];
		sprintf(filename, "Results_frame_%03d.txt", frame_count);

		FILE *pFile = fopen(filename, "w");
		for (int i = 1; i < NX + 1; i++){
			for (int j = 1; j < NY + 1; j++){
				int INDEX = i * (NY+2) + j;
				float X = (i - 0.5) * dx;
				float Y = (j - 0.5) * dy;
				fprintf(pFile, "%.4f\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", X, Y, h_p0[INDEX], h_p1[INDEX], h_p2[INDEX], h_p3[INDEX],
													       h_p4[INDEX], h_p5[INDEX], h_p6[INDEX], h_B[INDEX]);
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
	//Get the data from device
	Get_From_Device(&h_p0, &d_p0, N_CELLS);
	Get_From_Device(&h_p1, &d_p1, N_CELLS);
	Get_From_Device(&h_p2, &d_p2, N_CELLS);
	Get_From_Device(&h_p3, &d_p3, N_CELLS);
	Get_From_Device(&h_p4, &d_p4, N_CELLS);
	Get_From_Device(&h_p5, &d_p5, N_CELLS);
	Get_From_Device(&h_p6, &d_p6, N_CELLS);
	Get_From_Device(&h_B, &d_B, N_CELLS);

	FILE *pFile = fopen("Results_of_96x155_cells_10s.txt", "w");
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY+2) + j;
			float X = (i - 0.5) * dx;
			float Y = (j - 0.5) * dy;
			fprintf(pFile, "%.4f\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
					X, Y, h_p0[INDEX], h_p1[INDEX], h_p2[INDEX], h_p3[INDEX], h_p4[INDEX], h_p5[INDEX], h_p6[INDEX], h_B[INDEX]);
		}
	}
	fclose(pFile);
*/
	Free_memory(&h_p0, &h_p1, &h_p2, &h_p3, &h_p4, &h_p5, &h_p6, &h_B,
		    &d_p0, &d_p1, &d_p2, &d_p3, &d_p4, &d_p5, &d_p6, &d_B,
		    &d_mass, &d_momentum_X, &d_momentum_Y, &d_energy, &d_mass_fraction,
		    &d_interface_p, &d_flux_X, &d_flux_Y,
		    &MAX_Freq, &d_MAX_CFL);
return 0;
}
