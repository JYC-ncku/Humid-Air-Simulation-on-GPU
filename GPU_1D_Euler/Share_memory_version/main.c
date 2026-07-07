#include <stdlib.h>
#include <stdio.h>
#include "memory.h"
#include "Initial.h"
#include "GPU_calc_flux.h"
#include "Calc_conserved.h"
#include "Calc_variable.h"
#include "Boundary.h"

int main(){
	int N_CELLS = 200000;
	float *h_x, *h_p0, *h_p1, *h_p2, *h_p3, *h_mass, *h_momentum, *h_energy, *h_block_max,
	      *d_x, *d_p0, *d_p1, *d_p2, *d_p3, *d_mass, *d_momentum, *d_energy, *d_block_max,
	      *d_mass_flux, *d_momentum_flux, *d_energy_flux;
	      // p0 is density, p1 is velocity, p2 is temperature, p3 is pressure, h mean host, d mean device.
	float L = 1.0;
	float t = 0;
	float t_FINAL = 0.2;
	float R = 1.0;
	float GAMMA = 1.4;
	float CFL = 0.5;
	float dx = L/N_CELLS;
	float W_GLOBAL_MAX = 1e-10;
	int TPB = 128;
	int BPG=  (N_CELLS + TPB - 1) / TPB;
	Allocate_memory(&h_x, &h_p0, &h_p1, &h_p2, &h_p3, &h_mass, &h_momentum, &h_energy, &h_block_max,
			&d_p0, &d_p1, &d_p2, &d_p3, &d_mass, &d_momentum, &d_energy, &d_block_max, &d_mass_flux, &d_momentum_flux, &d_energy_flux, N_CELLS);

	Initial(h_x, h_p0, h_p1, h_p2, h_p3, h_mass, h_momentum, h_energy, GAMMA, dx, N_CELLS);
	Send_To_Device(&d_p0, &h_p0, N_CELLS);
	Send_To_Device(&d_p1, &h_p1, N_CELLS);
	Send_To_Device(&d_p2, &h_p2, N_CELLS);
	Send_To_Device(&d_p3, &h_p3, N_CELLS);
	Send_To_Device(&d_mass, &h_mass, N_CELLS);
	Send_To_Device(&d_momentum, &h_momentum, N_CELLS);
	Send_To_Device(&d_energy, &h_energy, N_CELLS);

	while (t < t_FINAL){
		Calc_flux(d_p0, d_p1, d_p2, d_p3, d_mass_flux, d_momentum_flux, d_energy_flux, d_block_max, GAMMA, R, N_CELLS);

		Get_From_Device_for_block(&h_block_max, &d_block_max, N_CELLS);
		for (int i = 0; i < BPG; i++){
			if (h_block_max[i] > W_GLOBAL_MAX){
				W_GLOBAL_MAX = h_block_max[i];
			}
		}
		float dt = CFL * (dx / W_GLOBAL_MAX);

		Boundary(d_mass_flux, d_momentum_flux, d_energy_flux, N_CELLS);
		Calc_conserved(d_mass, d_momentum, d_energy, d_mass_flux, d_momentum_flux, d_energy_flux, dt, dx, N_CELLS);
		Calc_variable(d_p0, d_p1, d_p2, d_p3, d_mass, d_momentum, d_energy, GAMMA, R, N_CELLS);

		t += dt;
	}

	Get_From_Device(&h_p0, &d_p0, N_CELLS);
	Get_From_Device(&h_p1, &d_p1, N_CELLS);
	Get_From_Device(&h_p2, &d_p2, N_CELLS);
	Get_From_Device(&h_p3, &d_p3, N_CELLS);

	FILE *pFile = fopen("Reuslts_of_200000_cells.txt", "w");
	for (int i = 0; i < N_CELLS; i++){
		fprintf(pFile, "%g\t%g\t%g\t%g\t%g\n", h_x[i], h_p0[i], h_p1[i], h_p2[i], h_p3[i]);
	}
	fclose(pFile);

	Free_memory(&h_x, &h_p0, &h_p1, &h_p2, &h_p3, &h_mass, &h_momentum, &h_energy, &h_block_max,
		    &d_p0, &d_p1, &d_p2, &d_p3, &d_mass, &d_momentum, &d_energy, &d_block_max, &d_mass_flux, &d_momentum_flux, &d_energy_flux);
}
