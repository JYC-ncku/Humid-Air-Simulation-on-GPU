#include <stdlib.h>
#include <stdio.h>
#include "memory.h"
#include "Initial.h"

int main(){
	int N_CELLS = 200;
	float *h_x, *h_p0, *h_p1, *h_p2, *h_p3, *h_mass, *h_momentum, *h_energy,
	      *d_x, *d_p0, *d_p1, *d_p2, *d_p3, *d_mass, *d_momentum, *d_energy,
	      *d_mass_flux, *d_momentum_flux, *d_energy_flux; // p0 is density, p1 is velocity, p2 is temperature, p3 is pressure, h mean host, d mean device.
	float L = 1.0;
	float t = 0;
	float t_FINAL = 0.2;
	float R = 1.0;
	float GAMMA = 1.4;
	float CFL = 0.5;
	float dx = L/N_CELLS;
	Allocate_memory(&h_x, &h_p0, &h_p1, &h_p2, &h_p3, &h_mass, &h_momentum, &h_energy,
			&d_x, &d_p0, &d_p1, &d_p2, &d_p3, &d_mass, &d_momentum, &d_energy, &d_mass_flux, &d_momentum_flux, &d_energy_flux, N_CELLS);
	Initial(h_x, h_p0, h_p1, h_p2, h_p3, h_mass, h_momentum, h_energy, GAMMA, dx, N_CELLS);
//	Send_To_Device();

	printf("Hello world!\n");
	Free_memory(&h_x, &h_p0, &h_p1, &h_p2, &h_p3, &h_mass, &h_momentum, &h_energy,
		    &d_x, &d_p0, &d_p1, &d_p2, &d_p3, &d_mass, &d_momentum, &d_energy, &d_mass_flux, &d_momentum_flux, &d_energy_flux);
}
