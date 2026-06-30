#include <stdlib.h>
#include <stdio.h>
#include "memory.h"
#include "Initial.h"

int main(){
	int N_CELLS = 200;
	float *x, *p0, *p1, *p2, *p3, *mass, *momentum, *energy, *mass_flux, *momentum_flux, *energy_flux;
	float L = 1.0;
	float t = 0;
	float t_FINAL = 0.2;
	float R = 1.0;
	float GAMMA = 1.4;
	float CFL = 0.25;
	float dx = L/N_CELLS;
	float W_GLOBAL_MAX;
	Allocate_memory(&x, &p0, &p1, &p2, &p3, &mass, &momentum, &energy, &mass_flux, &momentum_flux, &energy_flux, N_CELLS);
	Initial(x, p0, p1, p2, p3, mass, momentum, energy, dx, GAMMA, N_CELLS);
	for (int i = 0; i < N_CELLS; i++){
		printf("p0[%d] = %g\n", i, p0[i]);
	}
	printf("Hello world\n");
	Free_memory(&x, &p0, &p1, &p2, &p3, &mass, &momentum, &energy, &mass_flux, &momentum_flux, &energy_flux);
return 0;
}

