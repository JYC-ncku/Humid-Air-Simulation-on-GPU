#include <stdlib.h>
#include <stdio.h>
#include "memory.h"

int main(){
	int N_CELLS = 200;
	float *x, *p0, *p1, *p2, *p3, *mass, *momentum, *energy, *mass_flux, *momentum_flux, *energy_flux; // p0 is density, p1 is velocity, p2 is temperature, p3 is pressure
	float L = 1.0;
	float t = 0;
	float t_FINAL = 0.2;
	float R = 1.0;
	float GAMMA = 1.4;
	float CFL = 0.5;
	float dx = L/N_CELLS;
	Allocate_memory(&x, &p0, &p1, &p2, &p3, &mass, &momentum, &energy, &mass_flux, &momentum_flux, &energy_flux, N_CELLS);
	printf("Hello world!\n");
	Free_memory(&x, &p0, &p1, &p2, &p3, &mass, &momentum, &energy, &mass_flux, &momentum_flux, &energy_flux);
}
