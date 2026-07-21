#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "memory.h"

int main(){
	int NX = 1000;
	int NY = 5;
	int N_CELLS = (NX+2) * (NY+2);
	float L = 1.0;
	float H = 0.005;
	float dx = L/NX;
	float dy = H/NY;
	float *x, *y, *p0, *p1, *p2, *p3, *p4, *mass, *momentum_X, *momentum_Y, *energy, *mass_flux, *momentum_X_flux, *momentum_Y_flux, *energy_flux;
	float t = 0.0;
	float t_FINAL = 0.2;

	Allocate_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &mass, &momentum_X, &momentum_Y, &energy, &mass_flux, &momentum_X_flux, &momentum_Y_flux, &energy_flux, N_CELLS);
	printf("Hello world!\n");
	Free_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &mass, &momentum_X, &momentum_Y, &energy, &mass_flux, &momentum_X_flux, &momentum_Y_flux, &energy_flux);
return 0;
}
