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
	float R = 1.0;
	float GAMMA = 1.4;
	float CFL = 0.5;
	Allocate_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &mass, &momentum_X, &momentum_Y, &energy, &mass_flux, &momentum_X_flux, &momentum_Y_flux, &energy_flux, N_CELLS);
	//Initial condition (p0 is density, p1 is X-direction veloctiy, p2 is Y-direction veloctiy, p3 is temperature, p4 is pressure)
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY+2) + j;
			if (i < NX/2){
				p0[INDEX] = 10.0;
				p1[INDEX] = 0.0;
				p2[INDEX] = 0.0;
				p3[INDEX] = 1.0;
			} else {
				p0[INDEX] = 1.0;
				p1[INDEX] = 0.0;
				p2[INDEX] = 0.0;
				p3[INDEX] = 1.0;
			}
			p4[INDEX] = p0[INDEX] * R * p3[INDEX];
		}
	}
	FILE *pFile = fopen("Results_of_5000_cells", "w");
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY; j++){
			int INDEX = i *(NY+2) + j;
			float X = (i+0.5) * dx;
			float Y = (j+0.5) * dy;
			fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", X, Y, p0[INDEX], p1[INDEX], p2[INDEX], p3[INDEX], p4[INDEX]);
		}
	}
	fclose;
	Free_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &mass, &momentum_X, &momentum_Y, &energy, &mass_flux, &momentum_X_flux, &momentum_Y_flux, &energy_flux);
return 0;
}
