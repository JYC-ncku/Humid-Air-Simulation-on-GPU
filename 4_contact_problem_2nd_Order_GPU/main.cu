#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memory.h"
#include "Calc_rho_u_P_T.h"
#include "Boundary.h"
#include "Initial.h"
#include "Calc_new_properties.h"

/*
	  GHOST				  GHOST
	一一一一一一一一一一一一一一一一一一一一
	|	|	|       |	|	|
	|   0   |   1   |  ...  |   N   |  N+1  |
	|       |       |       |	|	|
	一一一一一一一一一一一一一一一一一一一一
		0	1      ...	N
*/

int main(){
	int NX = 1000;
	int NY = 1000;
	int N_CELLS = (NX+4) * (NY+4); //+2 for Ghost cells
	float L = 1.0;
	float H = 1.0;
	float dx = L/NX;
	float dy = H/NY;
	float t = 0;
	float t_FINAL = 0.2;
	float R = 1.0;
	float GAMMA = 1.4;
//	int wall_flag = 0;
//	float MAX_CFL = 1e-10;
	float *x, *y, *CFL,
	      *h_p0, *h_p1, *h_p2, *h_p3, *h_p4,
	      *d_p0, *d_p1, *d_p2, *d_p3, *d_p4,
	      *d_interface_p, *d_flux_X, *d_flux_Y; //p0 is density, p1 is x-dir velocity, p2 is y-dir veloctiy, p3 is temperature, p4 si pressure.
//	float flxnmn, flxpmn, flxqmn;
//	float CFL_t = 0.5;

	Allocate_memory(&x, &y, &CFL,
			&h_p0, &h_p1, &h_p2, &h_p3, &h_p4,
			&d_p0, &d_p1, &d_p2, &d_p3, &d_p4,
			&d_interface_p, &d_flux_X, &d_flux_Y, N_CELLS);
	//Send the data to device
	Send_To_Device(&d_p0, &h_p0, N_CELLS);
	Send_To_Device(&d_p1, &h_p1, N_CELLS);
	Send_To_Device(&d_p2, &h_p2, N_CELLS);
	Send_To_Device(&d_p3, &h_p3, N_CELLS);
	Send_To_Device(&d_p4, &h_p4, N_CELLS);
	//Initial condition
	Initial(CFL, d_p0, d_p1, d_p2, d_p3, d_p4, R, NX, NY, N_CELLS);

	while (t<t_FINAL){
		float dt = Compute_MAX_CFL(CFL, d_p0, d_p1, d_p2, d_p3,  dx,  dy,  NX,  NY,  N_CELLS);
		// Boundary condition for compute flux.
		Boundary(d_p0, d_p1, d_p2, d_p3, d_p4, NX, NY, N_CELLS);

//		float dt = CFL_t / MAX_CFL;

		Calc_flux(d_interface_p, d_flux_X, d_flux_Y, d_p0, d_p1, d_p2, d_p3, d_p4, R, GAMMA, dx, dy, NX, NY, N_CELLS);
		cudaDeviceSynchronize();
		Calc_new_properties(d_flux_X, d_flux_Y, d_p0, d_p1, d_p2, d_p3, d_p4, R, GAMMA, dt, dx, dy, NX, NY, N_CELLS);
		cudaDeviceSynchronize();
		t += dt;
	}
	//Get the data from device
	Get_From_Device(&h_p0, &d_p0, N_CELLS);
	Get_From_Device(&h_p1, &d_p1, N_CELLS);
	Get_From_Device(&h_p2, &d_p2, N_CELLS);
	Get_From_Device(&h_p3, &d_p3, N_CELLS);
	Get_From_Device(&h_p4, &d_p4, N_CELLS);

	FILE * pFile = fopen("Results_of_1000000_cells.txt","w");
	for (int i = 2; i < NX + 2; i++){
		for (int j = 2; j < NY + 2; j++){
			int INDEX = i * (NY + 4) + j;
			float X = (i - 1.5) * dx;
			float Y = (j - 1.5) * dy;
			fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.2f\n", X, Y, h_p0[INDEX], h_p1[INDEX], h_p2[INDEX], h_p3[INDEX], h_p4[INDEX]);
		}
	}
	fclose(pFile);

	Free_memory(&x, &y, &CFL, &h_p0, &h_p1, &h_p2, &h_p3, &h_p4, &d_p0, &d_p1, &d_p2, &d_p3, &d_p4, &d_interface_p, &d_flux_X, &d_flux_Y);
return 0;
}

