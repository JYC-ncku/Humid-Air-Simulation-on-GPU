#include <stdlib.h>
#include <math.h>

__global__ void GPU_Boundary(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6, float R_dry, float R_v, float dx, float dy, int NX, int NY, int N_CELLS){

	float R_mix, P_sat, P_v;
	float X_WALL = 0.344;
	float X_INLET_1 = 0.494; // INLET 1 hole width = 150 mm; 0.344+0.15=0.494
	float X_INLET_2 = 0.709; // INLET 2 hole width = 250 mm; 0.959-0.25=0.709
	float Y_LOWER_WALL = 0.662;
	float Y_UPPER_WALL = 0.982;
	float Y_LOWER_HOLE = 0.932;
	int i_wall = (int)round(X_WALL / dx);
	int i_inlet_1 = (int)round(X_INLET_1 / dx);
	int i_inlet_2 = (int)round(X_INLET_2 / dx);
	int j_lower_wall = (int)round(Y_LOWER_WALL / dy);
	int j_upper_wall = (int)round(Y_UPPER_WALL / dy);
	int j_lower_hole = (int)round(Y_LOWER_HOLE / dy);

	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+4);
	int j = (int)INDEX - i * (NY+4);

	if (INDEX < N_CELLS){
		// RIGHT (X = 959 mm)
		if (j >= 2 && j < NY + 2){
			int RIGHT_GHOST = (NX+2) * (NY+4) + j;
			int RIGHT_RIGHT_GHOST = (NX+3) * (NY+4) + j;;
			int RIGHT_INNER = (NX+1) * (NY+4) + j;
			int RIGHT_RIGHT_INNER = NX * (NY+4) + j;

			d_p0[RIGHT_GHOST] = d_p0[RIGHT_INNER];
			d_p1[RIGHT_GHOST] = -d_p1[RIGHT_INNER];
			d_p2[RIGHT_GHOST] = -d_p2[RIGHT_INNER];
			d_p3[RIGHT_GHOST] = d_p3[RIGHT_INNER];
			d_p4[RIGHT_GHOST] = d_p4[RIGHT_INNER];
			d_p5[RIGHT_GHOST] = d_p5[RIGHT_INNER];

			d_p0[RIGHT_RIHGT_GHOST] = d_p0[RIGHT_RIHGT_INNER];
			d_p1[RIGHT_RIHGT_GHOST] = -d_p1[RIGHT_RIHGT_INNER];
			d_p2[RIGHT_RIHGT_GHOST] = -d_p2[RIGHT_RIHGT_INNER];
			d_p3[RIGHT_RIHGT_GHOST] = d_p3[RIGHT_RIHGT_INNER];
			d_p4[RIGHT_RIHGT_GHOST] = d_p4[RIGHT_RIHGT_INNER];
			d_p5[RIGHT_RIHGT_GHOST] = d_p5[RIGHT_RIHGT_INNER];
		}

		// LEFT (X = 0 and Y = 662 ~982 mm)
		if (j >= 2 && j < NY + 2){
			int LEFT_GHOST = 0 * (NY+4) + j;
			int LEFT_INNER = 1 * (NY+4) + j;
			if (j >= j_lower_hole && j <= j_upper_wall){
				d_p1[LEFT_GHOST] = 1.0; // Inflow u = 1.0 m/s
				d_p2[LEFT_GHOST] = 0.0; // Inflow v = 0.0 m/s
				d_p3[LEFT_GHOST] = d_p3[LEFT_INNER];
				d_p4[LEFT_GHOST] = d_p4[LEFT_INNER];
				d_p5[LEFT_GHOST] = 0.0; // Inflow phi = 0;
			} else{
				d_p1[LEFT_GHOST] = -d_p1[LEFT_INNER];
				d_p2[LEFT_GHOST] = -d_p2[LEFT_INNER];
				d_p3[LEFT_GHOST] = d_p3[LEFT_INNER];
				d_p4[LEFT_GHOST] = d_p4[LEFT_INNER];
				d_p5[LEFT_GHOST] = d_p5[LEFT_INNER];
			}
		}

		// TOP (X = 344 ~ 959 mm and Y = 1550 mm)
		if (i >= 2 && i < NX + 2){
			int TOP_GHOST = i * (NY+4) + (NY+1);
			int TOP_INNER = i * (NY+4) + NY;
			// Inlet 1
			if (i >= i_wall + 1 && i <= i_inlet_1){
				d_p1[TOP_GHOST] = 0.0; // Inlet u = 0.0 m/s
				d_p2[TOP_GHOST] = -1.0; // Inlet v = -1.0 m/s
				d_p3[TOP_GHOST] = d_p3[TOP_INNER];
				d_p4[TOP_GHOST] = d_p4[TOP_INNER];
				d_p5[TOP_GHOST] = 0.0; // Inlet phi = 0;
			// Inlet 2
			} else if( i >= i_inlet_2 && i < NX + 1){
				d_p1[TOP_GHOST] = 0.0; // Inlet u = 0.0 m/s
				d_p2[TOP_GHOST] = -0.5; // Inlet v = -0.5 m/s
				d_p3[TOP_GHOST] = d_p3[TOP_INNER];
				d_p4[TOP_GHOST] = d_p4[TOP_INNER];
				d_p6[TOP_GHOST] = 0.3; // Inlet RH = 30%;
				P_sat = 611.0 * exp((17.27 * (d_p3[TOP_GHOST] - 273.15)) / ((d_p3[TOP_GHOST] - 273.15) + 237.3)); // Tetens equation
				P_v = 0.30 * P_sat;
				d_p5[TOP_GHOST] = (0.622 * P_v) / (d_p4[TOP_GHOST] - 0.378 * P_v);
			} else {
				d_p1[TOP_GHOST] = -d_p1[TOP_INNER];
				d_p2[TOP_GHOST] = -d_p2[TOP_INNER];
				d_p3[TOP_GHOST] = d_p3[TOP_INNER];
				d_p4[TOP_GHOST] = d_p4[TOP_INNER];
				d_p5[TOP_GHOST] = d_p5[TOP_INNER];
			}
		}
		//BOTTOM (X = 0 ~ 344 mm and Y = 662 mm)
		if (i >= 2 && i < NX + 2){
			int BOTTOM_GHOST = i * (NY+4) + 0;
			int BOTTOM_INNER = i * (NY+4) + 1;
			if (i >= i_wall + 1) {
				d_p1[BOTTOM_GHOST] = d_p1[BOTTOM_INNER];
				d_p2[BOTTOM_GHOST] = d_p2[BOTTOM_INNER];
				d_p3[BOTTOM_GHOST] = d_p3[BOTTOM_INNER];
				d_p4[BOTTOM_GHOST] = d_p4[BOTTOM_INNER];
				d_p5[BOTTOM_GHOST] = d_p5[BOTTOM_INNER];
			} else{
				//BOTTOM (X = 344 ~ 959 mm and Y = 0) Outlet
				d_p1[BOTTOM_GHOST] = -d_p1[BOTTOM_INNER];
				d_p2[BOTTOM_GHOST] = -d_p2[BOTTOM_INNER];
				d_p3[BOTTOM_GHOST] = d_p3[BOTTOM_INNER];
				d_p4[BOTTOM_GHOST] = d_p4[BOTTOM_INNER];
				d_p5[BOTTOM_GHOST] = d_p5[BOTTOM_INNER];
			}
		}
	}
}
void Boundary(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6, float R_dry, float R_v, float dx, float dy, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Boundary<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6, R_dry, R_v, dx, dy, NX, NY, N_CELLS);
}

