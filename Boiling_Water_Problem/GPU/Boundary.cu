#include <stdlib.h>
#include <math.h>

__global__ void GPU_Boundary(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6, float R_dry, float R_v, int NX, int NY, int N_CELLS){
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+2);
	int j = (int)INDEX - i * (NY+2);
	float R_mix_L, R_mix_B, P_sat, phi_max;
	if (INDEX < N_CELLS){
		// LEFT and RIGHT (INFLOW and OUTFLOW)
		if (j >= 1 && j < NY + 1 && i == 0){
			int LEFT_GHOST = 0 * (NY+2) + j;
			int RIGHT_GHOST = (NX+1) * (NY+2) + j;
			int RIGHT_INNER = NX * (NY+2) + j;
			d_p1[LEFT_GHOST] = 5.0; // u = 5 m/s
			d_p2[LEFT_GHOST] = 0.0; // v = 0 m/s
			d_p3[LEFT_GHOST] = 300.0; // T = 300 K
			d_p4[LEFT_GHOST] = 101325.0; // P = 1 atm
			d_p5[LEFT_GHOST] = 0.0; // phi = 0
			R_mix_L = R_dry * (1 - d_p5[LEFT_GHOST]) + R_v * d_p5[LEFT_GHOST];
			d_p0[LEFT_GHOST] = d_p4[LEFT_GHOST] / (R_mix_L * d_p3[LEFT_GHOST]);

			d_p0[RIGHT_GHOST] = d_p0[RIGHT_INNER];
			d_p1[RIGHT_GHOST] = d_p1[RIGHT_INNER];
			d_p2[RIGHT_GHOST] = d_p2[RIGHT_INNER];
			d_p3[RIGHT_GHOST] = d_p3[RIGHT_INNER];
			d_p4[RIGHT_GHOST] = d_p4[RIGHT_INNER];
			d_p5[RIGHT_GHOST] = d_p5[RIGHT_INNER];
		}
		//BOTTOM and TOP
		if (i >= 1 && i < NX + 1 && j == 0){
			int BOTTOM_GHOST = i * (NY+2) + 0;
			int TOP_GHOST = i * (NY+2) + (NY+1);
			int BOTTOM_INNER = i * (NY+2) + 1;
			int TOP_INNER = i * (NY+2) + NY;
			d_p1[BOTTOM_GHOST] = d_p1[BOTTOM_INNER];
			d_p2[BOTTOM_GHOST] = -d_p2[BOTTOM_INNER];
			d_p3[BOTTOM_GHOST] = 373.0; //T = 373 K (Boiling water)
			d_p4[BOTTOM_GHOST] = d_p4[BOTTOM_INNER]; // P = 1 atm
			P_sat = 611.0 * exp((17.27 * (d_p3[BOTTOM_GHOST] - 273.15)) / ((d_p3[BOTTOM_GHOST] - 273.15) + 237.3)); // Tetens equation
			if (d_p4[BOTTOM_GHOST] <= P_sat) {
				phi_max = 1.0;
			} else {
				phi_max = (P_sat / R_v) / (((d_p4[BOTTOM_GHOST] - P_sat) / R_dry) + (P_sat / R_v));
			}

			d_p5[BOTTOM_GHOST] = phi_max;
			R_mix_B = R_dry * (1 - d_p5[BOTTOM_GHOST]) + R_v * d_p5[BOTTOM_GHOST];
			d_p0[BOTTOM_GHOST] = d_p4[BOTTOM_GHOST] / (R_mix_B * d_p3[BOTTOM_GHOST]);

			d_p0[TOP_GHOST] = d_p0[TOP_INNER];
			d_p1[TOP_GHOST] = d_p1[TOP_INNER];
			d_p2[TOP_GHOST] = d_p2[TOP_INNER];
			d_p3[TOP_GHOST] = d_p3[TOP_INNER];
			d_p4[TOP_GHOST] = d_p4[TOP_INNER];
			d_p5[TOP_GHOST] = d_p5[TOP_INNER];
		}
	}
}

void Boundary(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6, float R_dry, float R_v, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Boundary<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6,  R_dry,  R_v,  NX,  NY,  N_CELLS);
}
