#include <stdlib.h>
#include <math.h>

__global__ void GPU_Boundary(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6, float R_dry, float R_v, float dx, float dy, int NX, int NY, int N_CELLS){
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+4);
	int j = (int)INDEX - i * (NY+4);
	float x_center = (i - 0.5) * dx;
	float y_center = (j - 0.5) * dy;
	float U_max = 5.0;
	float Delta = 0.1; // Boundary thickness
	float R_mix_L, R_mix_B, P_sat, phi_max;
	if (INDEX < N_CELLS){
		// LEFT and RIGHT (INFLOW and OUTFLOW)
		if (j >= 2 && j < NY + 1 && i == 0){
			int LEFT_LEFT_GHOST = 0 * (NY+4) + j;
			int LEFT_GHOST = 1 * (NY+4) + j;
			int RIGHT_RIGHT_GHOST = (NX+3) * (NY+4) + j;
			int RIGHT_GHOST = (NX+2) * (NY+4) + j;

			int LEFT_INNER = 2 * (NY+4) + j;
			int RIGHT_INNER = (NX+1) * (NY+4) + j;
			int RIGHT_INNER_INNER = NX * (NY+4) + j;
//			d_p1[LEFT_GHOST] = 5.0; // u = 5 m/s
			if (y_center < Delta) {
				float ratio = y_center / Delta;
				d_p1[LEFT_GHOST] = U_max * (2.0 * ratio - ratio * ratio);
			} else {
				d_p1[LEFT_GHOST] = U_max;
			}

			d_p2[LEFT_GHOST] = 0.0; // v = 0 m/s
			d_p3[LEFT_GHOST] = 300.15; // T = 300.15 K
			d_p4[LEFT_GHOST] = d_p4[LEFT_INNER]; // P = 1 atm
			d_p5[LEFT_GHOST] = 0.0; // phi = 0
			R_mix_L = R_dry * (1 - d_p5[LEFT_GHOST]) + R_v * d_p5[LEFT_GHOST];
			d_p0[LEFT_GHOST] = d_p4[LEFT_GHOST] / (R_mix_L * d_p3[LEFT_GHOST]);

			d_p0[LEFT_LEFT_GHOST] = d_p0[LEFT_GHOST];
			d_p1[LEFT_LEFT_GHOST] = d_p1[LEFT_GHOST];
			d_p2[LEFT_LEFT_GHOST] = d_p2[LEFT_GHOST];
			d_p3[LEFT_LEFT_GHOST] = d_p3[LEFT_GHOST];
			d_p4[LEFT_LEFT_GHOST] = d_p4[LEFT_GHOST];
			d_p5[LEFT_LEFT_GHOST] = d_p5[LEFT_GHOST];

			d_p0[RIGHT_GHOST] = d_p0[RIGHT_INNER];
			d_p1[RIGHT_GHOST] = d_p1[RIGHT_INNER];
			d_p2[RIGHT_GHOST] = d_p2[RIGHT_INNER];
			d_p3[RIGHT_GHOST] = d_p3[RIGHT_INNER];
			d_p4[RIGHT_GHOST] = d_p4[RIGHT_INNER];
			d_p5[RIGHT_GHOST] = d_p5[RIGHT_INNER];

			d_p0[RIGHT_RIGHT_GHOST] = d_p0[RIGHT_INNER_INNER];
			d_p1[RIGHT_RIGHT_GHOST] = d_p1[RIGHT_INNER_INNER];
			d_p2[RIGHT_RIGHT_GHOST] = d_p2[RIGHT_INNER_INNER];
			d_p3[RIGHT_RIGHT_GHOST] = d_p3[RIGHT_INNER_INNER];
			d_p4[RIGHT_RIGHT_GHOST] = d_p4[RIGHT_INNER_INNER];
			d_p5[RIGHT_RIGHT_GHOST] = d_p5[RIGHT_INNER_INNER];
		}
		//BOTTOM and TOP
		if (i >= 2 && i < NX + 1 && j == 0){
			int BOTTOM_BOTTOM_GHOST = i * (NY+4) + 0;
			int BOTTOM_GHOST = i * (NY+4) + 1;
			int TOP_TOP_GHOST = i * (NY+4) + (NY+3);
			int TOP_GHOST = i * (NY+4) + (NY+2);

			int BOTTOM_INNER = i * (NY+4) + 2;
			int BOTTOM_INNER_INNER = i * (NY+4) + 3;
			int TOP_INNER = i * (NY+4) + (NY+1);
			int TOP_INNER_INNER = i * (NY+4) + NY;

			d_p1[BOTTOM_GHOST] = -d_p1[BOTTOM_INNER];
			d_p2[BOTTOM_GHOST] = -d_p2[BOTTOM_INNER];
			d_p4[BOTTOM_GHOST] = d_p4[BOTTOM_INNER]; // P = 1 atm

			d_p1[BOTTOM_BOTTOM_GHOST] = -d_p1[BOTTOM_INNER_INNER];
			d_p2[BOTTOM_BOTTOM_GHOST] = -d_p2[BOTTOM_INNER_INNER];
			d_p4[BOTTOM_BOTTOM_GHOST] = d_p4[BOTTOM_INNER_INNER];


			if (x_center >= 0.2 && x_center <= 0.7){
				d_p3[BOTTOM_GHOST] = 373.15; //T = 373.15 K (Boiling water)
				P_sat = 611.0 * exp((17.27 * (d_p3[BOTTOM_GHOST] - 273.15)) / ((d_p3[BOTTOM_GHOST] - 273.15) + 237.3)); // Tetens equation
				if (d_p4[BOTTOM_GHOST] <= P_sat) {
					phi_max = 1.0;
				} else {
					phi_max = (P_sat / R_v) / (((d_p4[BOTTOM_GHOST] - P_sat) / R_dry) + (P_sat / R_v));
				}
				d_p5[BOTTOM_GHOST] = phi_max;
			} else {
				d_p3[BOTTOM_GHOST] = 300.15;
				d_p5[BOTTOM_GHOST] = 0.0;
			}

			R_mix_B = R_dry * (1 - d_p5[BOTTOM_GHOST]) + R_v * d_p5[BOTTOM_GHOST];
			d_p0[BOTTOM_GHOST] = d_p4[BOTTOM_GHOST] / (R_mix_B * d_p3[BOTTOM_GHOST]);


			d_p3[BOTTOM_BOTTOM_GHOST] = d_p3[BOTTOM_GHOST];
			d_p5[BOTTOM_BOTTOM_GHOST] = d_p5[BOTTOM_GHOST];
			d_p0[BOTTOM_BOTTOM_GHOST] = d_p0[BOTTOM_GHOST];

			d_p0[TOP_GHOST] = d_p0[TOP_INNER];
			d_p1[TOP_GHOST] = d_p1[TOP_INNER];
			d_p2[TOP_GHOST] = d_p2[TOP_INNER];
			d_p3[TOP_GHOST] = d_p3[TOP_INNER];
			d_p4[TOP_GHOST] = d_p4[TOP_INNER];
			d_p5[TOP_GHOST] = d_p5[TOP_INNER];

			d_p0[TOP_TOP_GHOST] = d_p0[TOP_INNER_INNER];
			d_p1[TOP_TOP_GHOST] = d_p1[TOP_INNER_INNER];
			d_p2[TOP_TOP_GHOST] = d_p2[TOP_INNER_INNER];
			d_p3[TOP_TOP_GHOST] = d_p3[TOP_INNER_INNER];
			d_p4[TOP_TOP_GHOST] = d_p4[TOP_INNER_INNER];
			d_p5[TOP_TOP_GHOST] = d_p5[TOP_INNER_INNER];

		}
	}
}

void Boundary(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6, float R_dry, float R_v, float dx, float dy, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Boundary<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6, R_dry, R_v, dx, dy, NX, NY, N_CELLS);
}
