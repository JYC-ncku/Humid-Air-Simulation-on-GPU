#include <stdlib.h>
#include <math.h>

void Boundary(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6, float R_dry, float R_v, float dx, float dy, int NX, int NY, int N_CELLS){
/*
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+4);
	int j = (int)INDEX - i * (NY+4);
*/
	float U_max = 0.1;
	float Delta = 0.1; // Boundary thickness
	float R_mix_L, R_mix_B, P_sat, phi_max;
//	if (INDEX < N_CELLS){
		// LEFT and RIGHT (INFLOW and OUTFLOW)
//		if (j >= 2 && j < NY + 2 && i == 0){
	for (int j = 1; j < NY + 1; j++){
//		int LEFT_LEFT_GHOST = 0 * (NY+4) + j;
		int LEFT_GHOST = 1 * (NY+4) + j;
//		int RIGHT_RIGHT_GHOST = (NX+3) * (NY+4) + j;
		int RIGHT_GHOST = (NX+2) * (NY+4) + j;

		int LEFT_INNER = 2 * (NY+4) + j;
		int RIGHT_INNER = (NX+1) * (NY+4) + j;
		float y_center = (j - 0.5) * dy;
//		p1[LEFT_GHOST] = 5.0; // u = 5 m/s
		if (y_center < Delta) {
			float ratio = y_center / Delta;
			p1[LEFT_GHOST] = U_max * (2.0 * ratio - ratio * ratio);
		} else {
			p1[LEFT_GHOST] = U_max;
		}

		p2[LEFT_GHOST] = 0.0; // v = 0 m/s
		p3[LEFT_GHOST] = 300.15; // T = 300.15 K
		p4[LEFT_GHOST] = p4[LEFT_INNER]; // P = 1 atm
		p5[LEFT_GHOST] = 0.0; // phi = 0
		R_mix_L = R_dry * (1 - p5[LEFT_GHOST]) + R_v * p5[LEFT_GHOST];
		p0[LEFT_GHOST] = p4[LEFT_GHOST] / (R_mix_L * p3[LEFT_GHOST]);
/*
		p0[LEFT_LEFT_GHOST] = p0[LEFT_GHOST];
		p1[LEFT_LEFT_GHOST] = p1[LEFT_GHOST];
		p2[LEFT_LEFT_GHOST] = p2[LEFT_GHOST];
		p3[LEFT_LEFT_GHOST] = p3[LEFT_GHOST];
		p4[LEFT_LEFT_GHOST] = p4[LEFT_GHOST];
		p5[LEFT_LEFT_GHOST] = p5[LEFT_GHOST];
*/
		p0[RIGHT_GHOST] = p0[RIGHT_INNER];
		p1[RIGHT_GHOST] = p1[RIGHT_INNER];
		p2[RIGHT_GHOST] = p2[RIGHT_INNER];
		p3[RIGHT_GHOST] = p3[RIGHT_INNER];
		p4[RIGHT_GHOST] = p4[RIGHT_INNER];
		p5[RIGHT_GHOST] = p5[RIGHT_INNER];
/*
		p0[RIGHT_RIGHT_GHOST] = p0[RIGHT_INNER];
		p1[RIGHT_RIGHT_GHOST] = p1[RIGHT_INNER];
		p2[RIGHT_RIGHT_GHOST] = p2[RIGHT_INNER];
		p3[RIGHT_RIGHT_GHOST] = p3[RIGHT_INNER];
		p4[RIGHT_RIGHT_GHOST] = p4[RIGHT_INNER];
		p5[RIGHT_RIGHT_GHOST] = p5[RIGHT_INNER];
*/
	}
	//BOTTOM and TOP
//		if (i >= 2 && i < NX + 2 && j == 0){
	for (int i = 2; i < NX + 2; i++){
//		int BOTTOM_BOTTOM_GHOST = i * (NY+4) + 0;
		int BOTTOM_GHOST = i * (NY+4) + 1;
//		int TOP_TOP_GHOST = i * (NY+4) + (NY+3);
		int TOP_GHOST = i * (NY+4) + (NY+2);

		int BOTTOM_INNER = i * (NY+4) + 2;
//		int BOTTOM_INNER_INNER = i * (NY+4) + 3;
		int TOP_INNER = i * (NY+4) + (NY+1);
		float x_center = (i - 0.5) * dx;

		p1[BOTTOM_GHOST] = -p1[BOTTOM_INNER];
		p2[BOTTOM_GHOST] = -p2[BOTTOM_INNER];
		p4[BOTTOM_GHOST] = p4[BOTTOM_INNER]; // P = 1 atm
/*
		p1[BOTTOM_BOTTOM_GHOST] = -p1[BOTTOM_INNER_INNER];
		p2[BOTTOM_BOTTOM_GHOST] = -p2[BOTTOM_INNER_INNER];
		p4[BOTTOM_BOTTOM_GHOST] = p4[BOTTOM_INNER_INNER];
*/

		if (x_center >= 0.2 && x_center <= 0.7){
			p3[BOTTOM_GHOST] = 373.15; //T = 373.15 K (Boiling water)
			P_sat = 611.0 * exp((17.27 * (p3[BOTTOM_GHOST] - 273.15)) / ((p3[BOTTOM_GHOST] - 273.15) + 237.3)); // Tetens equation
			if (p4[BOTTOM_GHOST] <= P_sat) {
				phi_max = 1.0;
			} else {
				phi_max = (P_sat / R_v) / (((p4[BOTTOM_GHOST] - P_sat) / R_dry) + (P_sat / R_v));
			}
			p5[BOTTOM_GHOST] = phi_max;
		} else {
			p3[BOTTOM_GHOST] = 300.15;
			p5[BOTTOM_GHOST] = 0.0;
		}

		R_mix_B = R_dry * (1 - p5[BOTTOM_GHOST]) + R_v * p5[BOTTOM_GHOST];
		p0[BOTTOM_GHOST] = p4[BOTTOM_GHOST] / (R_mix_B * p3[BOTTOM_GHOST]);

/*
		p3[BOTTOM_BOTTOM_GHOST] = p3[BOTTOM_GHOST];
		p5[BOTTOM_BOTTOM_GHOST] = p5[BOTTOM_GHOST];
		p0[BOTTOM_BOTTOM_GHOST] = p0[BOTTOM_GHOST];
*/
		p0[TOP_GHOST] = p0[TOP_INNER];
		p1[TOP_GHOST] = p1[TOP_INNER];
		p2[TOP_GHOST] = p2[TOP_INNER];
		p3[TOP_GHOST] = p3[TOP_INNER];
		p4[TOP_GHOST] = p4[TOP_INNER];
		p5[TOP_GHOST] = p5[TOP_INNER];
/*
		p0[TOP_TOP_GHOST] = p0[TOP_INNER];
		p1[TOP_TOP_GHOST] = p1[TOP_INNER];
		p2[TOP_TOP_GHOST] = p2[TOP_INNER];
		p3[TOP_TOP_GHOST] = p3[TOP_INNER];
		p4[TOP_TOP_GHOST] = p4[TOP_INNER];
		p5[TOP_TOP_GHOST] = p5[TOP_INNER];
*/
	}
}
/*
void Boundary(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6, float R_dry, float R_v, float dx, float dy, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Boundary<<<GPB, TPB>>>(p0, p1, p2, p3, p4, p5, p6, R_dry, R_v, dx, dy, NX, NY, N_CELLS);
}
*/
