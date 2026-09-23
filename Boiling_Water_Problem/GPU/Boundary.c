#include <stdlib.h>
#include <math.h>

void Boundary(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6, float R_dry, float R_v, int NX, int NY){
	float R_mix_L, R_mix_B, P_sat, phi_max;
	// LEFT and RIGHT (INFLOW and OUTFLOW)
	for (int j = 1 ; j < NY + 1; j++){
		int LEFT_GHOST = 0 * (NY+2) + j;
		int RIGHT_GHOST = (NX+1) * (NY+2) + j;
		int RIGHT_INNER = NX * (NY+2) + j;
		p1[LEFT_GHOST] = 5.0; // u = 5 m/s
		p2[LEFT_GHOST] = 0.0; // v = 0 m/s
		p3[LEFT_GHOST] = 300.0; // T = 300 K
		p4[LEFT_GHOST] = 101325.0; // P = 1 atm
		p5[LEFT_GHOST] = 0.0; // phi = 0
		R_mix_L = R_dry * (1 - p5[LEFT_GHOST]) + R_v * p5[LEFT_GHOST];
		p0[LEFT_GHOST] = p4[LEFT_GHOST] / (R_mix_L * p3[LEFT_GHOST]);

		p0[RIGHT_GHOST] = p0[RIGHT_INNER];
		p1[RIGHT_GHOST] = p1[RIGHT_INNER];
		p2[RIGHT_GHOST] = p2[RIGHT_INNER];
		p3[RIGHT_GHOST] = p3[RIGHT_INNER];
		p4[RIGHT_GHOST] = p4[RIGHT_INNER];
		p5[RIGHT_GHOST] = p5[RIGHT_INNER];
	}
	//BOTTOM and TOP
	for (int i = 1 ; i < NX + 1; i++){
		int BOTTOM_GHOST = i * (NY+2) + 0;
		int TOP_GHOST = i * (NY+2) + (NY+1);
		int BOTTOM_INNER = i * (NY+2) + 1;
		int TOP_INNER = i * (NY+2) + NY;
		p1[BOTTOM_GHOST] = p1[BOTTOM_INNER];
		p2[BOTTOM_GHOST] = -p2[BOTTOM_INNER];
		p3[BOTTOM_GHOST] = 373.0; //T = 373 K (Boiling water)
		p4[BOTTOM_GHOST] = p4[BOTTOM_INNER]; // P = 1 atm
		P_sat = 611.0 * exp((17.27 * (p3[BOTTOM_GHOST] - 273.15)) / ((p3[BOTTOM_GHOST] - 273.15) + 237.3)); // Tetens equation
		if (p4[BOTTOM_GHOST] <= P_sat) {
			phi_max = 1.0;
		} else {
			phi_max = (P_sat / R_v) / (((p4[BOTTOM_GHOST] - P_sat) / R_dry) + (P_sat / R_v));
		}

		p5[BOTTOM_GHOST] = phi_max;
		R_mix_B = R_dry * (1 - p5[BOTTOM_GHOST]) + R_v * p5[BOTTOM_GHOST];
		p0[BOTTOM_GHOST] = p4[BOTTOM_GHOST] / (R_mix_B * p3[BOTTOM_GHOST]);

		p0[TOP_GHOST] = p0[TOP_INNER];
		p1[TOP_GHOST] = p1[TOP_INNER];
		p2[TOP_GHOST] = p2[TOP_INNER];
		p3[TOP_GHOST] = p3[TOP_INNER];
		p4[TOP_GHOST] = p4[TOP_INNER];
		p5[TOP_GHOST] = p5[TOP_INNER];
	}
}
