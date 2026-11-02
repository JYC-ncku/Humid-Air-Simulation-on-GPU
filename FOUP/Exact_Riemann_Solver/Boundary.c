#include <stdlib.h>
#include <math.h>

void Boundary(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6, float R_dry, float R_v, float dx, float dy, int NX, int NY, int N_CELLS){

	float R_mix_L, R_mix_B, P_sat, phi_max;
	float X_WALL = 0.344;
	float Y_LOWER_WALL = 0.662;
	float Y_UPPER_WALL = 0.982;
	float Y_LOWER_HOLE = 0.932;
	int i_wall = (int)round(X_WALL / dx);
	int j_lower_wall = (int)round(Y_LOWER_WALL / dy);
	int j_upper_wall = (int)round(Y_UPPER_WALL / dy);
	int j_lower_hole = (int)round(Y_LOWER_HOLE / dy);

	// RIGHT (X = 959 mm)
	for (int j = 1; j < NY + 1; j++){
		int RIGHT_GHOST = (NX+1) * (NY+2) + j;
		int RIGHT_INNER = NX * (NY+2) + j;

		p0[RIGHT_GHOST] = p0[RIGHT_INNER];
		p1[RIGHT_GHOST] = -p1[RIGHT_INNER];
		p2[RIGHT_GHOST] = -p2[RIGHT_INNER];
		p3[RIGHT_GHOST] = p3[RIGHT_INNER];
		p4[RIGHT_GHOST] = p4[RIGHT_INNER];
		p5[RIGHT_GHOST] = p5[RIGHT_INNER];
	}

	// LEFT (X = 344.0 mm and Y = 0 ~ 662 mm)
	for (int j = 1; j < j_lower_wall; j++){
		int LEFT_GHOST = i_wall * (NY+2) + j;
		int LEFT_INNER = (i_wall+1) * (NY+2) + j;

		p0[LEFT_GHOST] = p0[LEFT_INNER];
		p1[LEFT_GHOST] = -p1[LEFT_INNER];
		p2[LEFT_GHOST] = -p2[LEFT_INNER];
		p3[LEFT_GHOST] = p3[LEFT_INNER];
		p4[LEFT_GHOST] = p4[LEFT_INNER];
		p5[LEFT_GHOST] = p5[LEFT_INNER];
	}

	// LEFT (X = 344.0 mm and Y = 982 ~ 1550 mm)
	for (int j = j_upper_wall + 1 ; j < NY + 1; j++){
		int LEFT_GHOST = i_wall * (NY+2) + j;
		int LEFT_INNER = (i_wall+1) * (NY+2) + j;

		p0[LEFT_GHOST] = p0[LEFT_INNER];
		p1[LEFT_GHOST] = -p1[LEFT_INNER];
		p2[LEFT_GHOST] = -p2[LEFT_INNER];
		p3[LEFT_GHOST] = p3[LEFT_INNER];
		p4[LEFT_GHOST] = p4[LEFT_INNER];
		p5[LEFT_GHOST] = p5[LEFT_INNER];
	}

	// LEFT (X = 0 and Y = 662 ~982 mm)
	for (int j = j_lower_wall + 1; j < j_upper_wall + 1; j++){
		int LEFT_GHOST = 0 * (NY+2) + j;
		int LEFT_INNER = 1 * (NY+2) + j;
		if (j >= j_lower_hole && j <= j_upper_wall){
			p1[LEFT_GHOST] = 1.0; // Inflow u = 1.0 m/s
			p2[LEFT_GHOST] = 0.0; // Inflow v = 0.0 m/s
			p3[LEFT_GHOST] = p3[LEFT_INNER];
			p4[LEFT_GHOST] = p4[LEFT_INNER];
			p5[LEFT_GHOST] = 0.0; // Inflow phi = 0;
		} else{
			p1[LEFT_GHOST] = -p1[LEFT_INNER];
			p2[LEFT_GHOST] = -p2[LEFT_INNER];
			p3[LEFT_GHOST] = p3[LEFT_INNER];
			p4[LEFT_GHOST] = p4[LEFT_INNER];
			p5[LEFT_GHOST] = p5[LEFT_INNER];
		}
	}


	for (int i = 1; i < NX + 1; i++){
//		int BOTTOM_BOTTOM_GHOST = i * (NY+4) + 0;
		int BOTTOM_GHOST = i * (NY+2) + 0;
//		int TOP_TOP_GHOST = i * (NY+4) + (NY+3);
		int TOP_GHOST = i * (NY+2) + (NY+1);

		int BOTTOM_INNER = i * (NY+2) + 1;
//		int BOTTOM_INNER_INNER = i * (NY+4) + 3;
		int TOP_INNER = i * (NY+2) + NY;
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
