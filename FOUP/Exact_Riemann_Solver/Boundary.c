#include <stdlib.h>
#include <math.h>

void Boundary(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6, float R_dry, float R_v, float dx, float dy, int NX, int NY, int N_CELLS){

	float R_mix_L, R_mix_B, P_sat, P_v, phi_max;
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

	// TOP (X = 344 ~ 959 mm and Y = 1550 mm)
	for (int i = i_wall + 1; i < NX + 1; i++){
		int TOP_GHOST = i * (NY+2) + (NY+1);
		int TOP_INNER = i * (NY+2) + NY;
		if (i >= i_wall + 1 && i <= i_inlet_1){
			p1[TOP_GHOST] = 0.0; // Inflow u = 0.0 m/s
			p2[TOP_GHOST] = -1.0; // Inflow v = -1.0 m/s
			p3[TOP_GHOST] = p3[TOP_INNER];
			p4[TOP_GHOST] = p4[TOP_INNER];
			p5[TOP_GHOST] = 0.0; // Inflow phi = 0;
		} else if( i >= i_inlet_2 && i < NX + 1){
			p1[TOP_GHOST] = 0.0; // Inflow u = 0.0 m/s
			p2[TOP_GHOST] = -0.5; // Inflow v = -0.5 m/s
			p3[TOP_GHOST] = p3[TOP_INNER];
			p4[TOP_GHOST] = p4[TOP_INNER];
			p6[TOP_GHOST] = 0.3; // Inflow RH = 30%;
			P_sat = 611.0 * exp((17.27 * (p3[TOP_GHOST] - 273.15)) / ((p3[TOP_GHOST] - 273.15) + 237.3)); // Tetens equation
			P_v = 0.30 * P_sat;
			p5[TOP_GHOST] = (0.622 * P_v) / (p4[TOP_GHOST] - 0.378 * P_v);
		} else {
			p1[TOP_GHOST] = -p1[TOP_INNER];
			p2[TOP_GHOST] = -p2[TOP_INNER];
			p3[TOP_GHOST] = p3[TOP_INNER];
			p4[TOP_GHOST] = p4[TOP_INNER];
			p5[TOP_GHOST] = p5[TOP_INNER];
		}
	}

	//TOP (X = 0 ~ 344 mm and Y = 982 mm)
	for (int i = 1; i <= i_wall; i++){
		int TOP_GHOST = i * (NY+2) + (j_upper_wall+1);
		int TOP_INNER = i * (NY+2) + (j_upper_wall);
		p1[TOP_GHOST] = -p1[TOP_INNER];
		p2[TOP_GHOST] = -p2[TOP_INNER];
		p3[TOP_GHOST] = p3[TOP_INNER];
		p4[TOP_GHOST] = p4[TOP_INNER];
		p5[TOP_GHOST] = p5[TOP_INNER];
	}

	//BOTTOM (X = 0 ~ 344 mm and Y = 662 mm)
	for (int i = 1; i <= i_wall; i++){
		int BOTTOM_GHOST = i * (NY+2) + (j_lower_wall);
		int BOTTOM_INNER = i * (NY+2) + (j_lower_wall+1);
		p1[BOTTOM_GHOST] = -p1[BOTTOM_INNER];
		p2[BOTTOM_GHOST] = -p2[BOTTOM_INNER];
		p3[BOTTOM_GHOST] = p3[BOTTOM_INNER];
		p4[BOTTOM_GHOST] = p4[BOTTOM_INNER];
		p5[BOTTOM_GHOST] = p5[BOTTOM_INNER];
	}

	//BOTTOM (X = 344 ~ 959 mm and Y = 0)
	for (int i = i_wall + 1; i < NX + 1 ; i++){
		int BOTTOM_GHOST = i * (NY+2) + 0;
		int BOTTOM_INNER = i * (NY+2) + 1;
		p1[BOTTOM_GHOST] = p1[BOTTOM_INNER];
		p2[BOTTOM_GHOST] = p2[BOTTOM_INNER];
		p3[BOTTOM_GHOST] = p3[BOTTOM_INNER];
		p4[BOTTOM_GHOST] = p4[BOTTOM_INNER];
		p5[BOTTOM_GHOST] = p5[BOTTOM_INNER];
	}
}
/*
void Boundary(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6, float R_dry, float R_v, float dx, float dy, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Boundary<<<GPB, TPB>>>(p0, p1, p2, p3, p4, p5, p6, R_dry, R_v, dx, dy, NX, NY, N_CELLS);
}
*/
