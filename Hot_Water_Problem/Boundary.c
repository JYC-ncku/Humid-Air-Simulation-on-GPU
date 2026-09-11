#include <stdlib.h>
#include <math.h>

void Boundary(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float Ru, float R_v, float R_dry, int NX, int NY){
	for (int j = 2 ; j <= NY+1; j++){
		int LEFT_LEFT_GHOST = 0 * (NY+4) + j;
		int LEFT_GHOST = 1 * (NY+4) + j;
		int RIGHT_RIGHT_GHOST = (NX+3) * (NY+4) + j;
		int RIGHT_GHOST = (NX+2) * (NY+4) + j;

		int LEFT_INNER = 2 * (NY+4) + j;
		int LEFT_LEFT_INNER = 3 * (NY+4) + j;
		int RIGHT_INNER = (NX+1) * (NY+4) + j;
		int RIGHT_RIGHT_INNER = NX * (NY+4) + j;

		p0[LEFT_GHOST] = p0[LEFT_INNER];
		p0[LEFT_LEFT_GHOST] = p0[LEFT_LEFT_INNER];
		p0[RIGHT_GHOST] = p0[RIGHT_INNER];
		p0[RIGHT_RIGHT_GHOST] = p0[RIGHT_RIGHT_INNER];

		// Velocity is reflective becuase left and right have the wall
		p1[LEFT_GHOST] = -p1[LEFT_INNER];
		p1[LEFT_LEFT_GHOST] = -p1[LEFT_LEFT_INNER];
		p1[RIGHT_GHOST] = -p1[RIGHT_INNER];
		p1[RIGHT_RIGHT_GHOST] = -p1[RIGHT_RIGHT_INNER];

		p2[LEFT_GHOST] = -p2[LEFT_INNER];
		p2[LEFT_LEFT_GHOST] = -p2[LEFT_LEFT_INNER];
		p2[RIGHT_GHOST] = -p2[RIGHT_INNER];
		p2[RIGHT_RIGHT_GHOST] = -p2[RIGHT_RIGHT_INNER];

		p3[LEFT_GHOST] = p3[LEFT_INNER];
		p3[LEFT_LEFT_GHOST] = p3[LEFT_LEFT_INNER];
		p3[RIGHT_GHOST] = p3[RIGHT_INNER];
		p3[RIGHT_RIGHT_GHOST] = p3[RIGHT_RIGHT_INNER];

		p4[LEFT_GHOST] = p4[LEFT_INNER];
		p4[LEFT_LEFT_GHOST] = p4[LEFT_LEFT_INNER];
		p4[RIGHT_GHOST] = p4[RIGHT_INNER];
		p4[RIGHT_RIGHT_GHOST] = p4[RIGHT_RIGHT_INNER];

		p5[LEFT_GHOST] = p5[LEFT_INNER];
		p5[LEFT_LEFT_GHOST] = p5[LEFT_LEFT_INNER];
		p5[RIGHT_GHOST] = p5[RIGHT_INNER];
		p5[RIGHT_RIGHT_GHOST] = p5[RIGHT_RIGHT_INNER];
	}
	//BOTTOM and TOP
	float T_base = 373.15; // The heat source in the BOTTOM is maintained at 100 C = 373.15 K
	float Y_base = 1.0; // The  mass fraction is maintained at 100% relative humidity in the BOTTOM
	float P_atm = 101325;
	for (int i = 2 ; i <= NX+1; i++){
		int BOTTOM_BOTTOM_GHOST = i * (NY+4) + 0;
		int BOTTOM_GHOST = i * (NY+4) + 1;
		int TOP_TOP_GHOST = i * (NY+4) + (NY+3);
		int TOP_GHOST = i * (NY+4) + (NY+2);

		int BOTTOM_INNER = i * (NY+4) + 2;
		int BOTTOM_BOTTOM_INNER = i * (NY+4) + 3;
		int TOP_INNER = i * (NY+4) + (NY+1);
		int TOP_TOP_INNER = i * (NY+4) + NY;

		// BOTTOM
		p1[BOTTOM_GHOST] = -p1[BOTTOM_INNER];
		p1[BOTTOM_BOTTOM_GHOST] = -p1[BOTTOM_BOTTOM_INNER];

		p2[BOTTOM_GHOST] = -p2[BOTTOM_INNER];
		p2[BOTTOM_BOTTOM_GHOST] = -p2[BOTTOM_BOTTOM_INNER];

		p3[BOTTOM_GHOST] = 2 * T_base - p3[BOTTOM_INNER];
		p3[BOTTOM_BOTTOM_GHOST] = 2 * T_base - p3[BOTTOM_BOTTOM_INNER];

		p4[BOTTOM_GHOST] = p4[BOTTOM_INNER];
		p4[BOTTOM_BOTTOM_GHOST] = p4[BOTTOM_BOTTOM_INNER];

		p5[BOTTOM_GHOST] = 2 * Y_base - p5[BOTTOM_INNER];
		p5[BOTTOM_BOTTOM_GHOST] = 2 * Y_base - p5[BOTTOM_BOTTOM_INNER];

		// rho = P / RT
		float R_mix_B = (1 - p5[BOTTOM_GHOST]) * R_dry + p5[BOTTOM_GHOST] * R_v;
		float R_mix_BB = (1 - p5[BOTTOM_BOTTOM_GHOST]) * R_dry + p5[BOTTOM_BOTTOM_GHOST] * R_v;
		p0[BOTTOM_GHOST] = p4[BOTTOM_GHOST] / (R_mix_B * p3[BOTTOM_GHOST]);
		p0[BOTTOM_BOTTOM_GHOST] = p4[BOTTOM_BOTTOM_GHOST] / (R_mix_BB * p3[BOTTOM_BOTTOM_GHOST]);

		// TOP
		p1[TOP_GHOST] = 2 * p1[TOP_INNER] - p1[TOP_TOP_INNER];
		p1[TOP_TOP_GHOST] = 3 * p1[TOP_INNER] - 2 * p1[TOP_TOP_INNER];

		p2[TOP_GHOST] = 2 * p2[TOP_INNER] - p2[TOP_TOP_INNER];
		p2[TOP_TOP_GHOST] = 3 * p2[TOP_INNER] - 2 * p2[TOP_TOP_INNER];

		if (p2[TOP_INNER] > 0.0){
			p3[TOP_GHOST] = 2 * p3[TOP_INNER] - p3[TOP_TOP_INNER];
			p3[TOP_TOP_GHOST] = 3 * p3[TOP_INNER] - 2 * p3[TOP_TOP_INNER];
		} else {
			p3[TOP_GHOST] = 298.15;
			p3[TOP_TOP_GHOST] = 298.15;
		}

		p4[TOP_GHOST] = 2 * P_atm - p4[TOP_INNER];
		p4[TOP_TOP_GHOST] = 2 * P_atm - p4[TOP_TOP_INNER];

		if (p2[TOP_INNER] > 0.0){
			p5[TOP_GHOST] = fmax(2 * p5[TOP_INNER] - p5[TOP_TOP_INNER], 0.0);
			p5[TOP_TOP_GHOST] = fmax(3 * p5[TOP_INNER] - 2 * p5[TOP_TOP_INNER], 0.0);
		} else {
			p5[TOP_GHOST] = 0.0;
			p5[TOP_TOP_GHOST] = 0.0;
		}

		// rho = P / RT
		float R_mix_T = (1 - p5[TOP_GHOST]) * R_dry + p5[TOP_GHOST] * R_v;
		float R_mix_TT = (1 - p5[TOP_TOP_GHOST]) * R_dry + p5[TOP_TOP_GHOST] * R_v;
		p0[TOP_GHOST] = p4[TOP_GHOST] / (R_mix_T * p3[TOP_GHOST]);
		p0[TOP_TOP_GHOST] = p4[TOP_TOP_GHOST] / (R_mix_TT * p3[TOP_TOP_GHOST]);
	}
}
