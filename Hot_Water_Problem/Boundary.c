#include <stdlib.h>
#include <math.h>

void Boundary(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float Ru, float R_v, float R_dry, int NX, int NY){
	float P_atm = 1013.25;  // 總環境壓力
	// LEFT and RIGHT
	for (int j = 2 ; j <= NY+1; j++){
		int LEFT_LEFT_GHOST = 0 * (NY+4) + j;
		int LEFT_GHOST = 1 * (NY+4) + j;
		int RIGHT_RIGHT_GHOST = (NX+3) * (NY+4) + j;
		int RIGHT_GHOST = (NX+2) * (NY+4) + j;

		int LEFT_INNER = 2 * (NY+4) + j;
//		int LEFT_LEFT_INNER = 3 * (NY+4) + j;
		int RIGHT_INNER = (NX+1) * (NY+4) + j;
//		int RIGHT_RIGHT_INNER = NX * (NY+4) + j;

		// Left hand side is Inflow and Right hand side is Outflow
		p1[LEFT_GHOST] = 10.0; // Inflow u = 10 m/s
		p1[LEFT_LEFT_GHOST] = p1[LEFT_GHOST];
		p1[RIGHT_GHOST] = p1[RIGHT_INNER];
		p1[RIGHT_RIGHT_GHOST] = p1[RIGHT_INNER];

		p2[LEFT_GHOST] = 0; // Inflow v = 0 m/s
		p2[LEFT_LEFT_GHOST] = p2[LEFT_GHOST];
		p2[RIGHT_GHOST] = p2[RIGHT_INNER];
		p2[RIGHT_RIGHT_GHOST] = p2[RIGHT_INNER];

		p3[LEFT_GHOST] = 300.0; // Inflow T = 300 K
		p3[LEFT_LEFT_GHOST] = p3[LEFT_GHOST];
		p3[RIGHT_GHOST] = p3[RIGHT_INNER];
		p3[RIGHT_RIGHT_GHOST] = p3[RIGHT_INNER];

		p4[LEFT_GHOST] = p4[LEFT_INNER]; // Inflow P = 1 atm = 101325 Pa (Kg/m-s^2)
		p4[LEFT_LEFT_GHOST] = p4[LEFT_INNER];
		p4[RIGHT_GHOST] = P_atm;
		p4[RIGHT_RIGHT_GHOST] = P_atm;

		p5[LEFT_GHOST] = 0.0;
		p5[LEFT_LEFT_GHOST] = p5[LEFT_GHOST];
		p5[RIGHT_GHOST] = p5[RIGHT_INNER];
		p5[RIGHT_RIGHT_GHOST] = p5[RIGHT_INNER];

		float R_mix_L = (1 - p5[LEFT_GHOST]) * R_dry + p5[LEFT_GHOST] * R_v;
		float R_mix_LL = (1 - p5[LEFT_LEFT_GHOST]) * R_dry + p5[LEFT_LEFT_GHOST] * R_v;
		float R_mix_R = (1 - p5[RIGHT_GHOST]) * R_dry + p5[RIGHT_GHOST] * R_v;
		float R_mix_RR = (1 - p5[RIGHT_RIGHT_GHOST]) * R_dry + p5[RIGHT_RIGHT_GHOST] * R_v;
		p0[LEFT_GHOST] = p4[LEFT_GHOST] / (R_mix_L * p3[LEFT_GHOST]);
		p0[LEFT_LEFT_GHOST] = p4[LEFT_LEFT_GHOST] / (R_mix_LL * p3[LEFT_LEFT_GHOST]);
		p0[RIGHT_GHOST] = p4[RIGHT_GHOST] / (R_mix_R * p3[RIGHT_GHOST]);
		p0[RIGHT_RIGHT_GHOST] = p4[RIGHT_RIGHT_GHOST] / (R_mix_RR * p3[RIGHT_RIGHT_GHOST]);
	}

	float T_base = 373.15;  // The bottom base is maintained at 373.15 K
	float Phi_base = 1.0;  // Relative humidtiy at bottom (100%)

	float T_C = T_base - 273.15; // Back to Celsius degree
	float P_sat = 610.78 * exp((17.27 * T_C) / (T_C + 237.3)); // Tetens equation
	float P_v = Phi_base * P_sat;
	if (P_v > P_atm) {
		P_v = P_atm;
	} // 防呆機制：水氣分壓不能超過總壓
	float Y_base = (P_v / R_v) / ((P_v / R_v) + ((P_atm - P_v) / R_dry));

	//BOTTOM and TOP
	for (int i = 2 ; i <= NX+1; i++){
		int BOTTOM_BOTTOM_GHOST = i * (NY+4) + 0;
		int BOTTOM_GHOST = i * (NY+4) + 1;
		int TOP_TOP_GHOST = i * (NY+4) + (NY+3);
		int TOP_GHOST = i * (NY+4) + (NY+2);

		int BOTTOM_INNER = i * (NY+4) + 2;
		int BOTTOM_BOTTOM_INNER = i * (NY+4) + 3;
		int TOP_INNER = i * (NY+4) + (NY+1);
//		int TOP_TOP_INNER = i * (NY+4) + NY;

		// BOTTOM
		p1[BOTTOM_GHOST] = -p1[BOTTOM_INNER];
		p1[BOTTOM_BOTTOM_GHOST] = -p1[BOTTOM_BOTTOM_INNER];

		p2[BOTTOM_GHOST] = -p2[BOTTOM_INNER];
		p2[BOTTOM_BOTTOM_GHOST] = -p2[BOTTOM_BOTTOM_INNER];

		p3[BOTTOM_GHOST] = 2 * T_base - p3[BOTTOM_INNER];
		p3[BOTTOM_BOTTOM_GHOST] = p3[BOTTOM_GHOST];

		p4[BOTTOM_GHOST] = p4[BOTTOM_INNER];
		p4[BOTTOM_BOTTOM_GHOST] = p4[BOTTOM_BOTTOM_INNER];

		p5[BOTTOM_GHOST] = fmax(0.0, fmin(2.0 * Y_base - p5[BOTTOM_INNER], 1.0));
		p5[BOTTOM_BOTTOM_GHOST] = p5[BOTTOM_GHOST];
//		p5[BOTTOM_GHOST] = 1.0; // The max relative humidity (phi_max) setting 1
//		p5[BOTTOM_BOTTOM_GHOST] = p5[BOTTOM_GHOST];

		// rho = P / RT
		float R_mix_B = (1 - p5[BOTTOM_GHOST]) * R_dry + p5[BOTTOM_GHOST] * R_v;
		float R_mix_BB = (1 - p5[BOTTOM_BOTTOM_GHOST]) * R_dry + p5[BOTTOM_BOTTOM_GHOST] * R_v;
		p0[BOTTOM_GHOST] = p4[BOTTOM_GHOST] / (R_mix_B * p3[BOTTOM_GHOST]);
		p0[BOTTOM_BOTTOM_GHOST] = p4[BOTTOM_BOTTOM_GHOST] / (R_mix_BB * p3[BOTTOM_BOTTOM_GHOST]);

		// TOP
		p1[TOP_GHOST] = p1[TOP_INNER];
		p1[TOP_TOP_GHOST] = p1[TOP_INNER];

		p2[TOP_GHOST] = p2[TOP_INNER];
		p2[TOP_TOP_GHOST] = p2[TOP_INNER];

		p3[TOP_GHOST] = p3[TOP_INNER];
		p3[TOP_TOP_GHOST] = p3[TOP_INNER];

		p4[TOP_GHOST] = P_atm;
		p4[TOP_TOP_GHOST] = P_atm;

		p5[TOP_GHOST] = p5[TOP_INNER];
		p5[TOP_TOP_GHOST] = p5[TOP_INNER];

		// rho = P / RT
		float R_mix_T = (1 - p5[TOP_GHOST]) * R_dry + p5[TOP_GHOST] * R_v;
		float R_mix_TT = (1 - p5[TOP_TOP_GHOST]) * R_dry + p5[TOP_TOP_GHOST] * R_v;
		p0[TOP_GHOST] = p4[TOP_GHOST] / (R_mix_T * p3[TOP_GHOST]);
		p0[TOP_TOP_GHOST] = p4[TOP_TOP_GHOST] / (R_mix_TT * p3[TOP_TOP_GHOST]);
	}
}
