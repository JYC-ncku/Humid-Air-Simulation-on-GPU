#include <stdlib.h>
#include <math.h>

float Compute_Cv(float T, float Y){
	float a1, a2, a3, a4, a5; // N2
	float b1, b2, b3, b4, b5; // O2
	float c1, c2, c3, c4, c5; // H2O
	float R_N2 = 296.8; // R of N2 (J/kg.K)
	float R_O2 = 259.84; // R of O2 (J/kg.k)
	float R_H2O = 461.5; // R of H2O (J/kg.k)
	if(T >= 1000.0){
		a1 = 2.926640;
		a2 = 0.14879768e-2;
		a3 = -0.05684760e-5;
		a4 = 0.10097038e-9;
		a5 = -0.06753351e-13;

		b1 = 3.28253784;
		b2 = 1.48308754e-03;
		b3 = -7.57966669e-07;
		b4 = 2.09470555e-10;
		b5 = -2.16717794e-14;

		c1 = 3.03399249;
		c2 = 2.17691804e-03;
		c3 = -1.64072518e-07;
		c4 = -9.70419870e-11;
		c5 = 1.68200992e-14;
	} else {
		a1 = 3.29877;
		a2 = 0.1408204e-2;
		a3 = -0.03963222e-4;
		a4 = 0.05641514e-7;
		a5 = -0.02444854e-10;

		b1 = 3.78245636;
		b2 = -2.99673416e-03;
		b3 = 9.84730201e-06;
		b4 = -9.68129509e-09;
		b5 = 3.24372837e-12;

		c1 = 4.19864056;
		c2 = -2.03643410e-03;
		c3 = 6.52040211e-06;
		c4 = -5.48797062e-09;
		c5 = 1.77197817e-12;
	}
	//Compute Cv
	float Cv_N2 = R_N2 * (a1 + a2*T + (a3 * pow(T,2)) + (a4 * pow(T,3)) + (a5 * pow(T,4))) - R_N2;
	float Cv_O2 = R_O2 * (b1 + b2*T + (b3 * pow(T,2)) + (b4 * pow(T,3)) + (b5 * pow(T,4))) - R_O2;
	float Cv_H2O = R_H2O * (c1 + c2*T + (c3 * pow(T,2)) + (c4 * pow(T,3)) + (c5 * pow(T,4))) - R_H2O;
	float Cv_dry = 0.79 * Cv_N2 + 0.21 * Cv_O2; // dry air = 79% N2 + 21% O2
	float Cv_mix = Cv_dry * (1 - Y) + Y * Cv_H2O;
	return Cv_mix;
}
