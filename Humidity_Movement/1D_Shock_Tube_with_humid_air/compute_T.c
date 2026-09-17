#include <stdlib.h>
#include <stdio.h>
#include <math.h>

float compute_T(float T, float e_target){
	float T_new, Cv, H, e, R; // R is Residual function.
	float a1, a2, a3, a4, a5, a6;
	float R_gas = 296.8; // R of N2 is 296.8 J/(kg*K)
	float Torlerance = 0.001;
	float error = 100.0;

	while(error > Torlerance){
		//Decide whether to use the high-temperature or low-temperature coefficient.
		if(T<=1000){
			a1 = 3.29877;
			a2 = 0.1408204e-2;
			a3 = -0.03963222e-4;
			a4 = 0.05641514e-7;
			a5 = -0.02444854e-10;
			a6 = -0.10208999e4;
		} else{
			a1 = 2.926640;
			a2 = 0.14879768e-2;
			a3 = -0.05684760e-5;
			a4 = 0.10097038e-9;
			a5 = -0.06753351e-13;
			a6 = -0.09227977e4;
		}
		//Compute T
		Cv = R_gas * (a1 + a2*T + (a3 * pow(T,2)) + (a4 * pow(T,3)) + (a5 * pow(T,4))) - R_gas;
		e = Cv * T;
		R = e - e_target;
		T_new = T - (R/Cv);

		//Compute error;
		error = fabs(T_new - T);
		T = T_new;
	}
return T;
}
