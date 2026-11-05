#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Compute_Cv.h"

__device__ float compute_T(float T, float Y, float e_target){
	float T_new, e, R; // R is Residual function.
	float Torlerance = 0.001;
	float error = 100.0;
	int iter = 0;
	while(error > Torlerance){
		float Cv_mix = Compute_Cv(T, Y);
		e = Cv_mix * T;
		R = e - e_target;
		T_new = T - (R/Cv_mix);

		//Compute error;
		error = fabs(T_new - T);
		T = T_new;
		iter++;
		if (iter > 1000) {
			printf("STOP!\n");
	        break;
		}
	}
return T;
}
