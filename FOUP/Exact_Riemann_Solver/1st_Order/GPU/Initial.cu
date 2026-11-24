#include <stdlib.h>
#include <math.h>
#include "Compute_Cv.h"

__global__ void GPU_Initial(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6, float *d_B,
			    float *d_mass, float *d_momentum_X, float *d_momentum_Y, float *d_energy, float *d_mass_fraction,
			    float *MAX_Freq, float *d_MAX_CFL, float R_dry, float R_v, float dx, float dy, int NX, int NY, int N_CELLS){
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+2);
	int j = (int)INDEX - i * (NY+2);
	float X = (i - 0.5) * dx;
	float Y = (j - 0.5) * dy;
	int is_wafer, num_wafers;
	float R_mix, wafer_thickness, total_span_Y, start_Y, pitch_Y, y_bottom, y_top;

	if (INDEX < N_CELLS){
		if (i >= 1 && i < NX + 1 && j >= 1 && j < NY + 1){
			d_p3[INDEX] = 300.0; // T = 300 K
			d_p4[INDEX] = 101325.0; // P = 1 atm

			// If B[INDEX] = 1 ==> walls, else B[INDEX] = 0 ==> Not wall
			if ( (X > 0.0 && X < 0.344 && ( (Y > 0.0 && Y < 0.662) || (Y > 0.982 && Y < 1.55) ) ) ){
				d_B[INDEX] = 1.0;
			} else {
				is_wafer = 0;
				if (X >= 0.039 && X <= 0.339){
					wafer_thickness = dy;
					num_wafers = 10;

					total_span_Y = 0.3; // Total hight of all wafers
					start_Y = 0.68; // Hight of first wafer
					pitch_Y = total_span_Y / num_wafers;
					for (int k = 0; k < num_wafers; k++) {
				        	y_bottom = start_Y + k * pitch_Y;
						y_top = y_bottom + wafer_thickness;
						if (Y >= y_bottom && Y <= y_top) {
							is_wafer = 1;
							break;
						}
					}
				}
				if (is_wafer){
					d_B[INDEX] = 1.0;
				} else{
					d_B[INDEX] = 0.0;
				}
			}

			if (X >= 0.494 && X <= 0.709 && Y == 1.55){
				d_B[INDEX] = 1.0;
			}

			if (d_B[INDEX] == 1.0){
				d_p1[INDEX] = 0.0;
				d_p2[INDEX] = 0.0;
				d_p5[INDEX] = 0.0;
				d_p6[INDEX] = 0.0;
			}

			if(d_B[INDEX] == 0.0){
				d_p1[INDEX] = 0.0; // u = 0 m/s
				d_p2[INDEX] = 0.0; // v = 0 m/s
				d_p5[INDEX] = 0.0; // phi = 0
				d_p6[INDEX] = 0.0; // RH = 0
			}

			R_mix = R_dry * (1 - d_p5[INDEX]) + R_v * d_p5[INDEX];
			d_p0[INDEX] = d_p4[INDEX] / (R_mix * d_p3[INDEX]);
			MAX_Freq[INDEX] = 0.0;
			d_MAX_CFL[INDEX] = 0.0;

			d_mass[INDEX] = d_p0[INDEX];
			d_momentum_X[INDEX] = d_p0[INDEX] * d_p1[INDEX];
			d_momentum_Y[INDEX] = d_p0[INDEX] * d_p2[INDEX];

			float Cv_mix_init = Compute_Cv(d_p3[INDEX], d_p5[INDEX]);
			float e_init = d_p3[INDEX] * Cv_mix_init;
			d_energy[INDEX] = d_mass[INDEX] * (e_init + 0.5 * (d_p1[INDEX] * d_p1[INDEX] + d_p2[INDEX] * d_p2[INDEX]));

			d_mass_fraction[INDEX] = d_mass[INDEX] * d_p5[INDEX];
		}
	}
}


void Initial(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float *d_p5, float *d_p6, float *d_B,
	     float *d_mass, float *d_momentum_X, float *d_momentum_Y, float *d_energy, float *d_mass_fraction,
	     float *MAX_Freq, float *d_MAX_CFL, float R_dry, float R_v, float dx, float dy, int NX, int NY, int N_CELLS){
	     int TPB = 128;
	     int GPB = (TPB + N_CELLS - 1) / TPB;
	     GPU_Initial<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_p4, d_p5, d_p6, d_B,
				       d_mass, d_momentum_X, d_momentum_Y, d_energy, d_mass_fraction,
				       MAX_Freq, d_MAX_CFL, R_dry, R_v, dx, dy, NX, NY, N_CELLS);
}

