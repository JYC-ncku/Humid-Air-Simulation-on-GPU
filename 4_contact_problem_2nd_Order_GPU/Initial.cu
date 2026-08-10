#include <stdlib.h>

__global__ void GPU_Initial(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float R, int NX, int NY, int N_CELLS){
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	if (INDEX < N_CELLS){
		int i = (int)INDEX / (NY+4);
		int j = (int)INDEX - i * (NY+4);
		//Quadrant I (region A)
		if (i >= ((NX/2) + 1) && j >= ((NY/2) + 1)){
			p0[INDEX] = 1.0;
			p1[INDEX] = 0.75;
			p2[INDEX] = -0.5;
			p4[INDEX] = 1.0;
		//Quadrant II (region B)
		} else if (i < ((NX/2) + 1) && j >= ((NY/2) + 1)){
			p0[INDEX] = 2.0;
			p1[INDEX] = 0.75;
			p2[INDEX] = 0.5;
			p4[INDEX] = 1.0;
		//Quadrant III (region C)
		} else if (i < ((NX/2) + 1) && j < ((NY/2) + 1)){
			p0[INDEX] = 1.0;
			p1[INDEX] = -0.75;
			p2[INDEX] = 0.5;
			p4[INDEX] = 1.0;
		//Quadrant IV (region D)
		} else if (i >= ((NX/2) + 1) && j < ((NY/2) + 1)){
			p0[INDEX] = 3.0;
			p1[INDEX] = -0.75;
			p2[INDEX] = -0.5;
			p4[INDEX] = 1.0;
		}
		p3[INDEX] = p4[INDEX] / (R * p0[INDEX]);
}

void Initial(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float R, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (N_CELLS + TPB - 1) / TPB;
	GPU_Initial<<<GPB,TPB>>>(d_p0, d_p1, d_p2, d_p3, d_p4, R, NX, NY, N_CELLS);
}
