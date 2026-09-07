#include <stdlib.h>

__global__ void GPU_Boundary(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, int NX, int NY, int N_CELLS){
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+4);
	int j = (int)INDEX - i * (NY+4);
	int LEFT_LEFT_GHOST = 0 * (NY+4) + j;
	int LEFT_GHOST = 1 * (NY+4) + j;
	int RIGHT_RIGHT_GHOST = (NX+3) * (NY+4) + j;
	int RIGHT_GHOST = (NX+2) * (NY+4) + j;
	int LEFT_INNER = 2 * (NY+4) + j;
	int RIGHT_INNER = (NX+1) * (NY+4) + j;

	int BOTTOM_BOTTOM_GHOST = i * (NY+4) + 0;
	int BOTTOM_GHOST = i * (NY+4) + 1;
	int TOP_TOP_GHOST = i * (NY+4) + (NY+3);
	int TOP_GHOST = i * (NY+4) + (NY+2);
	int BOTTOM_INNER = i * (NY+4) + 2;
	int TOP_INNER = i * (NY+4) + (NY+1);

	if (INDEX < N_CELLS){
		if (j >= 2 && j <= NY+1 && i == 0){
			d_p0[LEFT_GHOST] = d_p0[LEFT_INNER];
			d_p0[LEFT_LEFT_GHOST] = d_p0[LEFT_GHOST];
			d_p0[RIGHT_GHOST] = d_p0[RIGHT_INNER];
			d_p0[RIGHT_RIGHT_GHOST] = d_p0[RIGHT_GHOST];

			d_p1[LEFT_GHOST] = d_p1[LEFT_INNER];
			d_p1[LEFT_LEFT_GHOST] = d_p1[LEFT_GHOST];
			d_p1[RIGHT_GHOST] = d_p1[RIGHT_INNER];
			d_p1[RIGHT_RIGHT_GHOST] = d_p1[RIGHT_GHOST];

			d_p2[LEFT_GHOST] = d_p2[LEFT_INNER];
			d_p2[LEFT_LEFT_GHOST] = d_p2[LEFT_GHOST];
			d_p2[RIGHT_GHOST] = d_p2[RIGHT_INNER];
			d_p2[RIGHT_RIGHT_GHOST] = d_p2[RIGHT_GHOST];

			d_p3[LEFT_GHOST] = d_p3[LEFT_INNER];
			d_p3[LEFT_LEFT_GHOST] = d_p3[LEFT_GHOST];
			d_p3[RIGHT_GHOST] = d_p3[RIGHT_INNER];
			d_p3[RIGHT_RIGHT_GHOST] = d_p3[RIGHT_GHOST];

			d_p4[LEFT_GHOST] = d_p4[LEFT_INNER];
			d_p4[LEFT_LEFT_GHOST] = d_p4[LEFT_GHOST];
			d_p4[RIGHT_GHOST] = d_p4[RIGHT_INNER];
			d_p4[RIGHT_RIGHT_GHOST] = d_p4[RIGHT_GHOST];
		}
		//BOTTOM and TOP
		if (i >= 2 && i <= NX+1 && j == 0){
			d_p0[BOTTOM_GHOST] = d_p0[BOTTOM_INNER];
			d_p0[BOTTOM_BOTTOM_GHOST] = d_p0[BOTTOM_GHOST];
			d_p0[TOP_GHOST] = d_p0[TOP_INNER];
			d_p0[TOP_TOP_GHOST] = d_p0[TOP_GHOST];

			d_p1[BOTTOM_GHOST] = d_p1[BOTTOM_INNER];
			d_p1[BOTTOM_BOTTOM_GHOST] = d_p1[BOTTOM_GHOST];
			d_p1[TOP_GHOST] = d_p1[TOP_INNER];
			d_p1[TOP_TOP_GHOST] = d_p1[TOP_GHOST];

			d_p2[BOTTOM_GHOST] = d_p2[BOTTOM_INNER];
			d_p2[BOTTOM_BOTTOM_GHOST] = d_p2[BOTTOM_GHOST];
			d_p2[TOP_GHOST] = d_p2[TOP_INNER];
			d_p2[TOP_TOP_GHOST] = d_p2[TOP_GHOST];

			d_p3[BOTTOM_GHOST] = d_p3[BOTTOM_INNER];
			d_p3[BOTTOM_BOTTOM_GHOST] = d_p3[BOTTOM_GHOST];
			d_p3[TOP_GHOST] = d_p3[TOP_INNER];
			d_p3[TOP_TOP_GHOST] = d_p3[TOP_GHOST];

			d_p4[BOTTOM_GHOST] = d_p4[BOTTOM_INNER];
			d_p4[BOTTOM_BOTTOM_GHOST] = d_p4[BOTTOM_GHOST];
			d_p4[TOP_GHOST] = d_p4[TOP_INNER];
			d_p4[TOP_TOP_GHOST] = d_p4[TOP_GHOST];
			// Reflect boundary
//			d_p2[BOTTOM_GHOST] = -d_p2[BOTTOM_INNER];
//			d_p2[TOP_GHOST] = -d_p2[TOP_INNER];
		}
	}
}

void Boundary(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Boundary<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_p4, NX, NY, N_CELLS);
}
