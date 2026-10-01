#include <stdlib.h>

void Boundary(float *p0, float *p1, float *p2, float *p3, float *p4, int N_CELLS){
	p0[1] = p0[2];
	p1[1] = -p1[2];
	p2[1] = p2[2];
	p3[1] = p3[2];
	p4[1] = p4[2];

	p0[0] = p0[1];
	p1[0] = p1[1];
	p2[0] = p2[1];
	p3[0] = p3[1];
	p4[0] = p4[1];

	p0[N_CELLS+2] = p0[N_CELLS+1];
	p1[N_CELLS+2] = -p1[N_CELLS+1];
	p2[N_CELLS+2] = p2[N_CELLS+1];
	p3[N_CELLS+2] = p3[N_CELLS+1];
	p4[N_CELLS+2] = p4[N_CELLS+1];

	p0[N_CELLS+3] = p0[N_CELLS+2];
	p1[N_CELLS+3] = p1[N_CELLS+2];
	p2[N_CELLS+3] = p2[N_CELLS+2];
	p3[N_CELLS+3] = p3[N_CELLS+2];
	p4[N_CELLS+3] = p4[N_CELLS+2];
}
