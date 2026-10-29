#include <stdlib.h>

void Boundary(float *p0, float *p1, float *p2, float *p3, float *p4, int NX, int NY){
	for (int j = 1 ; j <= NY; j++){
		int LEFT_GHOST = 0 * (NY+2) + j;
		int RIGHT_GHOST = (NX+1) * (NY+2) + j;
		int LEFT_INNER = 1 * (NY+2) + j;
		int RIGHT_INNER = NX * (NY+2) + j;
		p0[LEFT_GHOST] = p0[LEFT_INNER];
		p0[RIGHT_GHOST] = p0[RIGHT_INNER];
		p1[LEFT_GHOST] = p1[LEFT_INNER];
		p1[RIGHT_GHOST] = p1[RIGHT_INNER];
		p2[LEFT_GHOST] = p2[LEFT_INNER];
		p2[RIGHT_GHOST] = p2[RIGHT_INNER];
		p3[LEFT_GHOST] = p3[LEFT_INNER];
		p3[RIGHT_GHOST] = p3[RIGHT_INNER];
		p4[LEFT_GHOST] = p4[LEFT_INNER];
		p4[RIGHT_GHOST] = p4[RIGHT_INNER];
	}
	//BOTTOM and TOP
	for (int i = 1 ; i <= NX; i++){
		int BOTTOM_GHOST = i * (NY+2) + 0;
		int TOP_GHOST = i * (NY+2) + (NY+1);
		int BOTTOM_INNER = i * (NY+2) + 1;
		int TOP_INNER = i * (NY+2) + NY;
		p0[BOTTOM_GHOST] = p0[BOTTOM_INNER];
		p0[TOP_GHOST] = p0[TOP_INNER];
		p1[BOTTOM_GHOST] = p1[BOTTOM_INNER];
		p1[TOP_GHOST] = p1[TOP_INNER];
		p2[BOTTOM_GHOST] = p2[BOTTOM_INNER];
		p2[TOP_GHOST] = p2[TOP_INNER];
		p3[BOTTOM_GHOST] = p3[BOTTOM_INNER];
		p3[TOP_GHOST] = p3[TOP_INNER];
		p4[BOTTOM_GHOST] = p4[BOTTOM_INNER];
		p4[TOP_GHOST] = p4[TOP_INNER];
		// Reflect boundary
		p2[BOTTOM_GHOST] = -p2[BOTTOM_INNER];
		p2[TOP_GHOST] = -p2[TOP_INNER];
	}
}
