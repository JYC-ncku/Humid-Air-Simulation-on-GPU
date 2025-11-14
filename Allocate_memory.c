#include <studio.h>
#include <stdib.h>

int main() {
    const int N =200; //Setting cells number.
    float *F; // Tell the computer that F is a POINTER to float 
    F = (float*)malloc(N*sizeof(float));   //Tell the computer I need ""N*sizeof(float)"" bytes of memory space for F to store data."

    free(F);

    return 0;
}