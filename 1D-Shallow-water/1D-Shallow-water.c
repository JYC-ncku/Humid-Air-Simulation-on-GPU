#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(double **array1, double **array2, double **array3, int N){
    *array1 = (double*)malloc(N * sizeof(double));
    *array2 = (double*)malloc(N * sizeof(double));
    *array3 = (double*)malloc(N * sizeof(double));
    if(*array1 == NULL || *array2 == NULL || *array3 == NULL){
        printf("Memory allocation faild!\n");
        exit(1);
    }
    printf("Memory allocation successfully for %d elements\n", N);
}

void Free_memory(double *array1, double *array2, double *array3){
    free(array1);
    free(array2);
    free(array3);
    printf("Memory freed successfully!\n");
}

//Setting parameters
    double L = 1.0;
    double g = 9.81;
    double no_steps = 15000;


void Shallow_Water(int N_CELLS, int no_steps, double *x, double *h, double *u, double *mass, double *momentum, double *mass_flux, double *momentum_flux, double *p0, double *p1){
    double dx = L/N_CELLS;
    double PHI = 0.8;
    double S_max = 0.0;

    // Set the initial co conditions
    for (int i=0; i<N_CELLS; i++){
        x[i] = (i+0.5) * dx;
        if (x[i] < 0.5*N_CELLS){
            p0[i] = 1; // Water high = 1m in 0 to 0.5m.
        } else {
            p0[i] = 0.1; // Water high = 1m in 0.5 to 1m.
        }
        p1[i] = 0; // Water speed = 0 m/s in initial condition.
    }

    for (int i=0; i<N_CELLS; i++){
        double S_local = fabs(u[i] + sqrt(g*h[i]));
        if (S_local>S_max){
            S_max = S_local;
        } 
    }
    double dt = PHI * (dx/S_max);

}

int main(){
    const int N_CELLS = 200; // Number of cells
    const int N_INTERFACES = N_CELLS+1; // Number of interfaces
    double *x, *h, *u, *mass, *momentum, *mass_flux, *momentum_flux, *p0, *p1;
}