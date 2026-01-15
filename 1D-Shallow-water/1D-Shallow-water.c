#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(double **array1, double **array2, double **array3, double **array4, double **array5, double **array6, double **array7, double **array8, double **array9, int N){
    *array1 = (double*)malloc(N * sizeof(double));
    *array2 = (double*)malloc(N * sizeof(double));
    *array3 = (double*)malloc(N * sizeof(double));
    *array4 = (double*)malloc(N * sizeof(double));
    *array5 = (double*)malloc(N * sizeof(double));
    *array6 = (double*)malloc(N * sizeof(double));
    *array7 = (double*)malloc(N * sizeof(double));
    *array8 = (double*)malloc(N * sizeof(double));
    *array9 = (double*)malloc(N * sizeof(double));
    if(*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL || *array7 == NULL || *array8 == NULL || *array9 == NULL){
        printf("Memory allocation faild!\n");
        exit(1);
    }
    printf("Memory allocation successfully for %d elements\n", N);
}

void Free_memory(double *array1, double *array2, double *array3, double *array4, double *array5, double *array6, double *array7, double *array8, double *array9){
    free(array1);
    free(array2);
    free(array3);
    printf("Memory freed successfully!\n");
}

//Setting parameters
    double L = 1.0;
    double g = 9.81;
    const int no_steps = 15000;


void Shallow_Water(int N_CELLS, int no_steps, double *x, double *h, double *u, double *mass, double *momentum, double *mass_flux, double *momentum_flux, double *p0, double *p1){
    double dx = L/N_CELLS;
    double CFL = 0.8;
    double S_max = 0.0;

    // Set the initial co conditions
    for (int i=0; i<N_CELLS; i++){
        x[i] = (i+0.5) * dx;
        if (i < N_CELLS/2){
            p0[i] = 1; // Initial water high = 1m in 0 to 0.5m.
        } else {
            p0[i] = 0.1; // Initial water high = 1m in 0.5 to 1m (Shallow Water).
        }
        p1[i] = 0; // Water speed = 0 m/s in initial condition.
        mass[i] = p0[i];
        momentum[i] = p0[i] * p1[i];
    }

    // Set CFL
    for (int i=0; i<N_CELLS; i++){
        double S_local = fabs(u[i] + sqrt(g*h[i]));
        if (S_local>S_max){
            S_max = S_local;
        } 
    }
    double dt = CFL * (dx/S_max);

}

int main(){
    const int N_CELLS = 200; // Number of cells
    const int N_INTERFACES = N_CELLS+1; // Number of interfaces
    double *x, *h, *u, *mass, *momentum, *mass_flux, *momentum_flux, *p0, *p1; // x is position, h is water high, u is water velocity, po is primitive values of h, p1 is primitive values of u
}