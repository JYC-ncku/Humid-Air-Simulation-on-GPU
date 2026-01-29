#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(double **array1, double **array2, double **array3, double **array4, double **array5, int N_CELLS){
    *array1 = (double*)malloc(N_CELLS * sizeof(double));
    *array2 = (double*)malloc(N_CELLS * sizeof(double));
    *array3 = (double*)malloc(N_CELLS * sizeof(double));
    *array4 = (double*)malloc(N_CELLS * sizeof(double));
    *array5 = (double*)malloc((N_CELLS+1) * sizeof(double));
    if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL ){
        printf("Memory allocation failed!\n");
        exit(1);
    }
        printf("Memory allocation successfully for %d elements!\n", N_CELLS);
}

void Free_memory(double *array1, double *array2, double *array3, double *array4, double *array5){
    free(array1);
    free(array2);
    free(array3);
    free(array4);
    free(array5);
    printf("Memory freed successfully!\n");
}

void Rusanov_1D_Euler(int N_CELLS, double *x,double *p0, double *p1, double *p2,double *flux){
    float L = 1.0;
    float R = 1.0;
    float GAMMA = 1.4;
    float t = 0;
    float t_FINAL = 0.2;
    double CFL = 0.5;
     double dx = L/N_CELLS;

    // Set initial condition
    for (int i = 0; i < N_CELLS; i++){
        x[i] = (i+0.5) * dx;
        if (i <= N_CELLS/2){
            p0[i] = 10;
            p1[i] = 0;
            p2[i] = 1;
        } else {
            p0[i] = 1;
            p1[i] = 0;
            p2[i] = 1;
        }
    }
}

int main(){
    int N_CELLS = 200;
    int N_INTERFACES = N_CELLS+1;
    double *x, *p0, *p1, *p2, *flux; // p0 is density, p1 is velocity, p2 is temperature
    Allocate_memory(&x, &p0, &p1, &p2, &flux, N_CELLS);
    Free_memory(x, p0, p1, p2, flux);
}