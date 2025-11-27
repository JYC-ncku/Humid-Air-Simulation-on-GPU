#include <stdio.h>
#include <stdlib.h>
#include <math.h>
    
void Allocate_memory(float **array1, float **array2, float **array3, int N){    
    *array1 = (float*)malloc(N * sizeof(float));
    *array2 = (float*)malloc(N * sizeof(float));
    *array3 = (float*)malloc(N * sizeof(float));

    if (*array1 == NULL || *array2 == NULL || *array3 == NULL){
        printf("Memory allocation failed!\n");
        exit(1);
    }
    printf("Memory allocation successfully for %d elements\n", N);
}

void Free_memory( float *array1, float *array2, float *array3){
        free(array1);
        free(array2);
        free(array3);
        printf("Memory freed successfully\n");
}

//Calculate time when T=325K at x = 5cm.
int main(){
    const int N = 200;
    const float k = 0.026; // Thermal conductivity of air is 0.026 W/m K
    const float R = 287; // Gas constant of air is 287 J/kg K
    const float rho = 1.23; // density of air is 1.23 kg/m^3
    const float gamma = 1.4; // Gamma value of air is 1.4
    const float Cv = R/(gamma-1); // specific heat at constant volume of air
    const float Cp = Cv+R; // specific heat at constant pressure of air
    
    float dt = ; 
    float dx = 0.1/(N-1); // Total length = 0.1 m, and N-1 is because N start from 0.

    float *T, *T_new, *x ;
    float time = 0.0;
    int i;

    Allocate_memory(&T, &T_new, &x, N);

    // Setting initial condition
    for (i = 0; i<N; i++){
        T[i] = 300;
        x[i] = i * dx ;
    }

    Free_memory(T, T_new, x);

    return(0);
}