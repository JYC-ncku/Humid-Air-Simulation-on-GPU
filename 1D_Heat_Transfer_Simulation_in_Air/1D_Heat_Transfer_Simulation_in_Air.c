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

void 

int main(){
    const int N = 200;
    const float alpha = 0.25; 
    const float dt = 0.0005; // CFL must be <= 0.5, and dt cannot be set too large to avoid divergence. 
    const float dx = 0.1/(N-1); // Total length = 0.1 m, and N-1 is because N start from 0.

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