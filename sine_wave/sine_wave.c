#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(float **array1, float **array2, int N){
    *array1 = (float*)malloc(N * sizeof(float));
    *array2 = (float*)malloc(N * sizeof(float));

    if (*array1 == NULL || *array2 == NULL){
        printf("Memory allocation faild!\n");
        exit(1);
    }
        printf("Memory allocation successfully for %d elements\n", N);
}

void Free_memory(float *array1, float *array2){
    free(array1);
    free(array2);
    printf("Memory freed successfully!\n");
}

int main(){
    const int N = 200;
    float *x, *y;

    Allocate_memory(&x, &y, N);
    for(int i=0; i<N; i++){
        x[i] = 5.0f * i / (N-1); 
        y[i] = sinf(x[i]);
    }
    Free_memory(x, y);
    return 0;
}