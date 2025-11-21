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

int main(){
    const int N = 200;
    float *T, *T_new, *x ;

    Allocate_memory(&T, &T_new, &x, N);

    Free_memory(T, T_new, x);

    return(0);
}