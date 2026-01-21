#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
void Allocate_memory(double **array1, double **array2, double **array3, double **array4, int N){
    *array1 = (double*)malloc(N * sizeof(double));
    *array2 = (double*)malloc(N * sizeof(double));
    *array3 = (double*)malloc(N * sizeof(double));
    *array4 = (double*)malloc(N * sizeof(double));
    if(*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL){
        printf("Memory allocation failed!\n");
    }
        printf("Memory allocation successfully for %d elements!\n", N);
}

void Free_memory(double *array1, double *array2, double *array3, double *array4){
    free(array1);
    free(array2);
    free(array3);
    free(array4);
    printf("Memory freed successfully!\n");
}