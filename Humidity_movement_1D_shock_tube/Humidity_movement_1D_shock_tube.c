#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(double **array1, double **array2, double **array3, double **array4, double **array5, double **array6,
                     double **array7, double **array8, double **array9, double **array10, int N_CELLS){
    *array1 = (double*)malloc(N_CELLS * sizeof(double));
    *array2 = (double*)malloc(N_CELLS * sizeof(double));
    *array3 = (double*)malloc(N_CELLS * sizeof(double));
    *array4 = (double*)malloc(N_CELLS * sizeof(double));
    *array5 = (double*)malloc(N_CELLS * sizeof(double));
    *array6 = (double*)malloc(N_CELLS * sizeof(double));
    *array7 = (double*)malloc(N_CELLS * sizeof(double));
    *array8 = (double*)malloc(N_CELLS * sizeof(double));
    *array9 = (double*)malloc(N_CELLS * sizeof(double));
    *array10 = (double*)malloc(N_CELLS * sizeof(double));
    if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL ||
        *array7 == NULL || *array8 == NULL || *array9 == NULL || *array10 == NULL){
        printf("Memory allocation failed!\n");
        exit(1);
    }
        printf("Memory allocation successfully for %d elements\n", N_CELLS);
}

void Free_memory(double *array1, double *array2, double *array3, double *array4, double *array5, double *array6,
                  double *array7, double *array8, double *array9, double *array10){
    free(array1);
    free(array2);
    free(array3);
    free(array4);
    free(array5);
    free(array6);
    free(array7);
    free(array8);
    free(array9);
    free(array10);
    printf("Memory freed successfully\n");
}

int main(){
    int N_CELLS = 200;
    
    Allocate_memory();
    Free_memory();
    return 0;
}