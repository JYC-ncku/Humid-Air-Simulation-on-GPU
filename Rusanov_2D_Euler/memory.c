#include <stdlib.h>
#include <stdio.h>

void Allocate_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6, float **array7, float **array8, float **array9,
		     float **array10, float **array11, float **array12, float **array13, float **array14, float **array15, int N_CELLS){
	*array1 = (float*)malloc(N_CELLS * sizeof(float));
	*array2 = (float*)malloc(N_CELLS * sizeof(float));
	*array3 = (float*)malloc(N_CELLS * sizeof(float));
	*array4 = (float*)malloc(N_CELLS * sizeof(float));
	*array5 = (float*)malloc(N_CELLS * sizeof(float));
	*array6 = (float*)malloc(N_CELLS * sizeof(float));
	*array7 = (float*)malloc(N_CELLS * sizeof(float));
	*array8 = (float*)malloc(N_CELLS * sizeof(float));
	*array9 = (float*)malloc(N_CELLS * sizeof(float));
	*array10 = (float*)malloc(N_CELLS * sizeof(float));
	*array11 = (float*)malloc(N_CELLS * sizeof(float));
	*array12 = (float*)malloc(N_CELLS * sizeof(float));
	*array13 = (float*)malloc(N_CELLS * sizeof(float));
	*array14 = (float*)malloc(N_CELLS * sizeof(float));
	*array15 = (float*)malloc(N_CELLS * sizeof(float));
	if(*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL || *array7 == NULL || *array8 == NULL ||
	   *array9 == NULL || *array10 == NULL || *array11 == NULL || *array12 == NULL || *array13 == NULL|| *array14 == NULL || *array15 == NULL){
		printf("Memory allocation failed!\n");
		exit(1);
	}
		printf("Memory allocation successfully for %d elements!\n", N_CELLS);
}

void Free_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6, float **array7, float **array8, float **array9,
		 float **array10, float **array11, float **array12, float **array13, float **array14, float **array15){
	free(*array1);
	free(*array2);
	free(*array3);
	free(*array4);
	free(*array5);
	free(*array6);
	free(*array7);
	free(*array8);
	free(*array9);
	free(*array10);
	free(*array11);
	free(*array12);
	free(*array13);
	free(*array14);
	free(*array15);
	printf("Memory freed successfully!\n");\
}
