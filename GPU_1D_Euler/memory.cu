#include <stdlib.h>
#include <stdio.h>

void Allocate_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6,
                     float **array7, float **array8, float **array9, float **array10, float **array11, int N_CELLS){
	//HOST
	*array1 = (float*)malloc(N_CELLS * sizeof(float));
	*array2 = (float*)malloc(N_CELLS * sizeof(float));
	*array3 = (float*)malloc(N_CELLS * sizeof(float));
	*array4 = (float*)malloc(N_CELLS * sizeof(float));
	*array5 = (float*)malloc(N_CELLS * sizeof(float));
	*array6 = (float*)malloc(N_CELLS * sizeof(float));
	*array7 = (float*)malloc(N_CELLS * sizeof(float));
	*array8 = (float*)malloc(N_CELLS * sizeof(float));
	if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL ||
	    *array7 == NULL || *array8 == NULL){
		printf("Memory allocation failed!\n");
		exit(1);
	}
	printf("Memory allocation successfully for %d elements!\n", N_CELLS);
	//Device
	cudaError_t Error;
	Error = cudaMalloc((void**)array9, (N_CELLS+1) * sizeof(float));
	printf("CUDA error (malloc array9) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array10, (N_CELLS+1) * sizeof(float));
	printf("CUDA error (malloc array10) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array11, (N_CELLS+1) * sizeof(float));
	printf("CUDA error (malloc array11) = %s\n", cudaGetErrorString(Error));
}

void Free_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6,
                 float **array7, float **array8, float **array9, float **array10, float **array11){
    free(*array1);
    free(*array2);
    free(*array3);
    free(*array4);
    free(*array5);
    free(*array6);
    free(*array7);
    free(*array8);
    cudaFree(*array9);
    cudaFree(*array10);
    cudaFree(*array11);
    printf("Memory freed successfully!\n");
}

void Send_To_Device(float **d_a, float **h_a, int N_CELLS){
	size_t size = N_CELLS * sizeof(float);
	cudaMemcpy(*d_a, *h_a, size, cudaMemcpyHostToDevice);
}

void Get_From_Device(float **h_a, float **d_a, int N_CELLS){
	size_t size = N_CELLS * sizeof(float);
	cudaMemcpy(*h_a, *d_a, size, cudaMemcpyDeviceToHost);
}
