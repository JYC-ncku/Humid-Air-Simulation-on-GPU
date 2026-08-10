#include <stdlib.h>
#include <stdio.h>

void Allocate_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6, float **array7, float **array8, float **array9,
		     float **array10, float **array11, float **array12, float **array13, float **array14, float **array15, int N_CELLS){
	float size = N_CELLS * sizeof(float);
	*array1 = (float*)malloc(size);
	*array2 = (float*)malloc(size);
	*array3 = (float*)malloc(size);
	*array4 = (float*)malloc(size);
	*array5 = (float*)malloc(size);
	*array6 = (float*)malloc(size);
	*array7 = (float*)malloc(size);
//	*array8 = (float*)malloc((N_CELLS) * 6 * sizeof(float)); // interface_p have 6 output
//	*array9 = (float*)malloc((N_CELLS) * 5 * sizeof(float)); // flux_X have 5 ouuput
//	*array10 = (float*)malloc((N_CELLS) * 5 * sizeof(float)); // flux_Y have 5 ouuput
	if(*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL || *array7 == NULL){
		printf("Memory allocation failed!\n");
		exit(1);
	}
	printf("Memory allocation successfully for %d elements!\n", N_CELLS);
	//Device
	cudaMalloc((void**)array8, size);
	cudaMalloc((void**)array9, size);
	cudaMalloc((void**)array10, size);
	cudaMalloc((void**)array11, size);
	cudaMalloc((void**)array12, size);
	cudaMalloc((void**)array13, (N_CELLS) * 6 * sizeof(float)); // interface_p have 6 output
	cudaMalloc((void**)array14, (N_CELLS) * 5 * sizeof(float)); // flux_X have 5 ouuput
	cudaMalloc((void**)array15, (N_CELLS) * 5 * sizeof(float)); // flux_Y have 5 ouuput
}

void Free_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6, float **array7, float **array8, float **array9, float **array10,
		 float **array11, float **array12, float **array13, float **array14, float **array15){
	free(*array1);
	free(*array2);
	free(*array3);
	free(*array4);
	free(*array5);
	free(*array6);
	free(*array7);
	cudaFree(*array8);
	cudaFree(*array9);
	cudaFree(*array10);
	cudaFree(*array11);
	cudaFree(*array12);
	cudaFree(*array13);
	cudaFree(*array14);
	cudaFree(*array15);
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
