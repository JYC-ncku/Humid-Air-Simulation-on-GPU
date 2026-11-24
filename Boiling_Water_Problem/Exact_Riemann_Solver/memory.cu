#include <stdio.h>
#include <stdlib.h>

void Allocate_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6,
		     float **array7, float **array8, float **array9, float **array10,float **array11, float **array12,
		     float **array13, float **array14, float **array15, float **array16, float **array17, float **array18,
		     float **array19, float **array20, float **array21, float **array22, float **array23, float **array24,
		     int N_CELLS){
	*array1 = (float*)malloc(N_CELLS * sizeof(float));
	*array2 = (float*)malloc(N_CELLS * sizeof(float));
	*array3 = (float*)malloc(N_CELLS * sizeof(float));
	*array4 = (float*)malloc(N_CELLS * sizeof(float));
	*array5 = (float*)malloc(N_CELLS * sizeof(float));
	*array6 = (float*)malloc(N_CELLS * sizeof(float));
	*array7 = (float*)malloc(N_CELLS * sizeof(float));
	if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL ||
	    *array7 == NULL){
		printf("Memory allocation failed!\n");
		exit(1);
	}
	printf("Memory allocation successfully for %d elements!\n", N_CELLS);
	//Device
	cudaError_t Error;
	Error = cudaMalloc((void**)array8, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array8) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array9, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array9) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array10, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array10) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array11, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array11) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array12, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array12) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array13, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array13) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array14, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array14) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array15, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array15) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array16, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array16) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array17, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array17) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array18, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array18) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array19, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array19) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array20, (N_CELLS) * 6 * sizeof(float));     // interface_p have 6 output
	printf("CUDA error (malloc array20) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array21, (N_CELLS) * 6 * sizeof(float));     // flux_X have 6 ouuput
	printf("CUDA error (malloc array21) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array22, (N_CELLS) * 6 * sizeof(float));     // flux_Y have 6 ouuput
	printf("CUDA error (malloc array22) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array23, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array23) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array24, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array24) = %s\n", cudaGetErrorString(Error));
}

void Free_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6,
		 float **array7, float **array8, float **array9, float **array10,float **array11, float **array12,
		 float **array13, float **array14, float **array15, float **array16, float **array17, float **array18,
		 float **array19, float **array20, float **array21, float **array22, float **array23, float **array24){
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
	cudaFree(*array16);
	cudaFree(*array17);
	cudaFree(*array18);
	cudaFree(*array19);
	cudaFree(*array20);
	cudaFree(*array21);
	cudaFree(*array22);
	cudaFree(*array23);
	cudaFree(*array24);
	printf("Memory freed successfully!\n");
}

void Get_From_Device(float **h_a, float **d_a, int N_CELLS){
	size_t size = N_CELLS * sizeof(float);
	cudaMemcpy(*h_a, *d_a, size, cudaMemcpyDeviceToHost);
}
