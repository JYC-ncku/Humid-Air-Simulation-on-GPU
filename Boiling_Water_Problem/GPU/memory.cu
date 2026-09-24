#include <stdio.h>
#include <stdlib.h>

void Allocate_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6,
		     float **array7, float **array8, float **array9, float **array10,float **array11, float **array12,
		     float **array13, float **array14, float **array15, float **array16, float **array17, float **array18,
		     float **array19, float **array20, float **array21, float **array22, float **array23, float **array24,
		     float **array25, float **array26, float **array27, float **array28, float **array29, float **array30,
		     float **array31, float **array32, float **array33, float **array34, float **array35, float **array36,
		     float **array37, float **array38, float **array39, float **array40, float **array41, float **array42,
		     float **array43, float **array44, float **array45, int N_CELLS){
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
	*array16 = (float*)malloc(N_CELLS * sizeof(float));
	*array17 = (float*)malloc(N_CELLS * sizeof(float));
	*array18 = (float*)malloc(N_CELLS * sizeof(float));
	*array19 = (float*)malloc(N_CELLS * sizeof(float));
	*array20 = (float*)malloc(N_CELLS * sizeof(float));
	*array21 = (float*)malloc(N_CELLS * sizeof(float));
	*array22 = (float*)malloc(N_CELLS * sizeof(float));
	if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL ||
	    *array7 == NULL || *array8 == NULL || *array9 == NULL || *array10 == NULL || *array11 == NULL || *array12 == NULL ||
	    *array13 == NULL || *array14 == NULL || *array15 == NULL || *array16 == NULL || *array17 == NULL || *array18 == NULL ||
	    *array19 == NULL || *array20 == NULL || *array21 == NULL || *array22 == NULL){
		printf("Memory allocation failed!\n");
		exit(1);
	}
	printf("Memory allocation successfully for %d elements!\n", N_CELLS);
	//Device
	cudaError_t Error;
	Error = cudaMalloc((void**)array23, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array23) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array24, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array24) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array25, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array25) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array26, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array26) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array27, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array27) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array28, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array28) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array29, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array29) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array30, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array30) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array31, (N_CELLS+1) * sizeof(float));
	printf("CUDA error (malloc array31) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array32, (N_CELLS+1) * sizeof(float));
	printf("CUDA error (malloc array32) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array33, (N_CELLS+1) * sizeof(float));
	printf("CUDA error (malloc array33) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array34, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array34) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array35, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array35) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array36, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array36) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array37, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array37) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array38, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array38) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array39, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array39) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array40, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array40) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array41, (N_CELLS) * sizeof(float));
	printf("CUDA error (malloc array41) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array42, (N_CELLS+1) * sizeof(float));
	printf("CUDA error (malloc array42) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array43, (N_CELLS+1) * sizeof(float));
	printf("CUDA error (malloc array43) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array44, (N_CELLS+1) * sizeof(float));
	printf("CUDA error (malloc array44) = %s\n", cudaGetErrorString(Error));
	Error = cudaMalloc((void**)array45, (N_CELLS+1) * sizeof(float));
	printf("CUDA error (malloc array45) = %s\n", cudaGetErrorString(Error));
}

void Free_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6,
		 float **array7, float **array8, float **array9, float **array10,float **array11, float **array12,
		 float **array13, float **array14, float **array15, float **array16, float **array17, float **array18,
		 float **array19, float **array20, float **array21, float **array22, float **array23, float **array24,
		 float **array25, float **array26, float **array27, float **array28, float **array29, float **array30,
		 float **array31, float **array32, float **array33, float **array34, float **array35, float **array36,
		 float **array37, float **array38, float **array39, float **array40, float **array41, float **array42,
		 float **array43, float **array44, float **array45){
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
	free(*array16);
	free(*array17);
	free(*array18);
	free(*array19);
	free(*array20);
	free(*array21);
	free(*array22);
	cudaFree(*array23);
	cudaFree(*array24);
	cudaFree(*array25);
	cudaFree(*array26);
	cudaFree(*array27);
	cudaFree(*array28);
	cudaFree(*array29);
	cudaFree(*array30);
	cudaFree(*array31);
	cudaFree(*array32);
	cudaFree(*array33);
	cudaFree(*array34);
	cudaFree(*array35);
	cudaFree(*array36);
	cudaFree(*array37);
	cudaFree(*array38);
	cudaFree(*array39);
	cudaFree(*array40);
	cudaFree(*array41);
	cudaFree(*array42);
	cudaFree(*array43);
	cudaFree(*array44);
	cudaFree(*array45);
	printf("Memory freed successfully!\n");
}

void Sent_To_Device(float **d_a, float **h_a, int N_CELLS){
	size_t size = N_CELLS * sizeof(float);
	cudaMemcpy(*d_a, *h_a, size, cudaMemcpyHostToDevice);
}

void Get_From_Device(float **h_a, float **d_a, int N_CELLS){
	size_t size = N_CELLS * sizeof(float);
	cudaMemcpy(*h_a, *d_a, size, cudaMemcpyDeviceToHost);
}
