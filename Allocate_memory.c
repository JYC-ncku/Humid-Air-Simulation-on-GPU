#include <stdio.h>
#include <stdlib.h>

void Allocate_memory(float **array1, float **array2, int N){
    *array1 = (float*)malloc(N * sizeof(float));
    *array2 = (float*)malloc(N * sizeof(float));

    //Check if memory allocation was successful.
    if (*array1 == NULL || *array2 == NULL) {
        printf("Memory allocation failed!\n");
        exit(1);
    }
    printf("Memory allocation successfully for %d elements\n", N);  //N輸入多少%d的位置就會顯示多少，％d的意思是<整數輸出，也就是只會顯示整數。
}

void Free_Memory( float *array1, float *array2){
        free(array1);
        free(array2);
        printf("Memory freed successfully\n");
}

int main() {
    const int N =200; //Setting cells number.
    float *F, *u; // Tell the computer that F and u are POINTER to float. 
 
    Allocate_memory(&F, &u, N); //Tell computer array1 and 2 gets &F and &u (address of F and u) to modify F and u directly, and they all need N * sizeof(float) space to store data.
    

    //Initialize array
    for (int i = 0; i<N; i++){
        F[i]=0.0f;
        u[i]=0.0f;
    }

    Free_Memory(F, u);
    
    return 0;
}