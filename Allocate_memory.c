#include <stdio.h>
#include <stdlib.h>

void Allocate_memory(float **array1, float **array2, int N){ //*就是go to address的意思，所以*array1就是go to arryay1 address的意思 **array1就是 address array1 go to address，例如"*我"去"*你家"做一件事情，這件事情只有我可以做，然後完成之後想要清掉我做的事情的痕跡，你不需要特別去知道我在你家做了什麼事情，你就是直接把整個你家所有更新的東西全部清調恢復原狀就好，所以只需要Free(*)，不需要Free(**)。
    *array1 = (float*)malloc(N * sizeof(float));
    *array2 = (float*)malloc(N * sizeof(float));

    //Check if memory allocation was successful.
    if (*array1 == NULL || *array2 == NULL) {
        printf("Memory allocation failed!\n");
        exit(1);
    }
    printf("Memory allocation successfully for %d elements\n", N);  //N輸入多少%d的位置就會顯示多少，％d的意思是"整數輸出"，也就是只會顯示整數。
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
    
    return 0; //return 0 表示沒有問題 1或1以上的數字通常代表有錯誤。
}