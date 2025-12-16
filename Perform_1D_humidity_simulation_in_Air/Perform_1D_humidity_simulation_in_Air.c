#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Setting parameters
const float L = 0.1
const float hum_init = 0;
const float hum_final = 1;
const float hum_target = 0.25;
const float Pi = M_PI;
const float D = 2.0e-5; //Diffusivity of air.
//Set a maximum number of steps
const int no_steps = 10000;

void Allocate_memory(float **array1, float **array2, float **array3, int N){
    *array1 = (float*)malloc(N * sizeof(float));
    *array2 = (float*)malloc(N * sizeof(float));
    *array3 = (float*)malloc(N * sizeof(float));
    if (*array1 == NULL || *array2 == NULL || *array3 == NULL){
        printf("Memory allocation failed!\n");
        exit(1);
    }
        printf("Memory allcoation successfully for %d elements\n", N);
}

void Free_memory(float *array1, float *array2, float *array3){
    free(array1);
    free(array2);
    free(array3);
    printf("Memory freed successfully!\n");
}

void Humidity_1D(int N, no_steps){
    float *x, *hum, *hum_new;
    float dx = L / (N-1);
    float PHI = 0.25;
    float dt = (PHI * (dx*dx))/D;
    float time = 0; //Iintialize the time.
    
    Allocate_memory(&x, &hum, &hum_new, N);
        for (int i=0; i<N; i++){
            x[i]=i * dx;
            hum[i]=hum_init;
        }
        hum[0]=hum_init;
        hum[N-1]=hum_final;

        while(time < T_final){

            for (int cell = 0; cell<N; cell++){
                hum_new[cell] = hum[cell]
            }
        }
    Free_memory(x, hum, hum_new);
}

int main(){
    const int N=200;
    const int NIF=N+1;
    Humidity_1D(N, no_steps);
    return(0);
}
