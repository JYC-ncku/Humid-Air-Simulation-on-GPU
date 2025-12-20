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

void Allocate_memory(float **array1, float **array2, float **array3, float **array4, int N){
    *array1 = (float*)malloc(N * sizeof(float));
    *array2 = (float*)malloc(N * sizeof(float));
    *array3 = (float*)malloc(N * sizeof(float));
    *array4 = (float*)malloc((N+1) * sizeof(float));
    if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL){
        printf("Memory allocation failed!\n");
        exit(1);
    }
        printf("Memory allcoation successfully for %d elements\n", N);
}

void Free_memory(float *array1, float *array2, float *array3, float *array4){
    free(array1);
    free(array2);
    free(array3);
    free(array4);
    printf("Memory freed successfully!\n");
}
/*
    ----------------------------
    |        |        |        |
    |   0    |  ...   |   N-1  |   ===>   Have total N cells and N+1 interface. (because N is from 0)
    |        |        |        | 
    ----------------------------
    0        1       N-1       N
*/
void Humidity_1D(int N, int no_steps){
    float *x, *hum, *hum_new, *F;
    float dx = L / N;
    float PHI = 0.25;
    float dt = PHI * ((dx*dx)/D);
    float time = 0; //Iintialize the time.
    
    Allocate_memory(&x, &hum, &hum_new, &F, N);
        //Initialization
        for (int i=0; i<N; i++){
            x[i]=(i+0.5) * dx; //計算每個cell中間的flux
            hum[i]=hum_init; //每個cell一開始的濕度都為0
        }

        //Calculate interface flux
        //中間邊界flux
    for (int step=0; step < no_steps; step++){
        for (int j = 1; j < N; j++) {
        F[j] = -D * (hum[j] - hum[j-1]) / dx;
        }
        
        F[0] = -D * (hum[0] - 0.0) / (0.5 * dx); //左邊界flux
        F[N] = -D * (1 - hum[N-1]) / (0.5 * dx); //右邊界flux

        /*    Solving
        dh/dt + dF/dx = 0
    So
        h* = h - dt*dF/dx
        */
        for (int cell = 0; cell<N; cell++){
                hum_new[cell] = hum[cell]-(dt/dx)*(F[cell+1]-F[cell]);
            }
        for (int i = 0; i < N; i++){
            hum[i]=hum_new[i];
        }
        time = time+dt
        }
    Free_memory(x, hum, hum_new, F);
}

int main(){
    const int N=200; //Number of cell
    Humidity_1D(N, no_steps);
    return(0);
}
