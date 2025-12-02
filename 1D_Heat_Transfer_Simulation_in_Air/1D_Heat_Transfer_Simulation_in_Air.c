#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
// Setting parameters   
    const float k = 0.026; // Thermal conductivity of air is 0.026 W/m K
    const float R = 287; // Gas constant of air is 287 J/kg K
    const float rho = 1.23; // density of air is 1.23 kg/m^3
    const float gamma = 1.4; // Gamma value of air is 1.4

// 程式控制常數
    const float L = 0.10;      // Total length 10cm = 0.10m
    const float Tt = 325.0; // 目標（Target)溫度 K
    const float Ti = 300.0; // 初始(initial)溫度 K
    const float Tf = 450.0; // 右邊界固定(fix)溫度 K
    const float dt = 0.001; // 固定時間步長(Time step)，若r值檢測大於0.5則須將此值調小。有兩種方式：一是固定r值，另外就是固定dt，由於這題我們想得到的是時間，所以固定dt所得到的答案較精確。

// Allocate and free memory
    void Allocate_memory(float **array1, float **array2, float **array3, int N){    
        *array1 = (float*)malloc(N * sizeof(float));
        *array2 = (float*)malloc(N * sizeof(float));
        *array3 = (float*)malloc(N * sizeof(float));

        if (*array1 == NULL || *array2 == NULL || *array3 == NULL){
            printf("Memory allocation failed!\n");
            exit(1);
        }
        printf("Memory allocation successfully for %d elements\n", N);
    }

    void Free_memory( float *array1, float *array2, float *array3){
            free(array1);
            free(array2);
            free(array3);
            printf("Memory freed successfully\n");
    }

// Using FTCS to calculate time when T=325K at x = 5cm.
    int main(){
        const int N = 200;

        float dt = ; 
        float dx = 0.1/(N-1); // Total length = 0.1 m, and N-1 is because N start from 0.
        float *T, *T_new, *x ;
        float time = 0.0;
        const float Cv = R/(gamma-1); // specific heat at constant volume of air
        const float Cp = Cv+R; // specific heat at constant pressure of air
        const float alpha = k/(rho*Cp); // thermal diffusivity
        const float r = alpha*(dt/(dx)^2);
        int i;

    Allocate_memory(&T, &T_new, &x, N);

    

    Free_memory(T, T_new, x);

    return(0);
}