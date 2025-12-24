#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Setting parameters
const double L = 0.1;
const double hum_init = 0;
const double hum_final = 1;
const double hum_target = 0.25;
const double Pi = M_PI;
const double D = 2.0e-5; //Diffusivity of air.
//Set a maximum number of steps
const int no_steps = 15000;

void Allocate_memory(double **array1, double **array2, double **array3, double **array4, double **array5, int N){
    *array1 = (double*)malloc(N * sizeof(double));
    *array2 = (double*)malloc(N * sizeof(double));
    *array3 = (double*)malloc(N * sizeof(double));
    *array4 = (double*)malloc((N+1) * sizeof(double));
    *array5 = (double*)malloc(N * sizeof(double));
    if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL){
        printf("Memory allocation failed!\n");
        exit(1);
    }
        printf("Memory allcoation successfully for %d elements\n", N);
}

void Free_memory(double *array1, double *array2, double *array3, double *array4, double *array5){
    free(array1);
    free(array2);
    free(array3);
    free(array4); 
    free(array5);   
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
    FILE *pFile;
    double *x, *hum, *hum_new, *F, *hum_exact;
    double dx = L / N;
    double PHI = 0.25;
    double dt = PHI * ((dx*dx)/D);
    double time = 0; //Iintialize the time.
    int reached = 0;
    double reach_time = 0;
    
    Allocate_memory(&x, &hum, &hum_new, &F, &hum_exact, N);
    pFile = fopen("results.txt","w");
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

        if (reached == 0 && hum_new[N/2] >= hum_target) {
                reach_time = time; // 記下現在的時間
                reached = 1;       // 改成 1，之後的 step 就不會再進來這個 if
                printf("\n At %.3f seconds, the humidity at x = 5 cm reaches 25%%.\n", reach_time);
            }

        for (int i = 0; i < N; i++){
            hum[i]=hum_new[i];
        }
        time = time+dt;
    }
    //Calculate exact solution
    for (int i = 0; i < N; i++) {
        double current_x = x[i];
        double sum = 0.0;
        for (int n = 1; n < 200; n++) { // n可以自己定
            double An = (2.0 * pow(-1.0, n)) / (n * M_PI);
            double term = An * sin((n * M_PI * current_x) / L) * exp(-(n * n * M_PI * M_PI * D * time) / (L * L));
            sum += term;
        }
        hum_exact[i] = (current_x / L) + sum;
    }



    for (int i =0; i<N; i++){
        double error = (hum[i]-hum_exact[i]) * (hum[i]-hum_exact[i]);
        fprintf(pFile,"%g\t%g\t%g\t%e\n", x[i], hum[i], hum_exact[i], error);
    }
    fclose(pFile);
    printf("\n Total time spent: %.3f seconds. \n", time);

    Free_memory(x, hum, hum_new, F, hum_exact);
}

int main(){
    const int N=200; //Number of cell
    Humidity_1D(N, no_steps);
    return(0);
}
