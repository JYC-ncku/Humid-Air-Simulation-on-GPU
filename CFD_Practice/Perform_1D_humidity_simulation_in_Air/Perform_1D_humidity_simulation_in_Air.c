#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Setting parameters
const double L = 0.10;
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
double Humidity_1D(int N, int no_steps, double *x,double *hum, double*hum_new, double*F){  //這裡不用void的原因是因為void不會回傳數據
    double dx = L / N;
    double PHI = 0.25;
    double dt = PHI * ((dx*dx)/D);
    double time = 0; //Iintialize the time. 
    int reached = 0;
    double reach_time = 0;
    
        //Initialization
        for (int i=0; i<N; i++){
            x[i]=(i+0.5) * dx; //計算每個cell中間的flux
            hum[i]=hum_init; //每個cell一開始的濕度都為0
        }

        //Calculate interface flux
        //中間邊界flux    FILE *pFile;
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
                reached = 1;       // 改成1，之後的 step 就不會再進來這個 if
                printf("\n At %.3f seconds, the humidity at x = 5 cm reaches 25%%.\n", reach_time);
            }

        for (int i = 0; i < N; i++){
            hum[i]=hum_new[i];
        }
        time = time+dt;
    }
   
    //write fluxes to file: one line per interface (index, flux) 
    FILE *pFile = fopen("fluxes.txt","w");
    for (int j=0; j<N+1; j++){
        fprintf(pFile, "%d %.15e\n", j, F[j]);
    }
    fclose(pFile);
    printf("\n Total time spent: %.3f seconds. \n", time);
    return time;
}

void Calculate_Humidity_Exact(int N, int N_TERMS, double *x_arr, double sim_time, double *hum_exact){  //x_arr跟x[i]是一樣的東西，只是在這裡需要更改寫法，不然電腦會看不懂x[i]是什麼
    //Calculate exact solution
    double C0 = 0.0;
    double Csat = 1.0;  
    double Pi = M_PI;

    // 2. 使用 for 迴圈計算級數累加 (n 從 0 到 N_TERMS)
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int n = 1; n < N_TERMS; n++) {
            // 每次迴圈都要重新計算當前的 An
            double An = ((2.0 * pow(-1.0, n)) / (n*Pi));

            // 計算級數的當前項 (term)
            double term = An * sin((n * Pi * x_arr[i]) / L) * exp(-( n * n * Pi * Pi * D * sim_time) / (L * L));
            
            // 累加到 sum
            sum += term;
        }
        double C_xt = (x_arr[i]/L)+sum;
        hum_exact[i] = C_xt;
    }
}

int main(){
    const int N=200; //Number of cell
    const int N_TERMS=50; //Number of exact_solution steps
    double *x, *hum, *hum_new, *F, *hum_exact;

    Allocate_memory(&x, &hum, &hum_new, &F, &hum_exact, N);
    //定義模擬時間= double Humidity_1D 回傳的時間
    double sim_time = Humidity_1D(N, no_steps, x, hum, hum_new, F); //這麼做可以確保Numerical solution 與 Exact solution跑得時間一樣，就不用另外定義一個精確解要跑的時間

    Calculate_Humidity_Exact(N, N_TERMS, x, sim_time, hum_exact);
    //寫入數據
    FILE *pFile = fopen("results.txt", "w");
    for (int i = 0; i < N; i++) {
        fprintf(pFile, "%.6f\t%.15e\t%.15e\n", x[i], hum[i], hum_exact[i]);
    }
    fclose(pFile);
    //計算誤差值
    double error = 0.0;
    for (int i =0; i<N; i++){
        error = error + (hum[i]-hum_exact[i]) * (hum[i]-hum_exact[i]);
    }
    error = (double)(error / N);
    printf("error = %e\n", error);

    Free_memory(x, hum, hum_new, F, hum_exact);

    return(0);
}
