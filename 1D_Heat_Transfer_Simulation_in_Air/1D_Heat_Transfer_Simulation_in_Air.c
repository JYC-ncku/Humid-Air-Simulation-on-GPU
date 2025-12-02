#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
// Setting parameters   
    const float k = 0.026; // Thermal conductivity of air is 0.026 W/m K
    const float R = 287; // Gas constant of air is 287 J/kg K
    const float rho = 1.23; // density of air is 1.23 kg/m^3
    const float gamma = 1.4; // Gamma value of air is 1.4

// 程式控制常數
    const float L = 0.10;      // Total length 0.10 m
    const float Tt = 325.0; // 目標（target)溫度 K
    const float Ti = 300.0; // 初始(initial)溫度 K
    const float Tf = 450.0; // 右邊界固定(fix)溫度 K

// Allocate and free memory
    void Allocate_memory(float **array1, float **array2, float **array3, int N){    
        *array1 = (float*)malloc(N * sizeof(float));
        *array2 = (float*)malloc(N * sizeof(float));
        *array3 = (float*)malloc(N * sizeof(float));

        if (*array1 == NULL || *array2 == NULL || *array3 == NULL){
            printf("Memory allocation failed!\n");
            exit(1); //終止整個程式
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
    void FRCS_1D(int N, int no_step){
        // 變數宣告
        float *x, *T, *T_new;      
        // malloc x, T, T_new
        Allocate_memory(&T, &T_new, &x, N);
        // 1.參數計算
        float Cv = R/(gamma-1);
        float Cp = Cv+R;
        float alpha = k/(rho*Cp);
        // 2. 網格參數計算
        float dx = L / (N - 1); // Because N is start from zero, so need to minus one.
        float dt = 0.001; // 固定時間步長(Time step)，若r值檢測大於0.5則須將此值調小。有兩種方式：一是固定r值，另外就是固定dt，由於這題我們想得到的是時間，所以固定dt所得到的答案較精確。
        // 3. 穩定參數 R (PHI) 計算
        float PHI = alpha * dt / (dx * dx); // PHI = alpha*(dt/dx^2)
        // 4. 穩定性檢查(stability analysis) (PHI < 0.5) 
        if (PHI >= 0.5) {
            printf("ERROR: PHI must be < 0.5, your PHI is %.3f.\n", PHI);
        //執行記憶體釋放(Free memory)
            Free_memory(T, T_new, x);
            return 1; // 會退出當前函式，回到呼叫處。
        }

        // 5. 初始條件設定(Setting initial condition),t=0       
        float time_to_reach = 0.0; //initialize the time,為了計算達到目標的時間，必須先將時間初始化。
        int target_reached = 0; //target_reached是一個旗標 (Flag) 變數，它用來表示某個布林（Boolean）狀態：目標是否達成，0表示未達成(not completed)、1表示已達成(complete)。
        // 6. 時間疊代迴圈
        for (int step = 1; step <= no_steps; step++) {        
            // a. 空間計算迴圈 
            for (int cell = 1; cell < N - 1; cell++) {
                    T_new[cell] = T[cell] + PHI * (T[cell-1] + T[cell+1] - 2*T[cell]); // Tnew(cell) = T(cell) + PHI*(LEFT + RIGHT - 2.0*T(cell));
            }

            // b. 邊界條件強制實施
            T_new[0] = T[0]; // 左邊界(LEFT boundary)：絕熱 
            T_new[N-1] = Tf; // 右邊界(RIGHT boundary)：固定溫度(fix temperature)
            while(time_to_reach < tmax){
                 // d. 陣列更新
                for (int i = 0; i < N; i++) { 
                    T[i] = Tnew[i];           
                }   
                if (target_reached == 1) {
                    break; 
                }   
            }              
        }

    // 7. 最終結果與記憶體釋放
    if (target_reached == 0) {
    printf("Simulation ended: The target was not achieved within the specified time.\n");
    } 
    else {
    printf("Target achieved! Time spent: %.4f seconds.\n", time_to_reach);
    }
    Free_memory(T, T_new, x); 
}
// Setting the number of grid points and the time step..
    int main(){
        const int N = 201; 
        const int tm = 500; //tmax:模擬時間設定最多500秒就結束，如果在500秒之前得到結果就會提早結束。
        const int time_step = (int)(tm/dt) //總時間步數(total time step)，total time step = tmax / dt.

    FTCS_1D(N, time_step);

    return(0);
    } 
