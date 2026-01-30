#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(double **array1, double **array2, double **array3, double **array4, double **array5, int N_CELLS){
    *array1 = (double*)malloc(N_CELLS * sizeof(double));
    *array2 = (double*)malloc(N_CELLS * sizeof(double));
    *array3 = (double*)malloc(N_CELLS * sizeof(double));
    *array4 = (double*)malloc(N_CELLS * sizeof(double));
    *array5 = (double*)malloc((N_CELLS+1) * sizeof(double));
    if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL ){
        printf("Memory allocation failed!\n");
        exit(1);
    }
        printf("Memory allocation successfully for %d elements!\n", N_CELLS);
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
    |   0    |  ...   |    N   |   ===>   Have total N cells and N+1 interface. (because N is from 0)
    |        |        |        | 
    ----------------------------
    0        1       N        N+1
*/

// In order to compute dt, need to find the Max wave speed first.
double MAX_Wave_Speed(double u_L, double u_R, double a_L, double a_R){
    double W_L = fabs(u_L) + a_L;
    double W_R = fabs(u_R) + a_R;
    if (W_L > W_R){
        return W_L;
    }else {
        return W_R;
    }
}


int main(){
    int N_CELLS = 200;
    int N_INTERFACES = N_CELLS+1;
    double *x, *p0, *p1, *p2, *flux; // p0 is density, p1 is velocity, p2 is temperature
    float L = 1.0;
    float t = 0;
    float t_FINAL = 0.2;
    double R = 1.0;
    double GAMMA = 1.4;
    double CFL = 0.5;
    double dx = L/N_CELLS;

    Allocate_memory(&x, &p0, &p1, &p2, &flux, N_CELLS);
    // Set initial condition
    for (int i = 0; i < N_CELLS; i++){
        x[i] = (i+0.5) * dx;
        if (i < N_CELLS/2){
            p0[i] = 10;
            p1[i] = 0;
            p2[i] = 1;
        } else {
            p0[i] = 1;
            p1[i] = 0;
            p2[i] = 1;
        }
    }
   while (t<t_FINAL){
        // In order to compute dt, need to find the Max wave speed first.
        double W_MAX = 1e-10;
        for (int i = 0; i < N_CELLS-1; i++){
			double rho_L = p0[i];
    		double rho_R = p0[i+1];
    		double u_L  = p1[i];
    		double u_R  = p1[i+1];
    		double T_L  = p2[i];
    		double T_R  = p2[i+1];      
            double a_L = sqrt(GAMMA * R * T_L); // Sound speed a = (R*T)^0.5
            double a_R = sqrt(GAMMA * R * T_R);
            double W_LOCAL = MAX_Wave_Speed(u_L, u_R, a_L, a_R);

            if (W_LOCAL>W_MAX){
                W_MAX = W_LOCAL;
            }
        }
        double dt = CFL * (dx/W_MAX);
        t += dt;
    }

    Free_memory(x, p0, p1, p2, flux);
    return 0;
}



/*
    while (t<t_FINAL){
        // In order to compute dt, need to find the Max wave speed first.
        double W_MAX, W_LOCAL, W_L, W_R;
        for (int i = 0; i < N_CELLS; i++){
			double QL_rho = p0[i];
    		double QR_rho = p0[i+1];
    		double QL_u  = p1[i];
    		double QR_u  = p1[i+1];
    		double QL_T   = p2[i];
    		double QR_T   = p2[i+1];      
            double a_L = sqrt(R * QL_T); // Sound speed a = (R*T)^0.5
            double a_R = sqrt(R * QR_T);
        }
            // Calculate MAX Wave Speed 
        for (int j = 1; j<N_CELLS; j++){
            W_L = fabs(p1[j-1]) + a_L;
            W_R = fabs(p1[j]) - a_R;
            }
        if (W_L > W_R){
            W_LOCAL = W_L;
        }else {
            W_LOCAL = W_R;
        }
        if (W_LOCAL>W_MAX){
            W_MAX = W_LOCAL;
        }

        double dt = CFL * (dx/W_MAX);
        t += dt;
    }
}
*/