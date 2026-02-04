#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(double **array1, double **array2, double **array3, double **array4, double **array5, double **array6,
                     double **array7, double **array8, double **array9, double **array10, double **array11, double **array12, 
                     double **array13, double **array14, double **array15, double **array16, int N_CELLS){
    *array1 = (double*)malloc(N_CELLS * sizeof(double));
    *array2 = (double*)malloc(N_CELLS * sizeof(double));
    *array3 = (double*)malloc(N_CELLS * sizeof(double));
    *array4 = (double*)malloc(N_CELLS * sizeof(double));
    *array5 = (double*)malloc(N_CELLS * sizeof(double));
    *array6 = (double*)malloc(N_CELLS * sizeof(double));
    *array7 = (double*)malloc(N_CELLS * sizeof(double));
    *array8 = (double*)malloc(N_CELLS * sizeof(double));
    *array9 = (double*)malloc(N_CELLS * sizeof(double));
    *array10 = (double*)malloc(N_CELLS * sizeof(double));
    *array11 = (double*)malloc(N_CELLS * sizeof(double));
    *array12 = (double*)malloc(N_CELLS * sizeof(double));
    *array13 = (double*)malloc(N_CELLS * sizeof(double));
    *array14 = (double*)malloc((N_CELLS+1) * sizeof(double));
    *array15 = (double*)malloc((N_CELLS+1) * sizeof(double));
    *array16 = (double*)malloc((N_CELLS+1) * sizeof(double));    
    if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL || 
        *array7 == NULL || *array8 == NULL || *array9 == NULL || *array10 == NULL || *array11 == NULL || *array12 == NULL ||
        *array13 == NULL || *array14 == NULL || *array15 == NULL || *array16 == NULL ){
        printf("Memory allocation failed!\n");
        exit(1);
    }
        printf("Memory allocation successfully for %d elements!\n", N_CELLS);
}

void Free_memory(double *array1, double *array2, double *array3, double *array4, double *array5, double *array6,
                 double *array7, double *array8, double *array9, double *array10, double *array11, double *array12,
                 double *array13, double *array14, double *array15, double *array16){
    free(array1);
    free(array2);
    free(array3);
    free(array4);
    free(array5);
    free(array6);
    free(array7);
    free(array8);
    free(array9);
    free(array10);
    free(array11);
    free(array13);
    free(array14);
    free(array15);
    free(array16);
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
    double W_LOCAL_MAX;
    if (W_L > W_R){
        W_LOCAL_MAX = W_L;
    }else {
        W_LOCAL_MAX = W_R;
    }
    return W_LOCAL_MAX;
}

void Calc_Rusanov_Flux(double rho_L, double rho_R, double u_L, double u_R, double T_L, double T_R, double p_L, double p_R, double e_L, double e_R, double W_LOCAL_MAX,
                       double *mass_flux, double *momentum_flux, double *energy_flux, int N_CELLS, int j){
    double mass_L = rho_L;
    double momentum_L = rho_L * u_L;
    double energy_L = e_L;
    double mass_R = rho_R;
    double momentum_R = rho_R * u_R;
    double energy_R = e_R;

    double mass_flux_L = rho_L * u_L;
    double mass_flux_R = rho_R * u_R;
    double momentum_flux_L = rho_L * u_L * u_L + p_L;
    double momentum_flux_R = rho_R * u_R * u_R + p_R;
    double energy_flux_L = (e_L + p_L) * u_L;
    double energy_flux_R = (e_R + p_R) * u_R;

    mass_flux[j] = 0.5 * (mass_flux_L + mass_flux_R) - 0.5 * W_LOCAL_MAX * (mass_R - mass_L);
    momentum_flux[j] = 0.5 * (momentum_flux_L + momentum_flux_R) - 0.5 * W_LOCAL_MAX * (momentum_R - momentum_L);
    energy_flux[j] = 0.5 * (energy_flux_L + energy_flux_R) - 0.5 * W_LOCAL_MAX * (energy_R - energy_L);
}
    

int main(){
    int N_CELLS = 200;
    double *x, *p0, *p1, *p2, *p3, *p4, *R, *CV, *GAMMA, *a, // p0 is density, p1 is velocity, p2 is temperature, p3 is pressure, p4 is humidity, a is sound speed
           *mass, *momentum, *energy, 
           *mass_flux, *momentum_flux, *energy_flux;
    float L = 1.0;
    float t = 0;
    float t_FINAL = 0.2;
    double CFL = 0.5;
    double dx = L/N_CELLS;
    double W_GLOBAL_MAX;

    Allocate_memory(&x, &p0, &p1, &p2, &p3, &p4, &R, &CV, &GAMMA, &a, &mass, &momentum, &energy, &mass_flux, &momentum_flux, &energy_flux, N_CELLS);
    // Set initial condition
    for (int i = 0; i < N_CELLS; i++){
        x[i] = (i+0.5) * dx;
        if (i < N_CELLS/2){
            p0[i] = 10;
            p1[i] = 0;
            p2[i] = 1;
            p3[i] = 10;
        } else {
            p0[i] = 1;
            p1[i] = 0;
            p2[i] = 1;
            p3[i] = 1;
        }
    }

    for (int i=0; i<N_CELLS; i++){
    mass[i] = p0[i];
    momentum[i] = p0[i] * p1[i];
    energy[i] = 0.5 * p0[i] * p1[i] * p1[i] + (p3[i] / (GAMMA - 1));
    }

   while (t<t_FINAL){
        // In order to compute dt, need to find the Max wave speed first.
        double W_GLOBAL_MAX = 1e-10;
        for (int j = 1; j < N_CELLS; j++){
			double rho_L = p0[j-1];
    		double rho_R = p0[j];
    		double u_L = p1[j-1];
    		double u_R = p1[j];
    		double T_L = p2[j-1];
    		double T_R = p2[j];
            double p_L = p3[j-1];
            double p_R = p3[j];
            double e_L = 0.5 * rho_L * u_L * u_L + (p_L / (GAMMA - 1));
            double e_R = 0.5 * rho_R * u_R * u_R + (p_R / (GAMMA - 1));    
            double a_L = sqrt(GAMMA * R * T_L); // Sound speed a = (R*T)^0.5
            double a_R = sqrt(GAMMA * R * T_R);
            double W_LOCAL_MAX = MAX_Wave_Speed(u_L, u_R, a_L, a_R);
            Calc_Rusanov_Flux(rho_L, rho_R, u_L, u_R, T_L, T_R, p_L, p_R, e_L, e_R, W_LOCAL_MAX,
                              mass_flux, momentum_flux, energy_flux, N_CELLS, j);
            if (W_LOCAL_MAX > W_GLOBAL_MAX){
                W_GLOBAL_MAX = W_LOCAL_MAX;
            }
        }

        double dt = CFL * (dx / W_GLOBAL_MAX);

        //Set boundary condition
        mass_flux[0] = mass_flux[1];
        momentum_flux[0] = momentum_flux[1];
        energy_flux[0] = energy_flux[1];
        mass_flux[N_CELLS] = mass_flux[N_CELLS - 1];
        momentum_flux[N_CELLS] = momentum_flux[N_CELLS - 1];
        energy_flux[N_CELLS] = energy_flux[N_CELLS - 1];
    
        // Use FVM to get new conservation values
        for (int i=0; i<N_CELLS; i++){
            mass[i] = mass[i] - (dt / dx)*(mass_flux[i+1] - mass_flux[i]);
            momentum[i] = momentum[i] - (dt / dx)*(momentum_flux[i+1] - momentum_flux[i]);
            energy[i] = energy[i] - (dt / dx)*(energy_flux[i+1] - energy_flux[i]);
        }

        for (int i = 0; i<N_CELLS; i++){
            p0[i] = mass[i];
            p1[i] = momentum[i] / mass[i];
            p3[i] = (GAMMA - 1) * (energy[i] - 0.5 * p0[i] * p1[i] * p1[i]);
            p2[i] = p3[i] / (p0[i] * R);
        }

        t += dt;        
    }

    FILE * pFile = fopen("Results_of_200_cells.txt","w");
    for (int i = 0; i<N_CELLS; i++){
        fprintf(pFile, "%.3f\t%.6f\t%.6f\t%.6f\t%.6f\n", x[i], p0[i], p1[i], p2[i], p3[i]);
    }
    fclose(pFile);

    Free_memory(x, p0, p1, p2, p3, p4, R, CV, GAMMA, a, mass, momentum, energy, mass_flux, momentum_flux, energy_flux);
    return 0;
}
