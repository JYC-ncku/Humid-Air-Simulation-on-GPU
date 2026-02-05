#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(double **array1, double **array2, double **array3, double **array4, double **array5, double **array6,
                     double **array7, double **array8, double **array9, double **array10, double **array11, double **array12, 
                     double **array13, double **array14, double **array15, double **array16, double **array17, double **array18, 
                     int N_CELLS){
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
    *array14 = (double*)malloc(N_CELLS * sizeof(double));
    *array15 = (double*)malloc((N_CELLS+1) * sizeof(double));
    *array16 = (double*)malloc((N_CELLS+1) * sizeof(double));
    *array17 = (double*)malloc((N_CELLS+1) * sizeof(double));
    *array18 = (double*)malloc((N_CELLS+1) * sizeof(double)); 
    if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL || 
        *array7 == NULL || *array8 == NULL || *array9 == NULL || *array10 == NULL || *array11 == NULL || *array12 == NULL ||
        *array13 == NULL || *array14 == NULL || *array15 == NULL || *array16 == NULL || *array17 == NULL || *array18 == NULL){
        printf("Memory allocation failed!\n");
        exit(1);
    }
        printf("Memory allocation successfully for %d elements!\n", N_CELLS);
}

void Free_memory(double *array1, double *array2, double *array3, double *array4, double *array5, double *array6,
                 double *array7, double *array8, double *array9, double *array10, double *array11, double *array12,
                 double *array13, double *array14, double *array15, double *array16, double *array17, double *array18){
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
    free(array17);
    free(array18);
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
                       double *mass_flux, double *momentum_flux, double *energy_flux, double *rhov_flux, double H_L, double H_R, int j){
    double mass_L = rho_L;
    double momentum_L = rho_L * u_L;
    double energy_L = e_L;
    double mass_R = rho_R;
    double momentum_R = rho_R * u_R;
    double energy_R = e_R;
    double rhov_L = H_L * rho_L; // phov is vapor density = rho * humidity
    double rhov_R = H_R * rho_R;

    double mass_flux_L = rho_L * u_L;
    double mass_flux_R = rho_R * u_R;
    double momentum_flux_L = rho_L * u_L * u_L + p_L;
    double momentum_flux_R = rho_R * u_R * u_R + p_R;
    double energy_flux_L = (e_L + p_L) * u_L;
    double energy_flux_R = (e_R + p_R) * u_R;
    double rhov_flux_L = rhov_L * u_L;
    double rhov_flux_R = rhov_R * u_R;

    mass_flux[j] = 0.5 * (mass_flux_L + mass_flux_R) - 0.5 * W_LOCAL_MAX * (mass_R - mass_L);
    momentum_flux[j] = 0.5 * (momentum_flux_L + momentum_flux_R) - 0.5 * W_LOCAL_MAX * (momentum_R - momentum_L);
    energy_flux[j] = 0.5 * (energy_flux_L + energy_flux_R) - 0.5 * W_LOCAL_MAX * (energy_R - energy_L);
    rhov_flux[j] = 0.5 * (rhov_flux_L + rhov_flux_R) - 0.5 * W_LOCAL_MAX * (rhov_R - rhov_L);
}
    

int main(){
    int N_CELLS = 800;
    double *x, *p0, *p1, *p2, *p3, *p4, *R, *CV, *GAMMA, *a, // p0 is density, p1 is velocity, p2 is temperature, p3 is pressure, p4 is humidity, a is sound speed
           *mass, *momentum, *energy, *rhov,
           *mass_flux, *momentum_flux, *energy_flux, *rhov_flux;
    float L = 1.0;
    float t = 0;
    float t_FINAL = 0.2;
    double R_bar = 8.315; // unti:J/(mol*K)
    double MW_H2O = 0.01802; // unit:kg/mol
    double MW_air = 0.02897; // unit:kg/mol
    double R_v = R_bar / MW_H2O; // unit:J/(kg*k)
    double R_dry = R_bar / MW_air;
    double CV_v = 1410; // unit:J/(kg*K)
    double CV_dry = 717.5;
    double CFL = 0.5;
    double dx = L/N_CELLS;
    double W_GLOBAL_MAX;

    Allocate_memory(&x, &p0, &p1, &p2, &p3, &p4, &R, &CV, &GAMMA, &a, &mass, &momentum, &energy, &rhov, &mass_flux, &momentum_flux, &energy_flux, &rhov_flux, N_CELLS);
    // Set initial condition
    for (int i = 0; i < N_CELLS; i++){
        x[i] = (i+0.5) * dx;
        if (i < N_CELLS/2){
            p0[i] = 10;
            p1[i] = 0;
            p2[i] = 1;
            p4[i] = 1;
        } else {
            p0[i] = 1;
            p1[i] = 0;
            p2[i] = 1;
            p4[i] = 0;
        }
    }

    for (int i = 0; i<N_CELLS; i++){
        R[i] = ((1 - p4[i]) * (R_dry / R_dry) + p4[i] * (R_v / R_dry)); // 所有參數都(密度、速度、溫度、壓力)都是用無因次化去做計算，所以R跟CV也要無因次化，通常以dry air為基準。
        CV[i] = (1 - p4[i]) * (CV_dry / CV_dry) + p4[i] * (CV_v / CV_dry);
        GAMMA[i] = 1 + (R[i] / CV[i]);
        p3[i] = p0[i] * R[i] * p2[i];
    }    

    for (int i=0; i<N_CELLS; i++){
        mass[i] = p0[i];
        momentum[i] = p0[i] * p1[i];
        energy[i] = 0.5 * p0[i] * p1[i] * p1[i] + (p3[i] / (GAMMA[i] - 1));
        rhov[i] = p0[i] * p4[i];
    }

   while (t<t_FINAL){

        // In order to compute dt, need to find the Max wave speed first.
        double W_GLOBAL_MAX = 1.0e-10;
        
        for (int j = 1; j < N_CELLS; j++){
			double rho_L = p0[j-1];
    		double rho_R = p0[j];
    		double u_L = p1[j-1];
    		double u_R = p1[j];
    		double T_L = p2[j-1];
    		double T_R = p2[j];
            double p_L = p3[j-1];
            double p_R = p3[j];
            double H_L = p4[j-1];
            double H_R = p4[j];
            double GAMMA_L = GAMMA[j-1];
            double GAMMA_R = GAMMA[j];
            double e_L = 0.5 * rho_L * u_L * u_L + (p_L / (GAMMA_L - 1));
            double e_R = 0.5 * rho_R * u_R * u_R + (p_R / (GAMMA_R - 1));
            double R_L = R[j-1];
            double R_R = R[j];
            double a_L = sqrt(GAMMA_L * R_L * T_L);
            double a_R = sqrt(GAMMA_R * R_R * T_R);
            double W_LOCAL_MAX = MAX_Wave_Speed(u_L, u_R, a_L, a_R);
            Calc_Rusanov_Flux(rho_L, rho_R, u_L, u_R, T_L, T_R, p_L, p_R, e_L, e_R, W_LOCAL_MAX,
                              mass_flux, momentum_flux, energy_flux, rhov_flux,  H_L, H_R, j);
            if (W_LOCAL_MAX > W_GLOBAL_MAX){
                W_GLOBAL_MAX = W_LOCAL_MAX;
            }
        }

        double dt = CFL * (dx / W_GLOBAL_MAX);

        //Set boundary condition
        mass_flux[0] = mass_flux[1];
        momentum_flux[0] = momentum_flux[1];
        energy_flux[0] = energy_flux[1];
        rhov_flux[0] = rhov_flux[1];
        mass_flux[N_CELLS] = mass_flux[N_CELLS - 1];
        momentum_flux[N_CELLS] = momentum_flux[N_CELLS - 1];
        energy_flux[N_CELLS] = energy_flux[N_CELLS - 1];
        rhov_flux[N_CELLS] = rhov_flux[N_CELLS - 1];
    
        // Use FVM to get new conservation values
        for (int i=0; i<N_CELLS; i++){
            mass[i] = mass[i] - (dt / dx) * (mass_flux[i+1] - mass_flux[i]);
            momentum[i] = momentum[i] - (dt / dx) * (momentum_flux[i+1] - momentum_flux[i]);
            energy[i] = energy[i] - (dt / dx) * (energy_flux[i+1] - energy_flux[i]);
            rhov[i] = rhov[i] - (dt/dx) * (rhov_flux[i+1] - rhov_flux[i]);
        }

        for (int i = 0; i<N_CELLS; i++){
            p0[i] = mass[i];
            p1[i] = momentum[i] / mass[i];
            p4[i] = rhov[i]/p0[i];
            p3[i] = (GAMMA[i] - 1) * (energy[i] - 0.5 * p0[i] * p1[i] * p1[i]);
            p2[i] = p3[i] / (p0[i] * R[i]);
        }

        t += dt;        
    }

    FILE * pFile = fopen("Results_of_800_cells.txt","w");
    for (int i = 0; i<N_CELLS; i++){
        fprintf(pFile, "%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", x[i], p0[i], p1[i], p2[i], p3[i], p4[i]);
    }
    fclose(pFile);

    Free_memory(x, p0, p1, p2, p3, p4, rhov, R, CV, GAMMA, a, mass, momentum, energy, mass_flux, momentum_flux, energy_flux, rhov_flux);
    return 0;
}
