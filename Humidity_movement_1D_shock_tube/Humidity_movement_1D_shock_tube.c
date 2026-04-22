#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6,
                     float **array7, float **array8, float **array9, float **array10, float **array11, float **array12,
                     float **array13, float **array14, float **array15, float **array16, float **array17, float **array18,
                     int N_CELLS){
    *array1 = (float*)malloc(N_CELLS * sizeof(float));
    *array2 = (float*)malloc(N_CELLS * sizeof(float));
    *array3 = (float*)malloc(N_CELLS * sizeof(float));
    *array4 = (float*)malloc(N_CELLS * sizeof(float));
    *array5 = (float*)malloc(N_CELLS * sizeof(float));
    *array6 = (float*)malloc(N_CELLS * sizeof(float));
    *array7 = (float*)malloc(N_CELLS * sizeof(float));
    *array8 = (float*)malloc(N_CELLS * sizeof(float));
    *array9 = (float*)malloc(N_CELLS * sizeof(float));
    *array10 = (float*)malloc(N_CELLS * sizeof(float));
    *array11 = (float*)malloc(N_CELLS * sizeof(float));
    *array12 = (float*)malloc(N_CELLS * sizeof(float));
    *array13 = (float*)malloc(N_CELLS * sizeof(float));
    *array14 = (float*)malloc(N_CELLS * sizeof(float));
    *array15 = (float*)malloc((N_CELLS+1) * sizeof(float));
    *array16 = (float*)malloc((N_CELLS+1) * sizeof(float));
    *array17 = (float*)malloc((N_CELLS+1) * sizeof(float));
    *array18 = (float*)malloc((N_CELLS+1) * sizeof(float));
    if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL ||
        *array7 == NULL || *array8 == NULL || *array9 == NULL || *array10 == NULL || *array11 == NULL || *array12 == NULL ||
        *array13 == NULL || *array14 == NULL || *array15 == NULL || *array16 == NULL || *array17 == NULL || *array18 == NULL){
        printf("Memory allocation failed!\n");
        exit(1);
    }
        printf("Memory allocation successfully for %d elements!\n", N_CELLS);
}

void Free_memory(float *array1, float *array2, float *array3, float *array4, float *array5, float *array6,
                 float *array7, float *array8, float *array9, float *array10, float *array11, float *array12,
                 float *array13, float *array14, float *array15, float *array16, float *array17, float *array18){
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
float MAX_Wave_Speed(float u_L, float u_R, float a_L, float a_R){
    float W_L = fabs(u_L) + a_L;
    float W_R = fabs(u_R) + a_R;
    float W_LOCAL_MAX;
    if (W_L > W_R){
        W_LOCAL_MAX = W_L;
    }else {
        W_LOCAL_MAX = W_R;
    }
    return W_LOCAL_MAX;
}

void Calc_Rusanov_Flux(float rho_L, float rho_R, float u_L, float u_R, float T_L, float T_R, float p_L, float p_R, float e_L, float e_R, float W_LOCAL_MAX,
                       float *mass_flux, float *momentum_flux, float *energy_flux, float *rhov_flux, float H_L, float H_R, int j){
    float mass_L = rho_L;
    float momentum_L = rho_L * u_L;
    float energy_L = e_L;
    float mass_R = rho_R;
    float momentum_R = rho_R * u_R;
    float energy_R = e_R;
    float rhov_L = H_L * rho_L; // phov is vapor density = rho * humidity
    float rhov_R = H_R * rho_R;

    float mass_flux_L = rho_L * u_L;
    float mass_flux_R = rho_R * u_R;
    float momentum_flux_L = rho_L * u_L * u_L + p_L;
    float momentum_flux_R = rho_R * u_R * u_R + p_R;
    float energy_flux_L = (e_L + p_L) * u_L;
    float energy_flux_R = (e_R + p_R) * u_R;
    float rhov_flux_L = rhov_L * u_L;
    float rhov_flux_R = rhov_R * u_R;

    mass_flux[j] = 0.5 * (mass_flux_L + mass_flux_R) - 0.5 * W_LOCAL_MAX * (mass_R - mass_L);
    momentum_flux[j] = 0.5 * (momentum_flux_L + momentum_flux_R) - 0.5 * W_LOCAL_MAX * (momentum_R - momentum_L);
    energy_flux[j] = 0.5 * (energy_flux_L + energy_flux_R) - 0.5 * W_LOCAL_MAX * (energy_R - energy_L);
    rhov_flux[j] = 0.5 * (rhov_flux_L + rhov_flux_R) - 0.5 * W_LOCAL_MAX * (rhov_R - rhov_L);
}

int main(){
    int N_CELLS = 200;
    float *x, *p0, *p1, *p2, *p3, *p4, *R, *CV, *GAMMA, *a, // p0 is density, p1 is velocity, p2 is temperature, p3 is pressure, p4 is humidity, a is sound speed
           *mass, *momentum, *energy, *rhov,
           *mass_flux, *momentum_flux, *energy_flux, *rhov_flux;
    float L = 1.0;
    float t = 0;
    float D = 1.837e-5; //Diffusivity of water vapor. unit:(m^2/s)
    float R_bar = 8.315; // unti:kJ/(mol*K)
    float MW_H2O = 18.02; // unit:kg/kmol
    float MW_air = 28.97; // unit:kg/kmol
    float R_v = R_bar / MW_H2O; // unit:J/(kg*k) R = R_bar / Molecular weight
    float R_dry = R_bar / MW_air;
    float CV_v = 4.18; // unit:kJ/(kg*K) // H2O的定容比熱
    float CV_dry = 0.718;
    float CFL = 0.5;
    float dx = L/N_CELLS;
    float W_GLOBAL_MAX;
    float t_target = L / (0.2 * sqrt((R_dry/1000) * 298.15));

    Allocate_memory(&x, &p0, &p1, &p2, &p3, &p4, &R, &CV, &GAMMA, &a, &mass, &momentum, &energy, &rhov, &mass_flux, &momentum_flux, &energy_flux, &rhov_flux, N_CELLS);
    // Set initial condition (Because R will change, so P_L and P_R can not equal 10 and 1 directly)
    for (int i = 0; i < N_CELLS; i++){
	x[i] = (i+0.5) * dx;
        if (i < N_CELLS/2){
            //p0[i] = 10.0;
            p0[i] = 12.55; //unti: kg/m^3
            p1[i] = 0.0; //unit: m/s
            //p2[i] = 1.0;
            p2[i] = 298.15; //unit: K
            p4[i] = 1.0;
        } else {
            //p0[i] = 1.0;
            p0[i] = 1.225; // density of air
            p1[i] = 0.0;
            //p2[i] = 1.0;
            p2[i] = 298.15; // room temperature
            p4[i] = 1.0;
        }
    }

    for (int i = 0; i<N_CELLS; i++){
//        R[i] = ((1 - p4[i]) * (R_dry/R_dry) + p4[i] * (R_v/R_dry)); // 所有參數(密度、速度、溫度、壓力)都是用無因次化去做計算，所以R跟CV也要無因次化，通常以dry air為基準。
//        CV[i] = (1 - p4[i]) * (CV_dry/R_dry) + p4[i] * (CV_v/R_dry);
        R[i] = (1 - p4[i]) * R_dry + p4[i] * R_v;
        CV[i] = (1 - p4[i]) * CV_dry + p4[i] * CV_v;
        GAMMA[i] = 1 + (R[i] / CV[i]);
        p3[i] = p0[i] * R[i] * p2[i]; // Pressure = rho * R * T
    }

    for (int i=0; i<N_CELLS; i++){
        mass[i] = p0[i];
        momentum[i] = p0[i] * p1[i];
        energy[i] = 0.5 * p0[i] * p1[i] * p1[i] + (p3[i] / (GAMMA[i] - 1));
        rhov[i] = p0[i] * p4[i];
    }

   while (t<t_target){

        // In order to compute dt, need to find the Max wave speed first.
        float W_GLOBAL_MAX = 1.0e-10;

        for (int j = 1; j < N_CELLS; j++){
		float rho_L = p0[j-1];
    		float rho_R = p0[j];
    		float u_L = p1[j-1];
    		float u_R = p1[j];
    		float T_L = p2[j-1];
    		float T_R = p2[j];
            	float p_L = p3[j-1];
            	float p_R = p3[j];
	        float H_L = p4[j-1];
            	float H_R = p4[j];
           	float GAMMA_L = GAMMA[j-1];
            	float GAMMA_R = GAMMA[j];
            	float e_L = 0.5 * rho_L * u_L * u_L + (p_L / (GAMMA_L - 1));
            	float e_R = 0.5 * rho_R * u_R * u_R + (p_R / (GAMMA_R - 1));
            	float R_L = R[j-1];
            	float R_R = R[j];
            	float a_L = sqrt(GAMMA_L * R_L * T_L);
            	float a_R = sqrt(GAMMA_R * R_R * T_R);
            	float W_LOCAL_MAX = MAX_Wave_Speed(u_L, u_R, a_L, a_R);
            Calc_Rusanov_Flux(rho_L, rho_R, u_L, u_R, T_L, T_R, p_L, p_R, e_L, e_R, W_LOCAL_MAX,
                              mass_flux, momentum_flux, energy_flux, rhov_flux,  H_L, H_R, j);
            if (W_LOCAL_MAX > W_GLOBAL_MAX){
                W_GLOBAL_MAX = W_LOCAL_MAX;
            }
        }

        float dt = CFL * (dx / W_GLOBAL_MAX);

        //Set boundary condition
        mass_flux[0] = mass_flux[1];
        momentum_flux[0] = momentum_flux[1];
        energy_flux[0] = energy_flux[1];
        rhov_flux[0] = rhov_flux[1];
	rhov[0] = rhov[1];
        mass_flux[N_CELLS] = mass_flux[N_CELLS - 1];
        momentum_flux[N_CELLS] = momentum_flux[N_CELLS - 1];
        energy_flux[N_CELLS] = energy_flux[N_CELLS - 1];
        rhov_flux[N_CELLS] = rhov_flux[N_CELLS - 1];
	rhov[N_CELLS] = rhov[N_CELLS - 1];


        // Use FVM to get new conservation values
        for (int i = 1; i<=N_CELLS; i++){
            mass[i] = mass[i] - (dt / dx) * (mass_flux[i+1] - mass_flux[i]);
            momentum[i] = momentum[i] - (dt / dx) * (momentum_flux[i+1] - momentum_flux[i]);
            energy[i] = energy[i] - (dt / dx) * (energy_flux[i+1] - energy_flux[i]);
            rhov[i] = rhov[i] - (dt/dx) * (rhov_flux[i+1] - rhov_flux[i]) + (dt/(dx*dx)) * D * (rhov[i+1] - 2 * rhov[i] + rhov[i-1]);
        }

        for (int i = 1; i<=N_CELLS; i++){
            p0[i] = mass[i];
            p1[i] = momentum[i] / mass[i];
            p4[i] = rhov[i]/p0[i];
            p3[i] = (GAMMA[i] - 1) * (energy[i] - 0.5 * p0[i] * p1[i] * p1[i]);
            p2[i] = p3[i] / (p0[i] * R[i]);
        }

        t += dt;
    }

    FILE * pFile = fopen("Results_of_200_cells.txt","w");
    for (int i = 0; i<N_CELLS; i++){
        fprintf(pFile, "%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", x[i], p0[i], p1[i], p2[i], p3[i], p4[i]);
    }
    fclose(pFile);

    Free_memory(x, p0, p1, p2, p3, p4, rhov, R, CV, GAMMA, a, mass, momentum, energy, mass_flux, momentum_flux, energy_flux, rhov_flux);
    return 0;
}
