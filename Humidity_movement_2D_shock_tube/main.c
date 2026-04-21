#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6,
                     float **array7, float **array8, float **array9, float **array10, float **array11, float **array12,
                     float **array13, float **array14, float **array15, float **array16, float **array17, float **array18,
                     int N){
	*array1 = (float*)malloc(N * sizeof(float));
	*array2 = (float*)malloc(N * sizeof(float));
	*array3 = (float*)malloc(N * sizeof(float));
	*array4 = (float*)malloc(N * sizeof(float));
	*array5 = (float*)malloc(N * sizeof(float));
	*array6 = (float*)malloc(N * sizeof(float));
	*array7 = (float*)malloc(N * sizeof(float));
	*array8 = (float*)malloc(N * sizeof(float));
	*array9 = (float*)malloc(N * sizeof(float));
	*array10 = (float*)malloc(N * sizeof(float));
	*array11 = (float*)malloc(N * sizeof(float));
	*array12 = (float*)malloc(N * sizeof(float));
	*array13 = (float*)malloc(N * sizeof(float));
	*array14 = (float*)malloc(N * sizeof(float));
	*array15 = (float*)malloc((N+1) * sizeof(float));
	*array16 = (float*)malloc((N+1) * sizeof(float));
	*array17 = (float*)malloc((N+1) * sizeof(float));
	*array18 = (float*)malloc((N+1) * sizeof(float));
	if (*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL ||
	    *array7 == NULL || *array8 == NULL || *array9 == NULL || *array10 == NULL || *array11 == NULL || *array12 == NULL ||
	    *array13 == NULL || *array14 == NULL || *array15 == NULL || *array16 == NULL || *array17 == NULL || *array18 == NULL){
	printf("Memory allocation failed!\n");
	exit(1);
	}
	printf("Memory allocation successfully for %d elements!\n", N);
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
float MAX_Wave_Speed(float u_L, float u_R, float u_T, float u_B, float a_L, float a_R, float a_T, float a_B){
	float W_L = fabs(u_L) + a_L;
	float W_R = fabs(u_R) + a_R;
	float W_T = fabs(u_T) + a_T;
	float W_B = fabs(u_B) + a_B;
	float W_LOCAL_MAX;
	if (W_L > W_R && W_L > W_T && W_L > W_B){
		W_LOCAL_MAX = W_L;
	}else if (W_R > W_L && W_R > W_T && W_R > W_B){
		W_LOCAL_MAX = W_R;
	}else if (W_T > W_L && W_T > W_R && W_T > W_B){
		W_LOCAL_MAX = W_P;
	}else if (W_B > W_L && W_B > W_R && W_B > W_T
		W_LOCAL_MAX = W_B;
	}
return W_LOCAL_MAX;
}

void Calc_Rusanov_Flux(float rho_L, float rho_R, float rho_T, float rho_B, float u_L, float u_R, float u_T, float u_B, float T_L, float T_R, float T_T, float T_B,
		       float p_L, float p_R, float p_T, float p_B, float e_L, float e_R, float e_T, float e_B, float W_LOCAL_MAX,
                       float *mass_flux, float *momentum_flux, float *energy_flux, float *rhov_flux, float H_L, float H_R, float H_T, float H_B){
	float mass_L, mass_R, mass_T, mass_B, momentum_L, momentum_R, momentum_T, momentum_B, energy_L, energy_R, energy_T, energy_B, rhov_L, rhov_R, rhov_T, rhov_B,
	      mass_flux_L, mass_flux_R, mass_flux_T, mass_flux_B, momentum_flux_L, momentum__flux_R, momentum_flux_T, momentum_flux_B,
	      energy_flux_L, energy_flux_R, energy_flux_T, energy_flux_B, rhov_L, rhov_R, rhov_T, rhov_B;
	//LEFT
	mass_L = rho_L;
	momentum_L = rho_L * u_L;
	energy_L = e_L;
	rhov_L = H_L * rho_L; // rhov is vapor density = rho * humidity
	//RIGHT
	mass_R = rho_R;
	momentum_R = rho_R * u_R;
	energy_R = e_R;
	rhov_R = H_R * rho_R;
	//TOP
	mass_T = rho_T;
	momentum_T = rho_T * u_T;
	energy_T = e_T;
	rhov_T = H_T * rho_T;
	//BOTTOM
	mass_B = rho_B;
	momentum_B = rho_B * u_B;
	energy_B = e_B;
	rhov_B = H_B * rho_B;

	//Calculate flux
	//LEFT
	mass_flux_L = rho_L * u_L;
	momentum_flux_L = rho_L * u_L * u_L + p_L;
	energy_flux_L = (e_L + p_L) * u_L;
	rhov_flux_L = rhov_L * u_L;
	//RIGHT
	mass_flux_R = rho_R * u_R;
	momentum_flux_R = rho_R * u_R * u_R + p_R;
	energy_flux_R = (e_R + p_R) * u_R;
    	rhov_flux_R = rhov_R * u_R;
	//TOP
	mass_flux_T = rho_T * u_T;
	momentum_flux_T = rho_T * u_T * u_T + p_T;
	energy_flux_T = (e_T + p_T) * u_T;
    	rhov_flux_T = rhov_T * u_T;
    	//BOTTOM
    	mass_flux_B = rho_B * u_B;
	momentum_flux_B = rho_B * u_B * u_B + p_B;
	energy_flux_B = (e_B + p_B) * u_B;
    	rhov_flux_B = rhov_B * u_B;

	//待修改
	mass_flux[j] = 0.5 * (mass_flux_L + mass_flux_R) - 0.5 * W_LOCAL_MAX * (mass_R - mass_L);
	momentum_flux[j] = 0.5 * (momentum_flux_L + momentum_flux_R) - 0.5 * W_LOCAL_MAX * (momentum_R - momentum_L);
	energy_flux[j] = 0.5 * (energy_flux_L + energy_flux_R) - 0.5 * W_LOCAL_MAX * (energy_R - energy_L);
	rhov_flux[j] = 0.5 * (rhov_flux_L + rhov_flux_R) - 0.5 * W_LOCAL_MAX * (rhov_R - rhov_L);
}

int main(){
	int L = 1;
	float H = 0.5;
	int NX = 300;
	int NY = 100;
	int N = NX * NY;
	float DX = (float)L/NX;
	float DY = (float)H/NY;

    	float *x, *p0, *p1, *p2, *p3, *p4, *R, *CV, *GAMMA, *a, // p0 is density, p1 is velocity, p2 is temperature, p3 is pressure, p4 is humidity, a is sound speed
	      *mass, *momentum, *energy, *rhov,
	      *mass_flux, *momentum_flux, *energy_flux, *rhov_flux;
	float t = 0;
	float D = 1.837e-5; //Diffusivity of water vapor, unit:(m^2/s)
	float t_FINAL = 0.2;
	float R_bar = 8.315; // Ideal gas constant, unti:kJ/(kmol*K)
	float MW_H2O = 18.02; // Molecular weigh, unit:kg/kmol
	float MW_air = 28.97; // Molecular weigh, unit:kg/kmol
	float R_v = R_bar / MW_H2O; // Specific gas constant of H2O, unit:kJ/(kg*k); R = R_bar / Molecular weight
	float R_dry = R_bar / MW_air;// Specific gas constand of air.
	float CV_v = 4.184; // unit:kJ/(kg*K) // H2O的定容比熱
	float CV_dry = 0.718;
	float CFL = 0.5;
	float dx = L/N_CELLS;
	float W_GLOBAL_MAX;

	Allocate_memory(&x, &p0, &p1, &p2, &p3, &p4, &R, &CV, &GAMMA, &a, &mass, &momentum, &energy, &rhov, &mass_flux, &momentum_flux, &energy_flux, &rhov_flux, N_CELLS);
	// Set initial condition (Because R will change, so P_L and P_R can not equal 10 and 1 directly)
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

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			int INDEX = (i * NY) + j;
			R[INDEX] = ((1 - p4[INDEX]) * (R_dry / R_dry) + p4[INDEX] * (R_v / R_dry)); // 所有參數(密度、速度、溫度、壓力)都是用無因次化去做計算，所以R跟CV也要無因次化，通常以dry air為基準。
			CV[INDEX] = (1 - p4[INDEX]) * (CV_dry / CV_dry) + p4[INDEX] * (CV_v / CV_dry);
			GAMMA[INDEX] = 1 + (R[INDEX] / CV[INDEX]);
			p3[INDEX] = p0[INDEX] * R[INDEX] * p2[INDEX]; // Pressure = rho * R * T
		}
	}

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			int INDEX = (i * NY) + j;
			mass[INDEX] = p0[INDEX];
			momentum[INDEX] = p0[INDEX] * p1[INDEX];
			energy[INDEX] = 0.5 * p0[INDEX] * p1[INDEX] * p1[INDEX] + (p3[INDEX] / (GAMMA[INDEX] - 1));
			rhov[INDEX] = p0[INDEX] * p4[INDEX];
		}
	}

	while (t<t_FINAL){
	// In order to compute dt, need to find the Max wave speed first.
	float W_GLOBAL_MAX = 1.0e-10;
		for (int i = 0; i < NX; i++){
	        	for (int j = 0; j < NY; j++){
	        		float rho_L, rho_R, rho_T, rho_B, u_L, u_R, u_T, u_B, T_L, T_R, T_P, T_B, P_L, P_R, P_T, P_B,
	        		      H_L, H_R, H_T, H_P, GAMMA_L, GAMMA_R, GAMMA_T, GAMMA_B, e_L, e_R, e_T, e_B,
	        		      R_L, R_R, R_T, R_B, a_L, a_R, a_T, a_B; // L = LEFT, R = RIGHT, T = TOP, B = BOTTOM
        			int INDEX = (i * NY) + j;
				rho_L = p0[INDEX - NY];
				rho_R = p0[INDEX + NY];
				rho_T = p0[INDEX + 1];
				rho_B = p0[INDEX - 1];
				u_L = p1[INDEX - NY];
				u_R = p1[INDEX + NY];
				u_T = p1[INDEX + 1;
				u_B = p1[INDEX - 1];
				T_L = p2[INDEX - NY];
				T_R = p2[INDEX + NY];
				T_T = p2[INDEX + 1;
				T_B = p2[INDEX - 1];
				p_L = p3[INDEX - NY];
				p_R = p3[INDEX + NY];
				p_T = p3[INDEX + 1;
				p_B = p3[INDEX - 1];
				H_L = p4[INDEX - NY];
				H_R = p4[INDEX + NY];
				H_T = p4[INDEX + 1;
				H_B = p4[INDEX - 1];
				GAMMA_L = GAMMA[INDEX - NY];
				GAMMA_R = GAMMA[INDEX + NY];
				GAMMA_T = GAMMA[INDEX + 1;
				GAMMA_B = GAMMA[INDEX - 1];
				R_L = R[INDEX + NY];
				R_R = R[INDEX - NY];
				R_T = R[INDEX + 1];
				R_B = R[INDEX - 1];
				e_L = 0.5 * rho_L * u_L * u_L + (p_L / (GAMMA_L - 1));
				e_R = 0.5 * rho_R * u_R * u_R + (p_R / (GAMMA_R - 1));
				e_T = 0.5 * rho_T * u_T * u_T + (p_T / (GAMMA_T - 1));
				e_B = 0.5 * rho_B * u_B * u_B + (p_B / (GAMMA_B - 1));
				a_L = sqrt(GAMMA_L * R_L * T_L);
				a_R = sqrt(GAMMA_R * R_R * T_R);
				a_T = sqrt(GAMMA_T * R_T * T_T);
				a_B = sqrt(GAMMA_B * R_B * T_B);
			float W_LOCAL_MAX = MAX_Wave_Speed(u_L, u_R, u_T, u_B, a_L, a_R, a_T, a_B);
			Calc_Rusanov_Flux(rho_L, rho_R, rho_T, rho_B, u_L, u_R, u_T, u_B, T_L, T_R, T_T, T_B, p_L, p_R, p_T, p_B, e_L, e_R, e_T, e_B, W_LOCAL_MAX,
					  mass_flux, momentum_flux, energy_flux, rhov_flux, H_L, H_R, H_T, H_B );
			if (W_LOCAL_MAX > W_GLOBAL_MAX){
				W_GLOBAL_MAX = W_LOCAL_MAX;
			}
		}
	float dt = CFL * (dx / W_GLOBAL_MAX);

	//Set boundary condition 待修改
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

        // Use FVM to get new conservation values 待修改
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
/*
FILE * pFile = fopen("Results_of_200_cells.txt","w");
	for (int i = 0; i<N_CELLS; i++){
		fprintf(pFile, "%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", x[i], p0[i], p1[i], p2[i], p3[i], p4[i]);
	}
fclose(pFile);
*/
Free_memory(x, p0, p1, p2, p3, p4, rhov, R, CV, GAMMA, a, mass, momentum, energy, mass_flux, momentum_flux, energy_flux, rhov_flux);
return 0;
}
