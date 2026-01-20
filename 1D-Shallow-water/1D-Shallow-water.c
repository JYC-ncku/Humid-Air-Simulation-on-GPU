#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(double **array1, double **array2, double **array3, double **array4, double **array5, double **array6, double **array7, int N_CELLS, int N_INTERFACES){
    *array1 = (double*)malloc(N_CELLS * sizeof(double));
    *array2 = (double*)malloc(N_CELLS * sizeof(double));
    *array3 = (double*)malloc(N_CELLS * sizeof(double));
    *array4 = (double*)malloc(N_INTERFACES * sizeof(double));
    *array5 = (double*)malloc(N_INTERFACES * sizeof(double));
    *array6 = (double*)malloc(N_CELLS * sizeof(double));
    *array7 = (double*)malloc(N_CELLS * sizeof(double));
    if(*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL || *array7 == NULL){
        printf("Memory allocation faild!\n");
        exit(1);
    }
    printf("Memory allocation successfully for %d elements\n", N_CELLS);
}

void Free_memory(double *array1, double *array2, double *array3, double *array4, double *array5, double *array6, double *array7){
    free(array1);
    free(array2);
    free(array3);
    free(array4);
    free(array5);
    free(array6);
    free(array7);
    printf("Memory freed successfully!\n");
}

//Setting parameters
    double L = 1.0;
    double g = 9.81;

/*
    ----------------------------
    |        |        |        |
    |   0    |  ...   |    N   |   ===>   Have total N cells and N+1 interface. (because N is from 0)
    |        |        |        | 
    ----------------------------
    0        1       N        N+1
*/

void Shallow_Water(int N_CELLS, int N_INTERFACES, double *x, double *mass, double *momentum, double *mass_flux, double *momentum_flux, double *p0, double *p1){
    double dx = L/N_CELLS;
    double CFL = 0.8;
    double t = 0; // Initialize time.
    double t_final = 0.5; // Simulation time.

    // Set the initial co conditions
    for (int i=0; i<N_CELLS; i++){
        x[i] = (i+0.5) * dx;
        if (i < N_CELLS/2){
            p0[i] = 1; // Initial water high = 1m in 0 to 0.5m.
        } else {
            p0[i] = 0.1; // Initial water high = 0.1m in 0.5 to 1m (Shallow Water).
        }
        p1[i] = 0; // Water speed = 0 m/s in initial condition.
    }    

    for (int i=0; i<N_CELLS; i++){
    mass[i] = p0[i];
    momentum[i] = p0[i] * p1[i];
    }

    while (t <= t_final){
        double S_max, S_local;
        // Set CFL
        for (int j=1; j<N_CELLS;j++){
        double S_L = fabs(p1[j-1]) + sqrt(g*p0[j-1]);
        double S_R = fabs(p1[j]) + sqrt(g*p0[j]);
            if (S_L>S_R){
                S_local = S_L;
            }else {
                S_local = S_R;
            }
            if (S_local>S_max){
                S_max = S_local;
            } 
        }
        double dt = CFL * (dx/S_max);
    // Compute the fluxes now using both left and right
    // left cell = interface -1, right cell = interface
    for (int j=1; j<N_CELLS; j++ ){
        double mass_flux_L =  p0[j-1]*p1[j-1]; // left mass flux
        double momentum_flux_L = p0[j-1]*p1[j-1]*p1[j-1] + 0.5*g*p0[j-1]*p0[j-1]; // left momentum flux
        double mass_flux_R = p0[j]*p1[j]; // right mass flux
        double momentum_flux_R = p0[j]*p1[j]*p1[j] + 0.5*g*p0[j]*p0[j]; //right momentum flux
    // Compute the Rusanov Flux
        mass_flux[j] = 0.5*(mass_flux_L + mass_flux_R)-0.5*S_max* (mass[j]-mass[j-1]);
        momentum_flux[j] = 0.5*(momentum_flux_L + momentum_flux_R)-0.5*S_max*(momentum[j]-momentum[j-1]);
    }

    //Set boundary condition
    mass_flux[0] = mass_flux[1];
    momentum_flux[0] = momentum_flux[1];
    mass_flux[N_CELLS] = mass_flux[N_CELLS-1]; //我們的網格是從0開始，最後一格為N-1。
    momentum_flux[N_CELLS] = momentum_flux[N_CELLS-1]; 
    // Compute fluxes
    for (int i=0; i<N_CELLS; i++){
        mass[i] = mass[i] - (dt/dx)*(mass_flux[i+1]-mass_flux[i]);
        momentum[i] = momentum[i] - (dt/dx)*(momentum_flux[i+1]-momentum_flux[i]);
    }
    for (int i=0; i<N_CELLS; i++){
        p0[i]=mass[i];
        p1[i]=momentum[i]/mass[i];
    }
    t += dt;
    }
}

int main(){
    const int N_CELLS = 800; // Number of cells
    const int N_INTERFACES = N_CELLS+1; // Number of interfaces
    double *x, *mass, *momentum, *mass_flux, *momentum_flux, *p0, *p1; // x is position, h is water high, u is water velocity, po is primitive values of h, p1 is primitive values of u

    Allocate_memory(&x, &mass, &momentum, &mass_flux, &momentum_flux, &p0, &p1, N_CELLS, N_INTERFACES);

    Shallow_Water(N_CELLS, N_INTERFACES, x, mass, momentum, mass_flux, momentum_flux, p0, p1);
    FILE *  pFile = fopen("Results_of_800_cells.txt","w");
    for (int i=0; i<N_CELLS; i++){  
        fprintf(pFile, "%.3f\t%.6f\t%.6f\n", x[i], p0[i], p1[i]);
    }
    fclose(pFile);

    Free_memory(x, mass, momentum, mass_flux, momentum_flux, p0, p1);
}