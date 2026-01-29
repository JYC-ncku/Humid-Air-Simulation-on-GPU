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

double MAX_CFL(double *p0, double *p1, double *p2, float dX, float dt, float N_CELLS) {
	float MAX_CFL = -1.0;
	for (int cell = 0; cell < N_CELLS; cell++) {
		// Extract rho, u and T out from the arrays
		double rho = p0[cell];
		double u = p1[cell];
		double T = p2[cell];
		// Check the temperature
		if (T < 0) {
			printf("Error: Negative temperature in cell %d: T = %f\n. Aborting.", cell, T);
			exit(1);
		}
		double a = sqrt(1.4 * 1.0 * T); // GAMMA = 1.4, R = 1.0
		double CFL = (fabs(u) + a) * (dt / dX);
		if (CFL > MAX_CFL) {
			MAX_CFL = CFL;
		}
	}

	// Check the return value
	if (MAX_CFL < 0) {
		printf("Error: MAX_CFL is negative!. Aborting. \n");
		exit(1);
	} else {
		return MAX_CFL;
	}
}



int main(){
    int N_CELLS = 200;
    float L = 1.0;
    float R = 1.0;
    float GAMMA = 1.4;
    float t = 0;
    float t_FINAL = 0.2;
    float dt = 0.0001;
    double dx = L/N_CELLS;
    double *x, *p0, *p1, *p2, *flux;
    Allocate_memory(&x, &p0, &p1, &p2, &flux, N_CELLS);
    Free_memory(x, p0, p1, p2, flux);
}