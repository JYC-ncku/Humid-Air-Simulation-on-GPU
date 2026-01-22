#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(double **array1, double **array2, double **array3, int N_CELLS){
    *array1 = (double*)malloc(N_CELLS * sizeof(double));
    *array2 = (double*)malloc(N_CELLS * sizeof(double));
    *array3 = (double*)malloc(N_CELLS * sizeof(double));
    if(*array1 == NULL || *array2 == NULL || *array3 == NULL){
        printf("Memory allocation failed!\n");
    }
        printf("Memory allocation successfully for %d elements!\n", N_CELLS);
}

void Free_memory(double *array1, double *array2, double *array3){
    free(array1);
    free(array2);
    free(array3);
    printf("Memory freed successfully!\n");
}

// Return the maximum CFL number across all cells
double CPU_Compute_MAX_CFL(double *p0, double *p1, double *p2, float dx, float dt, int N_CELLS){
    float MAX_CFL = -1.0;
    for (int cell=0; cell<N_CELLS; cell++){
        double rho = p0[cell];
        double u = p1[cell];
        double T = p2[cell];
        if (T<0){
            printf("Error: Negative temperature in cell %d: T = %f\n. Aborting.", cell, T);
            exit(1);
        }
        double a = sqrt(1.4 * 1.0 *T); // gamma = 1.4, R = 1.0
        double CFL = (fabs(u) + a) * (dt / dx);
        if (CFL>MAX_CFL){
            MAX_CFL = CFL;
        }
    }

    if(MAX_CFL<0){
        printf("Error: MAX_CFL is negative!. Aborting. \n");
        exit(1);
    } else {
        return MAX_CFL;
    }
}
	//     Purpose...
	//     -------
	//     Given the flow states either side of the interface,
	//     compute the interface flow state and the fluxes across
	//     the interface.
	// 
	//     Input...
	//     -----
	//     real(dbl), intent(in) :: QL_rho,QL_ux,QL_vy,QL_cRT
	//     QL_rho  : Left density
	//     QL_ux   : Left x-velocity
	//     QL_vy   : Left y-velocity
	//     QL_cRT  : Left thermal velocity sqrt(RT)
	//     QR_rho, QR_ux, QR_vy, QR_cRT : Right state
	//     real(dbl), intent(in) :: QR_rho,QR_ux,QR_vy,QR_cRT
	// 
	//     nx, ny  : unit normals for the interface
	//     px, py  : unit normals for the interface
	//     real(dbl), intent(in) :: nx,ny,px,py
	//     Output...
	//     ------
	//     flxmss  : mass flux across the interface
	//              (mass/unit-area/unit-time)
	//     flxmnx  : flux of x-momentum across the interface
	//     flxmny  : flux of y-momentum
	//     flxeng  : flux of energy

	// DECLARE CONSTANTS
	double QL_u, QL_v, QL_w, QR_u, QR_v, QR_w; // 一個界面會有來自x,y,z三個不同方向的速度分量進入，然後我們有左邊界面跟右邊界面。
	double QR_RT, QL_RT, QR_T, QL_T, QR_E, QL_E, QR_p, QL_p, QR_a, QL_a, QL_e, QR_e;
	double sqrR, sqrL, gm1, gp1, base,expon,pwr,z;
	double uLbar, uRbar,vacuum, term1, term2, F, geff;
	double ustar, pstar, eLstar, eRstar, rhoLstar, rhoRstar, TLstar, TRstar, aLstar, aRstar;
	double dFdpstar, delp, temporary, wspeedL, wspeedR;
	double VA, VB;
	double QI_rho, QI_u, QI_p, QI_e, QI_a, QI_T, frac, QI_v, QI_w, E_tot;
	double RHOMIN = 1.0e-15; 
	double EMIN = 1.0e-15;
	double AMIN = 1.0e-15;
	double BIGRAT = 1.5;
	double PMIN = EMIN*RHOMIN*(GAMMA - 1.0); 
	double CV = R/(GAMMA - 1.0);
	double TMIN = EMIN/CV;
    double mflx, pxflx, pyflx, pzflx, eflx; 
	int option;
	
	// Left hand normals
	QL_u = nx*QL_ux + ny*QL_vy + nz*QL_vz; 
	QL_v = px*QL_ux + py*QL_vy + pz*QL_vz;
	QL_w = qx*QL_ux + qy*QL_vy + qz*QL_vz;
	
	// Right hand normals
	QR_u = nx*QR_ux + ny*QR_vy + nz*QR_vz;
	QR_v = px*QR_ux + py*QR_vy + pz*QR_vz;
	QR_w = qx*QR_ux + qy*QR_vy + qz*QR_vz;

    // Switch normal components in case one is a wall (Reflective conditions) ==> 邊界條件
	if (wall_flag == 1) {
		QR_u = -1.0*QR_u;
	} else if (wall_flag == -1) {
		QL_u = -1.0*QL_u;
	}

	// Compute the initial state properties.
	QR_RT = QR_cRT*QR_cRT; //cRT is thermal velocity, cRT = sqrt(RT), cRT also mean Isothermal Sound Speed. So sqrt(RT) * sqrt(RT) = RT.
	QL_RT = QL_cRT*QL_cRT;
	QL_T = QL_RT/R; // RT/R = T
	QR_T = QR_RT/R;
	QL_e = QL_T*CV; // e = T * CV
	QR_e = QR_T*CV;
	QL_p = QL_rho*QL_RT; // p = rho * R * T (ideal gas equation)
	QR_p = QR_rho*QR_RT;
	QL_a = sqrt(GAMMA*QL_RT); // a is sound speed ==> a = sqrt(GAMMA * R * T) = sqrt(GAMMA) * cRT
	QR_a = sqrt(GAMMA*QR_RT);
	sqrL = sqrt(QL_rho);
	sqrR = sqrt(QR_rho);

	geff = GAMMA; // Effective GAMMA

	gm1 = geff - 1.0; // Pretty obvious   gm1 is mean gamma minus one ==> gamma-1
	gp1 = geff + 1.0; //                  gp1 is mean gamma plus one ==> gamma+1

    //                          計算p* and u*
	//     -----------------------------------------------------
	//     STAGE 1: Explicit solution using two isentropic waves.
	//              This gives pstar and ustar
	//     -----------------------------------------------------
	//
	//     Intermediate variable. 
    //利用等熵過程中壓力語音速之間的關係式：aL/aR =(PL/PR)^((1-gamma)/2*gamma)，z就是利用這個關係式來判斷這個波形為何種波形的一個指標，為Riemann solver中的一個方法。 
	base = QL_p/QR_p;
	expon = gm1/(2.0*geff);
	pwr = pow(base,expon);
	z = (QR_a/QL_a)*pwr;

	//   Riemann invariants. 
	uLbar = QL_u + (2.0/gm1*QL_a);
	uRbar = QR_u - (2.0/gm1*QR_a);

	vacuum = 0; // Assume no vacuum between states 
	if ((uLbar - uRbar) <= 0.0 ) {
	    // We have a situation in which a (near) vacuum is formed
	    //  between the left and right states.
	    vacuum = 1; // A near vacuum is present
	    ustar = 0.0;
	    pstar = PMIN; 
	    eLstar = EMIN;
	    eRstar = EMIN;
	    rhoLstar = RHOMIN;
	    rhoRstar = RHOMIN;
	    TLstar = TMIN;
	    TRstar = TMIN;
	    aLstar = AMIN;
	    aRstar = AMIN;
	}
	if (vacuum == 0) {
	    // Positive-pressure solution.
	    ustar = (uLbar*z + uRbar)/(1.0 + z);
	    base = (0.5*gm1*(uLbar - uRbar)/(QL_a*(1.0 + z)));
	    expon = 2.0*geff/gm1;
	    pwr = pow(base,expon);
	    pstar = QL_p*pwr;
	}

int main(){
    int N_CELLS = 200;
    float L = 1.0;
    double *p0, *p1, *p2; 
    float dx = L/N_CELLS; 
    float dt = 0.0001; //隨便設，如果CFL有error就調整dt值。
    Allocate_memory(&p0, &p1, &p2, N_CELLS);
    CPU_Compute_MAX_CFL(p0, p1, p2, dx, dt, N_CELLS);
    Free_memory(p0, p1, p2);
}