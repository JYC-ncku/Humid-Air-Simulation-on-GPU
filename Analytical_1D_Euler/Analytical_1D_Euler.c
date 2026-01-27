#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Allocate_memory(double **array1, double **array2, double **array3, double **array4, double **array5, double **array6, int N_CELLS){
    *array1 = (double*)malloc(N_CELLS * sizeof(double));
    *array2 = (double*)malloc(N_CELLS * sizeof(double));
    *array3 = (double*)malloc(N_CELLS * sizeof(double));
	*array4 = (double*)malloc(N_CELLS * sizeof(double));
	*array5 = (double*)malloc((N_CELLS+1) * 5 * sizeof(double));
	*array6 = (double*)malloc((N_CELLS+1) * 5 * sizeof(double));
    if(*array1 == NULL || *array2 == NULL || *array3 == NULL || *array4 == NULL || *array5 == NULL || *array6 == NULL){
        printf("Memory allocation failed!\n");
    }
        printf("Memory allocation successfully for %d elements!\n", N_CELLS);
}

void Free_memory(double *array1, double *array2, double *array3, double *array4, double *array5, double *array6){
    free(array1);
    free(array2);
    free(array3);
	free(array4);
	free(array5);
	free(array6);
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
        double a = sqrt(1.4 * 1.0 * T); // GAMMA = 1.4, R = 1.0
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

void CPU_Calc_rho_u_P_T(double *interface_p, double *flux,
    double QL_rho, double QL_ux, double QL_vy, double QL_vz, double QL_cRT,
    double QR_rho, double QR_ux, double QR_vy, double QR_vz, double QR_cRT, double R, double GAMMA,
    double nx, double ny, double nz,
    double px, double py, double pz,
    double qx, double qy, double qz, int wall_flag){
	// DECLARE CONSTANTS
	double QL_u, QL_v, QL_w, QR_u, QR_v, QR_w; // 一個界面會有來自x,y,z三個不同方向的速度分量進入，然後我們有左邊界面跟右邊界面。
	double QR_RT, QL_RT, QR_T, QL_T, QR_E, QL_E, QR_p, QL_p, QR_a, QL_a, QL_e, QR_e;
	double sqrR, sqrL, gm1, gp1, base,expon,pwr,z;
	double uLbar, uRbar,vacuum, term1, term2, F, geff;
	double ustar, pstar, eLstar, eRstar, rhoLstar, rhoRstar, TLstar, TRstar, aLstar, aRstar;
	double dFdpstar, delp, temporary, wspeedL, wspeedR;
	double VA, VB;
	double QI_rho, QI_u, QI_p, QI_e, QI_a, QI_T, frac, QI_v, QI_w, E_tot;
	double RHOMIN = 1.0e-15; //密度的最小值
	double EMIN = 1.0e-15; //能量的最小值
	double AMIN = 1.0e-15; // 音速的最小值
	double BIGRAT = 1.5; // 大比例(用在壓力比較)，當 STAGE 1 算出來的壓力比左右兩側的壓力還大超過這個比例時就須改用牛頓疊代法求p。
	double PMIN = EMIN*RHOMIN*(GAMMA - 1.0); // 壓力的最小值
	double CV = R/(GAMMA - 1.0);
	double TMIN = EMIN/CV; //溫度的最小值
	double flxnmn, flxpmn, flxqmn;
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

	gm1 = geff - 1.0; // Pretty obvious   gm1 is mean GAMMA minus one ==> GAMMA-1
	gp1 = geff + 1.0; //                  gp1 is mean GAMMA plus one ==> GAMMA+1

    //                          計算p* and u*
	//     -----------------------------------------------------
	//     STAGE 1: Explicit solution using two isentropic waves.
	//              This gives pstar and ustar
	//     -----------------------------------------------------
	//
	//     Intermediate variable. 
    //這一階段是在利用等熵關係來預測出當兩個網格撞在一起時，中間那塊區域的壓力跟速度。
     
    //利用等熵過程中壓力與音速之間的關係式：aL/aR =(PL/PR)^((1-GAMMA)/2*GAMMA)，z就是利用這個關係式來判斷這個波形為何種波形的一個指標，為Riemann solver中的一個判斷波形的方法。 
	base = QL_p/QR_p;
	expon = gm1/(2.0*geff);
	pwr = pow(base,expon);
	z = (QR_a/QL_a)*pwr; //z為阻抗(Impedance)，這個公式的意思就是右邊的音速相對於左邊的音速。

	//   Riemann invariants. 
	uLbar = QL_u + (2.0/gm1*QL_a); //左側黎曼不變量(u+a)
	uRbar = QR_u - (2.0/gm1*QR_a); //右側黎曼不變量(u-a)

	vacuum = 0; // Assume no vacuum between states 
	if ((uLbar - uRbar) <= 0.0 ) {
	    // We have a situation in which a (near) vacuum is formed
	    //  between the left and right states.
        //如果uLbar < uRbar，表示壓力波的速度追不上流體散開的速度，中間壓力變為0而形成真空。此時須將各個properites設定為最小值，以防止程式崩潰。 
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
        //如果沒有產生真空，就可以利用前面得到的z來求出u*跟p*。
	    // Positive-pressure solution.
	    ustar = (uLbar*z + uRbar)/(1.0 + z); //如果z=1，表示左右兩邊勢均力敵，此時u*=兩邊相加除以2。如果z很大，則左邊的係數會趨近於1，右邊則趨近於0，所以ustar=uLstar，反之亦然。
	    base = (0.5*gm1*(uLbar - uRbar)/(QL_a*(1.0 + z)));
	    expon = 2.0*geff/gm1;
	    pwr = pow(base,expon);
	    pstar = QL_p*pwr; //因這段公式是推導自等熵公式，故PL*=PR*，所以不管代QL_p還是QR_p結果都一樣。
	}

	//     ------------------------------------
	//     STAGE 2: The strong-shock correction 
	//     ------------------------------------
    // 前面 STAGE 1 是假設左右兩邊都是平滑的膨脹波(Rarefaction wave)，因此可假設等熵狀態來得到各個Properties，接下來要來考慮如果是激波(Shock wave)的情況
    // Shock wave 會開始產生能量耗損，所以無法繼續假設等熵狀態來得到各個Properties。
    // 簡單來說 STAGE 1 是在假設 Rarefaction wave 的情況下利用等熵狀態來猜測一個p*值，但如果我們猜出來發現p*值太大
    // 表示不是 左右兩邊不再是 Rarefaction wave 了，會變成其中一邊或是兩邊都是 Shock wave的情況，此時我們就必須改用牛頓疊代法來求出p* and u*。
	if (vacuum == 0)  {
	    if ((pstar > (BIGRAT*QL_p)) && (pstar > (BIGRAT*QR_p))) {  //如果左右兩邊都是 Shock wave。
		// Analytical solution to strong shock using shock relations
		// if both of the pressure jumps are large enough.
		// Take 4 Newton steps to get accurate estimates of pstar and ustar
		// --- Iteration 1 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p+gm1/(2.0*geff)); // check later
		term2 = sqrt(gp1/(2.0 * geff)*pstar/QR_p+gm1/(2.0*geff));
		F = QL_u-QL_a/geff*(pstar/QL_p-1.0)/term1- QR_u- QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    * ( pstar / QL_p + (3.0 * geff - 1.0) / gp1 ) 
		    / (term1 * term1 * term1)
		    - gp1 * QR_a / (4.0 * geff * geff * QR_p)
		    * ( pstar / QR_p + (3.0 * geff - 1.0) / gp1 )
		    / (term2 * term2 * term2);
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		
		//       --- Iteration 2 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p 
		    + gm1/(2.0*geff));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1)
		    /(term1*term1*term1)
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		pstar = pstar - F/dFdpstar;
		if (pstar < PMIN) {
		    pstar = PMIN;
		}
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//        --- Iteration 3 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p 
		    + gm1/(2.0*geff));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1)
		    /(term1*term1*term1)
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		pstar = pstar - F/dFdpstar;
		if (pstar < PMIN) {
		    pstar = PMIN;
		}
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1 * pstar;
		}
		//           --- Iteration 4 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p 
		    + gm1/(2.0*geff));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1)
		    /(term1*term1*term1)
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		pstar = pstar - F/dFdpstar;
		if (pstar < PMIN) {
		    pstar = PMIN;
		}
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		
        // 計算完壓力換算速度。
		//       --- Calculate Velocity ---
		ustar = QL_u-QL_a/geff*(pstar/QL_p-1.0)/term1;

	    } else if (pstar > (BIGRAT*QR_p)) {   //如果只有右邊是 Shock wave，左邊仍是 Rarefaction wave。
		//           Treat the right-moving wave as a shock, the
		//           left-moving wave as an isentropic wave, and take
		//           four Newton steps to improve the guess for pstar, ustar.
		//           --- Iteration 1 ---
		term1 = pow((pstar/QL_p),(gm1/(2.0*geff)));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p + gm1/(2.0*geff));
		F = QL_u
		    - 2.0*QL_a/gm1*(term1-1.0)
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -QL_a/(geff*pstar)*term1
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		delp = F / dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//           --- Iteration 2 ---
		term1 = pow((pstar/QL_p),(gm1/(2.0*geff)));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - 2.0*QL_a/gm1*(term1-1.0)
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -QL_a/(geff*pstar)*term1
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		delp = F / dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//          --- Iteration 3 ---
		term1 = pow((pstar/QL_p),(gm1/(2.0*geff)));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - 2.0*QL_a/gm1*(term1-1.0)
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -QL_a/(geff*pstar)*term1
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		delp = F / dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//           --- Iteration 4 ---
		term1 = pow((pstar/QL_p),(gm1/(2.0*geff)));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - 2.0*QL_a/gm1*(term1-1.0)
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -QL_a/(geff*pstar)*term1
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		delp = F / dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//           --- Velocity ---
		ustar = QR_u + QR_a/geff*(pstar/QR_p-1.0)/term2;

	    } else if (pstar > (BIGRAT*QL_p)) {  //如果只有左邊是 Shock wave，右邊仍是 Rarefaction wave。
		//           Treat the left-moving wave as a shock, the
		//           right-moving wave as an isentropic wave, and take
		//           four Newton steps to improve the guess for pstar, ustar.
		//           --- Iteration 1 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p
		    + gm1/(2.0*geff));
		term2 = pow((pstar/QR_p),(gm1/(2.0*geff)));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - 2.0*QR_a/gm1*(term2-1.0);
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1 )
		    /(term1*term1*term1)
		    -QR_a/(geff*pstar)*term2;
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//           --- Iteration 2 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p
		    + gm1/(2.0*geff));
		term2 = pow((pstar/QR_p),(gm1/(2.0*geff)));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - 2.0*QR_a/gm1*(term2-1.0);
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1 )
		    /(term1*term1*term1)
		    -QR_a/(geff*pstar)*term2;
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//           --- Iteration 3 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p
		    + gm1/(2.0*geff));
		term2 = pow((pstar/QR_p),(gm1/(2.0*geff)));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - 2.0*QR_a/gm1*(term2-1.0);
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1 )
		    /(term1*term1*term1)
		    -QR_a/(geff*pstar)*term2;
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//          --- Iteration 4 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p
		    + gm1/(2.0*geff));
		term2 = pow((pstar/QR_p),(gm1/(2.0*geff)));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - 2.0*QR_a/gm1*(term2-1.0);
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1 )
		    /(term1*term1*term1)
		    -QR_a/(geff*pstar)*term2;
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}

		//           --- Velocity ---
		ustar = QL_u-QL_a/geff*(pstar/QL_p-1.0)/term1;
	    }
	}
	if (vacuum == 0 ) {  //在非真空狀態下
	    if (pstar > QL_p) { //如果左邊是 Shock wave。
		//   Back out values using the shock relations.
		//   Density -- from the Rankine-Hugoniot relations

        // 密度須透過Rankine-Hugoniot relations來得到
		rhoLstar = QL_rho*
		    (gp1*pstar+gm1*QL_p)/
		    (gp1*QL_p+gm1*pstar);
		
		//    Specific energy -- from the Equation of state 
		eLstar = pstar/(gm1*rhoLstar); // eL* = p* / ((GAMMA-1) * rhoL*)
		
		//    Local speed of sound -- Perfect gas version.
		aLstar = sqrt(geff*gm1*eLstar); // aL* = (GAMMA * (GAMMA-1) * eL*)
	    } else {  //如果不是 Shock wave 是 Rarefaction。
		
		//        Use the isentropic-wave relations. 
		//        Local speed of sound -- Riemann invariants. 
		aLstar = (uLbar-ustar)*0.5*gm1;
		
		//        Specific energy -- sound speed 
		eLstar = aLstar*aLstar/(geff*gm1);
		
		//       Density -- equation of state
        // 密度可透過理想狀態方程式得到
		rhoLstar = pstar/(gm1*eLstar);
	    }
	    if (pstar > QR_p) {  //如果右邊是 Shock wave。
		//           Back out values using the shock relations.
		//           Density -- from the Rankine-Hugoniot relations 
		rhoRstar = QR_rho*
		    (gp1*pstar+gm1*QR_p)/
		    (gp1*QR_p+gm1*pstar);
		
		//           Specific energy -- from the Equation of state 
		eRstar = pstar/(gm1*rhoRstar);
		
		//           Local speed of sound -- Perfect gas version. 
		aRstar = sqrt(geff*gm1*eRstar);
	    } else {  //如果不是 Shock wave 是 Rarefaction。
		
		//           Use the isentropic-wave relations. 
		//           Local speed of sound -- Riemann invariants. 
		aRstar = (ustar-uRbar)*0.5*gm1;
		
		//           Specific energy -- sound speed
		eRstar = aRstar*aRstar/(geff*gm1);
		
		//           Density -- equation of state
		rhoRstar = pstar/(gm1*eRstar);
	    }
	    
	    //    Temperatures -- equation of state also.
	    TLstar = QL_T*(pstar*QL_rho)/(QL_p*rhoLstar);
	    TRstar = QR_T*(pstar*QR_rho)/(QR_p*rhoRstar);
	}

	//     ***********
	//     Wave speeds. 
	//     ***********

	if ((pstar > QL_p) && (vacuum == 0)) {
        // Wave speed^2 = (p2-p1)/(rho1(1-(rho1/rho2))); From second Rankine-Hugoniot relations
	    //     Left wave is a shock. 
	    temporary = 0.5*gp1*QL_p/QL_rho*(pstar/QL_p+gm1/gp1); 
	    wspeedL = QL_u - sqrt(temporary);
	} else {
	    //     Left wave is an expansion fan.
        // 沒有shock wave時的wave speed。
	    wspeedL = QL_u - QL_a;
	}
	if ((pstar > QR_p) && (vacuum == 0)) {
	    //     Right wave is a shock. 
	    temporary = 0.5*gp1*QR_p/QR_rho*(pstar/QR_p+gm1/gp1);
	    wspeedR = QR_u+sqrt(temporary);
	} else {
	    //    Right wave is an expansion fan. 
	    wspeedR = QR_u + QR_a;
	}

    // 接下來要開始找"Location"！有ExapnisonL, ExpansionR, Contact, Shock.
	option = 0;


	if (ustar > 0.0) {
	    //        The left wave and the contact discontinuity determine
	    //        the cell interface quantities. 
	    if ((pstar - QL_p) >= 0.0) {
		//           The left wave is a shock. 
		if (wspeedL >= 0.0) {
		    //              All waves have gone into the right cell.
		    //              The values are taken from the left cell only. 
		    option = 1;
		} else {
		    //              The values are taken from behind the left shock.
		    option = 2;
		}
	    } else {
		//           The left wave is an expansion.
		VA = QL_u - QL_a;
		VB = ustar - aLstar;
		if (VA >= 0.0) {
		    //              All waves have gone into the right cell.
		    //              The values are taken from the left cell only.
		    option = 3;
		} else if (VB > 0.0) {
		    //              The left rarefaction straddles the cell interface.
		    //              Interpolate velocity inside the left rarefaction. 
		    option = 4;
		} else {
		    //              The values come from behind the left rarefaction. 
		    option = 5;
		}
	    }
	} else {
	    //        The right wave and the contact discontinuity determine
	    //        the cell interface quantities. 
	    if ((pstar - QR_p) >= 0.0) {
		//          The right wave is a shock. 
		if (wspeedR < 0.0) {
		    //              All waves have gone into the left cell.
		    //              The values are taken from the right cell only.
		    option = 11;
		} else {
		    //              The values are taken from behind the right shock.
		    option = 12;
		}
	    } else {
		//           The right wave is an expansion. 
		VA = QR_u + QR_a;
		VB = ustar + aRstar;
		if (VA <= 0.0) {
		    //              All waves have gone into the left cell.
		    //              The values are taken from the right cell only. 
		    option = 13;
		} else if (VB < 0.0) {
		    //              The right rarefaction straddles the cell interface.
		    //              Interpolate velocity inside the right rarefaction. 
		    option = 14;
		} else {
		    //              The values come from behind the right rarefaction.
		    option = 15;
		}
	    }
	}

	//     ************************************************
	//     *  Now, copy or interpolate the relevant data. *
	//     ************************************************

	if ((option == 1) || (option == 3)) {
	    //        All waves have gone into the right cell.
	    //        The values are taken from the left cell only. 
	    QI_rho = QL_rho;
	    QI_u   = QL_u;
	    QI_p   = QL_p;
	    QI_e   = QL_e;
	    QI_a   = QL_a;
	    QI_T   = QL_T;
	}

	if ((option == 2) || (option == 5)) {
	    //        The values are taken from behind the left wave.
	    QI_rho = rhoLstar;
	    QI_u   = ustar;
	    QI_p   = pstar;
	    QI_e   = eLstar;
	    QI_a   = aLstar;
	    QI_T   = TLstar;
	}

	if (option == 4) {
	    //        The left wave is an expansion. 
	    VA = QL_u - QL_a;
	    VB = ustar - aLstar;
	    //        The left rarefaction straddles the cell interface.
	    //        Interpolate velocity inside the left rarefaction. 
	    frac = (-VA)/(VB-VA);
	    QI_u = QL_u - frac * (QL_u - ustar);
	    //        Take the easy way out and linearly interpolate. 
	    QI_a   = QL_a-frac*(QL_a-aLstar);
	    QI_rho = QL_rho-frac*(QL_rho-rhoLstar);
	    QI_p   = QL_p-frac*(QL_p-pstar);
	    QI_e   = QL_e-frac*(QL_e-eLstar);
	    QI_T   = QL_T-frac*(QL_T-TLstar);
	}

	//     For the options below ...
	//     The right wave and the contact discontinuity determine
	//     the cell interface quantities. 

	if ((option == 11) || (option == 13)) {
	    //        All waves have gone into the left cell.
	    //        The values are taken from the right cell only. 
	    QI_rho = QR_rho;
	    QI_u   = QR_u;
	    QI_p   = QR_p;
	    QI_e   = QR_e;
	    QI_a   = QR_a;
	    QI_T   = QR_T;
	}

	if ((option == 12) || (option == 15)) {
	    //        The values are taken from behind the right wave. 
	    QI_rho = rhoRstar;
	    QI_u   = ustar;
	    QI_p   = pstar;
	    QI_e   = eRstar;
	    QI_a   = aRstar;
	    QI_T   = TRstar;
	}

	//如果 Rarefaction wave 剛好在interface上
	if (option == 14) { 
	    //        The right wave is an expansion. 
	    VA = QR_u + QR_a;
	    VB = ustar + aRstar;
	    //        The right rarefaction straddles the cell interface.
	    //        Interpolate velocity inside the right rarefaction.
		// 計算當前位置距離波尾還有多遠，frac mean fraction(比例)，意思就是當前位置與整段波長(VA-VB)的比例，可得知目前我在幾％的位置。 
	    frac = (-VB)/(VA-VB);
	    QI_u = ustar+frac*(QR_u-ustar);
	    //        Take the easy way out and linearly interpolate. 
		//在膨脹波內部，壓力、密度、速度的變化其實是等熵且非線性的（通常是冪次定律）。
		//為了省計算量，直接用「連連看」的方式，在 ustar（波尾的值）與 QR_u（波頭的值）之間一條直線，按比例 frac 取值。
		//在膨脹波內部有一套非常漂亮的等熵公式，但那要用到更複雜的次方運算。這裡選擇用線性插值，是因為算的快，且在網格夠細時精確度通常已經夠用了。
	    QI_a   = aRstar   + frac * (QR_a - aRstar);
	    QI_rho = rhoRstar + frac * (QR_rho - rhoRstar);
	    QI_p   = pstar   + frac * (QR_p - pstar);
	    QI_e   = eRstar   + frac * (QR_e - eRstar);
	    QI_T   = TRstar   + frac * (QR_T - TRstar);
	}

	//     ******************
	//     Passive Quantities.
	//     ******************
	//
	//     We assume that the transverse velocity is unaffected by
	//     the normal interactions.  We only need to select the
	//     correct value.
    // 這裡使用upwind的邏輯去計算v跟w，如果u<0表示流體從右往左流，那v跟w就得用右邊的，反之亦然。
	if (QI_u < 0.0) {
	    QI_v = QR_v;
	    QI_w = QR_w;
	} else {
	    QI_v = QL_v;
	    QI_w = QL_w;
	}

	//     ******** **********
	//     Combine the fluxes
	//     ******************

	//    Mass/unit-area/unit-time
	mflx = QI_rho*QI_u; // Mss flux

	//    Normal momentum
	flxnmn = QI_rho*QI_u*QI_u + QI_p;

	//    Tangential momentums
	flxpmn = mflx*QI_v;
	flxqmn = mflx*QI_w;

	//    Energy Flux
	E_tot = QI_e + 0.5 * (QI_u*QI_u + QI_v*QI_v + QI_w*QI_w);
	eflx = QI_rho*E_tot*QI_u + QI_p*QI_u; // Energy flux

	//     Convert back to global axes
	pxflx = flxnmn*nx + flxpmn*px + flxqmn*qx; // Momentum fluxes
	pyflx = flxnmn*ny + flxpmn*py + flxqmn*qy;
	pzflx = flxnmn*nz + flxpmn*pz + flxqmn*qz;

	// Final Flux calculations
	flux[0] = mflx;
	flux[1] = pxflx;
	flux[2] = pyflx;
	flux[3] = pzflx;
	flux[4] = eflx;

	// States now
	interface_p[0] = QI_rho;
	interface_p[1] = QI_u;
	interface_p[2] = QI_v;
	interface_p[3] = QI_w;
	interface_p[4] = QI_T;	
}

int main(){
    int N_CELLS = 200;
	float L = 1.0;
    float dx = L/N_CELLS; 
    float dt = 0.0001; //隨便設，如果CFL有error就調整dt值。
	float t = 0;
	float t_FINAL = 0.2;
	double R = 1.0;
	double GAMMA = 1.4;
	int wall_flag = 0; 
    double *x, *p0, *p1, *p2, *interface_p, *flux;
	double flxnmn, flxpmn, flxqmn;

	Allocate_memory(&x, &p0, &p1, &p2, &interface_p, &flux, N_CELLS);
	//Initial condition
	for ( int i = 0; i<N_CELLS; i++){
		x[i] = (i+0.5) * dx;
		if (i<N_CELLS/2){
			p0[i] = 10.0; //rho_L = 10
			p1[i] = 0; //u_L = 0
			p2[i] = 1.0; // T_L = 1
		}else{
			p0[i] = 1.0; //rho_R = 1
			p1[i] = 0; //u_R = 0
			p2[i] = 1.0; //T_R = 0
		}
	}
	// Because this code only consider 1D, so let other two direction equal 0)
	double QL_vy = 0, QL_vz = 0;
	double QR_vy = 0, QR_vz = 0;
    double nx = 1.0, ny = 0.0, nz = 0.0;
    double px = 0.0, py = 1.0, pz = 0.0;
    double qx = 0.0, qy = 0.0, qz = 1.0; 

	while (t<t_FINAL){
    	double MAX_CFL = CPU_Compute_MAX_CFL(p0, p1, p2, dx, dt, N_CELLS);

		for (int i = 0; i < N_CELLS-1; i++){
			double QL_rho = p0[i];
    		double QR_rho = p0[i+1];
    		double QL_ux  = p1[i];
    		double QR_ux  = p1[i+1];
    		double QL_T   = p2[i];
    		double QR_T   = p2[i+1];
			double QL_cRT = sqrt(R * QL_T);
    		double QR_cRT = sqrt(R  * QR_T);
			CPU_Calc_rho_u_P_T(&interface_p[i*5], &flux[i*5], //因為flux跟interface_p都有5個物理量需要儲存，如果不加這行的話數據就會一直不斷被覆蓋，最後變成只有儲存到最後一格的資料。
    		QL_rho, QL_ux, QL_vy, QL_vz, QL_cRT,
    		QR_rho, QR_ux, QR_vy, QR_vz, QR_cRT, R, GAMMA,
     		nx, ny, nz,
     		px, py, pz,
     		qx, qy, qz, wall_flag);
		}

		for (int i = 1; i < N_CELLS - 1; i++) {
    		// 我們是 i*5，所以左界面是 (i-1)*5，右界面是 i*5
    		int L_interface = (i - 1) * 5;
    		int R_interface = i * 5;
			double CV = R / (GAMMA - 1.0);
			// 先將舊的值儲存起來
    		double rho_old = p0[i];
    		double u_old   = p1[i];
    		double T_old   = p2[i];
			// p1 存的是速度 u，我們要先算動量 rho*u 的變化再去除以rho得到u。
    		double Mom_old = rho_old * u_old; // p0[i] * p1[i]
			// 先從溫度算總能 E，更新完 E 再扣掉動能回算 T
    		double E_old   = rho_old * (CV * T_old + 0.5 * u_old * u_old); //p0[i] * (CV * p2[i] + 0.5 * p1[i] * p1[i])
    	

			// 使用FVM計算新的值，參考程式碼712～716行，interface_p[0]是密度、[1]是u、[2]是v、[3]是w、[4]是溫度。
    		double rho_new = rho_old + (dt / dx) * (flux[L_interface + 0] - flux[R_interface + 0]);
    		double Mom_new = Mom_old + (dt / dx) * (flux[L_interface + 1] - flux[R_interface + 1]);
    		double E_new   = E_old   + (dt / dx) * (flux[L_interface + 4] - flux[R_interface + 4]);

    		//更新密度 (p0)
    		p0[i] = rho_new;
    		// 更新動量並回推速度 (p1)
    		p1[i] = Mom_new / rho_new;
    		// 更新能量並回推溫度 (p2)
			double internal_e = (E_new / rho_new) - 0.5 * (p1[i] * p1[i]);
			p2[i] = internal_e / CV;
		}

		// Boundary condition for compute flux.	
		// 左邊界
		p0[0] = p0[1];
		p1[0] = p1[1];
		p2[0] = p2[1];

		// 右邊界
		p0[N_CELLS-1] = p0[N_CELLS-2];
		p1[N_CELLS-1] = p1[N_CELLS-2];
		p2[N_CELLS-1] = p2[N_CELLS-2];

		t += dt;
	}
	
	FILE * pFile = fopen("Results_of_200_cells.txt","w");
    for (int i=0; i<N_CELLS; i++){  
        fprintf(pFile, "%.3f\t%.6f\t%.6f\t%.6f\t%.2f\n", x[i], p0[i], p1[i], p2[i], t);
    }
    fclose(pFile);

    Free_memory(x, p0, p1, p2, interface_p, flux);
	return 0;
}