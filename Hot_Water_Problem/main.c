#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memory.h"
#include "Calc_rho_u_P_T.h"
#include "Boundary.h"

/*
	  GHOST				  GHOST
	一一一一一一一一一一一一一一一一一一一一
	|	|	|       |	|	|
	|   0   |   1   |  ...  |   N   |  N+1  |
	|       |       |       |	|	|
	一一一一一一一一一一一一一一一一一一一一
		0	1      ...	N
*/

float MINMOD(float QL_rho, float QC_rho, float QR_rho, float dx){
	float dU_dx;
	float Forward = (QR_rho - QC_rho) / dx;
	float Backward = (QC_rho - QL_rho) / dx;
	if (Backward * Forward < 0){
		dU_dx = 0;
	} else if ( fabs(Forward) < fabs(Backward) ){
			dU_dx = Forward;
		} else {
			dU_dx = Backward;
	}
	return dU_dx;
}

int main(){
	int NX = 1000;
	int NY = 1000;
	int N_CELLS = (NX+4) * (NY+4); //+2 for Ghost cells
	float L = 1.0;
	float H = 1.0;
	float dx = L/NX;
	float dy = H/NY;
	float t = 0;
	float t_FINAL = 0.2;

	float Ru = 8314.5; // unit: J/mole-K
	float R_dry = Ru / 28.97; // Gas constant of dry air (MW_air = 28.97 g/mole)
	float R_v = Ru / 18.02; // Gas constant of water vapor (MW_water = 18.02 g/mole)

	float Cp_dry = 1005; // unit: J/kg-K
	float Cp_v = 1864; // unit: J/kg-K (water vapor, not liquid! liquid is 4.179)
	float D = 2.42e-5;
//	float GAMMA = 1.4;
	int wall_flag = 0;
	float *x, *y, *p0, *p1, *p2, *p3, *p4, *p5, *p6, *interface_p, *flux_X, *flux_Y;
	//p0 is density, p1 is x-dir velocity, p2 is y-dir veloctiy, p3 is temperature, p4 is pressure, p5 is mass fraction of water vapor, p6 is relative humidity
	float flxnmn, flxpmn, flxqmn;
	float CFL = 0.5;

	Allocate_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &p5, &p6, &interface_p, &flux_X, &flux_Y, N_CELLS);
	//Initial condition
	for ( int i = 2; i < NX + 2; i++){
		for (int j = 2; j < NY + 2; j++){
			int INDEX = i * (NY + 4) + j;
			p1[INDEX] = 0.0;
			p2[INDEX] = 0.0;
			p3[INDEX] = 298.15; // Room temperature = 25 C = 273.15 + 25 = 298.15 K
			p4[INDEX] = 101325; // 1 atm = 101325 Pa
			p0[INDEX] = p4[INDEX] / (R_dry * p3[INDEX]); //Density: rho = P / (R * T)
			p5[INDEX] = 0.0; // Mass fraction of water vapor is 0 at every cells
		}
	}
/* For x-dir
	// Because this code only consider 2D, so let z direction equal 0)
	float QL_vy = 1.0, QL_vz = 0;
	float QR_vy = 1.0, QR_vz = 0;
	float nx = 1.0, ny = 0.0, nz = 0.0;
	float px = 0.0, py = 1.0, pz = 0.0;
	float qx = 0.0, qy = 0.0, qz = 1.0;
*/

/* For y-dir
	// Because this code only consider 2D, so let z direction equal 0)
	float QL_vy = 1.0, QL_vz = 0;
	float QR_vy = 1.0, QR_vz = 0;
	float nx = 0.0, ny = 1.0, nz = 0.0;
	float px = -1.0, py = 0.0, pz = 0.0;
	float qx = 0.0, qy = 0.0, qz = 1.0;
*/
	while (t<t_FINAL){
		// Boundary condition for compute flux.
		Boundary(p0, p1, p2, p3, p4, p5, Ru, R_v, R_dry, NX, NY);
	    	float MAX_CFL = CPU_Compute_MAX_CFL(p0, p1, p2, p3, dx, dy, NX, NY);
		float dt_advection = CFL / MAX_CFL;
		float dt_diffusion = 0.25 / (D * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
		float dt = fmin(dt_advection, dt_diffusion);
	    	//X-dir (flux_X)
		for (int i = 1; i < NX + 2; i++){		//N cells have N+1 interface
			for (int j = 2; j < NY + 2; j++){	//j = 2 ~ 101(NY+1) is real cells
				int INDEX = i * (NY + 4) + j;
				int INDEX_R = (i + 1) * (NY + 4) + j;
				int INDEX_RR = (i + 2) * (NY + 4) + j;
				int INDEX_L = (i - 1) * (NY + 4) + j;
//				float R_mix = (1 - p5[INDEX]) * R_dry + p5[INDEX] * R_v;
//				float Cp_mix =(1 - p5[INDEX]) * Cp_dry + p5[INDEX] * Cp_v;
//				float Cv_mix = Cp_mix - R_mix;
//				float GAMMA = Cp_mix / Cv_mix;

				float QL_rho = p0[INDEX_L];
				float QC_rho = p0[INDEX];
		    		float QR_rho = p0[INDEX_R];
		    		float QRR_rho = p0[INDEX_RR];
		    		float drho_dx_L = MINMOD(QL_rho, QC_rho, QR_rho, dx);
		    		float drho_dx_R = MINMOD(QC_rho, QR_rho, QRR_rho, dx);
		    		float QL_rho_star = QC_rho + 0.5 * dx * drho_dx_L;
		    		float QR_rho_star = QR_rho - 0.5 * dx * drho_dx_R;

    				float QL_ux = p1[INDEX_L];
    				float QC_ux = p1[INDEX];
    				float QR_ux = p1[INDEX_R];
    				float QRR_ux = p1[INDEX_RR];
		    		float du_dx_L = MINMOD(QL_ux, QC_ux, QR_ux, dx);
		    		float du_dx_R = MINMOD(QC_ux, QR_ux, QRR_ux, dx);
		    		float QL_ux_star = QC_ux + 0.5 * dx * du_dx_L;
		    		float QR_ux_star = QR_ux - 0.5 * dx * du_dx_R;

    				float QL_vy = p2[INDEX_L];
    				float QC_vy = p2[INDEX];
    				float QR_vy = p2[INDEX_R];
    				float QRR_vy = p2[INDEX_RR];
		    		float dv_dx_L = MINMOD(QL_vy, QC_vy, QR_vy, dx);
		    		float dv_dx_R = MINMOD(QC_vy, QR_vy, QRR_vy, dx);
		    		float QL_vy_star = QC_vy + 0.5 * dx * dv_dx_L;
		    		float QR_vy_star = QR_vy - 0.5 * dx * dv_dx_R;

    				float QL_vz  = 0.0;
    				float QR_vz  = 0.0;

		    		float QL_T = p3[INDEX_L];
		    		float QC_T = p3[INDEX];
    				float QR_T = p3[INDEX_R];
    				float QRR_T = p3[INDEX_RR];
		    		float dT_dx_L = MINMOD(QL_T, QC_T, QR_T, dx);
		    		float dT_dx_R = MINMOD(QC_T, QR_T, QRR_T, dx);
		    		float QL_T_star = QC_T + 0.5 * dx * dT_dx_L;
		    		float QR_T_star = QR_T - 0.5 * dx * dT_dx_R;

				// Calculate gas constant at interface for calculate Flux
				float R_L = (1 - p5[INDEX_L]) * R_dry + p5[INDEX_L] * R_v;
				float R_C = (1 - p5[INDEX]) * R_dry + p5[INDEX] * R_v;
				float R_R = (1 - p5[INDEX_R]) * R_dry + p5[INDEX_R] * R_v;
				float R_RR = (1 - p5[INDEX_RR]) * R_dry + p5[INDEX_RR] * R_v;

				float QL_cRT = sqrt(R_L * QL_T);
				float QC_cRT = sqrt(R_C * QC_T);
  				float QR_cRT = sqrt(R_R * QR_T);
				float QRR_cRT = sqrt(R_RR * QRR_T);
		    		float dcRT_dx_L = MINMOD(QL_cRT, QC_cRT, QR_cRT, dx);
		    		float dcRT_dx_R = MINMOD(QC_cRT, QR_cRT, QRR_cRT, dx);
		    		float QL_cRT_star = QC_cRT + 0.5 * dx * dcRT_dx_L;
		    		float QR_cRT_star = QR_cRT - 0.5 * dx * dcRT_dx_R;

		    		float QL_Y = p5[INDEX_L];
		    		float QC_Y = p5[INDEX];
    				float QR_Y = p5[INDEX_R];
    				float QRR_Y = p5[INDEX_RR];
		    		float dY_dx_L = MINMOD(QL_Y, QC_Y, QR_Y, dx);
		    		float dY_dx_R = MINMOD(QC_Y, QR_Y, QRR_Y, dx);
		    		float QL_Y_star = QC_Y + 0.5 * dx * dY_dx_L;
		    		float QR_Y_star = QR_Y - 0.5 * dx * dY_dx_R;

				float Y_face = 0.5 * (QL_Y_star + QR_Y_star);
				float R_mix_face = (1 - Y_face) * R_dry + Y_face * R_v;
				float Cp_mix_face = (1 - Y_face) * Cp_dry + Y_face * Cp_v;
				float GAMMA_face = Cp_mix_face / (Cp_mix_face - R_mix_face);

				CPU_Calc_rho_u_P_T(&interface_p[INDEX*6], &flux_X[INDEX*6], //因為flux跟interface_p都有6個物理量需要儲存，如果不加這行的話數據就會一直不斷被覆蓋，最後變成只有儲存到最後一格的資料。
						   QL_rho_star, QL_ux_star, QL_vy_star, QL_vz, QL_cRT_star, QL_Y_star,
						   QR_rho_star, QR_ux_star, QR_vy_star, QR_vz, QR_cRT_star, QR_Y_star, R_mix_face, GAMMA_face,
						   1.0, 0.0, 0.0,
						   0.0, 1.0, 0.0,
						   0.0, 0.0, 1.0, wall_flag);

				float rho_face_X = 0.5 * (QC_rho + QR_rho);
				float Diff_flux_X = rho_face_X * D * ((QR_Y - QC_Y) / dx); // Central difference
				flux_X[INDEX*6 + 5] -= Diff_flux_X; // Total flux = Advection flux - diffusion flux, flux_X is adveciton flux from FVM
			}
		}

		//Y-dir (flux_Y) 把X軸往逆時針轉90度看。
		for (int i = 2; i < NX + 2; i++){		//200 cells have 201 interface
			for (int j = 1; j < NY + 2; j++){
				int INDEX_B = i * (NY + 4) + (j - 1);
				int INDEX = i * (NY + 4) + j;
				int INDEX_T = i * (NY + 4) + (j + 1);
				int INDEX_TT = i * (NY + 4) + (j + 2);
//				float R_mix = (1 - p5[INDEX]) * R_dry + p5[INDEX] * R_v;
//				float Cp_mix =(1 - p5[INDEX]) * Cp_dry + p5[INDEX] * Cp_v;
//				float Cv_mix = Cp_mix - R_mix;
//				float GAMMA = Cp_mix / Cv_mix;

				float QB_rho = p0[INDEX_B];
				float QC_rho = p0[INDEX];
		    		float QT_rho = p0[INDEX_T];
		    		float QTT_rho = p0[INDEX_TT];
		    		float drho_dy_B = MINMOD(QB_rho, QC_rho, QT_rho, dy);
		    		float drho_dy_T = MINMOD(QC_rho, QT_rho, QTT_rho, dy);
		    		float QB_rho_star = QC_rho + 0.5 * dy * drho_dy_B;
		    		float QT_rho_star = QT_rho - 0.5 * dy * drho_dy_T;

    				float QB_ux  = p1[INDEX_B];
    				float QC_ux  = p1[INDEX];
    				float QT_ux  = p1[INDEX_T];
    				float QTT_ux  = p1[INDEX_TT];
		    		float du_dy_B = MINMOD(QB_ux, QC_ux, QT_ux, dy);
		    		float du_dy_T = MINMOD(QC_ux, QT_ux, QTT_ux, dy);
		    		float QB_ux_star = QC_ux + 0.5 * dy * du_dy_B;
		    		float QT_ux_star = QT_ux - 0.5 * dy * du_dy_T;

    				float QB_vy = p2[INDEX_B];
    				float QC_vy = p2[INDEX];
    				float QT_vy = p2[INDEX_T];
    				float QTT_vy = p2[INDEX_TT];
		    		float dv_dy_B = MINMOD(QB_vy, QC_vy, QT_vy, dy);
		    		float dv_dy_T = MINMOD(QC_vy, QT_vy, QTT_vy, dy);
		    		float QB_vy_star = QC_vy + 0.5 * dy * dv_dy_B;
		    		float QT_vy_star = QT_vy - 0.5 * dy * dv_dy_T;

    				float QB_vz  = 0.0;
    				float QT_vz  = 0.0;

		    		float QB_T = p3[INDEX_B];
		    		float QC_T = p3[INDEX];
    				float QT_T = p3[INDEX_T];
    				float QTT_T = p3[INDEX_TT];
		    		float dT_dy_B = MINMOD(QB_T, QC_T, QT_T, dy);
		    		float dT_dy_T = MINMOD(QC_T, QT_T, QTT_T, dy);
		    		float QB_T_star = QC_T + 0.5 * dy * dT_dy_B;
		    		float QT_T_star = QT_T - 0.5 * dy * dT_dy_T;

				// Calculate gas constant at interface for calculate Flux
				float R_B = (1 - p5[INDEX_B]) * R_dry + p5[INDEX_B] * R_v;
				float R_C = (1 - p5[INDEX]) * R_dry + p5[INDEX] * R_v;
				float R_T = (1 - p5[INDEX_T]) * R_dry + p5[INDEX_T] * R_v;
				float R_TT = (1 - p5[INDEX_TT]) * R_dry + p5[INDEX_TT] * R_v;

				float QB_cRT = sqrt(R_B * QB_T);
				float QC_cRT = sqrt(R_C * QC_T);
    				float QT_cRT = sqrt(R_T * QT_T);
    				float QTT_cRT = sqrt(R_TT * QTT_T);
		    		float dcRT_dy_B = MINMOD(QB_cRT, QC_cRT, QT_cRT, dy);
		    		float dcRT_dy_T = MINMOD(QC_cRT, QT_cRT, QTT_cRT, dy);
		    		float QB_cRT_star = QC_cRT + 0.5 * dy * dcRT_dy_B;
		    		float QT_cRT_star = QT_cRT - 0.5 * dy * dcRT_dy_T;

		    		float QB_Y = p5[INDEX_B];
		    		float QC_Y = p5[INDEX];
    				float QT_Y = p5[INDEX_T];
    				float QTT_Y = p5[INDEX_TT];
		    		float dY_dy_B = MINMOD(QB_Y, QC_Y, QT_Y, dy);
		    		float dY_dy_T = MINMOD(QC_Y, QT_Y, QTT_Y, dy);
		    		float QB_Y_star = QC_Y + 0.5 * dy * dY_dy_B;
		    		float QT_Y_star = QT_Y - 0.5 * dy * dY_dy_T;

				float Y_face = 0.5 * (QB_Y_star + QT_Y_star);
				float R_mix_face = (1 - Y_face) * R_dry + Y_face * R_v;
				float Cp_mix_face = (1 - Y_face) * Cp_dry + Y_face * Cp_v;
				float GAMMA_face = Cp_mix_face / (Cp_mix_face - R_mix_face);

				CPU_Calc_rho_u_P_T(&interface_p[INDEX*6], &flux_Y[INDEX*6], //因為flux跟interface_p都有6個物理量需要儲存，如果不加這行的話數據就會一直不斷被覆蓋，最後變成只有儲存到最後一格的資料。
						   QB_rho_star, QB_ux_star, QB_vy_star, QB_vz, QB_cRT_star, QB_Y_star,
						   QT_rho_star, QT_ux_star, QT_vy_star, QT_vz, QT_cRT_star, QT_Y_star, R_mix_face, GAMMA_face,
						   0.0, 1.0, 0.0,
						   -1.0, 0.0, 0.0,
						   0.0, 0.0, 1.0, wall_flag);
				float rho_face_Y = 0.5 * (QC_rho + QT_rho);
				float Diff_flux_Y = rho_face_Y * D * ((QT_Y - QC_Y) / dy);
				flux_Y[INDEX*6 + 5] -= Diff_flux_Y;
			}
		}

		for (int i = 2; i < NX + 2; i++){		//The ghost cells on the left and right are not include in calculation.
			for (int j = 2; j < NY + 2; j++){
				int INDEX = i * (NY + 4) + j;
	    			// 我們是 i*6，所以左界面是 (i-1)*6，右界面是 i*6
				int L_interface = ((i - 1) * (NY + 4) + j) * 6;
				int R_interface = INDEX * 6; // T_intewrface = R_interface
				int B_interface = (i * (NY + 4) + (j - 1)) * 6;
				int T_interface = INDEX * 6; // T_intewrface = R_interface

				float R_mix_old = (1 - p5[INDEX]) * R_dry + p5[INDEX] * R_v;
				float Cp_mix_old = (1 - p5[INDEX]) * Cp_dry + p5[INDEX] * Cp_v;
				float Cv_mix_old = Cp_mix_old - R_mix_old;

				// Save old data
				float rho_old = p0[INDEX];
				float u_old = p1[INDEX];
				float v_old = p2[INDEX];
				float T_old = p3[INDEX];
				float Y_old = p5[INDEX];

				float MomX_old = rho_old * u_old; // p0[INDEX] * p1[INDEX]
				float MomY_old = rho_old * v_old; // p0[INDEX] * p1[INDEX]
				float E_old = rho_old * (Cv_mix_old * T_old + 0.5 * (u_old * u_old + v_old * v_old)); //p0[INDEX] * (CV * p2[INDEX] + 0.5 * (p1[INDEX] * p1[INDEX] + p2[INDEX] * p2[INDEX]))
				float rho_Y_old = rho_old * Y_old;
				// Use FVM to get new primitive variable，interface_p[0] is density、[1] is u、[2] is v、[3] is w、[4] is temperature、[5] is mass fraction (Y)。
				float rho_new = rho_old + (dt / dx) * (flux_X[L_interface + 0] - flux_X[R_interface + 0])
							+ (dt / dy) * (flux_Y[B_interface + 0] - flux_Y[T_interface + 0]);
				float MomX_new = MomX_old + (dt / dx) * (flux_X[L_interface + 1] - flux_X[R_interface + 1])
							  + (dt / dy) * (flux_Y[B_interface + 1] - flux_Y[T_interface + 1]);
				float MomY_new = MomY_old + (dt / dx) * (flux_X[L_interface + 2] - flux_X[R_interface + 2])
							  + (dt / dy) * (flux_Y[B_interface + 2] - flux_Y[T_interface + 2]);
				float E_new = E_old + (dt / dx) * (flux_X[L_interface + 4] - flux_X[R_interface + 4])
						    + (dt / dy) * (flux_Y[B_interface + 4] - flux_Y[T_interface + 4]);
				float rho_Y_new = rho_Y_old + (dt / dx) * (flux_X[L_interface + 5] - flux_X[R_interface + 5])
							    + (dt / dy) * (flux_Y[B_interface + 5] - flux_Y[T_interface + 5]);

		    		p0[INDEX] = rho_new;
    				p1[INDEX] = MomX_new / rho_new;
    				p2[INDEX] = MomY_new / rho_new;
				p5[INDEX] = rho_Y_new / rho_new;
				float Y_new = p5[INDEX];
				float R_mix_new = (1 - Y_new) * R_dry + Y_new * R_v;
				float Cp_mix_new = (1 - Y_new) * Cp_dry + Y_new * Cp_v;
				float Cv_mix_new = Cp_mix_new - R_mix_new;
				float internal_e = (E_new / rho_new) - 0.5 * (p1[INDEX] * p1[INDEX] + p2[INDEX] * p2[INDEX]);
				p3[INDEX] = internal_e / Cv_mix_new;
				p4[INDEX] = p0[INDEX] * R_mix_new * p3[INDEX];
				if (isnan(p3[INDEX]) || isnan(p0[INDEX])){
					printf("抓到 NaN！時間 t=%f, 座標 i=%d, j=%d\n", t, i, j); exit(1);
				}
				float T_C = p3[INDEX] - 273.15;
				float P_sat = 610.78 * exp((17.27 * T_C) / (T_C + 237.3));
				float P_vapor = p0[INDEX] * p5[INDEX] * R_v * p3[INDEX]; // p0*p5*R_v*p3 就是vapor分壓
				p6[INDEX] = P_vapor / P_sat;
			}
		}
		t += dt;
	}
	FILE * pFile = fopen("Results_of_1000000_cells.txt","w");
	for (int i = 2; i < NX + 2; i++){
		for (int j = 2; j < NY + 2; j++){
			int INDEX = i * (NY + 4) + j;
			float X = (i - 1.5) * dx;
			float Y = (j - 1.5) * dy;
			fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", X, Y, p0[INDEX], p1[INDEX], p2[INDEX], p3[INDEX], p4[INDEX], p5[INDEX]);
		}
	}
	fclose(pFile);

	Free_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &p5, &p6, &interface_p, &flux_X, &flux_Y);
return 0;
}

