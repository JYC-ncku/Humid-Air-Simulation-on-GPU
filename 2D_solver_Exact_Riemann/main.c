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
	int NY = 5;
	int N_CELLS = (NX+2) * (NY+2); //+2 for Ghost cells
	float L = 1.0;
	float H = 0.005;
	float dx = L/NX;
	float dy = H/NY;
	float t = 0;
	float t_FINAL = 0.2;
	float R = 1.0;
	float GAMMA = 1.4;
	int wall_flag = 0;
	float *x, *y, *p0, *p1, *p2, *p3, *p4, *interface_p, *flux_X, *flux_Y; //p0 is density, p1 is x-dir velocity, p2 is y-dir veloctiy, p3 is temperature, p4 si pressure.
	float flxnmn, flxpmn, flxqmn;
	float CFL = 0.5;

	Allocate_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &interface_p, &flux_X, &flux_Y, N_CELLS);
	//Initial condition
	for ( int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			//int INDEX = i * (NY + 2) + j;
			int INDEX = i * (NY + 2) + j;
			if (i <= NX/2){
				p0[INDEX] = 10.0; //rho_L = 10
				p1[INDEX] = 0.0; //u_L = 0
				p2[INDEX] = 0.0; //v_L = 0
				p3[INDEX] = 1.0; // T_L = 1
			}else{
				p0[INDEX] = 1.0; //rho_R = 1
				p1[INDEX] = 0.0; //u_R = 0
				p2[INDEX] = 0.0; //v_R = 0
				p3[INDEX] = 1.0; //T_R = 0
			}
			p4[INDEX] = p0[INDEX] * R * p3[INDEX];
		}
	}
/* For x-dir
	// Because this code only consider 1D, so let other two direction equal 0)
	float QL_vy = 1.0, QL_vz = 0;
	float QR_vy = 1.0, QR_vz = 0;
	float nx = 1.0, ny = 0.0, nz = 0.0;
	float px = 0.0, py = 1.0, pz = 0.0;
	float qx = 0.0, qy = 0.0, qz = 1.0;
*/

/* For y-dir
	// Because this code only consider 1D, so let other two direction equal 0)
	float QL_vy = 1.0, QL_vz = 0;
	float QR_vy = 1.0, QR_vz = 0;
	float nx = 0.0, ny = 1.0, nz = 0.0;
	float px = -1.0, py = 0.0, pz = 0.0;
	float qx = 0.0, qy = 0.0, qz = 1.0;
*/
	while (t<t_FINAL){
		// Boundary condition for compute flux.
		Boundary(p0, p1, p2, p3, p4, NX, NY);

	    	float MAX_CFL = CPU_Compute_MAX_CFL(p0, p1, p2, p3, dx, dy, NX, NY);
		float dt = CFL / MAX_CFL;
	    	//X-dir (flux_X)
		for (int i = 0; i < NX + 1; i++){		//N cells have N+1 interface
			for (int j = 1; j < NY + 1; j++){
				int INDEX = i * (NY + 2) + j;
				int INDEX_R = (i + 1) * (NY + 2) + j;
				int INDEX_RR = (i + 2) * (NY + 2) + j;
				int INDEX_L = (i - 1) * (NY + 2) + j;
				float QL_rho = p0[INDEX_L];
				float QC_rho = p0[INDEX];
		    		float QR_rho = p0[INDEX_R];
		    		float QRR_rho = p0[INDEX_RR];
		    		float drho_dx_L = MINMOD(QL_rho, QC_rho, QR_rho, dx);
		    		float drho_dx_R = MINMOD(QC_rho, QR_rho, QRR_rho, dx);
		    		float QL_rho_star = QL_rho + 0.5 * dx * drho_dx_L;
		    		float QR_rho_star = QR_rho - 0.5 * dx * drho_dx_R;

    				float QL_ux = p1[INDEX_L];
    				float QC_ux = p1[INDEX];
    				float QR_ux = p1[INDEX_R];
    				float QRR_ux = p1[INDEX_RR];
		    		float du_dx_L = MINMOD(QL_ux, QC_ux, QR_ux, dx);
		    		float du_dx_R = MINMOD(QC_ux, QR_ux, QRR_ux, dx);
		    		float QL_ux_star = QL_rho + 0.5 * dx * du_dx_L;
		    		float QR_ux_star = QR_rho - 0.5 * dx * du_dx_R;

    				float QL_vy = p2[INDEX_L];
    				float QC_vy = p2[INDEX];
    				float QR_vy = p2[INDEX_R];
    				float QRR_vy = p2[INDEX_RR];
		    		float dv_dx_L = MINMOD(QL_vy, QC_vy, QR_vy, dx);
		    		float dv_dx_R = MINMOD(QC_vy, QR_vy, QRR_vy, dx);
		    		float QL_vy_star = QL_rho + 0.5 * dx * dv_dx_L;
		    		float QR_vy_star = QR_rho - 0.5 * dx * dv_dx_R;

    				float QL_vz  = 0.0;
    				float QR_vz  = 0.0;

		    		float QL_T = p3[INDEX_L];
		    		float QC_T = p3[INDEX];
    				float QR_T = p3[INDEX_R];
    				float QRR_T = p3[INDEX_RR];
		    		float dT_dx_L = MINMOD(QL_T, QC_T, QR_T, dx);
		    		float dT_dx_R = MINMOD(QC_T, QR_T, QRR_T, dx);
		    		float QL_T_star = QL_rho + 0.5 * dx * dT_dx_L;
		    		float QR_T_star = QR_rho - 0.5 * dx * dT_dx_R;

				float QL_cRT = sqrt(R * QL_T);
				float QC_cRT = sqrt(R * QC_T);
    				float QR_cRT = sqrt(R * QR_T);
    				float QRR_cRT = sqrt(R * QRR_T);
		    		float dcRT_dx_L = MINMOD(QL_cRT, QC_cRT, QR_cRT, dx);
		    		float dcRT_dx_R = MINMOD(QC_cRT, QR_cRT, QRR_cRT, dx);
		    		float QL_cRT_star = QL_rho + 0.5 * dx * dcRT_dx_L;
		    		float QR_cRT_star = QR_rho - 0.5 * dx * dcRT_dx_R;

				CPU_Calc_rho_u_P_T(&interface_p[INDEX*6], &flux_X[INDEX*5], //因為flux跟interface_p都有5個物理量需要儲存，如果不加這行的話數據就會一直不斷被覆蓋，最後變成只有儲存到最後一格的資料。
						   QL_rho_star, QL_ux_star, QL_vy_star, QL_vz, QL_cRT_star,
						   QR_rho_star, QR_ux_star, QR_vy_star, QR_vz, QR_cRT_star, R, GAMMA,
						   1.0, 0.0, 0.0,
						   0.0, 1.0, 0.0,
						   0.0, 0.0, 1.0, wall_flag);
			}
		}

		//Y-dir (flux_Y) 把X軸往逆時針轉90度看。
		for (int i = 1; i < NX + 1; i++){		//200 cells have 201 interface
			for (int j = 0; j < NY + 1; j++){
				int INDEX_B = i * (NY + 2) + (j - 1);
				int INDEX = i * (NY + 2) + j;
				int INDEX_T = i * (NY + 2) + (j + 1);
				int INDEX_TT = i * (NY + 2) + (j + 2);

				float QL_rho = p0[INDEX_B];
				float QC_rho = p0[INDEX];
		    		float QR_rho = p0[INDEX_T];
		    		float QRR_rho = p0[INDEX_TT];
		    		float drho_dy_L = MINMOD(QL_rho, QC_rho, QR_rho, dy);
		    		float drho_dy_R = MINMOD(QC_rho, QR_rho, QRR_rho, dy);
		    		float QL_rho_star = QL_rho + 0.5 * dy * drho_dy_L;
		    		float QR_rho_star = QR_rho - 0.5 * dy * drho_dy_R;

    				float QL_ux  = p1[INDEX_B];
    				float QC_ux  = p1[INDEX];
    				float QR_ux  = p1[INDEX_T];
    				float QRR_ux  = p1[INDEX_TT];
		    		float du_dy_L = MINMOD(QL_ux, QC_ux, QR_ux, dy);
		    		float du_dy_R = MINMOD(QC_ux, QR_ux, QRR_ux, dy);
		    		float QL_ux_star = QL_rho + 0.5 * dy * du_dy_L;
		    		float QR_ux_star = QR_rho - 0.5 * dy * du_dy_R;

    				float QL_vy = p2[INDEX_B];
    				float QC_vy = p2[INDEX];
    				float QR_vy = p2[INDEX_T];
    				float QRR_vy = p2[INDEX_TT];
		    		float dv_dy_L = MINMOD(QL_vy, QC_vy, QR_vy, dy);
		    		float dv_dy_R = MINMOD(QC_vy, QR_vy, QRR_vy, dy);
		    		float QL_vy_star = QL_rho + 0.5 * dy * dv_dy_L;
		    		float QR_vy_star = QR_rho - 0.5 * dy * dv_dy_R;

    				float QL_vz  = 0.0;
    				float QR_vz  = 0.0;

		    		float QL_T = p3[INDEX_B];
		    		float QC_T = p3[INDEX];
    				float QR_T = p3[INDEX_T];
    				float QRR_T = p3[INDEX_TT];
		    		float dT_dy_L = MINMOD(QL_T, QC_T, QR_T, dy);
		    		float dT_dy_R = MINMOD(QC_T, QR_T, QRR_T, dy);
		    		float QL_T_star = QL_rho + 0.5 * dy * dT_dy_L;
		    		float QR_T_star = QR_rho - 0.5 * dy * dT_dy_R;

				float QL_cRT = sqrt(R * QL_T);
				float QC_cRT = sqrt(R * QC_T);
    				float QR_cRT = sqrt(R * QR_T);
    				float QRR_cRT = sqrt(R * QRR_T);
		    		float dcRT_dy_L = MINMOD(QL_cRT, QC_cRT, QR_cRT, dy);
		    		float dcRT_dy_R = MINMOD(QC_cRT, QR_cRT, QRR_cRT, dy);
		    		float QL_cRT_star = QL_rho + 0.5 * dy * dcRT_dy_L;
		    		float QR_cRT_star = QR_rho - 0.5 * dy * dcRT_dy_R;

				CPU_Calc_rho_u_P_T(&interface_p[INDEX*6], &flux_Y[INDEX*5], //因為flux跟interface_p都有5個物理量需要儲存，如果不加這行的話數據就會一直不斷被覆蓋，最後變成只有儲存到最後一格的資料。
						   QL_rho_star, QL_ux_star, QL_vy_star, QL_vz, QL_cRT_star,
						   QR_rho_star, QR_ux_star, QR_vy_star, QR_vz, QR_cRT_star, R, GAMMA,
						   0.0, 1.0, 0.0,
						   -1.0, 0.0, 0.0,
						   0.0, 0.0, 1.0, wall_flag);
			}
		}

		for (int i = 1; i < NX + 1; i++){		//The ghost cells on the left and right are not include in calculation.
			for (int j = 1; j < NY + 1; j++){
				int INDEX = i * (NY + 2) + j;
	    			// 我們是 i*5，所以左界面是 (i-1)*5，右界面是 i*5
				int L_interface = ((i - 1) * (NY + 2) + j) * 5;
				int R_interface = INDEX * 5; // T_intewrface = R_interface
				int B_interface = (i * (NY + 2) + (j - 1)) * 5;
				int T_interface = INDEX * 5; // T_intewrface = R_interface
				float CV = R / (GAMMA - 1.0);
				// 先將舊的值儲存起來
				float rho_old = p0[INDEX];
				float u_old   = p1[INDEX];
				float v_old   = p2[INDEX];
				float T_old   = p3[INDEX];
				// p1 存的是速度 u，我們要先算動量 rho*u 的變化再去除以rho得到u。
				float MomX_old = rho_old * u_old; // p0[INDEX] * p1[INDEX]
				float MomY_old = rho_old * v_old; // p0[INDEX] * p1[INDEX]
				// 先從溫度算總能 E，更新完 E 再扣掉動能回算 T
				float E_old = rho_old * (CV * T_old + 0.5 * (u_old * u_old + v_old * v_old)); //p0[INDEX] * (CV * p2[INDEX] + 0.5 * (p1[INDEX] * p1[INDEX] + p2[INDEX] * p2[INDEX]))
				// 使用FVM計算新的值，interface_p[0]是密度、[1]是u、[2]是v、[3]是w、[4]是溫度。
				float rho_new = rho_old + (dt / dx) * (flux_X[L_interface + 0] - flux_X[R_interface + 0])
							+ (dt / dy) * (flux_Y[B_interface + 0] - flux_Y[T_interface + 0]);
				float MomX_new = MomX_old + (dt / dx) * (flux_X[L_interface + 1] - flux_X[R_interface + 1])
							  + (dt / dy) * (flux_Y[B_interface + 1] - flux_Y[T_interface + 1]);
				float MomY_new = MomY_old + (dt / dx) * (flux_X[L_interface + 2] - flux_X[R_interface + 2])
							  + (dt / dy) * (flux_Y[B_interface + 2] - flux_Y[T_interface + 2]);
				float E_new = E_old + (dt / dx) * (flux_X[L_interface + 4] - flux_X[R_interface + 4])
						    + (dt / dy) * (flux_Y[B_interface + 4] - flux_Y[T_interface + 4]);
				//更新密度 (p0)
		    		p0[INDEX] = rho_new;
    				// 更新動量並回推速度 (p1, p2)
    				p1[INDEX] = MomX_new / rho_new;
    				p2[INDEX] = MomY_new / rho_new;
		    		// 更新能量並回推溫度 (p3)
				float internal_e = (E_new / rho_new) - 0.5 * (p1[INDEX] * p1[INDEX] + p2[INDEX] * p2[INDEX]);
				p3[INDEX] = internal_e / CV;
				// 更新壓力 (p4)
				p4[INDEX] = p0[INDEX] * R * p3[INDEX];
			}
		}
		t += dt;
	}
	FILE * pFile = fopen("Results_of_5000_cells_x_dir.txt","w");
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY + 2) + j;
			float X = (i - 0.5) * dx;
			float Y = (j - 0.5) * dy;
			fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.2f\n", X, Y, p0[INDEX], p1[INDEX], p2[INDEX], p3[INDEX], p4[INDEX]);
		}
	}
	fclose(pFile);

	Free_memory(&x, &y, &p0, &p1, &p2, &p3, &p4, &interface_p, &flux_X, &flux_Y);
return 0;
}
