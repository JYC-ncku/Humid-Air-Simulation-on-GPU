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

int main(){
	int NX = 400;
	int NY = 400;
	int N_CELLS = (NX+2) * (NY+2); //+2 for Ghost cells
	float L = 1.0;
	float H = 1.0;
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
			//Quadrant I (area A)
			if (i >= NX/2 && j >= NY/2){
				p0[INDEX] = 1.0;
				p1[INDEX] = 0.75;
				p2[INDEX] = -0.5;
				p4[INDEX] = 1.0;
			//Quadrant II (area B)
			} else if (i < NX/2 && j >= NY/2){
				p0[INDEX] = 2.0;
				p1[INDEX] = 0.75;
				p2[INDEX] = 0.5;
				p4[INDEX] = 1.0;
			//Quadrant III (area C)
			} else if (i < NX/2 && j < NY/2){
				p0[INDEX] = 1.0;
				p1[INDEX] = 0.75;
				p2[INDEX] = 0.5;
				p4[INDEX] = 1.0;
			//Quadrant IV (area D)
			} else if (i >= NX/2 && j < NY/2){
				p0[INDEX] = 3.0;
				p1[INDEX] = -0.75;
				p2[INDEX] = -0.5;
				p4[INDEX] = 1.0;
			}
			p3[INDEX] = p4[INDEX] / (R * p0[INDEX]);
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
				float QL_rho = p0[INDEX];
		    		float QR_rho = p0[INDEX_R];
    				float QL_ux = p1[INDEX];
    				float QR_ux = p1[INDEX_R];
    				float QL_vy = p2[INDEX];
    				float QR_vy = p2[INDEX_R];
    				float QL_vz  = 0.0;
    				float QR_vz  = 0.0;
		    		float QL_T = p3[INDEX];
    				float QR_T = p3[INDEX_R];
				float QL_cRT = sqrt(R * QL_T);
    				float QR_cRT = sqrt(R * QR_T);
				CPU_Calc_rho_u_P_T(&interface_p[INDEX*6], &flux_X[INDEX*5], //因為flux跟interface_p都有5個物理量需要儲存，如果不加這行的話數據就會一直不斷被覆蓋，最後變成只有儲存到最後一格的資料。
						   QL_rho, QL_ux, QL_vy, QL_vz, QL_cRT,
						   QR_rho, QR_ux, QR_vy, QR_vz, QR_cRT, R, GAMMA,
						   1.0, 0.0, 0.0,
						   0.0, 1.0, 0.0,
						   0.0, 0.0, 1.0, wall_flag);
			}
		}

		//Y-dir (flux_Y) 把X軸往逆時針轉90度看。
		for (int i = 1; i < NX + 1; i++){		//200 cells have 201 interface
			for (int j = 0; j < NY + 1; j++){
				int INDEX = i * (NY + 2) + j;
				int INDEX_T = i * (NY + 2) + (j + 1);
				float QL_rho = p0[INDEX];
		    		float QR_rho = p0[INDEX_T];
    				float QL_ux  = p1[INDEX];
    				float QR_ux  = p1[INDEX_T];
    				float QL_vy = p2[INDEX];
    				float QR_vy = p2[INDEX_T];
    				float QL_vz  = 0.0;
    				float QR_vz  = 0.0;
		    		float QL_T   = p3[INDEX];
    				float QR_T   = p3[INDEX_T];
				float QL_cRT = sqrt(R * QL_T);
    				float QR_cRT = sqrt(R * QR_T);
				CPU_Calc_rho_u_P_T(&interface_p[INDEX*6], &flux_Y[INDEX*5], //因為flux跟interface_p都有5個物理量需要儲存，如果不加這行的話數據就會一直不斷被覆蓋，最後變成只有儲存到最後一格的資料。
						   QL_rho, QL_ux, QL_vy, QL_vz, QL_cRT,
						   QR_rho, QR_ux, QR_vy, QR_vz, QR_cRT, R, GAMMA,
						   0.0, 1.0, 0.0,
						   -1.0, 0.0, 0.0,
						   0.0, 0.0, 1.0, wall_flag);
			}
		}

		for (int i = 1; i < NX + 1; i++){		//200 cells, the ghost cells on the left and right are not include in calculation.
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
				float MomX_old = rho_old * u_old; // p0[INDEX] * p1[INDEX]
				float MomY_old = rho_old * v_old; // p0[INDEX] * p1[INDEX]
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
		    		p0[INDEX] = rho_new;
    				p1[INDEX] = MomX_new / rho_new;
    				p2[INDEX] = MomY_new / rho_new;
				float internal_e = (E_new / rho_new) - 0.5 * (p1[INDEX] * p1[INDEX] + p2[INDEX] * p2[INDEX]);
				p3[INDEX] = internal_e / CV;
				p4[INDEX] = p0[INDEX] * R * p3[INDEX];
			}
		}
		t += dt;
	}
	FILE * pFile = fopen("Results_of_160000_cells_x_dir.txt","w");
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
