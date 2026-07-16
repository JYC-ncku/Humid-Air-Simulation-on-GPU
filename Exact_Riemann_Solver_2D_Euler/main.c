#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memory.h"
#include "Calc_rho_u_P_T.h"

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
	int NX = 300;
	int NY = 100;
	int N_CELLS = (NX+2) * (NY+2); //+2 is for Ghost cells
	float L = 3.0;
	float H = 1.0;
	float dx = L/NX;
	float dy = H/NY;
	float dt = 0.001; //隨便設，如果CFL有error就調整dt值。
	float t = 0;
	float t_FINAL = 0.2;
	float R = 1.0;
	float GAMMA = 1.4;
	int wall_flag = 0;
	float *x, *p0, *p1, *p2, *p3, *p4, *interface_p, *flux_X, *flux_Y; //p0 is density, p1 is x-dir velocity, p2 is y-dir veloctiy, p3 is temperature, p4 si pressure.
	float flxnmn, flxpmn, flxqmn;

	Allocate_memory(&x, &p0, &p1, &p2, &p3, &p4, &interface_p, &flux_X, &flux_Y, N_CELLS);
	//Initial condition
	for ( int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY + 2) + j;
			x[INDEX] = (i - 0.5) * dx;
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
		// LEFT and RIGHT
		for (int j = 1 ; j <= NY; j++){
			int LEFT_GHOST = 0 * (NY+2) + j;
			int RIGHT_GHOST = (NX+1) * (NY+2) + j;
			int LEFT_INNER = 1 * (NY+2) + j;
			int RIGHT_INNER = NX * (NY+2) + j;
			p0[LEFT_GHOST] = p0[LEFT_INNER];
			p0[RIGHT_GHOST] = p0[RIGHT_INNER];
			p1[LEFT_GHOST] = p1[LEFT_INNER];
			p1[RIGHT_GHOST] = p1[RIGHT_INNER];
			p2[LEFT_GHOST] = p2[LEFT_INNER];
			p2[RIGHT_GHOST] = p2[RIGHT_INNER];
			p3[LEFT_GHOST] = p3[LEFT_INNER];
			p3[RIGHT_GHOST] = p3[RIGHT_INNER];
			p4[LEFT_GHOST] = p4[LEFT_INNER];
			p4[RIGHT_GHOST] = p4[RIGHT_INNER];
		}
		//BOTTOM and TOP
		for (int i = 1 ; i <= NX; i++){
			int BOTTOM_GHOST = i * (NY+2) + 0;
			int TOP_GHOST = i * (NY+2) + (NY+1);
			int BOTTOM_INNER = i * (NY+2) + 1;
			int TOP_INNER = i * (NY+2) + NY;
			p0[BOTTOM_GHOST] = p0[BOTTOM_INNER];
			p0[TOP_GHOST] = p0[TOP_INNER];
			p1[BOTTOM_GHOST] = p1[BOTTOM_INNER];
			p1[TOP_GHOST] = p1[TOP_INNER];
			p2[BOTTOM_GHOST] = p2[BOTTOM_INNER];
			p2[TOP_GHOST] = p2[TOP_INNER];
			p3[BOTTOM_GHOST] = p3[BOTTOM_INNER];
			p3[TOP_GHOST] = p3[TOP_INNER];
			p4[BOTTOM_GHOST] = p4[BOTTOM_INNER];
			p4[TOP_GHOST] = p4[TOP_INNER];
			// Reflect boundary
			p2[BOTTOM_GHOST] = -p2[BOTTOM_INNER];
			p2[TOP_GHOST] = -p2[TOP_INNER];
		}

	    	float MAX_CFL = CPU_Compute_MAX_CFL(p0, p1, p2, p3, dx, dy, dt, NX, NY);
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
	FILE * pFile = fopen("Results_of_30000_cells.txt","w");
	for (int i = 1; i < NX + 1; i++){
		for (int j = 1; j < NY + 1; j++){
			int INDEX = i * (NY + 2) + j;
			fprintf(pFile, "%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.2f\n", x[INDEX], p0[INDEX], p1[INDEX], p2[INDEX], p3[INDEX], t);
		}
	}
	fclose(pFile);

	Free_memory(&x, &p0, &p1, &p2, &p3, &p4, &interface_p, &flux_X, &flux_Y);
return 0;
}
