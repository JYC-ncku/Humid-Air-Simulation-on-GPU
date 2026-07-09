#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memory.h"
#include "Calc_rho_u_P_T.h"

int main(){
	int NX = 300;
	int NY = 100;
	int N_CELLS = NX * NY;
	float L = 3.0;
	float H = 1.0;
	float dx = L/NX;
	float dy = H/NY;
	float dt = 0.0001; //隨便設，如果CFL有error就調整dt值。
	float t = 0;
	float t_FINAL = 0.2;
	float R = 1.0;
	float GAMMA = 1.4;
	int wall_flag = 0;
	float *x, *p0, *p1, *p2, *p3, *p4, *interface_p, *flux; //p0 is density, p1 is x-dir velocity, p2 is y-dir veloctiy, p3 is temperature, p4 si pressure.
	float flxnmn, flxpmn, flxqmn;

	Allocate_memory(&x, &p0, &p1, &p2, &p3, &p4, &interface_p, &flux, N_CELLS);
	//Initial condition
	for ( int i = 0; i<N_CELLS; i++){
		x[i] = (i+0.5) * dx;
		if (i<N_CELLS/2){
			p0[i] = 10.0; //rho_L = 10
			p1[i] = 0; //u_L = 0
			p2[i] = 0; //v_L = 0
			p3[i] = 1.0; // T_L = 1
		}else{
			p0[i] = 1.0; //rho_R = 1
			p1[i] = 0; //u_R = 0
			p2[i] = 0; //v_R = 0
			p3[i] = 1.0; //T_R = 0
		}
		p4[i] = p0[i] * R * p3[i];
	}
	// Because this code only consider 1D, so let other two direction equal 0)
	float QL_vy = 1.0, QL_vz = 0;
	float QR_vy = 1.0, QR_vz = 0;
	float nx = 1.0, ny = 0.0, nz = 0.0;
	float px = 0.0, py = 1.0, pz = 0.0;
	float qx = 0.0, qy = 0.0, qz = 1.0;

	while (t<t_FINAL){
    	float MAX_CFL = CPU_Compute_MAX_CFL(p0, p1, p2, p3, dx, dy, dt, N_CELLS);

		for (int i = 0; i < N_CELLS-1; i++){
			float QL_rho = p0[i];
	    		float QR_rho = p0[i+1];
    			float QL_ux  = p1[i];
    			float QR_ux  = p1[i+1];
	    		float QL_T   = p2[i];
    			float QR_T   = p2[i+1];
			float QL_cRT = sqrt(R * QL_T);
    			float QR_cRT = sqrt(R  * QR_T);
			CPU_Calc_rho_u_P_T(&interface_p[i*6], &flux[i*5], //因為flux跟interface_p都有5個物理量需要儲存，如果不加這行的話數據就會一直不斷被覆蓋，最後變成只有儲存到最後一格的資料。
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
			float CV = R / (GAMMA - 1.0);
			// 先將舊的值儲存起來
			float rho_old = p0[i];
			float u_old   = p1[i];
			float T_old   = p2[i];
			// p1 存的是速度 u，我們要先算動量 rho*u 的變化再去除以rho得到u。
			float Mom_old = rho_old * u_old; // p0[i] * p1[i]
			// 先從溫度算總能 E，更新完 E 再扣掉動能回算 T
			float E_old   = rho_old * (CV * T_old + 0.5 * u_old * u_old); //p0[i] * (CV * p2[i] + 0.5 * p1[i] * p1[i])
			// 使用FVM計算新的值，參考程式碼712～716行，interface_p[0]是密度、[1]是u、[2]是v、[3]是w、[4]是溫度。
			float rho_new = rho_old + (dt / dx) * (flux[L_interface + 0] - flux[R_interface + 0]);
			float Mom_new = Mom_old + (dt / dx) * (flux[L_interface + 1] - flux[R_interface + 1]);
			float E_new   = E_old   + (dt / dx) * (flux[L_interface + 4] - flux[R_interface + 4]);
			//更新密度 (p0)
	    		p0[i] = rho_new;
    			// 更新動量並回推速度 (p1)
    			p1[i] = Mom_new / rho_new;
	    		// 更新能量並回推溫度 (p2)
			float internal_e = (E_new / rho_new) - 0.5 * (p1[i] * p1[i]);
			p2[i] = internal_e / CV;
			// 更新壓力 (p3)
			float P_new = p0[i] * R * p2[i];
			p3[i] = P_new;
		}

		// Boundary condition for compute flux.
		// 左邊界
		p0[0] = p0[1];
		p1[0] = p1[1];
		p2[0] = p2[1];
		p3[0] = p3[1];

		// 右邊界
		p0[N_CELLS-1] = p0[N_CELLS-2];
		p1[N_CELLS-1] = p1[N_CELLS-2];
		p2[N_CELLS-1] = p2[N_CELLS-2];
		p3[N_CELLS-1] = p3[N_CELLS-2];

		t += dt;
	}
	FILE * pFile = fopen("Results_of_200_cells.txt","w");
	for (int i=0; i<N_CELLS; i++){
		fprintf(pFile, "%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.2f\n", x[i], p0[i], p1[i], p2[i], p3[i], t);
	}
	fclose(pFile);

	Free_memory(&x, &p0, &p1, &p2, &p3, &p4, &interface_p, &flux);
return 0;
}
