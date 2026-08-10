#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memory.h"
#include "Calc_rho_u_P_T.h"
#include "Boundary.h"
#include "Initial.h"

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
	float R = 1.0;
	float GAMMA = 1.4;
//	int wall_flag = 0;
	float *x, *y,
	      *h_p0, *h_p1, *h_p2, *h_p3, *h_p4,
	      *d_p0, *d_p1, *d_p2, *d_p3, *d_p4,
	      *interface_p, *flux_X, *flux_Y; //p0 is density, p1 is x-dir velocity, p2 is y-dir veloctiy, p3 is temperature, p4 si pressure.
	float flxnmn, flxpmn, flxqmn;
	float CFL = 0.5;

	Allocate_memory(&x, &y,
			&h_p0, &h_p1, &h_p2, &h_p3, &h_p4,
			&d_p0, &d_p1, &d_p2, &d_p3, &d_p4,
			&interface_p, &flux_X, &flux_Y, N_CELLS);
	//Send the data to device
	Send_To_Device(d_p0, h_p0, N_CELLS);
	Send_To_Device(d_p1, h_p1, N_CELLS);
	Send_To_Device(d_p2, h_p2, N_CELLS);
	Send_To_Device(d_p3, h_p3, N_CELLS);
	Send_To_Device(d_p4, h_p4, N_CELLS);
	//Initial condition
	Initial(d_p0, d_p1, d_p2, d_p3, d_p4, R, NX, NY, N_CELLS);

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
		Boundary(d_p0, d_p1, d_p2, d_p3, d_p4, NX, NY);

	    	float MAX_CFL = CPU_Compute_MAX_CFL(p0, p1, p2, p3, dx, dy, NX, NY);
		float dt = CFL / MAX_CFL;
		Calc_flux_X(&interface_p[INDEX*6], &flux_Y[INDEX*5], d_p0, d_p1, dp_2, d_p3, d_p4, R, GAMMA, dx, NX, NY, N_CELLS);


		for (int i = 2; i < NX + 2; i++){		//The ghost cells on the left and right are not include in calculation.
			for (int j = 2; j < NY + 2; j++){
				int INDEX = i * (NY + 4) + j;
	    			// 我們是 i*5，所以左界面是 (i-1)*5，右界面是 i*5
				int L_interface = ((i - 1) * (NY + 4) + j) * 5;
				int R_interface = INDEX * 5; // T_intewrface = R_interface
				int B_interface = (i * (NY + 4) + (j - 1)) * 5;
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
	FILE * pFile = fopen("Results_of_1000000_cells.txt","w");
	for (int i = 2; i < NX + 2; i++){
		for (int j = 2; j < NY + 2; j++){
			int INDEX = i * (NY + 4) + j;
			float X = (i - 1.5) * dx;
			float Y = (j - 1.5) * dy;
			fprintf(pFile, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.2f\n", X, Y, h_p0[INDEX], h_p1[INDEX], h_p2[INDEX], h_p3[INDEX], h_p4[INDEX]);
		}
	}
	fclose(pFile);

	Free_memory(&x, &y, &h_p0, &h_p1, &h_p2, &h_p3, &h_p4, &d_p0, &d_p1, &d_p2, &d_p3, &d_p4, &interface_p, &flux_X, &flux_Y);
return 0;
}

