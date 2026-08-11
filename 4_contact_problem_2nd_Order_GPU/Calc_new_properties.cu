#include <stdlib.h>

__global__ void GPU_Calc_new_properties(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float R, float GAMMA, int NX, int NY, int N_CELLS){
	int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (int)INDEX / (NY+4);
	int j = (int)INDEX - i * (NY+4);
	if (i >= 2 && i < NX + 2){		//The ghost cells on the left and right are not include in calculation.
		if (j >= 2 && j < NY + 2){
    			// 我們是 i*5，所以左界面是 (i-1)*5，右界面是 i*5
			int L_interface = ((i - 1) * (NY + 4) + j) * 5;
			int R_interface = INDEX * 5; // T_intewrface = R_interface
			int B_interface = (i * (NY + 4) + (j - 1)) * 5;
			int T_interface = INDEX * 5; // T_intewrface = R_interface
			float CV = R / (GAMMA - 1.0);
			// 先將舊的值儲存起來
			float rho_old = d_p0[INDEX];
			float u_old   = d_p1[INDEX];
			float v_old   = d_p2[INDEX];
			float T_old   = d_p3[INDEX];
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
	    		d_p0[INDEX] = rho_new;
   			// 更新動量並回推速度 (p1, p2)
   			d_p1[INDEX] = MomX_new / rho_new;
   			d_p2[INDEX] = MomY_new / rho_new;
	    		// 更新能量並回推溫度 (p3)
			float internal_e = (E_new / rho_new) - 0.5 * (d_p1[INDEX] * d_p1[INDEX] + d_p2[INDEX] * d_p2[INDEX]);
			d_p3[INDEX] = internal_e / CV;
			// 更新壓力 (p4)
			d_p4[INDEX] = d_p0[INDEX] * R * d_p3[INDEX];
			}

	}
}

void Calc_new_properties(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float R, float GAMMA, int NX, int NY, int N_CELLS){
	int TPB = 128;
	int GPB = (TPB + N_CELLS - 1) / TPB;
	GPU_Calc_new_properties<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_p4, R, GAMMA, NX, NY, N_CELLS);
}
