float CPU_Compute_MAX_CFL(float *p0, float *p1, float *p2, float *p3, float dx, float dy, float dt, int N_CELLS);

void CPU_Calc_rho_u_P_T(float *interface_p, float *flux,
			float QL_rho, float QL_ux, float QL_vy, float QL_vz, float QL_cRT,
			float QR_rho, float QR_ux, float QR_vy, float QR_vz, float QR_cRT, float R, float GAMMA,
			float nx, float ny, float nz,
			float px, float py, float pz,
			float qx, float qy, float qz, int wall_flag);
