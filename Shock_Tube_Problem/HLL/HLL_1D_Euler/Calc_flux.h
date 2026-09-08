float Max_Wave_Speed(float u_L, float u_R, float a_L, float a_R);

void Calc_HLL_flux(float rho_L, float rho_R, float u_L, float u_R, float T_L, float T_R, float P_L, float P_R, float e_L, float e_R, float a_L, float a_R,
		   float *mass_flux, float *momentum_flux, float *energy_flux, int i);
