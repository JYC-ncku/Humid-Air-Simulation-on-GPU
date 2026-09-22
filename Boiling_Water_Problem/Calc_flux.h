float Max_Wave_Speed(float u_L, float u_R, float a_L, float a_R);

void Calc_HLL_X_flux(float rho_L, float rho_R, float u_L, float u_R, float v_L, float v_R, float T_L, float T_R, float P_L, float P_R, float Y_L, float Y_R, float e_L, float e_R, float a_L, float a_R,
		   float *mass_flux, float *momentum_X_flux, float *momenutm_Y_flux, float *energy_flux, float *mass_fraction_flux, int INDEX);

void Calc_HLL_Y_flux(float rho_B, float rho_T, float u_B, float u_T, float v_B, float v_T, float T_B, float T_T, float P_B, float P_T, float Y_B, float Y_T, float E_B, float E_T, float a_B, float a_T,
		     float *mass_flux, float *momentum_X_flux, float *momentum_Y_flux, float *energy_flux, float *mass_fraction_flux, int INDEX);
