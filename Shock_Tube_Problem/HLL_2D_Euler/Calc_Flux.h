float MAX_WAVE_SPEED(float u_L, float u_R, float a_L, float a_R);

void Calc_flux_X(float rho_L, float rho_R, float u_L, float u_R, float v_L, float v_R, float T_L, float T_R, float P_L, float P_R, float e_L, float e_R, float a_L, float a_R,
		 float *mass_flux, float *momentum_X_flux, float *momentum_Y_flux, float *energy_flux, int INDEX);

void Calc_flux_Y(float rho_B, float rho_T, float u_B, float u_T, float v_B, float v_T, float T_B, float T_T, float P_B, float P_T, float e_B, float e_T, float a_B, float a_T,
		 float *mass_flux_Y, float *momentum_X_flux_Y, float *momentum_Y_flux_Y, float *energy_flux_Y, int INDEX);
