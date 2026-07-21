float MAX_WAVE_SPEED(float u_L, float u_R, float a_L, float a_R);

void Calc_flux(float rho_L, float rho_R, float u_L, float u_R, float v_L, float v_R, float T_L, float T_R, float P_L, float P_R, float e_L, float e_R, float W_LOCAL_MAX,
	       float *mass_flux, float *momentum_X_flux, float *momentum_Y_flux, float *energy_flux, int INDEX);
