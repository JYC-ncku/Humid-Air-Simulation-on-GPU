void Calc_primitive_variable(float *p0, float *p1, float *p2, float *p3, float *p4, float *p5, float *p6,
			     float *mass, float *momentum_X, float *momentum_Y, float *energy, float *mass_fraction,
			     float *flux_X, float *flux_Y,
			     float R_dry, float R_v, float dx, float dy, float dt, int NX, int NY, int N_CELLS);

float Compute_dt(float *MAX_Freq,  float *p0, float *p1, float *p2, float *p3, float *p5, float dx, float dy, float R_dry, float R_v, float CFL, int NX, int NY, int N_CELLS);
