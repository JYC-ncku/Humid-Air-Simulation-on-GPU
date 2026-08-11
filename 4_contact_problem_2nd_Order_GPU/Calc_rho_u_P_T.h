void Compute_MAX_CFL(float *CFL, float *d_p0, float *d_p1, float *d_p2, float *d_p3, float dx, float dy, int NX, int NY, int N_CELLS);

void Calc_flux_X(float *d_interface_p, float *d_flux_X, float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float R, float GAMMA, float dx, int NX, int NY, int N_CELLS);
void Calc_flux_Y(float *d_interface_p, float *d_flux_Y, float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_p4, float R, float GAMMA, float dy, int NX, int NY, int N_CELLS);
