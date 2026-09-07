#include <stdlib.h>
#include <math.h>

__device__ float MAX_Wave_Speed(float u_L, float u_R, float a_L, float a_R){
    float W_L = fabs(u_L) + a_L;
    float W_R = fabs(u_R) + a_R;
    float W_LOCAL_MAX;
    if (W_L > W_R){
        W_LOCAL_MAX = W_L;
    }else {
        W_LOCAL_MAX = W_R;
    }
    return W_LOCAL_MAX;
}

__global__ void GPU_Calc_flux(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_mass_flux, float *d_momentum_flux, float *d_energy_flux, float *d_block_max,
			      float GAMMA, float R, int N_CELLS){
	__shared__ float W_GLOBAL_MAX[128];
	int tid = threadIdx.x;
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	float W_LOCAL_MAX = 0.0;
	if (i < N_CELLS){
		if (i >= 1 && i < N_CELLS){
			float rho_L = d_p0[i-1];
			float rho_R = d_p0[i];
			float u_L = d_p1[i-1];
			float u_R = d_p1[i];
			float T_L = d_p2[i-1];
			float T_R = d_p2[i];
			float P_L = d_p3[i-1];
			float P_R = d_p3[i];
			float e_L = 0.5 * rho_L * u_L * u_L + (P_L / (GAMMA - 1));
			float e_R = 0.5 * rho_R * u_R * u_R + (P_R / (GAMMA - 1));
			float a_L = sqrt(GAMMA * R * T_L);
			float a_R = sqrt(GAMMA * R * T_R);

			W_LOCAL_MAX = MAX_Wave_Speed(u_L, u_R, a_L, a_R);

			float mass_L = rho_L;
			float momentum_L = rho_L * u_L;
			float energy_L = e_L;
			float mass_R = rho_R;
			float momentum_R = rho_R * u_R;
			float energy_R = e_R;

			float mass_flux_L = rho_L * u_L;
			float mass_flux_R = rho_R * u_R;
			float momentum_flux_L = rho_L * u_L * u_L + P_L;
			float momentum_flux_R = rho_R * u_R * u_R + P_R;
			float energy_flux_L = (e_L + P_L) * u_L;
			float energy_flux_R = (e_R + P_R) * u_R;

			d_mass_flux[i] = 0.5 * (mass_flux_L + mass_flux_R) - 0.5 * W_LOCAL_MAX * (mass_R - mass_L);
			d_momentum_flux[i] = 0.5 * (momentum_flux_L + momentum_flux_R) - 0.5 * W_LOCAL_MAX * (momentum_R - momentum_L);
			d_energy_flux[i] = 0.5 * (energy_flux_L + energy_flux_R) - 0.5 * W_LOCAL_MAX * (energy_R - energy_L);
		}
		//把所有cells的W_LOCAL_MAX寫進share memory
		W_GLOBAL_MAX[tid] = W_LOCAL_MAX;
		__syncthreads();

		//Tree-Base Reduction (決出每個block最快的wave speed）
		for(int stride = blockDim.x / 2; stride > 0; stride = stride / 2) {
			if (tid < stride) {
				if (W_GLOBAL_MAX[tid + stride] > W_GLOBAL_MAX[tid]) {
					W_GLOBAL_MAX[tid] = W_GLOBAL_MAX[tid + stride];
				}
			}
			__syncthreads();
		}
		if (tid == 0){
			d_block_max[blockIdx.x] = W_GLOBAL_MAX[0];
		}
	}
}

void Calc_flux(float *d_p0, float *d_p1, float *d_p2, float *d_p3, float *d_mass_flux, float *d_momentum_flux, float *d_energy_flux, float *d_block_max,
	       float GAMMA, float R, int N_CELLS){
	int TPB = 128;
	int GPB = (N_CELLS + TPB - 1) / TPB;
	GPU_Calc_flux<<<GPB, TPB>>>(d_p0, d_p1, d_p2, d_p3, d_mass_flux, d_momentum_flux, d_energy_flux, d_block_max, GAMMA, R, N_CELLS);
}
