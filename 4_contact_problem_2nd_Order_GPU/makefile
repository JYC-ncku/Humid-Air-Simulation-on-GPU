all: GPU CPU
	nvcc --fmad=false main.o memory.o Initial.o Boundary.o Calc_rho_u_P_T.o Calc_new_properties.o -code=sm_89 -arch=compute_89 -o main.exe
GPU:
	nvcc --fmad=false Initial.cu memory.cu Boundary.cu Calc_rho_u_P_T.cu Calc_new_properties.cu -c -code=sm_89 -arch=compute_89
CPU:
	nvcc --fmad=false main.cu -c
clean:
	rm *.o main.exe
