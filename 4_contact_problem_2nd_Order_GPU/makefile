all: GPU CPU
	nvcc main.o memory.o Initial.o Boundary.o Calc_rho_u_P_T.o -code=sm_89 -arch=compute_89 -o main.exe
GPU:
	nvcc memory.cu Boundary.cu Calc_rho_u_P_T.cu -c -code=sm_89 -arch=compute_89
CPU:
	nvcc main.cu Initial.cu -c
clean:
	rm *.o main.exe
