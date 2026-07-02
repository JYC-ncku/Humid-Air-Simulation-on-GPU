all: GPU CPU
	nvcc main.o memory.o Initial.o GPU_calc_flux.o Calc_conserved.o -code=sm_89 -arch=compute_89 -o main.exe
GPU:
	nvcc memory.cu GPU_calc_flux.cu Calc_conserved.cu -c -code=sm_89 -arch=compute_89
CPU:
	g++ main.c Initial.c -c
clean:
	rm *.o main.exe
