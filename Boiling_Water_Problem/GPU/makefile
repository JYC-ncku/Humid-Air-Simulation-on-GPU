all: GPU CPU
	nvcc main.o memory.o Initial.o Calc_flux.o Primitive_variable.o compute_T.o Compute_Cv.o Boundary.o -code=sm_89 -arch=compute_89 -o main.exe
GPU:
	nvcc memory.cu Initial.cu Calc_flux.cu Primitive_variable.cu compute_T.cu Compute_Cv.cu Boundary.cu -dc -code=sm_89 -arch=compute_89
CPU:
	g++ main.c -c
clean:
	rm *.o main.exe
