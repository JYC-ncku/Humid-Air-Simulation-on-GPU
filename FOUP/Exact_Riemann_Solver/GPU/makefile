all: GPU CPU
	nvcc main.o memory.o Initial.o Compute_Cv.o -code=sm_89 -arch=compute_89 -o main.exe
GPU:
	nvcc memory.cu Initial.cu Compute_Cv.cu -dc -code=sm_89 -arch=compute_89
CPU:
	g++ main.c -c
clean:
	rm *.o
