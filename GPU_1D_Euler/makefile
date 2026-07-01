all: GPU CPU
	nvcc main.o memory.o -code=sm_89 -arch=compute_89 -o main.exe
GPU:
	nvcc memory.cu -c -code=sm_89 -arch=compute_89
CPU:
	g++ main.c -c
clean:
	rm *.o main.exe
