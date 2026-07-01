all: GPU CPU
	nvcc main.o memory.o Initial.o -code=sm_89 -arch=compute_89 -o main.exe
GPU:
	nvcc memory.cu -c -code=sm_89 -arch=compute_89
CPU:
	g++ main.c Initial.c -c
clean:
	rm *.o main.exe
