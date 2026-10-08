Compile this code using makefile
```
all: GPU CPU
	nvcc main.o memory.o Initial.o Calc_flux.o Primitive_variable.o compute_T.o Compute_Cv.o Boundary.o -code=sm_89 -arch=compute_89 -o main.exe
GPU:
	nvcc memory.cu Initial.cu Calc_flux.cu Primitive_variable.cu compute_T.cu Compute_Cv.cu Boundary.cu -dc -code=sm_89 -arch=compute_89
CPU:
	g++ main.c -c
clean:
	rm *.o
 ```


Final results for X-dir Veloctiy, Temperature, Mass fraction, Relative Humidity:

X-direction Velocity vs Location:
![X-dir_velocity_at_1s.png](./X-dir_velocity_at_1s.png)
![X-dir_velocity_at_2s.png](./X-dir_velocity_at_2s.png)
![X-dir_velocity_at_3s.png](./X-dir_velocity_at_3s.png)

Temperature vs Location:
![Temperature_at_1s.png](./Temperature_at_1s.png)
![Temperature_at_2s.png](./Temperature_at_2s.png)
![Temperature_at_3s.png](./Temperature_at_3s.png)

Mass fraction vs Location:
![Mass_Fraction_at_1s.png](./Mass_Fraction_at_1s.png)
![Mass_Fraction_at_2s.png](./Mass_Fraction_at_2s.png)
![Mass_Fraction_at_3s.png](./Mass_Fraction_at_3s.png)


Relative Humidity vs Location:
![Relative_Humidity_at_1s.png](./Relative_Humidity_at_1s.png)
![Relative_Humidity_at_2s.png](./Relative_Humidity_at_2s.png)
![Relative_Humidity_at_3s.png](./Relative_Humidity_at_3s.png)

As observed, the results at t = 1s, 2s, and 3s exhibit negligible variations, indicating that the simulation has already achieved a steady state by t = 1s.
