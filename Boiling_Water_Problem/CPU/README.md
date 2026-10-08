Compile this code using makefile
```
all:
	gcc main.c memory.c Initial.c Boundary.c Calc_flux.c Primitive_variable.c compute_T.c Compute_Cv.c -O3 -o main.exe -lm
```

Final results for X-dir Veloctiy, Temperature, Mass fraction, Relative Humidity:

X-direction Velocity vs Location:
![X-dir_velocity_at_1s.png](./X-dir_velocity_at_1s.png)

Temperature vs Location:
![Temperature_at_1s.png](./Temperature_at_1s.png)

Mass fraction vs Location:
![Mass_Fraction_at_1s.png](./Mass_Fraction_at_1s.png)

Relative Humidity vs Location:
![Relative_Humidity_at_1s.png](./Relative_Humidity_at_1s.png)
