This program applies the MINMOD limiter to perform 2nd-order spatial reconstruction for the 1D humid-air shock tube problem,
and benchmarks the results against the 1st-order scheme.

Initial condition:
```
u = velocity = 0 m/s (everywhere), T = temperature = 300 K (everywhere)
Pressure = 5 bar   (x < 0.5L)
	   0.5 bar (x >= 0.5L)
L = 0.01 m
Computed time = 7e-6 s
```
compile this code using :
```
	gcc main.c memory.c Initial.c Boundary.c Calc_flux.c Primitive_variable.c compute_T.c Compute_Cv.c -O3 -o main.exe -lm
```
Final results for each number of cells (200, 400, 800).

Density of each cells:
![Result_of_Density__200_cells.png](./Result_of_Density__200_cells.png)
![Result_of_Density__400_cells.png](./Result_of_Density__400_cells.png)
![Result_of_Density__800_cells.png](./Result_of_Density__800_cells.png)

Velocity of each cells:
![Result_of_X-dir_velocity__200_cells.png](./Result_of_X-velocity__200_cells.png)
![Result_of_X-dir_velocity__400_cells.png](./Result_of_X-velocity__400_cells.png)
![Result_of_X-dir_velocity__800_cells.png](./Result_of_X-velocity__800_cells.png)

Temperature of each cells:
![Result_of_Temperature__200_cells.png](./Result_of_Temperature__200_cells.png)
![Result_of_Temperature__400_cells.png](./Result_of_Temperature__400_cells.png)
![Result_of_Temperature__800_cells.png](./Result_of_Temperature__800_cells.png)

Pressure of each cells:
![Result_of_Pressure__200_cells.png](./Result_of_Pressure__200_cells.png)
![Result_of_Pressure__400_cells.png](./Result_of_Pressure__400_cells.png)
![Result_of_Pressure__800_cells.png](./Result_of_Pressure__800_cells.png)

Mass Fraction of each cells:
![Result_of_Mass_Fraction__200_cells.png](./Result_of_Mass_Fraction__200_cells.png)
![Result_of_Mass_Fraction__400_cells.png](./Result_of_Mass_Fraction__400_cells.png)
![Result_of_Mass_Fraction__800_cells.png](./Result_of_Mass_Fraction__800_cells.png)

Relative Humidity of each cells:
![Result_of_Relative_Humidity__200_cells.png](./Result_of_Relative_Humidity__200_cells.png)
![Result_of_Relative_Humidity__400_cells.png](./Result_of_Relative_Humidity__400_cells.png)
![Result_of_Relative_Humidity__800_cells.png](./Result_of_Relative_Humidity__800_cells.png)
