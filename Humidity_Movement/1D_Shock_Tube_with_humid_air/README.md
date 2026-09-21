This program employs the HLL Riemann solver to solve the 1D shock tube problem, incorporating moisture calculations with empirical formulas for Cv (Constant-volume heat capacity).

The advection terms are discretized using a 1st-order scheme, while the diffusion terms are handled via a 2nd-order Central Difference.

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
![Result_of_density_with_200_cells.png](./Result_of_density_with_200_cells.png)
![Result_of_density_with_400_cells.png](./Result_of_density_with_400_cells.png)
![Result_of_density_with_800_cells.png](./Result_of_density_with_800_cells.png)

Velocity of each cells:
![Result_of_X-Velocity_with_200_cells.png](./Result_of_X-Velocity_with_200_cells.png)
![Result_of_X-Velocity_with_400_cells.png](./Result_of_X-Velocity_with_400_cells.png)
![Result_of_X-Velocity_with_800_cells.png](./Result_of_X-Velocity_with_800_cells.png)

Temperature of each cells:
![Result_of_Temperature_with_200_cells.png](./Result_of_Temperature_with_200_cells.png)
![Result_of_Temperature_with_400_cells.png](./Result_of_Temperature_with_400_cells.png)
![Result_of_Temperature_with_800_cells.png](./Result_of_Temperature_with_800_cells.png)

Pressure of each cells:
![Result_of_Pressure_with_200_cells.png](./Result_of_Pressure_with_200_cells.png)
![Result_of_Pressure_with_400_cells.png](./Result_of_Pressure_with_400_cells.png)
![Result_of_Pressure_with_800_cells.png](./Result_of_Pressure_with_800_cells.png)

Mass Fraction of each cells:
![Result_of_Mass_Fraction_with_200_cells.png](./Result_of_Mass_Fraction_with_200_cells.png)
![Result_of_Mass_Fraction_with_400_cells.png](./Result_of_Mass_Fraction_with_400_cells.png)
![Result_of_Mass_Fraction_with_800_cells.png](./Result_of_Mass_Fraction_with_800_cells.png)

Relative Humidity of each cells:
![Result_of_Relative_Humidity_with_200_cells.png](./Result_of_Relative_Humidity_with_200_cells.png)
![Result_of_Relative_Humidity_with_400_cells.png](./Result_of_Relative_Humidity_with_400_cells.png)
![Result_of_Relative_Humidity_with_800_cells.png](./Result_of_Relative_Humidity_with_800_cells.png)
