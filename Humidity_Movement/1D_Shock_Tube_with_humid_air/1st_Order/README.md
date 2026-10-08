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
Final results for each number of cells.

Density of each cells:
![Result_of_density.png](./Result_of_density.png)

Velocity of each cells:
![Result_of_X-Velocity.png](./Result_of_X-Velocity.png)

Temperature of each cells:
![Result_of_Temperature.png](./Result_of_Temperature.png)

Pressure of each cells:
![Result_of_Pressure.png](./Result_of_Pressure.png)

Mass Fraction of each cells:
![Result_of_Mass_Fraction.png](./Result_of_Mass_Fraction.png)

Relative Humidity of each cells:
![Result_of_Relative_Humidity.png](./Result_of_Relative_Humidity.png)
