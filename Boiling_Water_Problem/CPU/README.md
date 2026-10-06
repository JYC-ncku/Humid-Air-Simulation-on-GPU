This program simulates a boiling water problem within a 1.0 m × 0.5 m rectangular domain.
A heated patch is designated along the bottom wall between x = 0.2 m and 0.7 m (T = 100 C, with the mass fraction set to its maximum value).
At the left boundary, a velocity profile with a peak of U_max = 5 m/s continuously injects flow into the domain, while the top and right boundaries are modeled as open boundaries.

In this program, the gas constant (R) and constant-volume specific heat capacity (Cv) are not treated as constants;
instead, they are evaluated based on a binary gas mixture of dry air and H2O.
Specifically, the values of Cv for dry air and water vapor are computed using the NASA polynomial formulations and their corresponding coefficient datasets.
(http://combustion.berkeley.edu/gri-mech/data/nasa_plnm.html)

Additionally, I used the Newton–Raphson method to iteratively obtain more precise temperature values.

The final results will present the full flow-field contours of x-direction velocity, temperature, mass fraction, and relative humidity at t = 1s of physical simulation time.

```
Initial condition:
u = X-direction velocity = 5 m/s (everywhere), v = Y-direction veloctiy = 0m/s (everywhere), T = temperature = 300.15 K (everywhere)
P = 1 atm (everywhere), phi = Mass fraction = 0 (everywhere)

```

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
