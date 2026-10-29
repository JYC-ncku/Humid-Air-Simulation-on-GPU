This program is a Finite Volume Method (FVM) solver for 2D Euler equations. 

This program extends the 1D domain to a 2D domain, still solving it with an exact Riemann solver.
To validate the results, I ran two test cases.
Theoretically, the two results must be identical, and the final simulated data matched our expectations perfectly.

Initial condition:
```
Case1:
u = X-direction velocity = 0 m/s (everywhere), v = Y-direction veloctiy = 0m/s (everywhere), T = temperature = 1 K (everywhere)

Density = 10m  (x <= 0.5L)
           1m  (x > 0.5L)

Case2:
u = X-direction velocity = 0 m/s (everywhere), v = Y-direction veloctiy = 0m/s (everywhere), T = temperature = 1 K (everywhere)

Density = 10m  (y <= 0.5L)
           1m  (y > 0.5L)

Ratio of specific heats = 1.4, R = 1, L = 1m. Computed time = 0.2s.
```

compile this code using :
```
COMPLIER := gcc
OPT_FLAGS := -O3 -lm
all:
	${COMPLIER} main.c memory.c Calc_rho_u_P_T.c Boundary.c ${OPT_FLAGS} -o main.exe
```

Final results for 2 cases:

Density vs location:
![Result_of_Density.png](./Result_of_Density.png)

X-direction Velocity vs location:
![Result_of_X-Dir_Velocity.png](./Result_of_X-Dir_Velocity.png)

Y-direction Velocity vs location:
![Result_of_Y-Dir_Velocity.png](./Result_of_Y-Dir_Velocity.png)

Temperature vs location:
![Result_of_Temperature.png](./Result_of_Temperature.png)

Pressure vs location:
![Result_of_Pressure.png](./Result_of_Pressure.png)
