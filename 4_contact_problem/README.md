This program is a Finite Volume Method (FVM) solver for 2D Euler equations. 

This program extends the 1D domain to a 2D domain, still solving it with an Exact-Riemann solver.
I used different initial conditions to validate my solver.
These results were obtained using the initial conditions from Configuration B, Fig. 9 of the paper 'NUMERICAL SOLUTION OF THE RIEMANN PROBLEM FOR TWO-DIMENSIONAL GAS DYNAMICS'.
The outcome matches the paper's results.

Initial condition:
```
Region A (Quadrant I):
0.5 <= x <= 1, 0.5 <= y <= 1.
Density = 1.0
u = X-direction velocity = 0.75
v = Y-direction veloctiy = -0.5
P = Pressure = 1.0

Region B (Quadrant II):
0 <= x < 0.5, 0.5 <= y <= 1.
Density = 2.0
u = X-direction velocity = 0.75
v = Y-direction veloctiy = 0.5
P = Pressure = 1.0

Region C (Quadrant III):
0 <= x < 0.5, 0 <= y < 0.5.
Density = 1.0
u = X-direction velocity = -0.75
v = Y-direction veloctiy = 0.5
P = Pressure = 1.0

Region D (Quadrant IV):
0.5 <= x <= 1, 0 <= y < 0.5.
Density = 3.0
u = X-direction velocity = -0.75
v = Y-direction veloctiy = -0.5
P = Pressure = 1.0

Ratio of specific heats = 1.4, R = 1, L = 1m. Computed time = 0.2s.
```

compile this code using :
```
COMPLIER := gcc
OPT_FLAGS := -O3 -lm
all:
	${COMPLIER} main.c memory.c Calc_rho_u_P_T.c Boundary.c Calc_Flux.c Calc_variable.c ${OPT_FLAGS} -o main.exe
```

Final results of density, X & Y-dir velocity, temperature, pressure:

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
