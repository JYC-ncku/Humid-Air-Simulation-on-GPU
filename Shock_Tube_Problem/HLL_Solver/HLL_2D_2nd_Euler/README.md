This program is a Finite Volume Method (FVM) solver for 2D Euler equations. 

This program applies the Minmod limiter to my 2D HLL Euler problem and compares the results with the first-order solution.

Initial condition:
```
u = X-direction velocity = 0 m/s (everywhere), v = Y-direction veloctiy = 0m/s (everywhere), T = temperature = 1 K (everywhere)

Density = 10m  (x <= 0.5L)
           1m  (x > 0.5L)
```

compile this code using :
```
COMPLIER := gcc
OPT_FLAGS := -O3 -lm
all:
	${COMPLIER} main.c memory.c Calc_rho_u_P_T.c Boundary.c ${OPT_FLAGS} -o main.exe
```

Final results:

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
