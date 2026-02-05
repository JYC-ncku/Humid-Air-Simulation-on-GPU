This program also simulates the 1D shock tube problem using the Rusanov method with identical initial and boundary conditions. 

The key difference is the inclusion of humidity effects, allowing for a comparison between the results with and without humidity.

Initial condition:
```
u = velocity = 0 m/s (everywhere),   T = temperature = 1 K (everywhere)

Density = 10m  (x < 0.5L)
           1m  (x >= 0.5L)

L = 1m. Computed time = 0.2s.
```

compile this code using :
```
gcc Humidity_movement_1D_shock_tube.c -O3 -lm -o Humidity_movement_1D_shock_tube.exe
```

Final results for each number of cells(The red dashed lines from left to right represent x_Expansion_L、x_Expansion_R、x_contact、x_shock):

Density of each cell:
![Humidity_movement_1D_Euler_200_cells_of_Density.png](./Humidity_movement_1D_Euler_200_cells_of_Density.png)
![Humidity_movement_1D_Euler_400_cells_of_Density.png](./Humidity_movement_1D_Euler_400_cells_of_Density.png)
![Humidity_movement_1D_Euler_800_cells_of_Density.png](./Humidity_movement_1D_Euler_800_cells_of_Density.png)

Velocity of each cell:
![Humidity_movement_1D_Euler_200_cells_of_Velocity.png](./Humidity_movement_1D_Euler_200_cells_of_Velocity.png)
![Humidity_movement_1D_Euler_400_cells_of_Velocity.png](./Humidity_movement_1D_Euler_400_cells_of_Velocity.png)
![Humidity_movement_1D_Euler_800_cells_of_Velocity.png](./Humidity_movement_1D_Euler_800_cells_of_Velocity.png)

Temperature of each cell:
![Humidity_movement_1D_Euler_200_cells_of_Temperature.png](./Humidity_movement_1D_Euler_200_cells_of_Temperature.png)
![Humidity_movement_1D_Euler_400_cells_of_Temperature.png](./Humidity_movement_1D_Euler_400_cells_of_Temperature.png)
![Humidity_movement_1D_Euler_800_cells_of_Temperature.png](./Humidity_movement_1D_Euler_800_cells_of_Temperature.png)

Pressure of each cell:
![Humidity_movement_1D_Euler_200_cells_of_Pressure.png](./Humidity_movement_1D_Euler_200_cells_of_Pressure.png)
![Humidity_movement_1D_Euler_400_cells_of_Pressure.png](./Humidity_movement_1D_Euler_400_cells_of_Pressure.png)
![Humidity_movement_1D_Euler_800_cells_of_Pressure.png](./Humidity_movement_1D_Euler_800_cells_of_Pressure.png)

Humidity of each cell:
![Humidity_movement_1D_Euler_200_cells_of_Humidity.png](./Humidity_movement_1D_Euler_200_cells_of_Humidity.png)
![Humidity_movement_1D_Euler_400_cells_of_Humidity.png](./Humidity_movement_1D_Euler_400_cells_of_Humidity.png)
![Humidity_movement_1D_Euler_800_cells_of_Humidity.png](./Humidity_movement_1D_Euler_800_cells_of_Humidity.png)