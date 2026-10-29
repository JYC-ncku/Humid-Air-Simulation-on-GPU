The purpose of this program is to simulate a 1D Euler Equation problem using the Rusanov method and compare the final results with those obtained from the Exact Riemann Solver.

Initial condition:
```
u = velocity = 0 m/s (everywhere),   T = temperature = 1 K (everywhere)

Density = 10m  (x < 0.5L)
           1m  (x >= 0.5L)

Ratio of specific heats = 1.4, R = 1, L = 1m. Computed time = 0.2s.
```

compile this code using :
```
gcc Analytical_1D_Euler.c -O3 -lm -o Analytical_1D_Euler.exe
```

Final results for each number of cells(The red dashed lines from left to right represent x_Expansion_L、x_Expansion_R、x_contact、x_shock):

Density of each cell:
![Rusanov_vs_ExactRiemann_for_200_cells_of_Density.png](./Rusanov_vs_ExactRiemann_for_200_cells_of_Density.png)
![Rusanov_vs_ExactRiemann_for_400_cells_of_Density.png](./Rusanov_vs_ExactRiemann_for_400_cells_of_Density.png)
![Rusanov_vs_ExactRiemann_for_800_cells_of_Density.png](./Rusanov_vs_ExactRiemann_for_800_cells_of_Density.png)

Velocity of each cell:
![Rusanov_vs_ExactRiemann_for_200_cells_of_Velocity.png](./Rusanov_vs_ExactRiemann_for_200_cells_of_Velocity.png)
![Rusanov_vs_ExactRiemann_for_400_cells_of_Velocity.png](./Rusanov_vs_ExactRiemann_for_400_cells_of_Velocity.png)
![Rusanov_vs_ExactRiemann_for_800_cells_of_Velocity.png](./Rusanov_vs_ExactRiemann_for_800_cells_of_Velocity.png)

Temperature of each cell:
![Rusanov_vs_ExactRiemann_for_200_cells_of_Temperature.png](./Rusanov_vs_ExactRiemann_for_200_cells_of_Temperature.png)
![Rusanov_vs_ExactRiemann_for_400_cells_of_Temperature.png](./Rusanov_vs_ExactRiemann_for_400_cells_of_Temperature.png)
![Rusanov_vs_ExactRiemann_for_800_cells_of_Temperature.png](./Rusanov_vs_ExactRiemann_for_800_cells_of_Temperature.png)

Pressure of each cell:
![Rusanov_vs_ExactRiemann_for_200_cells_of_Pressure.png](./Rusanov_vs_ExactRiemann_for_200_cells_of_Pressure.png)
![Rusanov_vs_ExactRiemann_for_400_cells_of_Pressure.png](./Rusanov_vs_ExactRiemann_for_400_cells_of_Pressure.png)
![Rusanov_vs_ExactRiemann_for_800_cells_of_Pressure.png](./Rusanov_vs_ExactRiemann_for_800_cells_of_Pressure.png)
