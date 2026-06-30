The purpose of this program is to simulate a 1D Euler Equation problem using the "HLL flux" method and compare the final results with those obtained from the Rusanov flux method.

Initial condition:
```
u = velocity = 0 m/s (everywhere),   T = temperature = 1 K (everywhere)
Density = 10m  (x < 0.5L)
	   1m  (x >= 0.5L)
Ratio of specific heats = 1.4, R = 1, L = 1m. Computed time = 0.2s.

```

compile this code using :
```
gcc main.c memory.c Initial.c Calc_flux.c -O3 -o main.exe -lm
```

Final results for each number of cells (200, 400, 800).

Density of each cells:
![HLL_vs_Rusanov_for_200_cells_of_Density.png](./HLL_vs_Rusanov_for_200_cells_of_Density.png)
![HLL_vs_Rusanov_for_400_cells_of_Density.png](./HLL_vs_Rusanov_for_400_cells_of_Density.png)
![HLL_vs_Rusanov_for_800_cells_of_Density.png](./HLL_vs_Rusanov_for_800_cells_of_Density.png)

Velocity of each cells:
![HLL_vs_Rusanov_for_200_cells_of_Velocity.png](./HLL_vs_Rusanov_for_200_cells_of_Velocity.png)
![HLL_vs_Rusanov_for_400_cells_of_Velocity.png](./HLL_vs_Rusanov_for_400_cells_of_Velocity.png)
![HLL_vs_Rusanov_for_800_cells_of_Velocity.png](./HLL_vs_Rusanov_for_800_cells_of_Velocity.png)

Temperature of each cells:
![HLL_vs_Rusanov_for_200_cells_of_Temperature.png](./HLL_vs_Rusanov_for_200_cells_of_Temperature.png)
![HLL_vs_Rusanov_for_400_cells_of_Temperature.png](./HLL_vs_Rusanov_for_400_cells_of_Temperature.png)
![HLL_vs_Rusanov_for_800_cells_of_Temperature.png](./HLL_vs_Rusanov_for_800_cells_of_Temperature.png)

Pressure of each cells:
![HLL_vs_Rusanov_for_200_cells_of_Pressure.png](./Rusanov_vs_ExactRiemann_for_200_cells_of_Pressure.png)
![HLL_vs_Rusanov_for_400_cells_of_Pressure.png](./Rusanov_vs_ExactRiemann_for_400_cells_of_Pressure.png)
![HLL_vs_Rusanovn_for_800_cells_of_Pressure.png](./Rusanov_vs_ExactRiemann_for_800_cells_of_Pressure.png)
