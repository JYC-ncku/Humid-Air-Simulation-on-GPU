This program is a Finite Volume Method (FVM) solver for 1D Euler equations. 

It utilizes an Approximate Riemann Solver to determine numerical fluxes between grid cells, aimed at simulating the classic physics of the Sod Shock Tube problem.

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
Use MATLAB to plot the results of 4 different cells :
```
Use the 200 cells result as the code example:

F1 = load('Results_of_200_cells.txt')
plot(F1(:,1),F1(:,2))
xlabel('Location(m)')
ylabel('Density(kg/m^3)')
title('Density vs Location for 200 cells')
xline(x0, '-k')
xline(x1, '--r')
xline(x2, '--r')
xline(x3, '--r')
xline(x4, '--r')
print('Results_of_200_cells_for_Density_vs_Location', '-dpng', '-r1000')
plot(F1(:,1),F1(:,3))
xlabel('Location(m)')
ylabel('Velocity(m/s)')
title('Velocity vs Location for 200 cells')
xline(x0, '-k')
xline(x1, '--r')
xline(x2, '--r')
xline(x3, '--r')
xline(x4, '--r')
print('Results_of_200_cells_for_Velocity_vs_Location', '-dpng', '-r1000')
plot(F1(:,1),F1(:,4))
xlabel('Location(m)')
ylabel('Temperature(K)')
title('Temperature vs Location for 200 cells')
xline(x0, '-k')
xline(x1, '--r')
xline(x2, '--r')
xline(x3, '--r')
xline(x4, '--r')
print('Results_of_200_cells_for_Temperature_vs_Location', '-dpng', '-r1000')
plot(F1(:,1),F1(:,5))
xlabel('Location(m)')
ylabel('Pressure(Pa)')
title('Pressure vs Location for 200 cells')
xline(x0, '-k')
xline(x1, '--r')
xline(x2, '--r')
xline(x3, '--r')
xline(x4, '--r')
print('Results_of_200_cells_for_Pressure_vs_Location', '-dpng', '-r1000')
```
Final results for each number of cells:

Density of each cell:
![Results_of_200_cells_for_Density_vs_Location.png](./Results_of_200_cells_for_Density_vs_Location.png)
![Results_of_400_cells_for_Density_vs_Location.png](./Results_of_400_cells_for_Density_vs_Location.png)
![Results_of_800_cells_for_Density_vs_Location.png](./Results_of_800_cells_for_Density_vs_Location.png)
![Results_of_4000_cells_for_Density_vs_Location.png](./Results_of_4000_cells_for_Density_vs_Location.png)

Velocity of each cell:
![Results_of_200_cells_for_Velocity_vs_Location.png](./Results_of_200_cells_for_Velocity_vs_Location.png)
![Results_of_400_cells_for_Velocity_vs_Location.png](./Results_of_400_cells_for_Velocity_vs_Location.png)
![Results_of_800_cells_for_Velocity_vs_Location.png](./Results_of_800_cells_for_Velocity_vs_Location.png)
![Results_of_4000_cells_for_Velocity_vs_Location.png](./Results_of_4000_cells_for_Velocity_vs_Location.png)

Temperature of each cell:
![Results_of_200_cells_for_Temperature_vs_Location.png](./Results_of_200_cells_for_Temperature_vs_Location.png)
![Results_of_400_cells_for_Temperature_vs_Location.png](./Results_of_400_cells_for_Temperature_vs_Location.png)
![Results_of_800_cells_for_Temperature_vs_Location.png](./Results_of_800_cells_for_Temperature_vs_Location.png)
![Results_of_4000_cells_for_Temperature_vs_Location.png](./Results_of_4000_cells_for_Temperature_vs_Location.png)

Pressure of each cell:
![Results_of_200_cells_for_Pressure_vs_Location.png](./Results_of_200_cells_for_Pressure_vs_Location.png)
![Results_of_400_cells_for_Pressure_vs_Location.png](./Results_of_400_cells_for_Pressure_vs_Location.png)
![Results_of_800_cells_for_Pressure_vs_Location.png](./Results_of_800_cells_for_Pressure_vs_Location.png)
![Results_of_4000_cells_for_Pressure_vs_Location.png](./Results_of_4000_cells_for_Pressure_vs_Location.png)