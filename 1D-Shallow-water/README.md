This code uses Rusanov Method to solve 1D Shallow water problem.

Initial conditions:
```
Left (x < 0.5m):

height = 1m

water speed = 0 m/s.

Right (x >= 0.5m):

height = 0.1 m

water speed = 0 m/s

Simulation time = 0.5s
```
compile this code using :
```
gcc 1D-Shallow-water.c -o 1D-Shallow-water.exe -lm
```
Use MATLAB to plot the results of 4 different cells :
```
yyaxis left
plot(F1(:,1),F1(:,2),'b-', 'DisplayName','Cell=200 (Depth)')
hold on
plot(F2(:,1),F2(:,2),'b--', 'DisplayName','Cell=400 (Depth)')
plot(F3(:,1),F3(:,2),'b-.', 'DisplayName','Cell=800 (Depth)')
plot(F4(:,1),F4(:,2),'b',   'DisplayName','Cell=1600 (Depth)')
ylabel('Water Depth(m)')
yyaxis right
plot(F1(:,1),F1(:,3),'r-', 'DisplayName','Cell=200 (Veloctiy)')
plot(F2(:,1),F2(:,3),'r--', 'DisplayName','Cell=400 (Velocity)')
plot(F3(:,1),F3(:,3),'r-.', 'DisplayName','Cell=800 (Velocity)')
plot(F4(:,1),F4(:,3),'r',   'DisplayName','Cell=1600 (Velocity)')
ylabel('Water Velocity(m/s)')
xlabel('Location(m)')
title('200, 400, 800, 16000 Cells')
legend('Location', 'northeast');
```

![Results_of_4_different_cells.png](./Results_of_4_different_cells.png)