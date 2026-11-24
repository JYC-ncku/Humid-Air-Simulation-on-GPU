# Exact Riemann Solver
The code here uses the Exact Riemann Solver to simulate the FOUP problem, comparing 1st-order and 2nd-order spatial Reconstruction.

For the 2nd-order, I use two different limiters — MC and MINMOD.

To better observe the flow field evolution, I have compiled the results into animations, which you can download to view.

## 1st Order

The code in this folder uses first-order spatial accuracy for the simulation and includes both CPU and GPU versions.

I first ran the CPU version to verify that the simulation converged without diverging, and then ported it to GPU for parallel acceleration.

This also allows me to check whether the GPU computation results match the CPU baseline and are free of bugs.

Results of CPU are in [here](./1st_Order/CPU/Results_of_10s_CPU)

Results of GPU are in [here](./1st_Order/GPU/Results_of_10s_384x620_Cells)

## 2nd Order MINMOD
The code here performs second-order spatial reconstruction using the MINMOD limiter.

Simultaneously, to prevent numerical dissipation, the time stepping is upgraded to second-order accuracy using the RK2 method (Second-order Runge-Kutta method).

In addition, I also generated results for both variable Cv and constant Cv to facilitate comparing their effects on the flow field.

Results of MINMOD：

Variable Cv are in [here](./2nd_Order_MINMOD/Variable_Cv/Results_of_75s)

Constant Cv are in [here](./2nd_Order_MINMOD/Constant_Cv)

## 2nd Order MC
The code here performs second-order spatial reconstruction using the MC limiter.

Results of MC：

Variable Cv are in [here](./2nd_Order_MC/Variable_Cv)

Constant Cv are in [here](./2nd_Order_MC/Constant_Cv)

