# Humid-Air-Simulation-on-GPU
This is my two-year plan, and I hope to complete it on time and graduate successfully.

## Allocate memory 
Using C code to allocate and free memory and test if it work!

## 1D Heat Transfer Simulation in Air

Information on the 1D Heat transfer simulation of air can be found [here](./1D_Heat_Transfer_Simulation_in_Air/README.md)

## Sine-wave
This is for practicing how to transfer calculation results from C code to MATLAB for plotting. 

The results can be viewed by clicking [here](./sine_wave/README.md).

## Perform_1D_humidity_simulation_in_Air
Information on the 1D Humidity simulation of air can be found [here](./Perform_1D_humidity_simulation_in_Air/READE.md)

## 1D_Shallow_Water
This program uses c language to simulation 1D Shallow water problem, and use Rusanov method to calculate flux.

More informations are in [here](./1D-Shallow-water/README.md)

## Exact_Riemann_Solver_1D_Euler
The flux computation part of this code was provided by Prof. Smith. Its purpose is to obtain the exact solution to the 1D Euler equation, allowing for cross-comparison with results from other solvers later on.

More informations are in [here](Exact_Riemann_Solver_1D_Euler)

## Exact_Riemann_Solver_2D_Euler
This program upgrades Exact_Riemann_Solver_1D_Euler to a 2D domain and verifies the accuracy of the computation by modifying the initial conditions.

More informations are in [here](Exact_Riemann_Solver_2D_Euler)

## 2D_solver_Exact_Riemann
This program employs second-order reconstruction to enhance the spatial accuracy of the Exact Riemann Solver results.

More informations are in [here](2D_solver_Exact_Riemann)

## Rusanov_1D_Euler
This program utilizes the Rusanov method to calculate flux and compares the results against the Exact Riemann solver.

More informations are in [here](Rusanov_1D_Euler)

## GPU_1D_Euler
This program parallelizes Rusanov_1D_Euler using CUDA acceleration and benchmarks its performance against the single-CPU-core execution. To address the data racing issue inherent to this problem, a shared memory implementation was also created to evaluate the speed difference between using and not using shared memory.

More informations are in [here](GPU_1D_Euler)

## Rusanov_2D_Euler
This program upgrades the original 1D Euler code to a 2D domain and verifies the accuracy of the computation by modifying the initial conditions.

More informations are in [here](Rusanov_2D_Euler)

## HLL_1D_Euler
This program solves the 1D Euler problem using the HLL method and compares the results with those obtained from the Rusanov method.

More informations are in [here](HLL_1D_Euler)

## HLL_2D_Euler
This program upgrades the 1D HLL solver to a 2D domain and verifies the computational accuracy by swapping the X and Y directions for symmetry validation.

More informations are in [here](HLL_2D_Euler)

## HLL_2D_2nd_Euler
This program reruns the 2D HLL Euler solver with second-order reconstruction using the Minmod limiter and compares the results with those obtained from the first-order scheme.

More informations are in [here](HLL_2D_2nd_Euler)

## 4_contact_problem
This program solves the classic 2D Riemann problem using the Exact Riemann Solver.

More informations are in [here](4_contact_problem)

## 4_contact_problem_2nd_Order
This program resolves the 4-contact problem using second-order reconstruction alongside the Exact Riemann Solver.

More informations are in [here](4_contact-problem_2nd_Order)

## 4_contact_problem_2nd_Order_GPU
This program applies CUDA acceleration to parallelize the 2nd-order 4-contact problem solver.

More informations are in [here](4_contact_problem_2nd_Order_GPU)


