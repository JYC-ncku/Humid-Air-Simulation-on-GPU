
# Riemann_Solver_1D_Euler
The flux computation part of this code was provided by Prof. Smith. Its purpose is to obtain the exact solution to the 1D Euler equation, allowing for cross-comparison with results from other solvers later>

More informations are in [here](./Peter_jacobs/Exact_Riemann_Solver_1D_Euler/README.md)

# Riemann_Solver_2D_Euler
This program upgrades Exact_Riemann_Solver_1D_Euler to a 2D domain and verifies the accuracy of the computation by modifying the initial conditions.

More informations are in [here](./Peter_Jacobs/Exact_Riemann_Solver_2D_Euler/README.md)

# Riemann_2nd_Solver_2D_Euler
This program employs second-order reconstruction to enhance the spatial accuracy of the Exact Riemann Solver results.

More informations are in [here](./Peter_Jacobs/Riemann_2nd_Solver_2D_Euler/README.md)

# Rusanov_1D_Euler
This program utilizes the Rusanov method to calculate flux and compares the results against the Exact Riemann solver.

More informations are in [here](./Rusanov/Rusanov_1D_Euler/README.md)

# Rusanov_1D_Euler_GPU
This program parallelizes Rusanov_1D_Euler using CUDA acceleration and benchmarks its performance against the single-CPU-core execution. To address the data racing issue inherent to this problem, a shared >

More informations are in [here](./Rusanov/Rusanov_1D_Euler_GPU/README.md)

# Rusanov_2D_Euler
This program upgrades the original 1D Euler code to a 2D domain and verifies the accuracy of the computation by modifying the initial conditions.

More informations are in [here](./Rusanov/Rusanov_2D_Euler/README.md)

# HLL_1D_Euler
This program solves the 1D Euler problem using the HLL method and compares the results with those obtained from the Rusanov method.

More informations are in [here](./HLL/HLL_1D_Euler/README.md)

# HLL_2D_Euler
This program upgrades the 1D HLL solver to a 2D domain and verifies the computational accuracy by swapping the X and Y directions for symmetry validation.

More informations are in [here](./HLL/HLL_2D_Euler/README.md)

# HLL_2D_2nd_Euler
This program reruns the 2D HLL Euler solver with second-order reconstruction using the Minmod limiter and compares the results with those obtained from the first-order scheme.

More informations are in [here](./HLL/HLL_2D_2nd_Euler/README.md)
