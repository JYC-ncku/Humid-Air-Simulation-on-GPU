import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('Results_of_1000000_cells.txt')
X_raw = data[:,0]
Y_raw = data[:,1]
rho_raw = data[:,2]

nx = 1000
ny = 1000
X = X_raw.reshape((nx,ny))
Y = Y_raw.reshape((nx,ny))
rho = rho_raw.reshape((nx,ny))

plt.figure(figsize=(12,5))

cp = plt.contourf(X, Y, rho, levels=100, cmap='jet')
plt.colorbar(cp, label='Density (rho)')

plt.title('Density with 1000000 cells (GPU)')
plt.xlabel('X')
plt.ylabel('Y')

plt.savefig('Results_of_density_GPU.png')
