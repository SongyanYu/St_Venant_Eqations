#---
# solve Saint Venant Equations with Lax-Friedichs method
# author: Songyan Yu
# date created: 22/10/2024
#---

import numpy as np
import matplotlib.pyplot as plt

# Constants
g = 9.81  # Gravitational acceleration (m/s^2)
n_manning = 0.03  # Manning's roughness coefficient
dx = 10  # Space step (m)
dt = 1    # Time step (s)
B = 1.5  # channel width (m)
L = 500 # Length of the channel (m)
T = 20  # Simulation time (s)
n_steps = int(T / dt)  # Number of time steps
n_grid = int(L / dx)   # Number of grid points
S0 = 0    # channel bed slope

h = np.ones(n_grid)
Q = np.zeros(n_grid)
A = h * B

# Boundary conditions
h_bc = np.ones(n_steps)
Q_bc = np.zeros(n_steps)
for i in range(n_steps):
    h_bc[i] = 5 - i*0.2
    Q_bc[i] = 3 - i*0.1

A_bc = h_bc*B

def fric_slope (Q, A, h):
    return n_manning**2*Q*abs(Q) / (A**2*h**(4/3))

for n in range(0, n_steps):

    h[0] = h_bc[n]
    Q[0] = Q_bc[n]
    A[0] = A_bc[n]

    h_new = np.copy(h)
    Q_new = np.copy(Q)
    A_new = np.copy(A)

    for i in range(1, n_grid-1):

        h_new[i] = (h[i+1] + h[i-1])/2 - (Q[i+1] - Q[i-1])/2/dx*dt/B  # Lax-Friedrichs method (central finite difference approximation)

#        print(f"h_new[{i}]: {h_new[i]}")
        
        #h[i] = h_new[i]
        A_new[i] = h_new[i] * B
        #A[i] = A_new[i]

        # momentum equation: dQ/dt + d(Q**2/A)/dx + g*A*dh/dx = g*A*(S0-Sf)
        Sf = fric_slope(Q[i], A[i], h[i])
        dQ_dx = (Q[i+1]**2/A[i+1] - Q[i-1]**2/A[i-1])/2/dx    
        dh_dx = (h[i+1] - h[i-1])/2/dx
        Q_new[i] = dt*(g*A[i]*(S0-Sf) - dQ_dx - g*A[i]*dh_dx) + (Q[i+1] + Q[i-1])/2    # Lax-Friedrichs method
        
        #Q[i] = Q_new[i]

    # Update for next time step
    h_new[-1] = h_new[-2]
    A_new[-1] = A_new[-2]
    Q_new[-1] = Q_new[-2]

    h = np.copy(h_new)
    Q = np.copy(Q_new)
    A = np.copy(A_new)
    
#    print(f"Sf for {n}: {Sf}")
#    print(f"h for {n}: {h_new}")
#    print(f"Q for {n}: {Q_new}")
#    print(f"A for {n}: {A_new}")

    # Plot at some time steps
    if n % 1 == 0:
        plt.plot(h, label=f'time={n*dt} s')

plt.xlabel('Position along the channel (m)')
plt.ylabel('Water depth (m)')
plt.legend()
plt.title('Water Surface Profile Over Time')
#plt.savefig('row_update.png')
plt.show()

