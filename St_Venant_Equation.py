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
L = 1000 # Length of the channel (m)
T = 100  # Simulation time (s)
n_steps = int(T / dt)  # Number of time steps
n_grid = int(L / dx)   # Number of grid points
S0 = 0    # channel bed slope

h = np.ones(n_grid)
Q = np.zeros(n_grid)
A = h * B

# Boundary conditions
h[0] = 5  # Upstream water depth (m)
Q[0] = 3   # Upstream discharge (m^3/s)
A[0] = h[0] *B

def fric_slope (Q, A, h):
    return n_manning**2*Q*abs(Q) / (A**2*h**(4/3))

for n in range(0, n_steps):

    h_new = np.copy(h)
    Q_new = np.copy(Q)
    A_new = np.copy(A)

    for i in range(0, n_grid-1):

        # Continuity equation: dA/dt + dQ/dx = 0
        if i == 0:
            h_new[i] = h[i] - dt*(Q[i+1] - Q[i])/dx/B   # forward finite-difference
#        else:
#            if i == (n_grid-1):
#                print('i==n_grid-1')
#                h_new[i] = h[i] - dt*(Q[i]-Q[i-1])/dx/B   # backward finite-difference
        else:
            h_new[i] = (h[i+1] + h[i-1])/2 - (Q[i+1] - Q[i-1])/2/dx*dt/B  # Lax-Friedrichs method (central finite difference approximation)
        
#        print(f"h_new[{i}]: {h_new[i]}")
        
        h[i] = h_new[i]
        A_new[i] = h_new[i] * B
        A[i] = A_new[i]

        # momentum equation: dQ/dt + d(Q**2/A)/dx + g*A*dh/dx = g*A*(S0-Sf)
        if i == 0:
            Sf = fric_slope(Q[i], A[i], h[i])
            dQ_dx = (Q[i+1]**2/A[i+1] - Q[i]**2/A[i])/dx
            dh_dx = (h[i+1]-h[i])/dx
            Q_new[i] = dt*(g*A[i]*(S0-Sf) - dQ_dx - g*A[i]*dh_dx) + Q[i]
        else:
#            if  i == n_grid-1:
#                Sf = fric_slope(Q[i], A[i], h[i])
#                dQ_dx = (Q[i]**2/A[i] - Q[i-1]**2/A[i-1])/dx
#                dh_dx = (h[i]-h[i-1])/dx
#                Q_new[i] = dt*(g*A[i]*(S0-Sf) - dQ_dx - g*A[i]*dh_dx) + Q[i]
#            else:
            Sf = fric_slope((Q[i+1]+Q[i-1])/2, (A[i+1]+A[i-1])/2, (h[i+1]+h[i-1])/2)
            dQ_dx = (Q[i+1]**2/A[i+1] - Q[i-1]**2/A[i-1])/2/dx    
            dh_dx = (h[i+1] - h[i-1])/2/dx
            Q_new[i] = dt*(g*(A[i+1]+A[i-1])/2*(S0-Sf) - dQ_dx - g*(A[i+1]+A[i-1])/2*dh_dx) + (Q[i+1] + Q[i-1])/2    # Lax-Friedrichs method
        
        Q[i] = Q_new[i]

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
    if n % 10 == 0:
        plt.plot(h, label=f'time={n*dt} s')

plt.xlabel('Position along the channel (m)')
plt.ylabel('Water depth (m)')
plt.legend()
plt.title('Water Surface Profile Over Time')
#plt.savefig('row_update.png')
plt.show()

