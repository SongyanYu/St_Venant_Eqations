#---
# model the confluence of two rivers joining into one downstream river.
# author: Songyan Yu
# date created: 11/11/2024
#---

import numpy as np
import matplotlib.pyplot as plt

# Genenral constants
g = 9.81  # Gravitational acceleration (m/s^2)
n_manning = 0.03  # Manning's roughness coefficient

# Parameters and discretization
# Define parameters for each river and downstream section
dx = 100.0  # Spatial step (m)
dt = 1.0    # Time step (s)
nx_river1 = 50  # Number of grids
nx_river2 = 50
nx_out = 100
B_river1 = 10  # river width (m)
B_river2 = 8
B_out = 15
S0 = 0    # channel bed slope

# Initialize variables
A_river1 = np.ones(nx_river1) * 20.0  # Cross-sectional area (m²)
Q_river1 = np.ones(nx_river1) * 5.0   # Discharge (m³/s)
A_river2 = np.ones(nx_river2) * 15.0
Q_river2 = np.ones(nx_river2) * 4.0
A_out = np.ones(nx_out) * 30.0
Q_out = np.zeros(nx_out)


def fric_slope (Q, A, h):
    return n_manning**2*Q*abs(Q) / (A**2*h**(4/3))

# Functions to update each river based on the Saint-Venant equations
def update_river(A, Q, h, dx, dt, nx, B):
    # Constants
    #L = dx*nx # Length of the channel (m)

    # Boundary conditions
 #   h[0] = 3  # fixed value for now, but can obtain thru loading input
 #   Q[0] = 8
 #   A[0] = h[0]*B

    h_new = np.copy(h)
    Q_new = np.copy(Q)
    A_new = np.copy(A)

    for i in range(1, nx-1):

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

    return(A, Q)

# Main time-stepping loop

fig, (ax1,ax2,ax3) = plt.subplots(1,3)

for t in range(3000):  # Run for 300 time steps or more
    # Update each river section
    h_river1 = A_river1 / B_river1 
    h_river2 = A_river2 / B_river2
    h_out = A_out / B_out 

    # Update both rivers separately
    h_river1[0] = 2 - t*0.0005
    A_river1, Q_river1 = update_river(A_river1, Q_river1, h_river1, dx, dt, nx_river1, B_river1)
    if t <=2:
        print(f'A_river1 : {A_river1}')

    h_river2[0] = 15/8 - t*0.0005
    A_river2, Q_river2 = update_river(A_river2, Q_river2, h_river2, dx, dt, nx_river2, B_river2)

    # Apply confluence/junction conditions
    Q_out[0] = Q_river1[-1] + Q_river2[-1]
    h_junction = (h_river1[-1] + h_river2[-1]) / 2  # Average water level at junction
    A_out[0] = h_junction * B_out  # Calculate cross-sectional area based on combined width

    # Update downstream (outflow) section
    A_out, Q_out = update_river(A_out, Q_out, h_out, dx, dt, nx_out, B_out)

    if t % 300 == 0:
        ax1.plot(A_river1/B_river1, label=f'time = {t} s')
        ax2.plot(A_river2/B_river2, label=f'time = {t} s')
        ax3.plot(A_out/B_out, label = f'time = {t} s') 

ax1.legend()
ax2.legend()
ax1.set_xlabel('Position along the river channel')
ax1.set_ylabel('Water surface profile (m)')
ax2.set_xlabel('Position along the river channel')
ax2.set_ylabel('Water surface profile (m)')
#plt.legend()
plt.show()
