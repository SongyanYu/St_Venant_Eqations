#---
# model the confluence of two rivers into one downstream river.
# author: Songyan Yu
# date created: 11/11/2024
#---

import numpy as np

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

# Initialize variables
A_river1 = np.ones(nx_river1) * 20.0  # Cross-sectional area (m²)
Q_river1 = np.ones(nx_river1) * 5.0   # Discharge (m³/s)
A_river2 = np.ones(nx_river2) * 15.0
Q_river2 = np.ones(nx_river2) * 4.0
A_out = np.ones(nx_out) * 30.0
Q_out = np.zeros(nx_out)

# Functions to update each river based on the Saint-Venant equations
def update_river(A, Q, h, dx, dt, nx):
    # Constants
    B = A/h  # channel width (m)
    L = dx*nx # Length of the channel (m)
    S0 = 0    # channel bed slope

    # Boundary conditions
    h_bc = np.ones(n_steps)
    Q_bc = np.zeros(n_steps)
    for i in range(n_steps):
        h_bc[i] = 5 - i*0.2
        Q_bc[i] = 3 - i*0.1

    A_bc = h_bc*B

    def fric_slope (Q, A, h):
        return n_manning**2*Q*abs(Q) / (A**2*h**(4/3))


    # Similar to update function from previous example, adapted for each river
    pass  # Complete update function

# Main time-stepping loop
for t in range(300):  # Run for 300 time steps or more
    # Update each river section
    h_river1 = A_river1 / 10.0  # Assume width of 10m
    h_river2 = A_river2 / 8.0   # Assume width of 8m
    h_out = A_out / 15.0        # Assume width of 15m

    # Update both rivers separately
    A_river1, Q_river1 = update_river(A_river1, Q_river1, h_river1, dx, dt, nx_river1)
    A_river2, Q_river2 = update_river(A_river2, Q_river2, h_river2, dx, dt, nx_river2)

    # Apply confluence/junction conditions
    Q_out[0] = Q_river1[-1] + Q_river2[-1]
    h_junction = (h_river1[-1] + h_river2[-1]) / 2  # Average water level at junction
    A_out[0] = h_junction * 15.0  # Calculate cross-sectional area based on combined width

    # Update downstream (outflow) section
    A_out, Q_out = update_river(A_out, Q_out, h_out, dx, dt, nx_out)
