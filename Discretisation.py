"""
Defines functions regarding the Discretisation of Field Equation (5)
and Energy Equation (6):

next_timestep
next_frame
energy

Defines the variables:
    
L, N, lamb, dx, dt, frame_space
"""
import numpy as np
import matplotlib.pyplot as plt


# system variables
L = 100                     # length of domain
N = 512                     # number of nodes
lamb = 0.5                  # potential coefficient λ

# units 
dx = L / N                  # spatial distance between nodes
dt = dx / 4                 # time-step between states
frame_space = 10            # time between measurements (timesteps)


# coefficients
C_2 = dt**2 / dx**2
C_3 = -lamb * dt**2
C_1 = 2 - 2 * C_2 - C_3


def next_timestep(f_old, f):
    """
    Given the previous and current field configurations, 'f_old' and 'f',
    updates the field configurations to the next timestep
    according to equation (13)
    
    In particular uses the Rolling Array equation (16)
    
    Returns
    -------
    f :     current field configuration
    f_new : next field configuration 
    """
    return f, (-f_old + C_1 * f +
            C_2 * (np.roll(f, 1) + np.roll(f, -1)) +
            C_3 * f ** 3)


def next_frame(f_old, f):
    """
    Given the previous and current field configurations, 'f_old' and 'f',
    updates the field configurations to the next frame
    using 'next_timestep'
    
    Returns
    -------
    f :     current field configuration
    f_new : next field configuration 
    """
    # until next frame
    for _ in range( frame_space ):
        # update the system by one timestep
        f_old, f = next_timestep(f_old, f)
        
    return f_old, f


def energy(f_old, f):
    """
    Given the previous and current field configurations, 'f_old' and 'f',
    calculates the Kinetic, Interaction, Potential and Total Energy 
    according to equations (17), (18), (19), (20)
    
    Returns
    -------
    K :  kinetic term
    I :  interaction term
    P :  potential term
    E :  total energy
    """
    K = np.sum( ( f - f_old )**2 ) / (2 * dt * dt)
    I = np.sum( f * ( f - np.roll( f, 2 ) ) ) / (4 * dx * dx)
    P = lamb * np.sum( ( f * f - 1 )**2 ) / 4
    E = K + I + P
    return K, I, P, E