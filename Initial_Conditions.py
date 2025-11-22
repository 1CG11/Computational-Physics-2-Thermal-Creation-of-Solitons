"""
Defines functions to implement the Metropolis Hastings Algorithm:

linear_coefficient
energy_diff
prob_accept
heat_bath_iteration
heat_bath
"""
from Discretisation import np, N, lamb, dx, dt, next_timestep


# Quadratic coefficient (lambda * p / 2) of the Node Energy Polynomial (29)
quad = 1 / (2 * dt * dt) + 1 / (4 * dx * dx) - 0.5*lamb


def linear_coefficient(f_old, f, k):
    """
    Linear coefficient (lambda * q) of the Node Energy Polynomial (29)
    """
    return - f_old[k] / (dt * dt) - (f[k-2] + f[(k+2)%N]) / (4 * dx * dx)


def energy_diff(f_old, f, k, z, y):
    """
    Calculates the Energy Difference according to Equation (30)
    using 'linear_coefficient'
    
    Behaves better factorised.
    """
    lin = linear_coefficient(f_old, f, k)
    D_E = (z-y) * (  (z+y) * ( lamb * (z*z+y*y) / 4 + quad )  + lin ) 
    return  D_E


def prob_accept(f_old, f, T, k, z, y):
    """
    Calculates the Acceptance Probability A(z|y) according to Equation (39)
    using 'energy_diff'
    """
    Diff_E = energy_diff(f_old, f, k, z, y)
    return np.exp(- Diff_E / T)


def heat_bath_iteration(f_old, f, T, sigma):
    """
    Updates the field configuration 'f' 
    according to one iteration of the Metropolis Hastings Algorithm
    using 'prob_accept'
    """    
    # Attach a Heat Bath
    # For each node, randomly ordered
    for k in np.random.permutation(N):
        y = f[k]
        
        # Propose a new value according to G(z|y)
        z = np.random.normal(y, sigma)
            
        # Accept change according to A(z|y)
        r = np.random.rand()
        if r < prob_accept(f_old, f, T, k, z, y):
            f[k] = z
    
    return f


def heat_bath(T, iter_max=100, sigma_factor=0.05):
    """
    Prepares a thermalised state of tmperature 'T'
    by applying 'iter_max' iterations of the Metropolis Hastings Algorithm
    using 'heat_bath_iteration'
    """
    # standard deviation
    sigma = sigma_factor * np.sqrt(T)
    
    # prepare the ground state
    f_old = -np.ones(N)
    f = -np.ones(N)

    # For a number of iterations
    # Evolve in contact with a Heat Bath
    for iter_num in range(iter_max):
        f = heat_bath_iteration(f_old, f, T, sigma)
        
        # evolve by a timestep
        f_old, f = next_timestep(f_old, f)

    return f_old, f