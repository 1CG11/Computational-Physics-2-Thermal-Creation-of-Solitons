"""
Defines functions used in the majority of plots:

heat_bath_T_test
zeros_and_wide_gaps_test
pairs_test
Gamma_and_tau_test
"""
from Discretisation import np, plt, N, next_timestep, next_frame, energy
from Initial_Conditions import heat_bath
from Kinks_and_Creations import zeros_and_wide_gaps, pairs, \
                            creation_rates, buff_frame

def heat_bath_T_test(T, iter_max, sigma_factor):
    """
    Prepares initial condition of temperature 'T'
    according to the Metropolis-Hastings Algorithm
    using 'heat_bath'
    Tracks Energy Distribution over the next 'iter_max' timesteps,
    returns averages
    using 'energy'
    
    Returns
    -------
    Kf_avg :  average kinetic fraction
    If_avg :  average interaction fraction
    Pf_avg :  average potential fraction
    E_avg :   average total energy
    alpha :   energy-temperature proportionality constant
    """
    # initialise state heat bath algorithm
    f_old, f = heat_bath(T, sigma_factor=sigma_factor)
    
    # Initialize arrays 
    Kf_array = np.zeros( iter_max )
    If_array = np.zeros( iter_max )
    Pf_array = np.zeros( iter_max )
    E_array = np.zeros( iter_max )
        
    # For a number of iterations
    # Evolve independent of Heat Bath
    # Measuring Energy
    for iter_num in range(iter_max):
        
        # Once per iteration, measure energies and store
        K, I, P, E = energy(f_old, f)
        Kf_array[iter_num] = K/E
        If_array[iter_num] = I/E
        Pf_array[iter_num] = P/E
        E_array[iter_num] = E
        
        # Evolve by one time-step
        f_old, f = next_timestep(f_old, f)
        
        
    #   calculate means since disconnection from heat bath
    Kf_avg = np.mean((Kf_array))
    If_avg = np.mean((If_array))
    Pf_avg = np.mean((Pf_array))
    E_avg = np.mean(E_array)
    
    #   calculate the constant of proportionality
    alpha = E_avg / (N*T)
    
    return Kf_avg, If_avg, Pf_avg, E_avg, alpha


def zeros_and_wide_gaps_test( T, e_tests, tmax_frame):
    """
    Average zeros and wide gaps at temperature
    
    Prepares initial condition of temperature 'T'
    according to the Metropolis-Hastings Algorithm
    using 'heat_bath'
    Evololved for 'tmax_frame' frames,
    using 'next_frame'
    Measures numbers of zero-crossings and wide gaps every frame,
    using 'zeros_and_wide_gaps'
    Measures energy 'e_tests' times throughout,
    using 'energy'
    
    
    Returns
    -------
    E_avg :   average total energy
    z_avg :   average number of zero-crossings
    g_avg :   average number of wide gaps
    """
    
    # intial conditions of this temperature
    f_old, f = heat_bath(T)
    
    # reset counters
    z = 0
    g = 0
    E = 0
    
    # for each frame until tmax
    for j in range( tmax_frame ):
            
        # evolve to next frame
        f_old, f = next_frame(f_old, f)
                
        # on frame, count zeros and wide gaps 
        z_new, g_new = zeros_and_wide_gaps( f )
        z += z_new
        g += g_new
            
        # rarely evaluate energy
        if j % (tmax_frame // e_tests) == 0:
            E += energy(f_old, f)[-1]
            
    # calculate average energy over all measurements
    E_avg = E / e_tests
    
    # calculate average average zeros and wide gaps over all frames
    z_avg = z / tmax_frame
    g_avg = g / tmax_frame
    
    return E_avg, z_avg, g_avg


def pairs_test( T, e_tests, tmax_frame):
    """
    Average pair number 'n' at temperature
    
    Prepares initial condition of temperature 'T'
    according to the Metropolis-Hastings Algorithm
    using 'heat_bath'
    Evololved for 'tmax_frame' frames,
    using 'next_frame'
    Measures numbers of pairs every frame,
    using 'pairs'
    Measures total energy 'e_tests' times throughout,
    using 'energy'
    
    Returns
    -------
    E_avg :   average total energy
    n_avg :   average number of pairs
    """
    # intial conditions of this temperature
    f_old, f = heat_bath(T)
    
    # reset counters
    n = 0
    E = 0
    
    # for each frame until tmax
    for j in range( tmax_frame ):
            
        # evolve to next frame
        f_old, f = next_frame(f_old, f)
                
        # on frame, count pairs 
        n += pairs( f )

            
        # rarely evaluate energy
        if j % (tmax_frame // e_tests) == 0:
            E += energy(f_old, f)[-1]
            
    # calculate average energy over all measurements
    E_avg = E / e_tests
    
    # calculate mean pair number over all frames
    n_avg = n / tmax_frame
    
    return E_avg, n_avg

def Gamma_and_tau_test( T, e_tests = 1000, tmax_frame=10**5 ):
    """
    Creation rate, creation time over 'tmax_frame' at temperature 'T'
    
    Prepares initial condition of temperature 'T'
    according to the Metropolis-Hastings Algorithm
    using 'heat_bath'
    Evololved for 'tmax_frame' frames,
    using 'next_frame'
    Measures numbers of pairs every frame,
    using 'pairs'
    Calculates creation rate 'gamma' and creation time 'tau',
    using 'creation_rates'
    Measures total energy 'e_tests' times throughout,
    using 'energy'
    
    Returns
    -------
    E_avg :   average total energy
    Gamma :   creation rate
    tau   :   creation time
    """
    # initial conditions of this temperature
    f_old, f = heat_bath(T)
        
    # reset counters
    k_array = np.zeros(tmax_frame + buff_frame)
    E = 0
        
    # for each frame until tmax_frame
    for j in range( tmax_frame + buff_frame ):
            
        # evolve to next frame
        f_old, f = next_frame(f_old, f)
                
        # on frame, count pairs 
        k_array[j] = pairs(f)
            
        # rarely evaluate energy
        if j % (tmax_frame // e_tests) == 0:
            E += energy(f_old, f)[-1]
            
    # calculate average energy over all measurements
    E_avg = E / e_tests
    
    # calculate the pair creation time
    Gamma, tau = creation_rates(k_array, tmax_frame)

    return E_avg, Gamma, tau
