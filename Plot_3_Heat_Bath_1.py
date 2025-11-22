"""
Produces Figure 3

Attachs ground state to heat bath
Evolves for short period of time
Returns plot of Energy Distribution over Time

Evolves for long time period
Returns plot of Energy Distribution over log of time
"""
from Discretisation import np, plt, N, next_timestep, energy
from Initial_Conditions import heat_bath_iteration

                            
T = 1.0                 # temperature of heat bath
sigma_factor = 0.05     # standard deviation factor
iter_max = 100          # iterations attached to heat bath
iter_max2 = 150         # timesteps of independent evolution
tmax_dt = 10**6         # timesteps to trace evolution over
num_tests = 1000        # number of energy evaluations in this time


# standard deviation
sigma = sigma_factor * np.sqrt(T)

# prepare the ground state
f_old = -np.ones(N)
f = -np.ones(N)

# initialise arrays
K_array = np.zeros( iter_max + iter_max2 )
I_array = np.zeros( iter_max + iter_max2 ) 
P_array = np.zeros( iter_max + iter_max2 )
E_array = np.zeros( iter_max + iter_max2 )

# For a number of iterations
# Evolve in contact with a Heat Bath
for iter_num in range(iter_max):
    f = heat_bath_iteration(f_old, f, T, sigma)
    
    # Once per iteration, measure energies, print and store
    K, I, P, E = energy(f_old, f)
    K_array[iter_num] = K
    I_array[iter_num] = I
    P_array[iter_num] = P
    E_array[iter_num] = E
    
    # evolve by a timestep
    f_old, f = next_timestep(f_old, f)
    
#   For a number of iterations
#   Evolve independent of Heat Bath
for iter_num in range( iter_max, iter_max + iter_max2 ):
    
    # Once per iteration, measure energies and store
    K, I, P, E = energy(f_old, f)
    K_array[iter_num] = K
    I_array[iter_num] = I
    P_array[iter_num] = P
    E_array[iter_num] = E
    
    # Evolve by one time-step
    f_old, f = next_timestep(f_old, f)


        # plot results
fig = plt.figure(figsize=(12, 6), dpi=80)
fig.suptitle(
    'Energy Distribution Acheived by Heat Bath Algorithm at Temperature '
              +r'$T=$'+str(T) )

ax1 = fig.add_subplot(1, 2, 1)
ax1.set_xlabel('Timesteps')
ax1.set_ylabel('Energy')
ax1.axvline( iter_max, linestyle='dashed', color = 'black')
ax1.plot(E_array, color = 'black')
ax1.plot(K_array)
ax1.plot(I_array)
ax1.plot(P_array)


        # long term energy conservation ( logarithmic )
time_array_dt = 10**np.linspace( 
                    np.log10(iter_max), np.log10(tmax_dt), num_tests)

# initialise arrays
K_array = np.zeros( num_tests )
I_array = np.zeros( num_tests )
P_array = np.zeros( num_tests )
E_array = np.zeros( num_tests )
counter = iter_max

# For each specified time
for i, time_dt in enumerate( time_array_dt ):
    
    # evolve to next time
    while counter < time_dt:
        
        # progress bar
        # inside loop since each step takes longer than last
        print( str(counter) +' out of ' + str(tmax_dt), 
              str(100*counter/tmax_dt)+'%')
        
        f_old, f = next_timestep( f_old, f )
        counter += 1
    
    # At time, measure energies and store
    K, I, P, E = energy(f_old, f)
    K_array[i] = K
    I_array[i] = I
    P_array[i] = P
    E_array[i] = E

    
        #   plot results
ax2 = fig.add_subplot(1, 2, 2, sharey = ax1)
ax2.set_xscale('log')
ax2.set_xlabel('Timesteps')
ax2.plot(time_array_dt, E_array, \
         label = 'Total Energy', color = 'black')
ax2.plot(time_array_dt, K_array, label = 'Kinetic Term')
ax2.plot(time_array_dt, I_array, label = 'Interaction Term')
ax2.plot(time_array_dt, P_array, label = 'Potential term')
ax2.legend(loc='center left')