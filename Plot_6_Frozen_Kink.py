"""
Produces Figure 6

Attachs ground state to hot heat bath
Evolves for time
Attachs to cold heat bath

Returns plot of energy of time
Returns plot of 'frozen' kink
"""
from Discretisation import np, plt, L, N, energy, next_timestep
from Initial_Conditions import heat_bath_iteration

                            
T1 = 1.0                # temperature of hot heat bath
iter_max = 100          # iterations attached to hot heat bath
iter_max2 = 300         # timesteps of evolution between heat baths
T2 = 0.01               # temperature of cold heat bath
iter_max3 = 100         # iterations attached to cold heat bath
sigma_factor = 0.05     # standard deviation factor


# standard deviation
sigma1 = sigma_factor * np.sqrt(T1)
sigma2 = sigma_factor * np.sqrt(T2)

# prepare the ground state
f_old = -np.ones(N)
f = -np.ones(N)

# initialise arrays
K_array = np.zeros( iter_max + iter_max2 + iter_max3 )
I_array = np.zeros( iter_max + iter_max2 + iter_max3 ) 
P_array = np.zeros( iter_max + iter_max2 + iter_max3 )
E_array = np.zeros( iter_max + iter_max2 + iter_max3 )

# For a number of iterations
# Evolve in contact with a hot Heat Bath
for iter_num in range(iter_max):
    f = heat_bath_iteration(f_old, f, T1, sigma1)
    
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


fig1, ax1 = plt.subplots()
E = int( energy(f_old, f)[-1] )
axis = np.linspace(0,L,N)
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$f$')
ax1.set_ylim(-1.5, 1.5)
ax1.plot(axis, np.ones(N), color='black', \
     linestyle = 'dashed', alpha = 0.7)
ax1.plot(axis, np.zeros(N), color='black', \
                 linestyle = 'dashed', alpha = 0.7)
ax1.plot(axis, -np.ones(N), color='black', \
     linestyle = 'dashed', alpha = 0.7)
ax1.plot(axis, f, color='black')

# For a number of iterations
# Evolve in contact with a cold Heat Bath
for iter_num in range( iter_max + iter_max2, \
                      iter_max + iter_max2 + iter_max3 ):
    f = heat_bath_iteration(f_old, f, T2, sigma2)
    
    # Once per iteration, measure energies, print and store
    K, I, P, E = energy(f_old, f)
    K_array[iter_num] = K
    I_array[iter_num] = I
    P_array[iter_num] = P
    E_array[iter_num] = E
    
    # evolve by a timestep
    f_old, f = next_timestep(f_old, f)

        # plot results
fig2, ax2 = plt.subplots()
ax2.set_xlabel('Timesteps')
ax2.set_ylabel('Energy')
ax2.axvline( iter_max, linestyle='dashed', color = 'black')
ax2.axvline( iter_max + iter_max2, linestyle='dashed', color = 'black')
ax2.plot(E_array, label = 'Total Energy', color = 'black')
ax2.plot(K_array, label = 'Kinetic Term')
ax2.plot(I_array, label = 'Interaction Term')
ax2.plot(P_array, label = 'Potential term')
ax2.legend()


        # plot field
fig3, ax3 = plt.subplots()
E = int( energy(f_old, f)[-1] )
axis = np.linspace(0,L,N)
ax3.set_xlabel(r'$x$')
ax3.set_ylabel(r'$f$')
ax3.set_ylim(-1.5, 1.5)
ax3.plot(axis, np.ones(N), color='black', \
     linestyle = 'dashed', alpha = 0.7)
ax3.plot(axis, np.zeros(N), color='black', \
                 linestyle = 'dashed', alpha = 0.7)
ax3.plot(axis, -np.ones(N), color='black', \
     linestyle = 'dashed', alpha = 0.7)
ax3.plot(axis, f, color='black')
