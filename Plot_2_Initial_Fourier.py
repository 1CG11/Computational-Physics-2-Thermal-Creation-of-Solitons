"""
Produces Figure 2
"""
from Discretisation import np, plt, L, N, dt, next_timestep, energy


def initial_fourier(scale):
    """
    Mimics the thermal spectrum of a thermalised state
    Uses a Fourier transform to create an approximately thermalised state 
    Will serve as experimental initial conditions
    """
    # random phases
    phases_1 = 2*np.pi*np.random.rand(N)
    phases_2 = 2*np.pi*np.random.rand(N)
    
    # a range of amplitudes
    amplitudes_1 = 1 / np.concatenate(
                                      [100*np.ones(5), 
                                      np.arange(6, (N//2)+1), 
                                      np.arange((N//2), 0, -1)]
                                      )
    amplitudes_2 = 2 * np.pi / L
    
    # fourier transforms of these
    z1 = np.fft.fft( amplitudes_1 * np.exp( 1j * phases_1 ) )
    z2 = np.fft.fft( amplitudes_2 * np.exp( 1j * phases_2 ) )

    # build initial conditions from these
    f_old = -1 + scale * np.real(z1)
    f =  f_old + scale * np.real(z2) * dt

    return f_old, f
                            
                            
scale = 0.8             # fourier 'scale' parameter
iter_max = 250          # short term evolution 
tmax_dt = 10**6         # long term evolution
num_tests = 1000        # number of measurements to be made in long term


# Fourier Initial Conditions
f_old, f = initial_fourier( scale )

        # plot the initial conditions
E = int( energy(f_old, f)[-1] )
plt.figure()
plt.xlabel(r'$x$')
plt.ylabel(r'$f$')
axis = np.linspace(0,L,N)
plt.title('Fourier Initial Conditions with \'scale\' = ' \
              +str(scale) 
              +'\n Current Energy '+r'$E=$'+str(E))
plt.plot(axis, np.ones(N), color='black', \
         linestyle = 'dashed', alpha = 0.7)
plt.plot(axis, np.zeros(N), color='black', \
         linestyle = 'dashed', alpha = 0.7)
plt.plot(axis, -np.ones(N), color='black', \
         linestyle = 'dashed', alpha = 0.7)
plt.plot(axis, f, color='black')


        # short term energy conservation

# initialise arrays
K_array = np.zeros( iter_max )
I_array = np.zeros( iter_max )
P_array = np.zeros( iter_max )
E_array = np.zeros( iter_max )

# For a number of iterations
for iter_num in range( iter_max ):
    
    # Once per iteration, measure energies and store
    K, I, P, E = energy(f_old, f)
    K_array[iter_num] = K
    I_array[iter_num] = I
    P_array[iter_num] = P
    E_array[iter_num] = E
    
    # Evolve by one time-step
    f_old, f = next_timestep(f_old, f)
    

        #   plot results
fig = plt.figure(figsize=(12, 6), dpi=80)

ax1 = fig.add_subplot(1, 2, 1)
fig.suptitle('Fourier Initial Conditions with \'scale\' = ' \
              +str(scale) )
ax1.set_xlabel('Timesteps')
ax1.set_ylabel('Energy')
ax1.plot(E_array, label = 'Total Energy', color = 'black')
ax1.plot(K_array, label = 'Kinetic Term')
ax1.plot(I_array, label = 'Interaction Term')
ax1.plot(P_array, label = 'Potential term')


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
    
    # progress bar
    # inside loop since each step takes longer than last
    print( str(counter) +' out of ' + str(tmax_dt), 
          str(100*counter/tmax_dt)+'%')
    
    # evolve to next time
    while counter < time_dt:
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