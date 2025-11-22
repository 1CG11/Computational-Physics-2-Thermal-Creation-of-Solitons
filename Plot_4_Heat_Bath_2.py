"""
Produces Figure 3
"""
from Discretisation import np, plt
from Test_Functions import heat_bath_T_test


num_tests = 100         # number of simulations to run
T_min = 10**-4          # minimum temperature
T_max = 10**4           # maximum temperature
iter_max = 100          # evolution time
sigma_factor = 0.05     # standard deviation factor


# initialise arrays
T_array = 10**np.linspace( np.log10(T_min), np.log10(T_max), num_tests )
Kf_avg_array = np.zeros( num_tests )
If_avg_array = np.zeros( num_tests )
Pf_avg_array = np.zeros( num_tests )
E_avg_array = np.zeros( num_tests )
alpha_array = np.zeros( num_tests )

# for each temperature
# calculate energy distribution acheived
for i,T in enumerate(T_array):
    
    # progress bar
    print( str(i+1) + ' out of ' + str(num_tests) )

    Kf_avg, If_avg, Pf_avg, E_avg, alpha = \
        heat_bath_T_test(T, iter_max, sigma_factor)
    Kf_avg_array[i] = Kf_avg 
    If_avg_array[i] = If_avg 
    Pf_avg_array[i] = Pf_avg 
    E_avg_array[i] = E_avg 
    alpha_array[i] = alpha
    
    
        #   plot results
fig, (ax1, ax2) = plt.subplots(2, sharex=True)

ax1.set_ylabel('Energy '+r'$E$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.plot(T_array, E_avg_array, color='black')

ax2.set_xlabel('Temperature ' +r'$T$')
ax2.set_ylabel(r'$\alpha=\dfrac{E}{NT}$')
ax2.set_xscale('log')
ax2.scatter(T_array, alpha_array, color='black', marker ='+')

fig3, ax3 = plt.subplots()
ax3.set_xlabel('Temperature ' +r'$T$')
ax3.set_ylabel('Fraction of Total Energy')
ax3.set_ylim(0.0, 1.0)
ax3.set_xscale('log')
ax3.plot(T_array, Kf_avg_array, label = 'Kinetic Term')
ax3.plot(T_array, If_avg_array, label = 'Interaction Term')
ax3.plot(T_array, Pf_avg_array, label = 'Potential Term')
ax3.legend()