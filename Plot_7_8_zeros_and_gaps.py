"""
Produces Figures 7 and 8

Energy Dependence of Average Numbers of Zero-Crossings and Wide Gaps
"""
from Discretisation import np, plt, N, frame_space
from Initial_Conditions import heat_bath
from Test_Functions import zeros_and_wide_gaps_test


num_tests = 25          # number of tests
T_min = 10**-1          # minimum temperature
T_max = 10**3           # maximum Temperature
e_tests = 1000          # nuber of energy evaluations
tmax_frame = 10**5      # number of frames for evolution to be traced over


# initialise arrays
T_array = 10**np.linspace( np.log10(T_min), np.log10(T_max), num_tests )
E_array = np.zeros( num_tests )
zeros_array = np.zeros( num_tests )
gaps_array = np.zeros( num_tests )


# for each temperature
for i, T in enumerate(T_array):
    
    # progress bar
    print(str(i+1) + ' out of '+str(num_tests))
    
    # intial conditions of this temperature
    f_old, f = heat_bath(T)
    
    # calculate average energy, number of zeros and gaps
    E, z, g = zeros_and_wide_gaps_test(T, e_tests, tmax_frame)
    E_array[i] = E
    zeros_array[i] = z
    gaps_array[i] = g
    
    print('Energy: ' + str(int(E)) +' zeros: '+str(z) +' wide gaps: '+str(g))
    
        # plot results
power = np.log10( tmax_frame * frame_space)
fig1, ax1 = plt.subplots()
ax1.set_xlabel('Energy  '+r'$E$')
ax1.set_ylabel('Average Number of Zero Crossings')
ax1.set_xscale('log')
ax1.scatter( E_array, zeros_array, color ='black', marker ='+')

fig2, ax2 = plt.subplots()
ax2.set_xlabel('Energy  '+r'$E$')
ax2.set_ylabel('Average Number of Big Gaps')
ax2.set_xscale('log')
ax2.scatter( E_array, gaps_array, color ='black', marker ='+')

fig3, ax3 = plt.subplots()
ax3.set_xlabel('Temperature  '+r'$T$')
ax3.set_ylabel(r'$\alpha=\dfrac{E}{NT}$')
ax3.set_xscale('log')
ax3.scatter( T_array, E_array / (N*T_array), color ='black', marker ='+' )

