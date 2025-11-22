"""
Produces Figure 15

Energy Dependence of creation time 
"""
from Discretisation import np, N, plt
from Test_Functions import Gamma_and_tau_test


num_tests = 25          # number of tests
T_min = 0.5             # minimum temperature
T_max = 1.0             # maximum Temperature
e_tests = 1000          # nuber of energy evaluations
tmax_frame = 10**5      # number of frames for evolution to be traced over
            
T_array = np.linspace( T_min, T_max, num_tests )
            

#   initialise arrays
E_array = np.zeros( num_tests )
tau_array = np.zeros( num_tests )

for i, T in enumerate(T_array):
    
    #   progress bar
    print(str(i+1) + ' out of '+str(num_tests))
    
    E, gamma, tau = Gamma_and_tau_test(T, e_tests, tmax_frame)
    E_array[i] = E
    tau_array[i] = tau

    print( 'Energy: ' + str(int(E)) +' Creation Time: '+str(tau) )

        #   plot results
fig1, ax1 = plt.subplots()
ax1.set_xlabel('Energy  '+r'$E$')
ax1.set_ylabel('Creation Time ' + r'$\tau$')
ax1.scatter( E_array, tau_array, color='black', marker ='+')


fig2, ax2 = plt.subplots()
ax2.set_xlabel('Temperature  '+r'$T$')
ax2.set_ylabel(r'$\alpha=\dfrac{E}{NT}$')
ax2.scatter( T_array, E_array / (N*T_array), color='black', marker ='+' )