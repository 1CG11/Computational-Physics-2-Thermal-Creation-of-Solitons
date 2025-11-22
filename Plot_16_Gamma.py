"""
Produces Figure 16

Energy Dependence of creation rate gamma
"""
from Discretisation import np, plt, N
from Test_Functions import Gamma_and_tau_test


num_tests = 25          # number of tests
e_tests = 1000          # number of energy evaluations
tmax_frame = 10**5      # number of frames for evolution to be traced over
                        

#   initialise arrays
inverse_T_array = np.linspace( 0, 2, num_tests)
G_array = np.zeros( num_tests )
inverse_E_array = np.zeros( num_tests )

for i, inverse_T in enumerate(inverse_T_array):
    
    #   progress bar
    print(str(i+1) + ' out of '+str(num_tests))
    
    E, Gamma = Gamma_and_tau_test( 1/inverse_T, e_tests, tmax_frame )[:-1]
    
    inverse_E_array[i] = 1/E
    G_array[i] = Gamma



#   plot results
fig1, ax1 = plt.subplots()
ax1.set_xlabel('Inverse Energy  '+r'$1/E$')
ax1.set_ylabel('Creation Rate ' + r'$\Gamma$')
ax1.set_yscale('log')
ax1.scatter( inverse_E_array, G_array, color='black', marker ='+')


fig2, ax2 = plt.subplots()
ax2.set_xlabel('Temperature  '+r'$T$')
ax2.set_ylabel(r'$\alpha=\dfrac{E}{NT}$')
ax2.scatter( inverse_T_array, inverse_T_array / (N*inverse_E_array) , 
            color='black', marker ='+' )