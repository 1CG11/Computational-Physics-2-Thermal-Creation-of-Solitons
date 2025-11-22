"""
Produces Figures 1 and 5

Prepares a field configuration at temperature 'T'
Evolves for 'tmax' timesteps
Plots 'num_frames' field configurations over this evolution
"""

from Discretisation import np, plt, L, N, next_timestep, energy
from Initial_Conditions import heat_bath


T = 1.0             # temperature
tmax = 10**4        # max timesteps
num_frames = 10     # number of plots


space = tmax // num_frames  # timesteps between plots
axis = np.linspace(0,L,N)   # to plot field over
f_old, f = heat_bath(T)     # intitial conditions

for i in range(num_frames):
    
    # progress bar
    print( str(i+1) + ' out of ' + str(num_frames))
    
    # evolve over time
    for _ in range(space):
        f_old, f = next_timestep(f_old, f)
    # measure energy
    E = int( energy(f_old, f)[-1] )
    
    # plot
    plt.figure()
    plt.xlabel(r'$x$')
    plt.ylabel(r'$f$')
    plt.ylim(-1.5, 1.5)
    plt.plot(axis, np.ones(N), color='black', \
         linestyle = 'dashed', alpha = 0.7)
    plt.plot(axis, np.zeros(N), color='black', \
                     linestyle = 'dashed', alpha = 0.7)
    plt.plot(axis, -np.ones(N), color='black', \
         linestyle = 'dashed', alpha = 0.7)
    plt.plot(axis, f, color='black')