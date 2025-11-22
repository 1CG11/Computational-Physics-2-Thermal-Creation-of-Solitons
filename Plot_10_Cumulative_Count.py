"""
Produces Figure 10

Plots cumulative count of kink number over time
"""
from Discretisation import np, plt, next_frame, energy, frame_space
from Initial_Conditions import heat_bath
from Kinks_and_Creations import pairs

T = 1.0
target = 10**4

f_old, f = heat_bath(T)
E = int( energy(f_old, f)[-1] )
cum_sum_pairs = np.array( [pairs(f)] )
time_array_dt = np.array( [0] )

while cum_sum_pairs[-1] < target:
    
    # evolve to next frame
    f_old, f = next_frame(f_old, f)
    time_array_dt = np.append( time_array_dt, \
                              time_array_dt[-1] + frame_space)
    cum_sum_pairs = np.append( cum_sum_pairs, \
                              cum_sum_pairs[-1] + pairs(f))
    print( str(cum_sum_pairs[-1]) +' out of ' +str(target))

plt.figure()
plt.title('Field prepared at Temperature ' +r'$T=$'+str(T) \
              +'\n Initial Energy '+r'$E=$'+str(E))
plt.xlabel('Timesteps')
plt.ylabel('Cumulative Sum of Pair Detections')
plt.plot(time_array_dt, cum_sum_pairs, drawstyle='steps-post' )