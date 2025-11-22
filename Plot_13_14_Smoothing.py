"""
Produces Figures 13 and 14

Plots pair number over time, unsmoothed and smoothed
"""
from Discretisation import np, plt, frame_space, next_frame, energy
from Initial_Conditions import heat_bath
from Kinks_and_Creations import pairs, smooth,buff_frame

T = 1.0
tmax_frame = 200

time_array_dt = np.arange(0, (buff_frame+tmax_frame)*frame_space, frame_space)
num_frames = len(time_array_dt)
pairs_array = np.zeros( num_frames)
f_old, f = heat_bath(T)
E = 0

for i in range(num_frames):
    f_old, f = next_frame(f_old, f)
    
    pairs_array[i] = pairs(f)
    E += energy(f_old, f)[-1]
    
print( 'Average Energy: '+str( E / num_frames ) )
fig1, ax1 = plt.subplots()
ax1.set_xlim(0, tmax_frame*frame_space)
ax1.set_ylim(-0.1, 2.1)
ax1.set_xlabel('Timesteps')
ax1.set_ylabel('Pair Number '+r'$n$')
ax1.plot(time_array_dt, pairs_array, color='black', drawstyle='steps-post')

smoothed_array = smooth(pairs_array, tmax_frame)
fig2, ax2 = plt.subplots()
ax2.set_xlim(0, tmax_frame*frame_space)
ax2.set_ylim(-0.1, 2.1)
ax2.set_xlabel('Timesteps')
ax2.set_ylabel('Pair Number '+r'$n$')
ax2.plot(time_array_dt[:tmax_frame], smoothed_array, 
         color='black', drawstyle='steps-post')