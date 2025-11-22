"""
Animates field prepared at temperature 'T',
plotting once every 'ani_frame_space' timesteps

Not used in project directly, 
we found it useful to see how the field behaved
"""
from Discretisation import np, plt, L, N, next_timestep, energy
from Initial_Conditions import heat_bath
from Kinks_and_Creations import pairs
from matplotlib.animation import FuncAnimation


T = float( input("Enter Temperature: " ) )
ani_frame_space = int( input("Enter Frame Spacing: ") )

def animate(i):
    global f_old, f, TimeCounter
    if not pause:
        
        for _ in range(ani_frame_space):
            f_old, f = next_timestep(f_old, f)
            TimeCounter+=1
        E = int(energy(f_old, f)[-1])
        p = pairs(f)

            
    line.set_ydata(f)
    title.set_text("\n Time Elapsed = "+str(TimeCounter) +' timesteps'
                       +"\n Energy = "+str(E)
                       +"\n Pairs = "+str(p))
    return line,


# Pause function
def toggle_pause(event):
    global pause
    pause = not pause
    if pause:
        ani.event_source.stop()
    else:
        ani.event_source.start()
        

# Initialize 
TimeCounter=0
f_old, f = heat_bath(T)
fig, ax = plt.subplots()
ax.set_ylim(-1.5, 1.5)
x = np.linspace(0, L, N)
title = ax.text(0.5,0.85,"",bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")
line, = ax.plot(x, f, color='black')
line3, = ax.plot(np.zeros(L), alpha=0.5, \
                 linestyle='dashed', color = 'black')
line4, = ax.plot(1.0*np.ones(L), alpha=0.5, \
                 linestyle='dashed', color = 'black')
line5, = ax.plot(-1.0*np.ones(L), alpha=0.5, \
                 linestyle='dashed', color = 'black')
pause = False
plt.connect('key_press_event', toggle_pause)

ani = FuncAnimation(fig, animate, frames=None, interval=100, 
                    cache_frame_data=False)
plt.show()