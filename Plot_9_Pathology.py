"""
Produces Figure 9
Pathological Field Configuration (artificially)
"""

import numpy as np
import matplotlib.pyplot as plt
L = 100
N = 512
axis = np.linspace(0,L,N)
array=np.sin(2*np.pi*axis/L)


v = np.array( [0.2, 0.5, 0.8, 1.0, 1.2, 1.2, 1.0, 0.8, 0.5, 0.2] )
u = np.zeros(90)

array4 = np.concatenate( (u, -v,  u, -v, np.zeros(112), v, u, v, u ) )

f = np.cbrt(array) + array4 + np.random.random(N)*0.1

fig1, ax1 = plt.subplots()
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