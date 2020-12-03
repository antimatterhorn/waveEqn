# Looking at solving the 1D wave eqn as a precursor to an acoustics question
# TUE December  1, 2020

N = 1000
finalCycle = 1000000
tgoal = 1.0

# d^2 phi / dt^2 = c^2 d^2 phi / dx^2

# phi = sin(k*x - w*t) works
# phi'' = -w^2 * sin(k*x - w*t)
# phi,, = -k^2 * sin(x-w*t)
# ergo
# phi'' = (k/w)^2 phi,,

# IIRC, which solution you get all comes down to initial conditions
# and boundary conditions.

# RHS is second spatial derivative.  Easy.  From Taylor,
# phi_{i+1} = phi_i + dx * phi,_i + dx^2/2 phi,,_i + ...
# phi_{i-1} = phi_i - dx * phi,_i + dx^2/2 phi,,_i + ...
# --------- = -----------------------------------
#           = 2 * phi_i           + 2 * dx^2/2 phi,,_i + ...

#  phi+1 + phi-1 = 2phi + dx^2 * phi,,
#  phi,, = (phi+1 + phi-1 - 2phi)/dx^2

# Time derivative gives us many (sigh, many many) options.  Does this
# need to be implicit?  I dunno.  Let's start explicit and see what
# the issues are.

# split time derivative into two steps
# d phi / dt = xi
# d^2 phi / dt^2 = d xi / dt
# evolve xi and then use xi to evolve phi

# Let's go with forward Euler (omg yuck) for the first version!  Super
# easy!

# xi^{n+1}_i - xi^{n}_i = dt * c^2 * phi,,_i
# phi^{n+1}_i - phi^{n}_i = dt * xi^{n+1}_i

# phi^{n+1}_i = phi^{n}_i + dt * xi^{n+1}_i
# phi^{n+1}_i = phi^{n}_i + dt * xi^{n}_i + dt^2 * c^2 * phi,,_i
# phi^{n+1}_i = phi^{n}_i + dt * xi^{n}_i + dt^2 * c^2 * 0.5*(phi_{i-1} - 2*phi_i + phi_{i+1})/dx^2
# phi^{n+1}_i = phi^{n}_i + dt * xi^{n}_i + (dt*c/dx)^2 * 0.5*(phi_{i-1} - 2*phi_i + phi_{i+1})

# boundary conditions... reflecting for right hand boundary seems fine.
# Left hand can be a forcing function if we want.

import numpy as np
from math import pi, sin
plot = True
printIt = False
from matplotlib import pyplot as plt

L = 1.0
c = 1.0
# c = 1.0/0.7
# c = 1.35
c2 = c*c
dx = L/N
x = np.arange(0,1+dx,dx)
print (x)
dt = .1*dx/c # this is a guess --- I have not done any analysis.  It
             # could be (dx/c)**2 or something even worse.
phi = np.zeros(N+1)
xi = np.zeros(N+1)

# initial conditions
t = 0.0
cycle = 0

# traveling sine wave?
# k = 2*pi/(L/.25)
# phi = sin(x)

# # little tent on the left
# phi[1] = 0.5
# phi[2] = 1.0
# phi[3] = 0.5

# source on left
phi[0] = 1.0
def source(t):
    T = (L/3.0) / c  # period
    w = 2 * pi / T
    return sin(w*t)

d2phidx2 = np.zeros(N+1)

while cycle < finalCycle and t < tgoal:
    # compute phi,,
    # phi,,_i = 0.5*(phi_{i-1} - 2*phi_i + phi_{i+1})/dx^2
    d2phidx2[1:N] = 1.0*(phi[0:N-1] - 2 * phi[1:N] + phi[2:N+1])/dx**2 # 1.0 for doug's posterity
    
    # compute xi
    xi = xi + dt*c2*d2phidx2
    
    # compute phi
    phi = phi + dt*xi

    # reflecting wall boundary conditions on right
    phi[N] = phi[N-1]
    # source BC on left
    phi[0] = source(t)
    
    # update time and cycle
    t += dt
    cycle += 1

    # print state
    # print('phi = ',)
    if printIt:
        for i in range(N+1):
            print('%9.2e' % phi[i], end='')
        print()

# plot if we're plotting
if plot:
    plt.plot(x,phi, '-o')#,marker='o',color='blue',linewidth=4,markersize=8)
    plt.xlabel('cm')
    plt.ylabel('$\phi$')
    plt.text(.8, .8, 't = %4.2e' % t)
    plt.text(.8, .9, '%4d' % cycle)
    plt.show()

print('...all done...')

    
