# Looking at solving the 2D wave eqn as a precursor to an acoustics question
# WED December  2, 2020

N = 100  # number of points in each direction
finalCycle = 5000
tgoal = 2.0

plot = False
animate = True
printIt = False

# d^2 phi / dt^2 = c^2 del^2 phi

# IIRC, which solution you get all comes down to initial conditions
# and boundary conditions.

# RHS is second spatial derivative.  Easy.  From Taylor,
# phi_{i+1,j} = phi_ij + dx * phi,x_ij + dx^2/2 phi,,x_ij + ... # where ,,x means d^2/dx^2
# phi_{i-1,j} = phi_ij - dx * phi,x_ij + dx^2/2 phi,,x_ij + ...
# phi_{i,j+1} = phi_ij + dy * phi,y_ij + dy^2/2 phi,,y_ij + ... # where ,,y means d^2/dy^2
# phi_{i,j-1} = phi_ij - dy * phi,y_ij + dy^2/2 phi,,y_ij + ...
# --------- = -----------------------------------
#           = 4 * phi_ij          + 2 * dx^2/2 phi,,x_ij + 2 * dy^2/2 phi,,y_ij + ...
#           = 4 * phi_ij          + dx^2 phi,,x_ij + dy^2 phi,,y_ij + ...
#           = 4 * phi_ij          + dx^2 (phi,,x_ij + phi,,y_ij) + ... # assume dx=dy
#           = 4 * phi_ij          + dx^2 del^2(phi_ij) + ... # assume dx=dy

# del^2 phi_ij = (-4*phi_{i,j} + phi_{i+1,j} + phi_{i-1,j} phi_{i,j+1} + phi_{i,j-1})/dx^2

# Time derivative gives us many (sigh, many many) options.  Does this
# need to be implicit?  I dunno.  Let's start explicit and see what
# the issues are.

# split time derivative into two steps
# d phi / dt = xi
# d^2 phi / dt^2 = d xi / dt
# d xi / dt = c^2 del^2 phi
# evolve xi and then use xi to evolve phi

# Let's go with forward Euler (omg yuck) for the first version!  Super
# easy!

# xi^{n+1}_i - xi^{n}_i = dt * c^2 * del^2 phi^n_ij
# phi^{n+1}_i - phi^{n}_i = dt * xi^{n+1}_ij


# boundary conditions... reflecting for all boundaries seems fine to start with.
# We can put a forcing function in the middle to start.

import numpy as np
from math import pi, sin, sqrt
from matplotlib import pyplot as plt

L = 1.0  # an L x L square in physical space
c = 1.0  # sound speed
c2 = c*c
dx = L/N
y, x = np.meshgrid(np.linspace(0, L, N+1), np.linspace(0, L, N+1))
# x = np.arange(0,1+dx,dx)
# y = np.arange(0,1+dx,dx)
dt = .1*dx/c # this is a guess --- I have not done any analysis.  It
             # could be (dx/c)**2 or something even worse.
phi = np.zeros([N+1, N+1])
bounds = np.zeros([N+1,N+1])
c2  = np.ones([N+1, N+1])*c2

xi = np.zeros([N+1, N+1])

# initial conditions
t = 0.0
cycle = 0

def init():
    radius = L*0.25
    center = 0.5
    for i in range(0,N+1):
        for j in range(0,N+1):
            if sqrt((x[i,j]-center)**2+(y[i,j]-center)**2) > radius:
                c2[i,j]*=5.0

# oscillating source
def source(t):
    T = (L/5.0) / c  # period is such that wavelength is L/5
    w = 2 * pi / T   # angular frequency 
    return sin(w*t)

class box:
    def __init__(self,xmin,xmax):
        self.xmin = xmin
        self.xmax = xmax
    def addBounds(self):
        for i in range(N+1):
            for j in range(N+1):
                xmin = self.xmin
                xmax = self.xmax
                if x[i,j] > xmin[0] and x[i,j] < xmax[0] and y[i,j] > xmin[1] and y[i,j] < xmax[1]:
                    bounds[i,j] = 1.0
        return
    def removeBounds(self):
        for i in range(N+1):
            for j in range(N+1):
                xmin = self.xmin
                xmax = self.xmax
                if x[i,j] > xmin[0] and x[i,j] < xmax[0] and y[i,j] > xmin[1] and y[i,j] < xmax[1]:
                    bounds[i,j] = 0.0
        return

class circle:
    def __init__(self,center,radius):
        self.radius = radius
        self.center = center
    def addBounds(self):
        for i in range(N+1):
            for j in range(N+1):
                r2 = (x[i,j]-self.center[0])**2+(y[i,j]-self.center[1])**2
                if r2 < self.radius**2:
                    bounds[i,j] = 1.0
        return
    def removeBounds(self):
        for i in range(N+1):
            for j in range(N+1):
                r2 = (x[i,j]-self.center[0])**2+(y[i,j]-self.center[1])**2
                if r2 < self.radius**2:
                    bounds[i,j] = 0.0
        return

del2phi = np.zeros([N+1,N+1])

if animate:
    from matplotlibBlitManager import BlitManager
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    plt.xlabel('cm')
    plt.ylabel('cm')
    plt.title('$\phi$')
    # ax.set_ylim(-2.1, 2.1)
    cm = ax.pcolormesh(x, y, phi[:-1,:-1], vmin=0, vmax=1, animated = True)
    cyc = ax.annotate('0', (0.5,1.0), animated=True)
    fig.colorbar(cm)
    bm = BlitManager(fig.canvas, [cm, cyc])
    # make sure window is on the screen and drawn
    plt.show(block=False)
    plt.pause(0.2)

def applyBoundaries():
    # reflecting wall boundary conditions everywhere
    phi[0:N+1,0]    = phi[0:N+1,1]
    phi[0:N+1,N]    = phi[0:N+1,N-1]
    phi[0,0:N+1]    = phi[1,0:N+1]
    phi[N,0:N+1]    = phi[N-1,0:N+1] 

    for i in range(1,N):
        for j in range(1,N):
            if bounds[i,j] == 1.0:
                # get average phi from all adjacent cells that are not bounds
                lb = (1.0-bounds[i-1,j])
                rb = (1.0-bounds[i+1,j])
                ub = (1.0-bounds[i,j+1])
                db = (1.0-bounds[i,j-1])
                left    = lb*phi[i-1,j]
                right   = rb*phi[i+1,j]
                up      = ub*phi[i,j+1]
                down    = db*phi[i,j-1]
                nums = np.asarray([left,right,up,down])
                phimax = max(nums.min(),nums.max(),key=abs)
                phi[i,j] = phimax

box1 = box((0.2,0.2),(0.8,0.8))
box1.addBounds()
circle1 = circle((0.5,0.5),0.2)
circle1.removeBounds()
box3 = box((0.2,0.45),(0.5,0.55))
box3.removeBounds()

#init()
    
while cycle < finalCycle and t < tgoal:
    # compute del^2 phi
    del2phi[1:N,1:N] =(-4*phi[1:N,1:N] + phi[2:N+1,1:N] + phi[0:N-1,1:N] + phi[1:N,2:N+1] + phi[1:N,0:N-1])/dx**2
    
    # compute xi
    xi = xi + dt*np.multiply(c2,del2phi)
    
    # compute phi
    phi = phi + dt*xi

    # apply any boundaries you have
    applyBoundaries()

    # source in the middle
    phi[N//2, N//2] = source(t)
    
    # update time and cycle
    t += dt
    cycle += 1

    # print state
    # print('phi = ',)
    if printIt:
        for i in range(N+1):
            print('%9.2e' % phi[i], end='')
        print()

    # plot
    if animate:
        cm.set_array(phi[:-1, :-1].flatten())
        cyc.set_text('cycle: {}, time = {:4.2f}'.format(cycle,t))
        # tell the blitting manager to do its thing
        bm.update()

    
# plot if we're plotting
if plot:
    fig, ax = plt.subplots()
    cm = plt.pcolormesh(x,y,phi)

    # cm = ax.pcolormesh(x,y,np.zeros([N+1, N+1]))
    xxx = phi[:-1, :-1]
    q = np.zeros([N+1, N+1])
    for i in range(N):
        for j in range(N):
            q[i,j] = 1.0
    cm.set_array(q.ravel())
    cm.set_array(phi.ravel())
    cm.set_array(xxx.ravel())
    # plt.draw()

    ax.set_aspect('equal')
    plt.xlabel('cm')
    plt.ylabel('cm')
    plt.title('$\phi$')
    plt.text(.8, .8, 't = %4.2e' % t)
    plt.text(.8, .9, '%4d' % cycle)
    plt.show()

if animate:
    plt.show()
    
print('...all done...')

    
