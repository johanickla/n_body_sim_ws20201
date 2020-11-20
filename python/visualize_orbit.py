# Import the rebound module
import rebound
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# Create Simulation object
def setup():
    sim = rebound.Simulation()
    sim.add( m=1. )                         # sun
    sim.add( m=1e-3, a=1., e=0.1, inc = 45 )          # Planet 1 (Jupiter)
    sim.add( m=1e-3, a=1.001, e=0.1, inc = 0 )           # lesser Planet 2 (fictional)
    # sim.add( m=1e-6, a=1., e=0.2, inc = 0 )         # Planet 3 (fictional)
    return sim

def PAFintegrate(arbsim ,tmax, Ntimes):
    # setting time scale
    times = np.linspace(0,tmax,Ntimes)
    zerros = np.zeros((len(times),arbsim.N-1))
    # create equal arrays for each observable
    x = np.zeros((len(times),arbsim.N-1))
    y = np.zeros((len(times),arbsim.N-1))
    z = np.zeros((len(times),arbsim.N-1))
    # set integrator and dt
    arbsim.integrator = "whfast"
    arbsim.dt = 1e-3
    for i, t in enumerate(times):
        arbsim.integrate(t)
        particles = arbsim.particles
        for k in range(1,arbsim.N,1):
            x[i,k-1] = particles[k].x
            y[i,k-1] = particles[k].y
    return arbsim, times, x, y, z

sim1 = setup()
tmax = 10 * 2*np.pi
Ntimes = 1000
sim1, times, x, y, z = PAFintegrate(sim1, tmax, Ntimes)
#-------------------- test plotting ----------------------------
# fig1 = plt.figure(figsize=(5,5))
# plt.plot(times2,x2)
# fig1.savefig("test.png")
#--------------------   plotting    ----------------------------
def plotabsolutes(x,y,z):
    axes = ()
    fig1, (ax1,ax2,ax3) = plt.subplots(3, 1)

    ax1 = plt.subplot(311)
    ax1.set( ylabel= 'x')
    ax1.plot(times, x)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.grid()
    #---------
    ax2 = plt.subplot(312, sharex=ax1)
    ax2.set(ylabel = 'y')
    ax2.plot(times,y)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.grid()
    #---------
    ax3 = plt.subplot(313, sharex=ax1)
    ax3.set(xlabel = 'time', ylabel = 'z')
    ax3.plot(times,z)
    ax3.grid()

    fig1.savefig('real1.png')

def plotphase(x,y):
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='y', xlabel = 'x')
    ax1.plot(x,y)
    ax1.grid()
    fig.savefig('real2.png')



plotabsolutes(x,y,z)
plotphase(x,y)
# ax4 = plt.subplot(412)
# ax4.set(xlabel = 'x', ylabel = 'y')
# plt.plot(x,y)
# ax4.grid()
