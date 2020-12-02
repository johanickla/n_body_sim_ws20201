# Import the rebound module
import rebound
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#--------------- Create Simulation object

#h= semimajoraxis Helga/Jupiter
def setup(option,h = 1):
    if option=='Helga':
        sim = rebound.Simulation()
        sim.add( m=1. )                         # sun
        sim.add( m=0.000954, a=5.204, M=0.600, omega=0.257, e=0.0489,inc = 0.02276) # Planet 1 (Jupiter) Dtaen von Wiki
        sim.add(m=3e-13, a=5.204*h, M=38.41426796877275, omega=0.257, e=.08634521111588543 ,inc=0.0768 )  # Helga NASA
        #Masse unbekannt?
        #sim.add( "NAME=Pluto")
    elif option=='fictional':
        sim = rebound.Simulation()
        sim.add( m=1. )                         # sun
        sim.add( m=1e-3, a=1., e=0.1, inc = 45 )          # Planet 1 (Jupiter)
        sim.add(m=1e-3, x=1., y=0., z= -1., vx=0, vy=1, vz=0 )           # lesser Planet 2 (fictional)
        sim.add( m=1e-3, x=1., y=0., z= 1., vx=0, vy=1, vz=0 )         # Planet 3 (fictional)
    return sim


#------------- calculate -------------

def PAFintegrate(arbsim ,tmax, Ntimes,delta_t):
    # setting time scale
    times = np.linspace(0,tmax,Ntimes)
    #zerros = np.zeros((len(times),arbsim.N-1))
    # create equal arrays for each observable
    x = np.zeros((len(times),arbsim.N-1))
    y = np.zeros((len(times),arbsim.N-1))
    z = np.zeros((len(times),arbsim.N-1))
    # set integrator and dt
    arbsim.integrator = "whfast"
    arbsim.dt = delta_t #sets the timestep
    arbsim.move_to_com() #center of mass is at the origin and does not move
    for i, t in enumerate(times):
        arbsim.integrate(t)
        particles = arbsim.particles
        for k in range(1,arbsim.N,1):
            x[i,k-1] = particles[k].x
            y[i,k-1] = particles[k].y
            z[i,k-1] = particles[k].z
    return arbsim, times, x, y, z

#-------------------- test plotting ----------------------------
# fig1 = plt.figure(figsize=(5,5))
# plt.plot(times2,x2)
# fig1.savefig("test.png")
#--------------------   plotting xyz   ----------------------------
def plotabsolutes(x,y,z,n):
    axes = ()
    fig1, (ax1,ax2,ax3) = plt.subplots(3, 1)
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312, sharex=ax1)
    ax3 = plt.subplot(313, sharex=ax1)
    for k in range(n):
        ax1.set( ylabel= 'x')
        ax1.plot(times, x[:,k], label = 'Planet %d' %(k+1))
        plt.setp(ax1.get_xticklabels(), visible = False)
        ax1.grid()
        ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=n, mode="expand", borderaxespad=0.)
        #---------
        ax2.set(ylabel = 'y')
        ax2.plot(times,y[:,k])
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.grid()
        #---------
        ax3.set(xlabel = 'time', ylabel = 'z')
        ax3.plot(times,z[:,k])
        ax3.grid()
    fig1.savefig('orbit_in_cartesians.png')
#------------------ plot phase diagram (only useful for two coordinate)
def plotphase(x1,x2):
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='$x_2$', xlabel = '$x_1$')
    ax1.plot(x1,x2,'o-')
    ax1.grid()
    fig.savefig('orbit_as_phase_plot.png')

def PAForbitplot(sim):
    fig, ax = rebound.OrbitPlot(sim, fancy = True, color=True, lw=1)
    fig.savefig('orbitplot.png')

#---------------- the real deal ---------------
h=0.696
sim1 = setup('Helga',h)
# sim1.status()
PAForbitplot(sim1)
# sim1.status()
a = 5
tmax = a * 2*np.pi
Ntimes = a*10
delta_t=1e-3
sim1, times, x, y, z = PAFintegrate(sim1, tmax, Ntimes, delta_t)
#sim1.status()
plotabsolutes(x,y,z,sim1.N-1)
plotphase(x,y)
#sim1.status()
