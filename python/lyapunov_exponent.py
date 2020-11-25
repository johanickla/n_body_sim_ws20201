import rebound
import numpy as np
import multiprocessing
import warnings
import matplotlib
import matplotlib.pyplot as plt
import visualize_orbit


sim = visualize_orbit.setup()
# visualize_orbit.PAForbitplot(sim1)
# a = 1
# tmax = a * 2*np.pi*1e3
# Ntimes = a*10
# delta_t=1.
# sim1.init_megno()
# lyapunov = sim1.calculate_lyapunov()/(2.*np.pi) #in years
# sim1, times, x, y, z = visualize_orbit.PAFintegrate(sim1, tmax, Ntimes, delta_t)
# #visualize_orbit.plotabsolutes(x,y,z,sim1.N-1)
# #self, seed=None
# print(lyapunov)

def simu(t):
    sim.integrator = "whfast"
    #sim.min_dt = 5.
    sim.dt = 1.

    sim.move_to_com() #center of mass is at the origin and does not move
    sim.init_megno()
    # Hide warning messages (WHFast timestep too large)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        sim.integrate(t*2.*np.pi)
        return [sim.calculate_megno(),(sim.calculate_lyapunov()/(2.*np.pi))] # returns MEGNO and Lypunov exp in 1/years

def plotlyapunov(l,T): #lyapunov und zeit
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'Zeit')
    ax1.plot(T,l,'o-')
    ax1.grid()
    fig.savefig('lyapunov.png')

n=30
step=1e3
T=np.arange(1e3,1e3+n*step,step)
megno=np.zeros(n)
lyapunov=np.zeros(n)
for i,t in enumerate(T): #in years
    m,l=simu(t)
    megno[i]=m
    lyapunov[i]=l

print(lyapunov)
plotlyapunov(lyapunov,T)
