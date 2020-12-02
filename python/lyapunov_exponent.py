import rebound
import numpy as np
import multiprocessing
import warnings
import matplotlib
import matplotlib.pyplot as plt
import visualize_orbit



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

def simu(sim,t):
    sim.integrator = "whfast"
    #sim.min_dt = 5.
    sim.dt = 1.
    sim.move_to_com() #center of mass is at the origin and does not move
    sim.init_megno()
    # Hide warning messages (WHFast timestep too large)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        sim.integrate(t*2.*np.pi)
    return  sim, sim.calculate_megno(),(sim.calculate_lyapunov()/(2.*np.pi)) # returns MEGNO and Lypunov exp in 1/years

def plotlyapunov_t(fig, ax1, l,T): #lyapunov und zeit
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'Zeit')
    ax1.plot(T,l,'o-')
    ax1.grid()
    return fig

def plotlyapunov_a(l,h): #lyapunov und zeit
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'a/a_Jupiter')
    ax1.plot(h,l,'o-')
    ax1.grid()
    fig.savefig('lyapunov_a.png')


def lyapunov_t(n,step,start):
    T=np.arange(start,start+n*step,step)
    megno=np.zeros(n)
    lyapunov=np.zeros(n)
    fig, ax1 = plt.subplots(1,1)
    #h= semimajoraxis Helga/Jupiter
    for i,t in enumerate(T): #in years
        sim, m, l=simu(t)
        megno[i]=m
        lyapunov[i]=l

        print(lyapunov)

        plotlyapunov_t(fig, ax1, lyapunov, T)
    fig.savefig('lyapunov_t.png')

def lyapunov_a(n,step,start,T):
    #step and start beziehen sich auf das Verhältnis der großen halbachse a/Halbachse von a_Jupiter
    #T ist der betrachtete Zeitpunkt
    lyapunov=np.zeros(n)
    H=np.zeros(n)
    for i,h in enumerate(np.arange(start,start+n*step,step)):
        sim = visualize_orbit.setup('Helga',h) #jedes mal neu initialisiert
        sim, m, l=simu(sim,T)
        lyapunov[i]=l
        H[i]=h
    plotlyapunov_a(lyapunov,H)

#lyapunov zu verschiedenen Zeitpunkten
sim = visualize_orbit.setup('Helga',0.969)
#lyapunov_t(10,2e3,1e3)

#lyapunov für einen Zeitpunkt für verschiedene Halbachsen
lyapunov_a(31,0.0001,0.696,1e5)
