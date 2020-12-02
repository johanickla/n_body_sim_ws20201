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

def simu(t):
    sim.integrator = "whfast"
    #sim.min_dt = 5.
    sim.dt = 1.
    sim.move_to_com() # center of mass is at the origin and does not move
    sim.init_megno()
    # Hide warning messages (WHFast timestep too large)
    with warnings.catch_warnings(record=True) as w:
         warnings.simplefilter("always")
         sim.integrate(t*2.*np.pi)
    return  sim, sim.calculate_megno(),(sim.calculate_lyapunov()/(2.*np.pi)) # returns MEGNO and Lypunov exp in 1/years
#--------------- plotlyapunov_t --------------------
def plotlyapunov_t(fig, ax1, lyapunovs, times): #lyapunov und zeit
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'Zeit')
    ax1.plot(times,lyapunovs,'o-')
    ax1.grid()
    return fig
#-------------- plotlyapunov_a ---------------------

def plotlyapunov_a(l,a): #lyapunov und große Halbachse
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel = 'Lyapunov-Exponent', xlabel = 'a/a_Jupiter')
    ax1.plot(a,l,'o-')
    ax1.grid()
    fig.savefig('lyapunov_a.png')

# -------------- plot lyapunov_t_multiple ------------
# creates a plot of lypanov-exponents over time for n different simulations
# starting from given start time 'start' and stepping in time with 'step'
def plotlyapunov_t_multiple(n, start, end, step):
    #  r: number of steps
    r = (start-end)/step
    times = np.arange(start,end,step)
    megno = np.zeros(start)
    lyapunov = np.zeros(r)
    fig, ax1 = plt.subplots(1,1)
    # h =  ratio of semimajor axis (Helga/Jupiter)
    h=0.696
    sim = visualize_orbit.setup('Helga',h)
    # for each simulation calculate the lyapnov exponent, megno and
    # return the simulation
    return the simulation
    for k in range(0,n,1):
        for i,t in enumerate(times):
            #i is index of years
            sim, m, l = simu(t)
            megno[i] = m
            lyapunov[i] = l
        # print(lyapunov)
        plotlyapunov_t(fig, ax1, lyapunov, times)
    fig.savefig('lyapunov_mit_return_sim.png')
#------------------- plot lyapunov_a_multiple ------------
# creates a plot of lyapunov-exponents over the semimajor axis
#
def plotlyapunov_a_multiple(start,end,step):
    lyapunov_a = np.zeros(31)
    H = np.zeros(31)
    for i,h in enumerate(np.arange(start,end,step)):
        sim = visualize_orbit.setup('Helga',h) #jedes mal neu initialisiert
        sim, m, l = simu(start+n*step)
        lyapunov_a[i] = l
        H[i] = h
        plotlyapunov_a(lyapunov_a,H)
    return lyapunov_a, H
#--------------------------------------
# set parameters for functions
n_1=1
step_1=2e3
start_1=1e6

start = 0.696
end = 0.699
step = 0.0001

if __name__ == "__main__":
