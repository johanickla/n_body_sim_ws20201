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

def simu(sim, t):
    sim.integrator = "whfast"
    #sim.min_dt = 5.
    sim.dt = 1.
    # center of mass is at the origin and does not move
    sim.move_to_com()
    sim.init_megno()
    # Hide warning messages (WHFast timestep too large)
    with warnings.catch_warnings(record=True) as w:
         warnings.simplefilter("always")
         sim.integrate(t*2.*np.pi)
    return  sim, sim.calculate_megno(),(sim.calculate_lyapunov()/(2.*np.pi)) # returns MEGNO and Lypunov exp in 1/years
#--------------- plotlyapunov_t --------------------
def plotlyapunov_t(fig, ax1, lyapunovs, times): #lyapunov und zeit
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'Zeit $t$')
    ax1.plot(times,lyapunovs,'o-')
    ax1.grid()
    return fig
#-------------- plotlyapunov_a ---------------------
# lyapunov-exponent und große Halbachse
def plotlyapunov_a(l,a):
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel = 'Lyapunov-Exponent', xlabel = '$a$/$a_{Jupiter}$')
    ax1.plot(a,l,'o-')
    ax1.grid()
    # fig.savefig('lyapunov_a.png')
    return fig
# --------------- plotlyapunov_m ---------------------
def plotlyapunov_m(l,m):
    fig, ax1 = plt.subplots(1,1)
    ax1.set_xscale('log')
    ax1.set(ylabel = 'Lyapunov-Exponent', xlabel = '$m$/$M_{Helga}$')
    ax1.plot(m,l,'o-')
    ax1.grid()
    return fig
# -------------- plot lyapunov_t_multiple ------------
# creates a plot of lypanov-exponents over time for n EQUAL simulations of Helga
# starting from given 'start' time and stepping in time with 'stepsize'
def lyapunov_t_multiple(n, start, end, stepsize):
    #  r: number of steps
    r = int((end-start)/stepsize) + 1
    print('integrating', n, 'simulation(s)', 'with' ,r, 'steps')
    # initialize arrays for data to be stored
    times = np.arange(start,end,stepsize)
    # print(np.size(times))
    megno = np.zeros(r)
    lyapunov = np.zeros(r)
    fig, ax1 = plt.subplots(1,1)
    # h: ratio of semimajor axis (Helga/Jupiter)
    h = 0.696
    sim = visualize_orbit.setup('Helga',h)
    # for each simulation calculate the lyapnov exponent and megno
    for k in range(0,n,1):
        sim = visualize_orbit.setup('Helga',h)
        for i,t in enumerate(times):
            #i is index of years
            sim, m, l = simu(sim,t)
            megno[i] = m
            lyapunov[i] = l
        plotlyapunov_t(fig, ax1, lyapunov, times)
    fig.savefig('lyapunov_exp_t_mulitple.png')
#------------------- plot lyapunov_a_multiple ------------
# creates a plot of lyapunov-exponents over semi-major axis
#
def lyapunov_a_multiple(start,end,stepsize,t):
    # steps: number of setps
    steps = int((end-start)/stepsize) + 1
    # initialize arrays for data to be stored
    # H: Halbachsen
    # lyapunov: Lyapunov-Exponenten
    H = np.arange(start,end,stepsize)
    lyapunov = np.zeros(steps)

    for i,h in enumerate(H):
        sim = visualize_orbit.setup('Helga',h) #jedes mal neu initialisiert
        sim, m, l = simu(sim, t)
        lyapunov[i] = l
    fig = plotlyapunov_a(lyapunov,H)
    plt.title('Auswertung bei $t = %3d  2 \pi $' %t)
    fig.savefig('lyapunov_exp_a_variation.png')
# -----------------------------------------------------------------
def lyapunov_helga_m_stoerung(start, end, steps, t):
    # initialize arrays for data to be stored
    # lyapunov: Lyapunov-Exponenten
    lyapunov = np.zeros(steps)
    # masses: logaritmisch ansteigende Massen
    masses = np.logspace(start, end, steps)
    print(masses)
    m_Helga = 3e-13
    a_Helga = 5.204
    h = 0.696
    a_stoer = (a_Helga+0.001)*h
    for k in range(0,steps,1):
        sim = visualize_orbit.setup('Helga',h)
        m_stoer = masses[k]*m_Helga
        sim.add(m = m_stoer, a=a_stoer, M=38.41426796877275, omega=0.257, e=.08634521111588543 ,inc=4.4 )  # Helga NASA
        sim, m, l = simu(sim, t)
        lyapunov[k] = l
    print(lyapunov)
    fig = plotlyapunov_m(lyapunov, masses)
    plt.title('$t = %3d  2 \pi $, $a_{Stoer}$/$a_{Jupiter}$ = %5.4f ' %(t,a_stoer))
    fig.savefig('lyapunov_exp_m_variation.png')

#----------------------------------------------------------------
# set parameters for functions
t_n = 3
t_stepsize = 2e3
t_start = 1e3
t_end = 5e4

a_t = 1e4   # a_t: Auswertungszeitpunkt
a_start = 0.696
a_end = 0.699
a_stepsize = 0.00005

m_start = -1 # m_start: starting exponent base 10
m_end = 5 # m_end: ending exponent base 10
m_steps = 50
m_t = 1e5

if __name__ == "__main__":
    # lyapunov_t_multiple(t_n, t_start, t_end, t_stepsize)
    # lyapunov_a_multiple( a_start, a_end, a_stepsize, a_t)
    lyapunov_helga_m_stoerung(m_start,m_end,m_steps,m_t)
