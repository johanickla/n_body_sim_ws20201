import rebound
import numpy as np
import warnings
import matplotlib.pyplot as plt
from matplotlib import cm
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
    jov_yr = 11.86
    # Hide warning messages (WHFast timestep too large)
    with warnings.catch_warnings(record=True) as w:
         warnings.simplefilter("always")
         sim.integrate(t)
    return  sim,(sim.calculate_lyapunov()/(2.*np.pi)) # returns Lypunov exp in 1/years
#--------------- plotlyapunov_t --------------------
# lyapunov und zeit
def plotlyapunov_t(lyapunov, times, k):
    fig, ax1 = plt.subplots(1,1)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'Zeit $t$')
    ax1.plot(times,lyapunov,'o-', label = 'run %d' %(k+1))
    ax1.legend()
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
    ax1.set(ylabel = 'Lyap.-Expo.', xlabel = '$m$/$M_{Helga}$')
    if isinstance(l,list):
        n_l = len(l)
        for i in range(n_l):
            steps = len(m[i])
            ax1.plot(m[i],l[i],'o-', label = ' run %d with %d points' %(i,steps))
    else: ax1.plot(m,l,'o-', label = 'a fixed')
    ax1.legend()
    ax1.grid()
    return fig
#---------------- plotlyapunov_mt_surface ---------------------
def plotlyapunov_ma_surface(m,a,l):
    m = np.log10(m)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(10,8))
    # print('shape m:', np.shape(m))
    # print('shape a:', np.shape(a))
    m, a = np.meshgrid(m, a)
    print('shape m:', np.shape(m))
    print('shape a:', np.shape(a))
    ax.set(zlabel = 'Lyap.-Expo.', xlabel = '$\log_{10}$($m$/$M_{Helga})$',
                ylabel= '$a$/$a_{Jupiter}$')
    # ax.set_xlim(1e-13, 1e-8)
    # ax.set_ylim(3.62,3.64 )
    # transpose lyapunovs
    # print('shape l:', np.shape(l))
    # l = np.log10(l)
    l = np.transpose(l)
    print('shape l:', np.shape(l))
# Plot the surface.
    surf = ax.plot_surface(m, a, l, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    return fig
#---------------- plotlyapunov_max ---------------------
def plotlyapunov_max(masses,a_max_lyapunov,t):
    fig, ax1 = plt.subplots(1,1)
    ax1.set_xscale('log')
    ax1.plot(masses,a_max_lyapunov, label= '')
    ax1.set(xlabel = '$\log_{10}$($m$/$M_{Helga})$',
                ylabel= '$a$/$a_{Jupiter}$')
    plt.title('max. Lyapunovexponent nach $t = %3d$ years %(t)')

    ax1.grid()
    return fig
#------------------------------------------------------
# -------------- lyapunov_t_multiple ------------
# creates a plot of lypanov-exponents over time for n EQUAL simulations of Helga
# starting from given 'start' time and stepping in time with 'stepsize'
def lyapunov_t_multiple(n, start, end, steps):
    # initialize arrays for data to be stored
    times = np.logspace(start, end, steps)
    print('size(times):',np.size(times))
    #  r: number of steps
    r = np.size(times)
    print('integrating', n, 'simulation(s)', 'with' ,r, 'steps')
    megno = np.zeros(r)
    lyapunov = np.zeros(r)
    fig, ax1 = plt.subplots(1,1)
    # h: ratio of semimajor axis (Helga/Jupiter)
    h = 0.696
    sim = visualize_orbit.setup('Helga', h)
    # for each simulation calculate the lyapnov exponent and megno
    for k in range(0,n,1):
        sim = visualize_orbit.setup('Helga', h)
        # lyapunov = man_lyapunov(sim, times)
        for i,t in enumerate(times):
            #i is index of years
            sim, l = simu(sim,t)
            # megno[i] = m
            lyapunov[i] = l
            # if t == times[r-1]:
            #     sim.status()
    fig = plotlyapunov_t(lyapunov, times, k)
    fig.savefig('lyapunov_exp_t_mulitple.png')
#------------------- plot lyapunov_a_multiple ------------
# calculates lyap. expo. for a perturbation body of fixed mass and
# varying distance to Helga
# creates a plot of lyapunov-exponents over semi-major axis
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
        # lyapunov = man_lyapunov(sim, [t])
        sim, l = simu(sim, t)
        lyapunov[i] = l
    fig = plotlyapunov_a(lyapunov,H)
    plt.title('Auswertung bei $t = %3d  2 \pi $' %t)
    fig.savefig('lyapunov_exp_a_multiple.png')
# -----------------------------------------------------------------
def lyapunov_helga_m_stoerung(start, end, steps, t):
    # initialize arrays for data to be stored
    # lyapunov: Lyapunov-Exponenten
    lyapunov = np.zeros(steps)
    # masses: logaritmisch ansteigende Massen
    masses = np.logspace(start, end, steps)
    # print(masses)
    m_Helga = 3e-13
    a_Jupiter = 5.204
    h = 0.696
    a_stoer = (a_Jupiter+0.001)*h
    for k in range(0,steps,1):
        sim = visualize_orbit.setup('Helga',h)
        m_stoer = masses[k]*m_Helga
        sim.add(m = m_stoer, a=a_stoer, M=38.41426796877275, omega=0.257, e=.08634521111588543 ,inc=0.0768 )  # Helga NASA
        # lyapunov = man_lyapunov(sim, [t])
        sim, l = simu(sim, t)
        lyapunov[k] = l
    # print(lyapunov)
    fig = plotlyapunov_m(lyapunov, masses)
    plt.title('$t = %3d  2 \pi $, $a_{Stoer}$ = %5.4f ' %(t,a_stoer))
    fig.savefig('lyapunov_exp_m_stoerung.png')

#------------------ lyapunov_m_multiple_stoerung
def lyapunov_m_multiple_stoerung(m_start,m_end,m_steps,m_t,m_runs):
    lyapunovs_array = []
    masses_array = []
    m_Helga = 3e-13
    a_Jupiter = 5.204
    h = 0.696
    a_stoer = (a_Jupiter+0.001)*h
    for r in range(m_runs):
        masses = (np.logspace(m_start, m_end, (2**r)*m_steps))
        # print('Massen: ',masses, 'length: ', len(masses))
        lyapunovs = [0]*((2**r)*m_steps)
        for k in range((2**r)*m_steps):
            m_stoer = masses[k] * m_Helga
            sim = visualize_orbit.setup('Helga',h)
            sim.add(m=m_stoer, a=a_stoer, M=38.41426796877275, omega=0.257, e=.08634521111588543 ,inc=0.0768 )  # Helga NASA
            sim, l = simu(sim, m_t)
            print(sim.N)
            lyapunovs[k] = l
        # print('lyapunovs: ' , lyapunovs, 'length: ', len(lyapunovs))
        masses_array.append(masses)
        lyapunovs_array.append(lyapunovs)
    fig = plotlyapunov_m(lyapunovs_array, masses_array)
    plt.title('$t = %3d  2 \pi $, $a_{Stoer}$ = %5.4f ' %(m_t,a_stoer))
    fig.savefig('lyapunov_exp_m_multiple_stoerung.png')

#----------------------------------------------------------------
def lyapunov_helga_ma_stoerung(m_start,m_end,m_steps,a_start,a_end,a_steps,t):
    m_Helga = 3e-13
    a_Jupiter = 5.204
    h = 0.696
    # lyapunovs: O-array
    lyapunovs = np.zeros((m_steps, a_steps))
    a_max_lyapunov = np.zeros((m_steps))

    # print(lyapunovs)
    # masses: logaritmisch ansteigende Massenfaktoren
    masses = np.logspace(m_start, m_end, m_steps)
    masses2 = masses * m_Helga
    print('Massen: ',masses2)
    # a: linear ansteigende Halbachsenfaktoren
    a = np.linspace(a_start, a_end, a_steps)
    a2 = a*a_Jupiter
    print('Halbachsen: ',a2)
    for k in range(0,m_steps,1):
        for j in range(0,a_steps,1):
            sim = visualize_orbit.setup('Helga',h)
            m_stoer = masses2[k]
            a_stoer = a2[j]
            sim.add(m=m_stoer, a=a_stoer, M=38.41426796877275, omega=0.257, e=.08634521111588543 ,inc=0.0768 )  # Helga NASA
            # lyapunov = man_lyapunov(sim, [t])
            sim, l = simu(sim, t)
            if l < 0:
                l = 1e-8
            # print(l)
            lyapunovs[k,j] = l
        a_max_lyapunov[k] = a[ np.where( lyapunovs[k,:] == max(lyapunovs[k,:]) ) ]
    # print(lyapunovs)
    fig1 = plotlyapunov_ma_surface(masses, a, lyapunovs)
    plt.title('$t = %3d$ years' %(t))
    fig1.savefig('lyapunov_exp_ma_stoerung.png')
    fig2 = plotlyapunov_max(masses,a_max_lyapunov,t)
    fig2.savefig('maximaler Lyapunovexponent')

#----------------------------------------------------------------
#-------------------- set parameters for functions---------------
t_n = 1
t_steps = 30
t_start = 1
t_end = 7

a_t = 1e5   # a_t: Auswertungszeitpunkt
a_start_1 = 0.696
a_end_1 = 0.699
a_stepsize_1 = 0.00005

m_start = 1 # m_start: starting exponent base 10
m_end = 10 # m_end: ending exponent base 10
m_steps = 30 # m_steps: Mindestanzahl an Schritten
m_t = 1e5
m_runs = 4 # m_runs: wie oft soll die Schrittzahl verdoppelt werden

m_start_2 = 1 # m_start: starting exponent base 10
m_end_2 = 3 # m_end: ending exponent base 10
m_steps_2 = 50
a_start_2 = 0.695
a_end_2 = 0.701
a_steps_2 = 30
t = 1e5

if __name__ == "__main__":
    # lyapunov_t_multiple(t_n, t_start, t_end, t_steps)
    # lyapunov_a_multiple( a_start_1, a_end_1, a_stepsize_1, a_t)
    # lyapunov_helga_m_stoerung(m_start,m_end,m_steps,m_t)
    # lyapunov_m_multiple_stoerung(m_start,m_end,m_steps,m_t,m_runs)
    lyapunov_helga_ma_stoerung(m_start_2,m_end_2,m_steps_2,a_start_2,a_end_2,a_steps_2,t)
