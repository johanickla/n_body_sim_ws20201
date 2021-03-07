import rebound
import numpy as np
import warnings

import visualize_orbit
# from plotting_programs import *


#------------- calculate the hill sphere radius of an object with mass m, half achse a
#----------------- surrounding an object with size M at an excentricity e
def hill_sphere_radius(a,m,M,e):
    r = a*(1-e) * ( (m/(3*M))**(1/3) )
    return r
#------------- calculate h for a given ratio of orbit periods using Keplers Law ----------
def orbit_ratio(n,m):
    h = (n/m)**(2/3)
    return h

def simu(sim, t):
    #set integrator and timestep
    sim.integrator = "whfast"
    sim.dt = 1.
    # center of mass is at the origin and does not move
    sim.move_to_com()
    #initializes the MEGNO particles and makes integration possible
    sim.init_megno()
    # Hide warning messages (WHFast timestep too large)
    with warnings.catch_warnings(record=True) as w:
         warnings.simplefilter("always")
         sim.integrate(t)
    return  sim, sim.calculate_lyapunov() # returns Lypunov exp in 1/years

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
    # for each simulation calculate the lyapnov exponent
    for k in range(0,n,1):
        sim = visualize_orbit.setup('Helga', h)
        # lyapunov = man_lyapunov(sim, times)
        for i,t in enumerate(times):
            sim, l = simu(sim,t)
            lyapunov[i] = l
            # if t == times[r-1]:
            #     sim.status()
    fig = plotlyapunov_t(lyapunov, times, k)
    fig.savefig('plots/lyapunov_exp_t_mulitple.png')
#------------------- plot lyapunov_a_multiple ------------
# calculates lyap. expo. for varying position of Helga
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
    plt.title('Auswertung bei $t = %3d $' %t)
    fig.savefig('resonance_width_plots/lyapunov_exp_helga_a_variation.png')
#------------------- plot lyapunov_e_multiple ------------
# calculates lyap. expo. for varying excentricity of Helga
# creates a plot of lyapunov-exponents over e
def lyapunov_e_multiple(e_start,e_end,e_steps,e_t):
    e_stepsize = (e_end - e_start) / e_steps
    # initialize arrays for data to be stored
    # E: Exzentrizitäten
    # lyapunov: Lyapunov-Exponenten
    h = 0.696
    E = np.arange(e_start,e_end,e_stepsize)
    lyapunov = np.zeros(e_steps)
    for i,exc in enumerate(E):
        sim = visualize_orbit.setup('resonant',h,exc) #jedes mal neu initialisiert
        # lyapunov = man_lyapunov(sim, [t])
        sim, l = simu(sim, e_t)
        lyapunov[i] = l
    fig = plotlyapunov_e(lyapunov,E)
    plt.title('$t = %3d $' %e_t)
    fig.savefig('plots/lyap_exp_helga_e_variation.png')
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
    fig.savefig('plots/lyapunov_exp_m_stoerung.png')

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
    fig.savefig('plots/lyapunov_exp_m_multiple_stoerung.png')

#----------------------------------------------------------------
def lyapunov_helga_ma_stoerung(m_start,m_end,m_steps,a_start,a_end,a_steps,t,h):
    m_Helga = 3e-13
    a_Jupiter = 5.204
    # lyapunovs: O-array
    lyapunovs = np.zeros((m_steps, a_steps))
    a_max_lyapunov = np.zeros((m_steps))
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
            sim.add(m=m_stoer, a=a_stoer, M=38.414268, omega=0.257, e=.0863452 ,inc=0.0768 )  # Helga NASA
            # lyapunov = man_lyapunov(sim, [t])
            sim, l = simu(sim, t)
            if l < 0:
                l = 1e-8
            # print(l)
            lyapunovs[k,j] = l
        a_max_lyapunov[k] = a[ np.argmax(lyapunovs[k,:]) ]
    # print(lyapunovs)
    fig1 = plotlyapunov_ma_surface(masses, a, lyapunovs,h)
    plt.title('$t = %3d$ years' %(t))
    fig1.savefig('plots/lyapunov_exp_ma_stoerung.png')
    # fig2 = plotlyapunov_max(masses,a_max_lyapunov,t)
    # fig2.savefig('maximaler Lyapunovexponent')

#----------------------------------------------------------------
#-------------------- set parameters for functions---------------
t_n = 1
t_steps = 30
t_start = 1
t_end = 7

a_t = 1e5  # a_t: Auswertungszeitpunkt
a_start_1 = 0.696
a_end_1 = 0.699
a_stepsize_1 = 0.00005

e_start = .07
e_end = .09
e_steps = 50
e_t = 1e2

m_start = 1 # m_start: starting exponent base 10
m_end = 10 # m_end: ending exponent base 10
m_steps = 30 # m_steps: Mindestanzahl an Schritten
m_t = 1e5
m_runs = 4 # m_runs: wie oft soll die Schrittzahl verdoppelt werden

m_start_2 = 0 # m_start: starting exponent base 10
m_end_2 = 6 # m_end: ending exponent base 10
m_steps_2 = 100
h = 0.696    # orbit_ratio(11,12)
a_start_2 = h*0.998      # 0.696
a_end_2 = h*1.002        # 0.699
a_steps_2 = 100

t = 1e2

if __name__ == "__main__":
    # lyapunov_t_multiple(t_n, t_start, t_end, t_steps)
    lyapunov_a_multiple( a_start_1, a_end_1, a_stepsize_1, a_t)
    # lyapunov_e_multiple( e_start, e_end, e_steps, e_t)
    # lyapunov_helga_m_stoerung(m_start,m_end,m_steps,m_t)
    # lyapunov_m_multiple_stoerung(m_start,m_end,m_steps,m_t,m_runs)
    # lyapunov_helga_ma_stoerung(m_start_2,m_end_2,m_steps_2,a_start_2,a_end_2,a_steps_2,t,h)
