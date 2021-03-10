#-------------- resonance_width.py
# Import the rebound module
import rebound
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from plotting_programs import plotlyapunov_a
from scipy.optimize import curve_fit
from lyapunov_exponent import hill_sphere_radius,orbit_ratio
#---------------------------------
a_resonants = np.array([0.697956,0.697910,0.697700,0.697498,0.697155,0.696579])
#--------------- functions for fitting --------------------------------
def gaussian_2(x,sigma=1,mu=1,a=1):  # sigma: std dev., mu: mean, a: height
    return a * np.exp(-((x-mu)/(2*sigma))**2) # *1/np.sqrt(2*np.pi*sigma**2)
def gaussian_4(x,sigma=1,mu=1,a=1):    # sigma: std dev., mu: mean, a: height
    return a * np.exp(-((x-mu)/(2*sigma))**4) # *1/np.sqrt(2*np.pi*sigma**2)
def power2_3(x,alpha):
    return  alpha*x**2/3
def gen_power(x,alpha,beta):
    return alpha*x**beta
#--------------- Create Simulation object ----------------------------
def setup(option, h=1, exc=1, mu=1):
    sim = rebound.Simulation()
    if option=='Helga':
        sim.add( m=1. )                         # sun
        sim.add( m=0.000954, a=5.204, M=0.600, omega=0.257, e=0.0489,inc = 0.02276) # Planet 1 (Jupiter) Dtaen von Wiki
        sim.add( m=3e-13, a=5.204*h, M=38.41426796877275, omega=0.257, e=.08634521111588543 ,inc=0.0768 )  # Helga NASA
        #sim.add( "NAME=Pluto")
    elif option == 'resonant':
        sim.add( m=1. )                         # sun
        sim.add( m=mu, a=5.204, M=0.600, omega=0.257, e=0.0489,inc = 0.02276) # Planet 1 (generalized Jupiter)
        sim.add( m=3e-13, a=5.204*h, M=38.414268, omega=0.257, e= exc,inc=0.0768 )  # Helga NASA
    elif option == 'general_resonant':
        sim.add( m=1. )                         # sun
        sim.add( m=mu, a=5.204, M=0.600, omega=0.257, e=0.0489,inc = 0.02276) # Planet 1 (Jupiter) Dtaen von Wiki
        sim.add( m=1e-13, a=5.2044*h, M=38.414268, omega=0.257, e= 0.086 ,inc=0.0768 )  # Helga NASA
    elif option=='fictional':
        sim.add( m=1. )                                         # sun
        sim.add( m=1e-3, a=1., e=0.1, inc = 0 )                 # Planet 1 (Earth)
        sim.add( m=1e-4, a=1., e=0.1, inc = np.pi/2 )            # lesser Planet 2 (fictional)
        # sim.add( m=1e-3, x=1., y=0., z= 1., vx=0, vy=1, vz=0 )         # Planet 3 (fictional)
    return sim
#--------- returning an array of suitable size for resonance observation ------
def a_resonant(mu):
    b = 0.6981
    a = -0.0012/0.001
    a0 = 5.204
    c2 = 0.06
    # erg = b*(1-c2*mu**(1/2))
    erg = b+a*mu
    return erg
def presumed_res_width(m):
    return 0.15*(m**(2/3))
def semimajors_surrounding(mu,a_steps):
    a_res = orbit_ratio(3,4)
    a_start = a_res*0.995
    a_end = a_res*1.0
    A = np.linspace(a_start,a_end,a_steps)
    # a_res = a_resonant(mu)
    # c = presumed_res_width(mu)
    # C = np.linspace(a_res-c, a_res, int(a_steps/2)-1, endpoint = False)
    # D = np.linspace(a_res, a_res+c, int(a_steps/2)+1, endpoint = True)
    # A = np.concatenate((C,D), axis=0)
    return A
#----------------- fitting of all calc. resonant widths ------------------------
def multi_resonance_width_fit(rw,Mu):
    alpha = rw[-1]/(Mu[-1]**(2/3)) # for initial alpha assume power 2/3 and use last data point
    p0 = [alpha,2/3]
    params, params_conv = curve_fit(gen_power, Mu, rw, p0 )
    fitted_rw = gen_power(Mu, *params)
    return fitted_rw, params
#------------------ calculate a single lyap. exponent -------------------------
def lyapunov_calculator(sim, t):
    sim.integrator = "whfast"
    sim.dt = 1.         # sim.min_dt = 5.
    sim.move_to_com()   # center of mass is at the origin and does not move
    sim.init_megno()
    # jov_yr = 11.86
    # with warnings.catch_warnings(record=True) as w:
    #      warnings.simplefilter("always")
    sim.integrate(t)
    return  sim, sim.calculate_lyapunov()/(2*np.pi) # returns Lypunov exp in 1/years
#------------ calc. mulitple lyap. exponents for Helga at varying a -----------
def lyapunov_a_multiple(A,t,m_Jupiter):
    L = np.zeros(len(A))                     # L: Lyapunov exponents
    for i,a in enumerate(A):
        sim_helga = setup('general_resonant', h=a, mu=m_Jupiter) #jedes mal neu initialisiert
        sim, l = lyapunov_calculator(sim_helga, t)
        L[i] = l
    return L
# -------------- get the resonance width by fitting -----------------------
def single_resonance_width(L,A,mu):
    lyap_max = np.amax(L)
    maxindex = np.where(L == lyap_max)[0][0]
    print('maxindex: ', maxindex, '\nlyap_max: ', lyap_max)
    a_lyap_max = A[maxindex]
    c = (A[-1]-A[0])/2
    p_start = [c, A[len(A)//2], lyap_max]
    parameters,parameters_cov = curve_fit(gaussian_4, A, L, p0 = p_start
            # bounds=([0, A[0], lyap_max/4],[c, A[-1], 2*lyap_max])
            )
    print('Parameters (sigma,mu,a): ', parameters)
    return  parameters
#--------------------- main function: -----------------------------------------
#------ initialize, perform calculations, plot single res(mu) ----------------
def multi_resonance_width(a_steps,Mu,t):
    A = np.zeros( (len(Mu), a_steps) )
    L = np.zeros( (len(Mu), a_steps) )
    Res_Widths = np.zeros_like(Mu)

    for i in range(len(Mu)):
        A[i,:] = semimajors_surrounding(Mu[i],a_steps)
        L[i,:] = lyapunov_a_multiple(A[i,:], t, Mu[i])
        params = single_resonance_width(L[i,:], A[i,:],Mu[i])
        Res_Widths[i] = abs(params[0])
        # params =[0.001, 0.698, 1e-6]
        plotsingles(L[i,:], A[i,:], params, Mu[i])
    return Res_Widths
#------------- plot all resonant widths with fit ---------------------
def plot_res_widths(Res_Widths,Mu):
    # SMALL_SIZE = 8
    # matplotlib.rc('font', size=11)
    RWfit, params = multi_resonance_width_fit(Res_Widths,Mu)
    fig, ax1 = plt.subplots()
    alphastr = r'$\alpha$'
    betastr = r'$\beta$'
    ax1.set(xscale = 'log',yscale = 'log',
            title = '7/12-resonance'
            )
    plt.ylabel('resonance width ($a_{Jupiter}$)', fontsize=11.5)
    plt.xlabel('mass ratio $\mu = m_{planet}/m_{star}$', fontsize=11.5)
    ax1.plot(Mu,Res_Widths,'o-',label='calculated')
    ax1.plot(Mu, RWfit,'r-',label=('fit: a$\mu^b$' ) )
    ax1.grid()
    leg = ax1.legend(fontsize = 12)
    leg.set_title('fitting parameters: \n a: %.2e \n b: %.2e' %(params[0],params[1]),prop={'size':12})
    # ax1.ticklabel_format(axis="y", style="sci", scilimits = (-7,-4))
    return fig
#------------- plot single resonant with fit ---------------------
def plotsingles(L,A,params,mu):
    fig = plotlyapunov_a(L,A)
    plt.title('$t = %d\cdot2\pi$, mass ratio $\mu = %.5f$' %(t,mu))
    fit = np.zeros(len(A))
    for i in range(len(A)):
         fit[i] = gaussian_4(A[i], *params)
    plt.plot(A, fit, 'r-',label='fit: $b\cdot\exp (-[(a-a_0)/(2\sigma)]^4)$')
    leg = plt.legend(fontsize = 11.5)
    leg.set_title('fitting parameters: \n $\sigma$: %.4e \n $a_0$: %.4e \n $b$: %.4e' %(params[0],params[1],params[2]))
    fig.savefig('resonance_width_plots/single_lyapunovs_over_a/lyap_exp_a_variation_%.6f.png' %(mu))
#--------------------------------------------------------------


if __name__ == '__main__':
    t = 1e6  # t: Auswertungszeitpunkt
    n = 3
    m = 4
    a_start = orbit_ratio(n,m)*0.997
    a_end = orbit_ratio(n,m)*1.0
    a_steps = 50
    a_stepsize = (a_end - a_start)/a_steps
    # A = np.linspace(a_start, a_end, a_steps)
    # A: array großer Halbachsen
    mass_ratio_start = -5
    mass_ratio_end = -4
    mass_ratio_steps = 10
    Mu = np.logspace(mass_ratio_start,mass_ratio_end,mass_ratio_steps)
    # Mu: array von Massenverhältnissen
    print('Masses of largest Planet (Jupiter):\n', Mu)
    reswidth = multi_resonance_width(a_steps, Mu, t)
    print('Resonance Width: \n',reswidth)
    fig1 = plot_res_widths(reswidth, Mu)
    fig1.savefig('resonance_width_plots/res_width_over_mu')
