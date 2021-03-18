#-------------- resonance_width.py
# Import the rebound module
import rebound
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from plotting_programs import plotlyapunov_a
from scipy.optimize import curve_fit
from lyapunov_exponent import orbit_ratio
#---------------------------------
a_resonants = np.array([0.697956,0.697910,0.697700,0.697498,0.697155,0.696579])
#--------------- functions for fitting --------------------------------
def gaussian_2(x,sigma=1,mu=1,a=1):  # sigma: std dev., mu: mean, a: height
    return a * np.exp(-((x-mu)/(2*sigma))**2) # *1/np.sqrt(2*np.pi*sigma**2)
def gaussian_4(x,sigma=1,mu=1,a=1):    # sigma: std dev., mu: mean, a: height
    return a * np.exp(-((x-mu)/(2*sigma))**4) # *1/np.sqrt(2*np.pi*sigma**2)
# 2/3 power law without offset used for fitting    
def power2_3(x,alpha):
    return  alpha*x**2/3
# general power law without offset used for fitting
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

# presumed resonance width
def pres_res_width(m):
    return 0.15*(m**(2/3))

# to find the surrounding semimajor axes, one needs to find an individual fit!
def semimajors_surrounding(mu,a_steps):
    a_res = orbit_ratio(3,4)
    a_start = 0.8230#a_res*0.996
    a_end = 0.8245 #a_res
    A = np.linspace(a_start,a_end,a_steps)
    # a_res = a_resonant(mu)
    # c = pres_res_width(mu)
    # C = np.linspace(a_res-c, a_res, int(a_steps/2)-1, endpoint = False)
    # D = np.linspace(a_res, a_res+c, int(a_steps/2)+1, endpoint = True)
    # A = np.concatenate((C,D), axis=0)
    return A
#----------------- multi_resonance_width_fit ------------------------
#-------- fitting of all calc. resonant widths (to a power law)
def multi_resonance_width_fit(rw,Mu):
    # for initial alpha assume power 2/3 and use last data point
    alpha = rw[-1]/(Mu[-1]**(2/3))
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
#------------ lyapunov_a_multiple -----------
#------- initialize L array
#------- calculate mulitple lyap. exponents for Helga at varying semimajors a
#-------------- by seting up the simulation from anew each time
def lyapunov_a_multiple(A,t,m_Jupiter):
    L = np.zeros(len(A))                     # L: Lyapunov exponents
    for i,a in enumerate(A):
        sim_helga = setup('general_resonant', h=a, mu=m_Jupiter) #jedes mal neu initialisiert
        sim, l = lyapunov_calculator(sim_helga, t)
        L[i] = l
    return L
# -------------- single_resonance_width -----------------------
#------ get the resonance width by fitting a gauss to L(A)
#------ for starting gauss parameters assume
#------------- expectation value: middle of A (gauss centered)
#------------- sigma: half width of A
#------------- max. value: max. of L
#------ print result
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
#--------------------- multi_resonance_width (main function) ------------------
#------ initialize the arrays A, L and Res_Widths
#------ fill A, so that elements include the resonant regio for each mass ratio
#------ perform calculations for each mass ratio in Mu
#------ fit a augmented gauss to each resonance lyap(a)
#------ plot single resonance lyap(a) for each mass ratio
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
    RWfit, params = multi_resonance_width_fit(Res_Widths,Mu)
    fig, ax1 = plt.subplots()
    alphastr = r'$\alpha$'
    betastr = r'$\beta$'
    ax1.set(xscale = 'log',yscale = 'log',
            ylabel = 'Resonanzbreite (a.u.)',
            xlabel = 'Massenverhältnis $\mu = m_{Planet}/m_{Star}$',
            title = '7/12-Resonanz'
            )
    ax1.plot(Mu,Res_Widths,'o-',label='calculated')
    ax1.plot(Mu, RWfit,'r-',label=('fit: a$\mu^b$' ) )
    ax1.grid()
    ax1.legend(title = ('fitting parameters: \n a: %.3e \n b: %.3e' %(params[0],params[1])))
    # ax1.ticklabel_format(axis="y", style="sci", scilimits = (-7,-4))
    return fig
#------------- plot single resonant with fit ---------------------
def plotsingles(L,A,params,mu):
    fig = plotlyapunov_a(L,A)
    plt.title('$t = %3d \cdot  2\pi $ , $\mu = %.5f$' %(t,mu))
    fit = np.zeros(len(A))
    for i in range(len(A)):
         fit[i] = gaussian_4(A[i], *params)
    plt.plot(A, fit, 'r-',label='fit')
    plt.legend(title = 'fitting parameters: \n $\sigma$: %.4e \n $\mu$: %.4e \n $a$: %.4e' %(params[0],params[1],params[2]))
    fig.savefig('resonance_width_plots/single_lyapunovs_over_a/lyap_exp_a_variation_%.6f.png' %(mu))
#--------------------------------------------------------------


if __name__ == '__main__':
    t = 1e6  # t: time of evaluation
    # a_start = 0.6955
    # a_end = 0.6985
    # a_stepsize = (a_end - a_start)/a_steps

    a_steps = 50
    # # A: array of semimajor axises
    # A = np.linspace(a_start, a_end, a_steps)
    mass_ratio_start = -5
    mass_ratio_end = -3
    mass_ratio_steps = 20
    # Mu: array of mass ratios
    Mu = np.logspace(mass_ratio_start,mass_ratio_end,mass_ratio_steps)

    print('Masses of largest Planet (Jupiter):\n', Mu)

    reswidth = multi_resonance_width(a_steps, Mu, t)
    fig1 = plot_res_widths(reswidth, Mu)
    fig1.savefig('resonance_width_plots/res_width_over_mu')

    print('Resonance Width: \n',reswidth)
