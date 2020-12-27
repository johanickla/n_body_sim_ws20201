import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm

import numpy as np
from lyapunov_exponent import hill_sphere_radius
#--------------- plotlyapunov_t --------------------
# lyapunov und zeit
def plotlyapunov_t(lyapunov, times, k):
    times_j=times/11.863 #Zeit in Jupiterjahren, 1Jupiterjahr sind 11,863 erdjahre
    lyapunov_j=lyapunov*11.863 #in 1/jupiterjahrn
    fig, ax1 = plt.subplots(1,1)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'Zeit $t$ in Jupiterjahren')
    ax1.plot(times_j,lyapunov,'o-', label = 'run %d' %(k+1))
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
    # fig.savefig('plots/lyapunov_a.png')
    return fig
#-------------- plotlyapunov_e ---------------------
# lyapunov-exponent und Exzentrizitäten
def plotlyapunov_e(l,e):
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel = 'Lyapunov-Exponent', xlabel = 'Exzentrizität $e$')
    ax1.plot(e,l,'o-')
    ax1.grid()
    # fig.savefig('plots/lyapunov_a.png')
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
def plotlyapunov_ma_surface(m,a,l,h = None):
    # m = np.log10(m)
    fig, ax1 = plt.subplots(1, figsize=(8,6))
# caculate and plot Hill sphere radius
    if h != None:
        m_Helga = 3e-13
        e_Helga = .086345
        a_Jupiter = 5.204
        a_Helga = h*5.204

        r_Hill_plus = np.zeros(len(m))
        r_Hill_minus = np.zeros(len(m))
        for i in range(len(m)):
            r_Hill_plus[i] = (hill_sphere_radius(a_Helga, m[i]*m_Helga, 1,e_Helga) + a_Helga)/a_Jupiter
            r_Hill_minus[i] = (-hill_sphere_radius(a_Helga, m[i]*m_Helga, 1,e_Helga) + a_Helga)/a_Jupiter
        ax1.plot(m,r_Hill_plus,'-',linewidth=2, color = 'white',label = 'Hill sphere radius')
        ax1.plot(m,r_Hill_minus,'-',linewidth=2, color = 'white')
        ax1.legend()
# reformat m,a,l
    ax1.set_ylim(a[0],a[len(a)-1])
    m, a = np.meshgrid(m, a)
    l = np.transpose(l)
    # print('shape m:', np.shape(m))
    # print('shape a:', np.shape(a))
    # print('shape l:', np.shape(l))
# Plot the surface.
    # mylevels = np.logspace(-10,2,10)
    # c = ax1.contourf(m, a, l,mylevels, norm = LogNorm(vmin=l.min(), vmax=l.max()),
    #                 cmap='plasma')
    c = ax1.pcolor(m, a, l, norm = LogNorm(vmin=l.min(), vmax=l.max()),
                     shading = 'auto', cmap='plasma')
    ax1.set(xlabel = '$m$/$M_{Helga}$', xscale = 'log',
                ylabel= '$a$/$a_{Jupiter}$')
    fig.colorbar(c, ax=ax1)
    ax1.set(xscale = 'log')
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
