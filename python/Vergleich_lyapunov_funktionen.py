import lyapunov_manuell_skript
import rebound
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import visualize_orbit
import lyapunov_exponent
from scipy import stats
import os
import time
from statistics import mean


#Vergleich von calculate und manuell
##calculate mit PAF

h=0.6967
sim = visualize_orbit.setup('Helga',h)
tmax=1e4
Ntimes=20
###calculate_PAF
###lyapunov_PAF[i]=(arbsim.calculate_lyapunov()/(2.*np.pi))
arbsim, times_PAF, x, y, z, a, lyapunov_PAF=visualize_orbit.PAFintegrate(sim ,tmax, Ntimes,1.)
#print(lyapunov_PAF)
#print(times_PAF)


#calculate mit simu
def calculate_simu(h,sim,tmax,Ntimes):
    times_simu = np.linspace(0,tmax,Ntimes)
    sim = visualize_orbit.setup('Helga',0.696)
    lyapunov_simu = np.zeros(len(times_simu))
    for i in range(Ntimes):
        t=times_simu[i]
        sim = visualize_orbit.setup('Helga',h)
        sim,lyapunov_simu[i]=lyapunov_exponent.simu(sim, t)
    return lyapunov_simu,times_simu
#print(lyapunov_simu)


##manuell
def manuell(h,sim,tmax,Ntimes):
    abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2 = lyapunov_manuell_skript.Abstand_ln('Helga', h, 1e-10,tmax, Ntimes, 1.)
    lyapunov_manu,times,N = lyapunov_manuell_skript.lyapunov_manuell(ln_abstand, times1, x1, y1, z1, h, 1e-10, tmax, 1.)
    return lyapunov_manu,times
# print(lyapunov_manu)
# print(times)

def plotlyapunov_t(lyapunov, times, k,fig, ax1,form):
    times_j=times/11.863 #Zeit in Jupiterjahren, 1Jupiterjahr sind 11,863 erdjahre
    lyapunov_j=lyapunov*11.863 #in 1/jupiterjahrn
    ax1.plot(times_j,lyapunov,form, label = 'run %d' %(k+1))
    return fig

h=0.6967
sim = visualize_orbit.setup('Helga',h)
tmax=1e4
Ntimes=20

fig, ax1 = plt.subplots(1,1)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'Zeit $t$ in Jupiterjahren')
ax1.grid()
for k in range(7):
    sim = visualize_orbit.setup('Helga',h)
    h=0.6967
    lyapunov_manu,times=manuell(h,sim,tmax,Ntimes)
    print(lyapunov_manu)
    fig=plotlyapunov_t(lyapunov_manu, times, k,fig, ax1, 'o-')
    #fig.savefig('lyapunov_t_manuell_0.696.png')
    sim = visualize_orbit.setup('Helga',h)
    arbsim, times_PAF, x, y, z, a, lyapunov_PAF=visualize_orbit.PAFintegrate(sim ,tmax, Ntimes,1.)
    print(lyapunov_PAF)
    fig=plotlyapunov_t(lyapunov_PAF, times_PAF, k,fig, ax1,'v-')
    #fig.savefig('lyapunov_t_PAF_0.696.png')
    sim = visualize_orbit.setup('Helga',h)
    lyapunov_simu,times_simu=calculate_simu(h,sim,tmax,Ntimes)
    print(lyapunov_simu)
    fig=plotlyapunov_t(lyapunov_simu, times_simu, k,fig, ax1,'^-')
    #fig.savefig('lyapunov_t_simu_0.696.png')
    #lyapunov_exponent.lyapunov_t_multiple(10, 0, tmax, Ntimes)
ax1.legend()
plt.legend(loc='upper right')
plt.subplots_adjust(right=0.85)
fig.savefig('lyapunov_t_Vergleich_0.6967_tmax=1e6.png')

# lyapunov_manuell_skript.h_variieren(31 ,0.696,0.0001,'Helga',0.697, 1e-10, tmax, Ntimes, 1.)
# lyapunov_exponent.lyapunov_a_multiple(0.696,0.699,0.0001,tmax)
