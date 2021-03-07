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
#Compare different methods of calculating the lyapunov exponent

#---------calculate with simu-------------
def calculate_simu(h,tmax,Ntimes,t_start):
    #setting time scale
    times_simu = np.linspace(0,tmax,Ntimes)
    lyapunov_simu = np.zeros(len(times_simu))
    #current semimajoraxis
    a=np.zeros(len(times_simu))
    for i in range(Ntimes):
        t=times_simu[i] #time
        sim=visualize_orbit.setup('Helga',h) #always new setup
        if t_start>0: #if another startingtime is wanted
            sim.integrate(t_start)
        #calculate lyapunov with rebound function calculate_lyapunov
        sim,lyapunov_simu[i]=lyapunov_exponent.simu(sim, t+t_start)
        a[i]=sim.particles[2].a
    return lyapunov_simu,times_simu,a

#-------------lyapunov aus megno------
def calc_megno(h,tmax,Ntimes,t_start):
    #for comparing the calculate_lyapunov function to an manual from the megno value calculated lyapunov
    #set timescale
    times_simu = np.linspace(0,tmax,Ntimes)
    lyapunov_simu = np.zeros(len(times_simu))
    l = np.zeros(len(times_simu))#lyapunovexponent
    megno = np.zeros(len(times_simu))#megno value
    a=np.zeros(len(times_simu))#semimajoraxis
    for i in range(Ntimes):
        t=times_simu[i]
        sim=visualize_orbit.setup('Helga',h)
        if t_start>0:
            sim.integrate(t_start)
        sim,lyapunov_simu[i]=lyapunov_exponent.simu(sim, t+t_start)
        #calcalte megno value with rebound function
        megno[i]=sim.calculate_megno()
        a[i]=sim.particles[2].a
    return lyapunov_simu,times_simu,a,megno

def lyapunov_aus_megno(h,tmax,Ntimes,t_start):
    megno = np.zeros(Ntimes)
    #calculate megno values
    lyapunov_simu,times_simu,a,megno=calc_megno(h,tmax,Ntimes,t_start)
    mean_megno = np.zeros(len(times_simu))
    sum_megno = np.zeros(len(times_simu))
    slope = np.zeros(len(times_simu))
    intercept= np.zeros(len(times_simu))
    r_value= np.zeros(len(times_simu))
    p_value= np.zeros(len(times_simu))
    std_err= np.zeros(len(times_simu))
    lyapunov_megno= np.zeros(len(times_simu))
    #slope of megno value corresponds to lyapunov exponent
    for i in range(len(times_simu)-1):
        slope[i+1], intercept[i+1], r_value[i+1], p_value[i+1], std_err[i+1] = stats.linregress(times_simu[0:i+1], megno[0:i+1])
    lyapunov_megno=slope
    # lyapunov_megno=slope*2
    return lyapunov_megno, times_simu,megno, mean_megno, lyapunov_simu,sum_megno

##------------manuell--------------
def manuell_gesamt(h,sim,tmax,Ntimes):
    #direct method for calculating lyapunov, further information in lyapunov_manuell_skript.py
    abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1 = lyapunov_manuell_skript.Abstand_ln('Helga', h, 1e-10,tmax, Ntimes, 1.)
    lyapunov_manu,m,o, slope, intercept = lyapunov_manuell_skript.lyapunov_manuell_gesamt(ln_abstand, times1, x1, y1, z1, h, 1e-10, tmax, 1.)
    times=times1[m:o]
    return lyapunov_manu,times

def plotlyapunov_t(lyapunov, times, k,fig, ax1,form):
    #plot
    #k: index of calculation
    #form: plotform
    times_j=times/11.863 #time in Jovian years, 1 Jupiter year corresponds to 11,863 earth years
    lyapunov_j=lyapunov*11.863 #in 1/Jovian year
    ax1.plot(times_j,lyapunov,form, label = 'run %d' %(k+1))
    return fig

h=0.6967
tmax=0.8e6
Ntimes=[5000,10000,50000,100000]

fig, ax1 = plt.subplots(1,1)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set(ylabel ='lyapunov exponent', xlabel = 'time $t$ in Jovian year')
ax1.grid()
for k in range(1):
    #---------------manuell_alt---------------
    sim = visualize_orbit.setup('Helga',h)
    #start = time.time()
    lyapunov_manu,times=manuell_gesamt(h,sim,tmax,Ntimes[k])
    #ende = time.time()
    #print('manu-Berechnungszeit:','{:5.3f}s'.format(ende-start))
    #print('lyapunov_manu=',lyapunov_manu[-1])
    fig=plotlyapunov_t(lyapunov_manu, times, k,fig, ax1, '-')
    fig.savefig('lyapunov_t_manuell_0.696.png')

    #-------------simu------------------
    #start = time.time()
    lyapunov_simu,times_simu,a=calculate_simu(h,tmax,Ntimes[k],0)
    #ende = time.time()
    #print('simu-Berechnungszeit:','{:5.3f}s'.format(ende-start))
    #print('lyapunov_simu=', lyapunov_simu[-1])
    fig=plotlyapunov_t(lyapunov_simu, times_simu, k,fig, ax1,'-')
    fig.savefig('lyapunov_t_simu_0.696.png')

    #-------------simu späterer Startzeitpunkt------------------
    # sim = visualize_orbit.setup('Helga',h)
    # lyapunov_simu,times_simu,a=calculate_simu(h,sim,tmax,Ntimes,t_start)
    # print('lyapunov_simu_spaeter=', lyapunov_simu)
    # fig=plotlyapunov_t(lyapunov_simu, times_simu, k,fig, ax1,'o-')

plt.legend(['aus ln(Abstand) gesamt', 'aus zeitl.gem. megno'],loc='upper right')
fig.savefig('Zlyapunov_t_Vergleich_h='+str(h)+'_tmax='+str(tmax)+'_Ntimes='+str(Ntimes[k])+ '.png')

#----------------------lyapunov from megno-plot---------------------------
# fig, ax1 = plt.subplots(1,1)
# ax1.set_xscale('log')
# ax1.set_yscale('log')
# ax1.set(ylabel ='lyapunov exponent', xlabel = 'time $t$ in Jovian years')
# ax1.grid()
# lyapunov_megno, times_simu,megno,mean_megno,lyapunov_simu,sum_megno=lyapunov_aus_megno(h,tmax,Ntimes,0)
# fig=plotlyapunov_t(lyapunov_megno[1:-1], times_simu[1:-1], 0,fig, ax1,'^-')
# fig=plotlyapunov_t(lyapunov_simu[1:-1], times_simu[1:-1], 0,fig, ax1,'o-')
# plt.legend(['aus megno','aus zeitl. gemitt. megno'],loc='upper right')
# fig.savefig('lya_megno_und_simu_h='+str(h)+'_tmax='+str(tmax)+'.png')

#----------------------------megno plot-------------------------
# fig, ax1 = plt.subplots(1,1)
# ax1.set_xscale('log')
# ax1.set(ylabel ='MEGNO', xlabel = 'Zeit $t$ in Jupiterjahren')
# ax1.grid()
# -------nonresonant orbit
# lyapunov_megno, times_simu,megno,mean_megno,lyapunov_simu,sum_megno=lyapunov_aus_megno(0.696,tmax,Ntimes,0)
# fig=plotlyapunov_t(megno[1:-1], times_simu[1:-1], 0,fig, ax1,'^-')
#-------resonant orbit
# lyapunov_megno, times_simu,megno,mean_megno,lyapunov_simu,sum_megno=lyapunov_aus_megno(0.6967,tmax,Ntimes,0)
# fig=plotlyapunov_t(megno[1:-1], times_simu[1:-1], 0,fig, ax1,'^-')
# plt.legend(['h=0.696','h=0.6967'],loc='upper left')
# fig.savefig('megno_h=beides_tmax='+str(tmax)+'.png')
