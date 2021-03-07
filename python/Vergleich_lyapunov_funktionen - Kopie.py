import lyapunov_manuell_skript
import lyapunov_manuell_neu
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
###----------------calculate_PAF--------------------------
###lyapunov_PAF[i]=(arbsim.calculate_lyapunov()/(2.*np.pi))
arbsim, times_PAF, x, y, z,vx, vy, vz, a, lyapunov_PAF=visualize_orbit.PAFintegrate(sim ,tmax, Ntimes,1.)
#print(lyapunov_PAF)
#print(times_PAF)


#---------calculate mit simu-------------
def calculate_simu(h,tmax,Ntimes,t_start):
    times_simu = np.linspace(0,tmax,Ntimes)
    #sim = visualize_orbit.setup('Helga',h)
    lyapunov_simu = np.zeros(len(times_simu))
    a=np.zeros(len(times_simu))
    for i in range(Ntimes):
        t=times_simu[i]
        # sim = sim_setup # jedes mal neues setup
        sim=visualize_orbit.setup('Helga',h)
        #sim.integrate(t)
        if t_start>0:
            # sim,l=lyapunov_exponent.simu(sim,t_start)
            sim.integrate(t_start)
        sim,lyapunov_simu[i]=lyapunov_exponent.simu(sim, t+t_start)
        # sim.status()
        #sim_setup.status()
        a[i]=sim.particles[2].a
    return lyapunov_simu,times_simu,a
#print(lyapunov_simu)

#-------------lyapunov aus megno------
def calc_megno(h,tmax,Ntimes,t_start):
    times_simu = np.linspace(0,tmax,Ntimes)
    #sim = visualize_orbit.setup('Helga',h)
    lyapunov_simu = np.zeros(len(times_simu))
    l = np.zeros(len(times_simu))#test lya
    megno = np.zeros(len(times_simu))
    a=np.zeros(len(times_simu))
    for i in range(Ntimes):
        t=times_simu[i]
        # sim = sim_setup # jedes mal neues setup
        sim=visualize_orbit.setup('Helga',h)
        #sim.integrate(t)
        if t_start>0:
            # sim,l=lyapunov_exponent.simu(sim,t_start)
            sim.integrate(t_start)
        sim,lyapunov_simu[i]=lyapunov_exponent.simu(sim, t+t_start)
        megno[i]=sim.calculate_megno()
        # sim.status()
        #sim_setup.status()
        a[i]=sim.particles[2].a
    return lyapunov_simu,times_simu,a,megno

def lyapunov_aus_megno(h,tmax,Ntimes,t_start):
    megno = np.zeros(Ntimes)
    lyapunov_simu,times_simu,a,megno=calc_megno(h,tmax,Ntimes,t_start)
    mean_megno = np.zeros(len(times_simu))
    sum_megno = np.zeros(len(times_simu))
    # sum_megno[0]=megno[0]
    # mean_megno[0]=megno[0]
    slope = np.zeros(len(times_simu))
    intercept= np.zeros(len(times_simu))
    r_value= np.zeros(len(times_simu))
    p_value= np.zeros(len(times_simu))
    std_err= np.zeros(len(times_simu))
    lyapunov_megno= np.zeros(len(times_simu))
    # for i in range(len(times_simu)-1):
    #     sum_megno[i+1]=sum_megno[i]+megno[i+1]
    #     mean_megno[i+1]=sum_megno[i+1]/times_simu[i+1]
    #     #slope[i+1], intercept[i+1], r_value[i+1], p_value[i+1], std_err[i+1] = stats.linregress(times_simu[0:i+1], mean_megno[0:i+1])
    # slope, intercept, r_value, p_value, std_err = stats.linregress(times_simu, mean_megno)
    for i in range(len(times_simu)-1):
        slope[i+1], intercept[i+1], r_value[i+1], p_value[i+1], std_err[i+1] = stats.linregress(times_simu[0:i+1], megno[0:i+1])
    lyapunov_megno=slope
    # lyapunov_megno=slope*2
    return lyapunov_megno, times_simu,megno, mean_megno, lyapunov_simu,sum_megno

##------------manuell--------------
def manuell_gesamt(h,sim,tmax,Ntimes):
    abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2 = lyapunov_manuell_skript.Abstand_ln('Helga', h, 1e-10,tmax, Ntimes, 1.)
    lyapunov_manu,times, slope, intercept = lyapunov_manuell_skript.lyapunov_manuell_gesamt(ln_abstand, times1, x1, y1, z1, h, 1e-10, tmax, 1.)
    return lyapunov_manu,times
# print(lyapunov_manu)
# print(times)

#---------manuell neu--------------
def manuell_neu(h,sim,tmax,Ntimes,t_berechnung):
    delta_h=1e-10
    delta_t=1.
    N=20 #Anzahl der zeitpunkte aus der lyapunov berechnet wird
     #t_berechnung zeitpunkt bis zu dem lyapunov berechnet wird
    times_manu = np.linspace(0,tmax,Ntimes)
    lyapunov_manu = np.zeros(len(times_manu))
    for i in range(Ntimes):
        t=times_manu[i]
        sim = visualize_orbit.setup('Helga',h)
        abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2=lyapunov_manuell_neu.Abstand_ln(h,delta_h, N, delta_t,t,t_berechnung,sim)
        lyapunov_manu[i]=lyapunov_manuell_neu.lyapunov_manuell(ln_abstand, times1, x1, y1, z1,h,delta_h,t_berechnung, delta_t)
    return lyapunov_manu, times_manu

def plotlyapunov_t(lyapunov, times, k,fig, ax1,form):
    times_j=times/11.863 #Zeit in Jupiterjahren, 1Jupiterjahr sind 11,863 erdjahre
    lyapunov_j=lyapunov*11.863 #in 1/jupiterjahrn
    ax1.plot(times_j,lyapunov,form, label = 'run %d' %(k+1))
    return fig

h=0.6967
#sim = visualize_orbit.setup('Helga',h)
tmax=1e6
Ntimes=50
# t_start=tmax
# t_berechnung=1000

fig, ax1 = plt.subplots(1,1)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'Zeit $t$ in Jupiterjahren')
ax1.grid()
for k in range(1):
    #---------------manuell_alt---------------
    sim = visualize_orbit.setup('Helga',h)
    lyapunov_manu,times=manuell_gesamt(h,sim,tmax,Ntimes)
    #print(lyapunov_manu)
    fig=plotlyapunov_t(lyapunov_manu, times, k,fig, ax1, 'o-')
    #fig.savefig('lyapunov_t_manuell_0.696.png')
    #------------manuell_neu----------------
    # sim = visualize_orbit.setup('Helga',h[k])
    # lyapunov_manu, times_manu=manuell_neu(h[k],sim,tmax,Ntimes,t_berechnung)
    # print('lyapunov_manu=', lyapunov_manu)
    # fig=plotlyapunov_t(lyapunov_manu, times_manu, k,fig, ax1, 'o-')
    #-------------PAF---------------------
    # sim = visualize_orbit.setup('Helga',h)
    # arbsim, times_PAF, x, y, z,vx, vy, vz, a, lyapunov_PAF=visualize_orbit.PAFintegrate(sim ,tmax, Ntimes,1.)
    # print(lyapunov_PAF)
    # fig=plotlyapunov_t(lyapunov_PAF, times_PAF, k,fig, ax1,'v-')
    #fig.savefig('lyapunov_t_PAF_0.696.png')
    #-------------simu------------------
    lyapunov_simu,times_simu,a=calculate_simu(h,tmax,Ntimes,0)
    print('lyapunov_simu=', lyapunov_simu)
    fig=plotlyapunov_t(lyapunov_simu, times_simu, k,fig, ax1,'^-')
    #fig.savefig('lyapunov_t_simu_0.696.png')

    #-------------simu späterer Startzeitpunkt------------------
    # sim = visualize_orbit.setup('Helga',h)
    # lyapunov_simu,times_simu,a=calculate_simu(h,sim,tmax,Ntimes,t_start)
    # print('lyapunov_simu_spaeter=', lyapunov_simu)
    # fig=plotlyapunov_t(lyapunov_simu, times_simu, k,fig, ax1,'o-')

plt.legend(['aus ln(Abstand) gesamt', 'aus zeitl.gem. megno'],loc='upper right')
fig.savefig('Zlyapunov_t_Vergleich_h='+str(h)+'_tmax='+str(tmax)+ '.png')
# fig.savefig('lyapunov_t_Vergleich_neu_h='+str(h[0])+str(h[1])+'_tmax='+str(tmax)+'_t_berechnung='+str(t_berechnung)+ '.png')
# fig.savefig('lyapunov_t_simu_mit verschiedenen_Startzeiten_tmax='+str(tmax)+'_'+str(h)+ '.png')

# lyapunov_manuell_skript.h_variieren(31 ,0.696,0.0001,'Helga',0.697, 1e-10, tmax, Ntimes, 1.)
# lyapunov_exponent.lyapunov_a_multiple(0.696,0.699,0.0001,tmax)
#
#
# lyapunov_simu,times_simu,a=calculate_simu(h,sim,tmax,Ntimes,0)
#
# fig, ax1 = plt.subplots(1,1)
# ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'a aktuell')
# ax1.plot(a,lyapunov_simu,"o")
# plt.title('Lyapunov über a mit tmax='+str(tmax))
# fig.savefig('lyapunov_simu_ueber_aktuellem_a_sim=sim_setup_'+str(h)+'.png')
# h=0.6967
# fig, ax1 = plt.subplots(1,1)
# ax1.set_xscale('log')
# ax1.set_yscale('log')
# ax1.set(ylabel ='Lyapunovexponent', xlabel = 'Zeit $t$ in Jupiterjahren')
# ax1.grid()
# lyapunov_megno, times_simu,megno,mean_megno,lyapunov_simu,sum_megno=lyapunov_aus_megno(h,tmax,Ntimes,0)
# # print('mean_megno',mean_megno)
# # print('lyapunov_megno',lyapunov_megno)
# # print('lyapunov_simu',lyapunov_simu)
# fig=plotlyapunov_t(lyapunov_megno[1:-1], times_simu[1:-1], 0,fig, ax1,'^-')
# fig=plotlyapunov_t(lyapunov_simu[1:-1], times_simu[1:-1], 0,fig, ax1,'o-')
# # lyapunov_simu,times_simu,a,megno=megno(0.6967,tmax,Ntimes,0)
# # fig=plotlyapunov_t(megno, times_simu, 0,fig, ax1,'^-')
# plt.legend(['aus megno','aus zeitl. gemitt. megno'],loc='upper right')
#
# fig.savefig('lya_megno_und_simu_h='+str(h)+'_tmax='+str(tmax)+'.png')
#
# fig, ax1 = plt.subplots(1,1)
# ax1.set_xscale('log')
# ax1.set(ylabel ='MEGNO', xlabel = 'Zeit $t$ in Jupiterjahren')
# ax1.grid()
# lyapunov_megno, times_simu,megno,mean_megno,lyapunov_simu,sum_megno=lyapunov_aus_megno(0.696,tmax,Ntimes,0)
# fig=plotlyapunov_t(megno[1:-1], times_simu[1:-1], 0,fig, ax1,'^-')
# lyapunov_megno, times_simu,megno,mean_megno,lyapunov_simu,sum_megno=lyapunov_aus_megno(0.6967,tmax,Ntimes,0)
# fig=plotlyapunov_t(megno[1:-1], times_simu[1:-1], 0,fig, ax1,'^-')
# plt.legend(['h=0.696','h=0.6967'],loc='upper left')
# fig.savefig('megno_h=beides_tmax='+str(tmax)+'.png')
