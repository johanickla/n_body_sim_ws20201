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


#lyapunov_manuell mit virtuelle Helga für jeden Zeitpunkt neu berechnen
def Abstand_ln(h,delta_h, N, delta_t,t,t_berechnung, sim1):
    #sim.status()
    #sim1, lyapunov_simu=lyapunov_exponent.simu(sim1 ,t)
    #Virtuelle Helga
    #sim1.status()
    #sim2 = rebound.Simulation()
    #sim2=sim1
    sim1=visualize_orbit.setup('Helga',h)
    sim2=visualize_orbit.setup('Helga',h)
    #sim2.add( sim1.particles[0] )                         # sun
    #sim2.add( sim1.particles[1]) # Planet 1 (Jupiter) Dtaen von Wiki
    # sim2.add(m=sim1.particles[2].m, a=sim1.particles[2].a+delta_h*5.204, M=sim1.particles[2].M, omega=sim1.particles[2].omega, e=sim1.particles[2].e, inc=sim1.particles[2].inc)
    #sim2.add(m=sim1.particles[2].m, a=sim1.particles[2].a, M=sim1.particles[2].M, omega=sim1.particles[2].omega, e=sim1.particles[2].e, inc=sim1.particles[2].inc,Omega=sim1.particles[2].Omega)
    ###################################################################################################################primary,m,a,anom,e,omega,inv,Omega,MEAN
    #sim2.add( sim1.particles[2])
    #sim2.status()
    #+delta_h*5.204

    sim1, times1, x1, y1, z1, a1, lyapunov_PAF = visualize_orbit.PAFintegrate(sim1, t_berechnung, N, delta_t)
    sim2, times2, x2, y2, z2, a2, lyapunov_PAF = visualize_orbit.PAFintegrate(sim2, t_berechnung, N, delta_t)
    #x,y,z haben so viele Zeilen wie zusätzliche Objekte!! hier 2
    # abstand=np.sqrt((x1[:,0:2]-x2[:,0:2])**2+(y1[:,0:2]-y2[:,0:2])**2+(z1[:,0:2]-z2[:,0:2])**2)
    abstand=np.sqrt((x1[:,2]-x2[:,2])**2+(y1[:,2]-y2[:,2])**2+(z1[:,2]-z2[:,2])**2)
    print(x1.shape)
    print(abstand.shape)
    ln_abstand=np.log(abstand)
    #ln_abstand=abstand
    return abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2

def lyapunov_manuell(ln_abstand, times1, x1, y1, z1,h,delta_h,t_berechnung, delta_t):
    slope=0
    intercept=0
    r_value=0
    p_value=0
    std_err=0
    #N=shape(times1)[]

    #slope nur für Helga!!
    slope, intercept, r_value, p_value, std_err = stats.linregress(times1, ln_abstand[:,1])###
    lyapunov_manu=slope
    return lyapunov_manu/(2.*np.pi)

def manuell_t(h,sim,tmax,Ntimes):
    delta_h=1e-10
    delta_t=1.
    N=20 #Anzahl der zeitpunkte aus der lyapunov berechnet wird
    t_berechnung=50 #zeitpunkt bis zu dem lyapunov berechnet wird
    times_manu = np.linspace(0,tmax,Ntimes)
    lyapunov_manu = np.zeros(len(times_manu))
    for i in range(Ntimes):
        t=times_manu[i]
        sim = visualize_orbit.setup('Helga',h)
        abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2=lyapunov_manuell_neu.Abstand_ln(h,delta_h, N, delta_t,t,t_berechnung,sim)
        lyapunov_manu[i]=lyapunov_manuell_neu.lyapunov_manuell(ln_abstand, times1, x1, y1, z1,h,delta_h,t_berechnung, delta_t)
    return lyapunov_manu, times_manu

h=0.696
t_berechnung=1e3 #beeinflusst zeitbereich an den die Steigung gelegt wird
t=0 #zeitpunkt a la simu
N=20 #Anzahl der zeitpunkte aus der lyapunov berechnet wird
delta_h=1e-10

sim = visualize_orbit.setup('Helga',h)
abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2=Abstand_ln(h,delta_h, N,1.,t,t_berechnung,sim)
# lyapunov_manu=lyapunov_manuell(ln_abstand, times1, x1, y1, z1,h,delta_h,t_berechnung, 1.)
# print(lyapunov_manu)

fig, ax1 = plt.subplots(1,1)
ax1.set(ylabel ='ln(Abstand)', xlabel = 'Zeit t')
ax1.plot(times1,ln_abstand,'o-')
ax1.grid()
fig.savefig('lyapunov_manuell_neu_lnAbstand_'+str(h)+'_'+str(t)+'.png')
