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

h=0.6967
sim = visualize_orbit.setup('Helga',h)

tmax=1e3 #beeinflusst zeitbereich an den die Steigung gelegt wird
t=0 #zeitpunkt a la simu
Ntimes=20
delta_h=1e-10

#lyapunov_manuell mit virtuelle Helga für jeden Zeitpunkt neu berechnen
#simu
sim1, times, x, y, z, a, lyapunov_PAF=visualize_orbit.PAFintegrate(sim ,t, 1,1.)
sim1.status()
sim2=sim1
particles = sim2.particles
a_neu=sim2.particles[1].a+2
#sim2.particles[1].a=a_neu

#particles[2].a = particles[2].a+delta_h
sim2 = rebound.Simulation()
sim2.add( sim1.particles[0] )                         # sun
sim2.add( sim1.particles[1]) # Planet 1 (Jupiter) Dtaen von Wiki
#sim2.add(m=sim1.particles[2].m, a=a_neu, M=38.41426796877275, omega=0.257, e=.08634521111588543 ,inc=0.0768 )  # Helga NASA
sim2.add(m=sim1.particles[2].m, a=a_neu,)
print(sim1.particles[2].a)
print(sim2.particles[2].a)

def Abstand_ln(option,h,delta_h,tmax, Ntimes, delta_t):
    sim1 = arbsim1
    sim2 = arbsim2
    #sim1.status()
    sim1, times1, x1, y1, z1, a1, lyapunov_PAF = visualize_orbit.PAFintegrate(sim1, tmax, Ntimes,delta_t)
    sim2, times2, x2, y2, z2, a2, lyapunov_PAF = visualize_orbit.PAFintegrate(sim2, tmax, Ntimes,delta_t)
    #x,y,z haben so viele Zeilen wie zusätzliche Objekte!! hier 2
    a=sim1.particles[1].a
    abstand=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    index=x1.shape[0]-1
    if np.any(abstand[:,1]>a*0.3):
        index=next(i for i,v in enumerate(abstand[:,1]) if abstand[i,1]>a*0.5)
        print(index)
        #print(x1.shape)
        tmax=times1[index]
        sim1 = visualize_orbit.setup(option,h)
        sim2 = visualize_orbit.setup(option,h+delta_h)
        sim1, times1, x1, y1, z1, a1, lyapunov_PAF = visualize_orbit.PAFintegrate(sim1, tmax, Ntimes,delta_t)
        sim2, times2, x2, y2, z2, a2, lyapunov_PAF = visualize_orbit.PAFintegrate(sim2, tmax, Ntimes,delta_t)
        abstand=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    #sim1.status()
    #print(abstand)
    ln_abstand=np.log(abstand)
    #ln_abstand=abstand
    return abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2
