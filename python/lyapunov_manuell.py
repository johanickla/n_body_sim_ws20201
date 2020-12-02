import rebound
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import visualize_orbit
from scipy import stats
import os
import time

def lyapunov_manuell(option,h,delta_h,tmax, Ntimes,delta_t):
    sim1 = visualize_orbit.setup(option,h)
    sim2 = visualize_orbit.setup(option,h+delta_h)
    sim1, times1, x1, y1, z1 = visualize_orbit.PAFintegrate(sim1, tmax, Ntimes,delta_t)
    sim2, times2, x2, y2, z2 = visualize_orbit.PAFintegrate(sim2, tmax, Ntimes,delta_t)
    #x,y,z haben so viele Zeilen wie zusätzliche Objekte!! hier 2

    abstand=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    ln_abstand=np.log(abstand)

    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='ln(Abstand)', xlabel = 'Zeit t')
    ax1.plot(times1,ln_abstand,'o-')

    slope=np.zeros(x1.shape[1])
    intercept=np.zeros(x1.shape[1])
    r_value=np.zeros(x1.shape[1])
    p_value=np.zeros(x1.shape[1])
    std_err=np.zeros(x1.shape[1])
    #fit=np.zeros(x1.shape[1])

    for n in range(0,x1.shape[1],1):
        slope[n], intercept[n], r_value[n], p_value[n], std_err[n] = stats.linregress(times1[int(Ntimes/2):int(Ntimes-1)], ln_abstand[int(Ntimes/2):int(Ntimes-1),n])###
        print('Objekt',n+1, slope[n], intercept[n], r_value[n], p_value[n], std_err[n])
        fit=slope[n]*times1+intercept[n]
        ax1.plot(times1,fit,'--')

    ax1.grid()
    #fügt datum in dateiname ein, damit ich nichts versehentlich überspeicher...
    #fig.savefig(os.path.join('Abstand_ln_' + time.strftime('%H%M%S_%d%m%Y') + '.png'))
    fig.savefig(os.path.join('Abstand_ln_' + str(option) +'_'+ str(h) +'_'+ str(delta_h) +'_'+ str(tmax) +'_' + str(Ntimes) +'_'+ str(delta_t) +'_' + '.png'))
    lyapunov_manu=slope
    return lyapunov_manu

#delta_h erstmal nur zufällig gewählter wert
def delta_h_variieren():
    delta_h=1e-11
    steps=5
    lyapunov_manu=np.zeros((2,steps))
    delta_h_n=np.zeros(steps)
    for n in range(0,steps,1):
        delta_h=delta_h*10
        lyapunov_manu[:,n]=lyapunov_manuell('Helga',0.696, delta_h, 1e5, 100, 1.)
        delta_h_n[n]=delta_h
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'delta_h')

    ax1.plot(np.transpose(delta_h_n),lyapunov_manu[0,:],'o-')
    ax1.plot(np.transpose(delta_h_n),lyapunov_manu[1,:],'o-')
    fig.savefig('lyapunov_delta_h.png')

def h_variieren(n,start,step):
    delta_h=1e-10
    lyapunov_manu=np.zeros((2,n))
    H=np.zeros(n)
    for i,h in enumerate(np.arange(start,start+n*step,step)):
        lyapunov_manu[:,i]=lyapunov_manuell('Helga', h, delta_h, 1e7, 100, 1.)
        H[i]=h
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'h')
    ax1.plot(np.transpose(H),lyapunov_manu[0,:],'o-')
    ax1.plot(np.transpose(H),lyapunov_manu[1,:],'o-')
    fig.savefig('lyapunov_a_manuell.png')

#lyapunov_manuell('Helga',0.696, 1e-10, 1e7, 100, 1.)
h_variieren(31,0.696,0.0001)
