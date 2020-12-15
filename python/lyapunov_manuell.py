import rebound
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import visualize_orbit
import lyapunov_exponent
from scipy import stats
import os
import time

#a=große Halbachse Helga
def Abstand_ln(option,h,delta_h,tmax, Ntimes,delta_t):
    sim1 = visualize_orbit.setup(option,h)
    sim2 = visualize_orbit.setup(option,h+delta_h)
    #sim1.status()
    a=5.204*h
    sim1, times1, x1, y1, z1 = visualize_orbit.PAFintegrate(sim1, tmax, Ntimes,delta_t)
    sim2, times2, x2, y2, z2 = visualize_orbit.PAFintegrate(sim2, tmax, Ntimes,delta_t)
    #x,y,z haben so viele Zeilen wie zusätzliche Objekte!! hier 2

    abstand=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    index=x1.shape[0]-1
    if np.any(abstand[:,1]>a*0.3):
        index=next(i for i,v in enumerate(abstand[:,1]) if abstand[i,1]>a*0.5)
        print(index)
        #print(x1.shape)
        tmax=times1[index]
        sim1 = visualize_orbit.setup(option,h)
        sim2 = visualize_orbit.setup(option,h+delta_h)
        sim1, times1, x1, y1, z1 = visualize_orbit.PAFintegrate(sim1, tmax, Ntimes,delta_t)
        sim2, times2, x2, y2, z2 = visualize_orbit.PAFintegrate(sim2, tmax, Ntimes,delta_t)
        abstand=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    # abstand=abstand[0:index,:]
    # times1=times1[0:index]
    # x1=x1[0:index,:]
    # y1=y1[0:index,:]
    # z1=z1[0:index,:]
    ln_abstand=np.log(abstand)
    #print(abstand)
    #print(n)
    #print(x1.shape)
    return abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1

def lyapunov_manuell(ln_abstand, times1, x1, y1, z1,h,delta_h,tmax, Ntimes,delta_t):
    # fig, ax1 = plt.subplots(1,1)
    # ax1.set(ylabel ='ln(Abstand)', xlabel = 'Zeit t')
    # ax1.plot(times1,ln_abstand[:,1],'o-')
    # ax1.grid()

    slope=np.zeros(x1.shape[0])
    intercept=np.zeros(x1.shape[0])
    r_value=np.zeros(x1.shape[0])
    p_value=np.zeros(x1.shape[0])
    std_err=np.zeros(x1.shape[0])
    Ntimes=x1.shape[0]
    #fit=np.zeros(x1.shape[1])

    #slope nur für Helga!!
    #fit für immer einen Punkt mehr
    for n in range(int(Ntimes/5)+1,x1.shape[0],1):
        slope[n], intercept[n], r_value[n], p_value[n], std_err[n] = stats.linregress(times1[n-int(Ntimes/10):n], ln_abstand[n-int(Ntimes/10):n,1])###
        #print('Zeit=',times1[n], slope[n], intercept[n], r_value[n], p_value[n], std_err[n])
        #fit=slope[n]*times1[n-int(Ntimes/10):n]+intercept[n]
        #ax1.plot(times1[n-int(Ntimes/10):n],fit,'--')

    #fügt datum in dateiname ein, damit ich nichts versehentlich überspeicher...
    #fig.savefig(os.path.join('Abstand_ln_' + time.strftime('%H%M%S_%d%m%Y') + '.png'))
    #fig.savefig(os.path.join('Abstand_ln_' +'_'+ str(h) +'_'+ str(delta_h) +'_'+ str(tmax) +'_' + str(Ntimes) +'_'+ str(delta_t) +'_' + '.png'))
    lyapunov_manu=slope[int(Ntimes/10):x1.shape[0]]
    times1=times1[int(Ntimes/10):x1.shape[0]]
    #print(lyapunov_manu)
    return lyapunov_manu, times1


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
    lyapunov_manu=np.zeros(n)
    H=np.zeros(n)
    for i,h in enumerate(np.arange(start,start+n*step,step)):
        abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1 = Abstand_ln('Helga',h, 1e-10,0.2e7, 100, 1.)
        lyapunov,times1=lyapunov_manuell(ln_abstand, times1, x1, y1, z1, h, 1e-10, 0.2e7, 100, 1.)
        lyapunov_manu[i]=max(lyapunov[:])
        H[i]=h
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'h')
    ax1.plot(np.transpose(H),lyapunov_manu[:],'o-')
    #ax1.plot(np.transpose(H),lyapunov_manu[1,:],'o-')
    fig.savefig('lyapunov_a_manuell.png')

# abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1 = Abstand_ln('Helga',0.697, 1e-10,0.2e7, 100, 1.)
# lyapunov_manu,times1 = lyapunov_manuell(ln_abstand, times1, x1, y1, z1,0.697, 1e-10, 0.2e7, 100, 1.)
# fig=lyapunov_exponent.plotlyapunov_t(lyapunov_manu, times1, 0)
# fig.savefig('lyapunov_t_manuell.png')
h_variieren(31,0.696,0.0001)
