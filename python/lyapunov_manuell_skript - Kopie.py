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
import types

#a=große Halbachse Helga
def Abstand_ln(option,h,delta_h,tmax, Ntimes, delta_t):
    sim1 = visualize_orbit.setup(option,h)
    sim2 = visualize_orbit.setup(option,h+delta_h)
    sim1, times1, x1, y1, z1, vx1, vy1, vz1, a1, lyapunov_PAF = visualize_orbit.PAFintegrate(sim1, tmax, Ntimes,delta_t)
    sim2, times2, x2, y2, z2, vx2, vy2, vz2, a2, lyapunov_PAF = visualize_orbit.PAFintegrate(sim2, tmax, Ntimes,delta_t)
    #x,y,z haben so viele Zeilen wie zusätzliche Objekte!! hier 2 init megno wird jedes mal aufgerufen!
    a=sim1.particles[1].a
    abstand=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    ln_abstand=np.log(abstand)
    index=x1.shape[0]-1
    if np.any(ln_abstand[:,1]>1.5):
        index=next(i[0] for i in enumerate(ln_abstand[:,1]) if i[1]>1.5)
        tmax=times1[index]-tmax/(Ntimes*4)
        print(index)
        print(tmax)
        abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2=Abstand_ln(option,h,delta_h,tmax, Ntimes, delta_t)
    return abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2

# def Abstand_ln_aus_a(option,h,delta_h,tmax, Ntimes, delta_t):
#     sim1 = visualize_orbit.setup(option,h)
#     sim2 = visualize_orbit.setup(option,h+delta_h)
#     sim1, times1, x1, y1, z1, vx1, vy1, vz1, a1, lyapunov_PAF = visualize_orbit.PAFintegrate(sim1, tmax, Ntimes,delta_t)
#     sim2, times2, x2, y2, z2, vx2, vy2, vz2, a2, lyapunov_PAF = visualize_orbit.PAFintegrate(sim2, tmax, Ntimes,delta_t)
#     abstand_a=np.abs(a2-a1)
#     index=x1.shape[0]-1
#     a=5.204*h
#     if np.any(abstand_a[:,1]>a*0.3):
#         index=next(i for i,v in enumerate(abstand[:,1]) if abstand[i,1]>a*0.5)
#         tmax=times1[index]
#         sim1 = visualize_orbit.setup(option,h)
#         sim2 = visualize_orbit.setup(option,h+delta_h)
#         sim1, times1, x1, y1, z1, vx1, vy1, vz1, a1, lyapunov_PAF = visualize_orbit.PAFintegrate(sim1, tmax, Ntimes,delta_t)
#         sim2, times2, x2, y2, z2, vx2, vy2, vz2, a2, lyapunov_PAF = visualize_orbit.PAFintegrate(sim2, tmax, Ntimes,delta_t)
#         abstand_a=np.abs(a2-a1)
#     ln_abstand_a=np.log(abstand_a)
#     return abstand_a, ln_abstand_a, sim1, sim2, times1, x1, y1, z1, a1, a2

def lyapunov_manuell_kleiner_Zeitbereich(ln_abstand, times1, x1, y1, z1,h,delta_h,tmax, delta_t):
    slope=np.zeros(x1.shape[0])
    intercept=np.zeros(x1.shape[0])
    r_value=np.zeros(x1.shape[0])
    p_value=np.zeros(x1.shape[0])
    std_err=np.zeros(x1.shape[0])
    Ntimes=x1.shape[0]
    #fit=np.zeros(x1.shape[1])

    #slope nur für Helga!!
    #fit wird über die Kurve geschoben
    #erster punkt wird übersprungen!!!
    N=int(Ntimes/10) #zahl der gemittelten punkteF

    for n in range(N+1,x1.shape[0],1):
        slope[n], intercept[n], r_value[n], p_value[n], std_err[n] = stats.linregress(times1[n-N:n], ln_abstand[n-N:n,1])###

    lyapunov_manu=slope[N:x1.shape[0]] #die ersten einträge bleiben leer
    # times=np.zeros(x1.shape[0]-N)
    # for i,t in enumerate(times1[1:x1.shape[0]-N+1]):
    #     times[i]=mean(times1[i:i+N])
    times1=times1[N:x1.shape[0]] #slope wird letzem Zeitpunkt zugeordnet
    times=times1
    return lyapunov_manu/(2.*np.pi), times, N, slope, intercept

def lyapunov_manuell_gesamt(ln_abstand, times1, x1, y1, z1,h,delta_h,tmax, delta_t):
    slope=np.zeros(x1.shape[0])
    intercept=np.zeros(x1.shape[0])
    r_value=np.zeros(x1.shape[0])
    p_value=np.zeros(x1.shape[0])
    std_err=np.zeros(x1.shape[0])
    Ntimes=x1.shape[0]
    #slope nur für Helga #fit für immer einen Punkt mehr #erster punkt wird übersprungen
    for n in range(x1.shape[0]-2):
        slope[n], intercept[n], r_value[n], p_value[n], std_err[n] = stats.linregress(times1[1:n+2], ln_abstand[1:n+2,1])###
    slope=slope[0:x1.shape[0]-2]
    intercept=intercept[0:x1.shape[0]-2]
    lyapunov_manu=slope #die ersten einträge bleiben leer
    times=times1[1:x1.shape[0]-1] #slope wird letzem Zeitpunkt zugeordnet
    return lyapunov_manu/(2.*np.pi), times, slope, intercept

def lyapunov_manuell_gesamt_Endwert(ln_abstand, times1):
    #slope nur für Helga #fit für alle Punkte #erster punkt wird übersprungen
    slope, intercept, r_value, p_value, std_err = stats.linregress(times1[1:-1], ln_abstand[1:-1,1])###
    lyapunov_manu=slope #die ersten einträge bleiben leer
    return lyapunov_manu/(2.*np.pi)


#delta_h erstmal nur zufällig gewählter wert
def delta_h_variieren(delta_h,steps):
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
    plt.title('Variation Abstand Helga und virtuelle Helga')
    fig.savefig('lyapunov_delta_h.png')

def h_variieren(n,start,step, option ,h, delta_h, tmax, Ntimes, delta_t):
    lyapunov_manu=np.zeros(n)
    H=np.zeros(n)
    for i,h in enumerate(np.arange(start,start+n*step,step)):
        abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2 = Abstand_ln(option ,h, delta_h, tmax, Ntimes, delta_t)
        # lyapunov,times1, N=lyapunov_manuell(ln_abstand, times1, x1, y1, z1, h, delta_h, tmax, delta_t)
        # lyapunov_manu[i]=max(lyapunov[:])
        lyapunov_manu[i]=lyapunov_manuell_gesamt_Endwert(ln_abstand, times1)
        # lyapunov_manu[i]=lyapunov[-1]
        H[i]=h
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = '$a$/$a_{Jupiter}$')
    ax1.plot(np.transpose(H),lyapunov_manu[:],'o-')
    ax1.grid()
    print(H)
    print(lyapunov_manu)
    #ax1.plot(np.transpose(H),lyapunov_manu[1,:],'o-')
    plt.title('Lyapunovexponent bei unterschiedlichen Halbachsen')
    fig.savefig('Zlyapunov_a_manuell.png')

def lyapunov_ueber_aktuellem_a(a,lyapunov,N):
    A=np.zeros(a.shape[0]-N)
    #halbachse wird gemittelt über die anzahl der punkte wo die steigung ermittelt wird
    for i,t in enumerate(a[1:a.shape[0]-N+1,1]):
        A[i]=mean(a[i:i+N,1])

    #print(A)
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'a aktuell')
    ax1.plot(A,lyapunov,"o")
    #ax1.plot(a[0:a.shape[0]-N,1],lyapunov,"o")
    #ax1.plot(a[N:a.shape[0],1],lyapunov,"o")
    plt.title('Lyapunov über gemitteltem a')
    fig.savefig('lyapunov_manuell_ueber_aktuellem_a.png')

def plot_ln_absstand(ln_abstand,times,slope,intercept,h,tmax,Ntimes):
        fig, ax1 = plt.subplots(1,1)
        ax1.set(ylabel ='ln(Abstand)', xlabel = 'Zeit t in Jupiterjahren')
        times=times/11.863 #Zeit in Jupiterjahren
        ax1.plot(times,ln_abstand[:,1],'s-')
        ax1.grid()
        if np.any(slope!=0) or np.any(intercept!=0):
            m=1
            if type(slope) is tuple:
                m=slope.shape[0]
            for i in range(m):
                fit=np.zeros(len(times))
                fit=slope[i]*times*11.863+intercept[i]
                ax1.plot(times,fit,'-')
        plt.title('t_max = '+str(tmax)+', Ntimes = '+str(Ntimes))
        fig.savefig(os.path.join('ZAbstand_ln_' +'_'+ str(h) +'_'+ str(tmax) +'_' + str(Ntimes)+ '.png'))


# h=0.696
tmax=1e6
Ntimes=50
#
# abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2 = Abstand_ln('Helga', h, 1e-10,tmax, 100, 1.)
# # lyapunov_manu,times, slope, intercept = lyapunov_manuell_gesamt(ln_abstand, times1, x1, y1, z1, h, 1e-10, tmax, 1.)
# lyapunov_manu,times,N, slope, intercept = lyapunov_manuell_kleiner_Zeitbereich(ln_abstand, times1, x1, y1, z1, h, 1e-10, tmax, 1.)
#
# plot_ln_absstand(ln_abstand,times1, [slope[len(intercept)-1]], [intercept[len(intercept)-1]],h,tmax,Ntimes)

h_variieren(31,0.696,0.0001,'Helga',0.697, 1e-10, tmax, Ntimes, 1.)
# lyapunov_manu,times,N = lyapunov_manuell(ln_abstand, times1, x1, y1, z1, h, 1e-10, tmax, 1.)
# fig=lyapunov_exponent.plotlyapunov_t(lyapunov_manu, times, 0)
# abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2 = Abstand_ln_aus_a('Helga', h, 1e-10,tmax, 100, 1.)
# lyapunov_manu,times,N = lyapunov_manuell(ln_abstand, times1, x1, y1, z1, h, 1e-10, tmax, 1.)
# fig=lyapunov_exponent.plotlyapunov_t(lyapunov_manu, times, 0)
#fig.savefig('lyapunov_t_manuell_a_0.6967.png')
#



# abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1, a1, a2 = Abstand_ln('Helga', h, 1e-10,tmax, Ntimes, 1.)
#
