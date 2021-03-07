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


def Abstand_ln(option,h,delta_h,tmax, Ntimes, delta_t):
    #option: for setup, for example 'Helga', 'resonant' or 'fictional'
    #h: ratio of orbits of Helga and Jupiter, for example 0.6967
    #delta_h: initial ratio difference of Helga and its shadow particle
    #tmax: integraion time maximum in earth years
    #Ntimes: number of time steps
    #delta_t: time steps size for the integrator

    #set initial condition for Helga and its shadow particle
    sim1 = visualize_orbit.setup(option,h)
    sim2 = visualize_orbit.setup(option,h+delta_h)
    #integrates both developments for equal times and time steps, returns its locations
    sim1, times1, x1, y1, z1 = visualize_orbit.PAFintegrate(sim1, tmax, Ntimes,delta_t)
    sim2, times2, x2, y2, z2 = visualize_orbit.PAFintegrate(sim2, tmax, Ntimes,delta_t)
    #abstand = distance between all particles and its shadow particles (without sun)
    abstand=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    ln_abstand=np.log(abstand)
    return abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1

def lyapunov_manuell_gesamt(ln_abstand, times1, x1, y1, z1,h,delta_h,tmax, delta_t):
    #derive slope from ln_abstand
    m=next(i[0] for i in enumerate(times1) if i[1]>1e4) #Initial area is skipped
    o=len(ln_abstand)-1 #last fitted time
    #if saturation occurs, the last fitted time step is the first time when ln_abstand is larger than a threshhold
    if np.any(ln_abstand[:,1]>1.1):
        o=next(i[0] for i in enumerate(ln_abstand[:,1]) if i[1]>1.1)

    slope=np.zeros(o-m)
    intercept=np.zeros(o-m)
    r_value=np.zeros(o-m)
    p_value=np.zeros(o-m)
    std_err=np.zeros(o-m)
    #slope only for Helga
    #fit for one time step more each iteration
    for n in range(o-m):
        slope[n], intercept[n], r_value[n], p_value[n], std_err[n] = stats.linregress(times1[m:n+m+1], ln_abstand[m:n+m+1,1])
    lyapunov_manu=slope
    return lyapunov_manu/(2.*np.pi), m,o , slope, intercept

def lyapunov_manuell_gesamt_Endwert(ln_abstand, times1):
    #slope only for Helga if only the total fit is wanted
    m=next(i[0] for i in enumerate(times1) if i[1]>1e4) #Initial area is skipped
    o=len(ln_abstand)-1
    #fit for all points in time until the saturation
    if np.any(ln_abstand[:,1]>1.1):
        o=next(i[0] for i in enumerate(ln_abstand[:,1]) if i[1]>1.1)
    slope, intercept, r_value, p_value, std_err = stats.linregress(times1[m:o], ln_abstand[m:o,1])###
    lyapunov_manu=slope
    return lyapunov_manu/(2.*np.pi)

def lyapunov_manuell_kleiner_Zeitbereich(ln_abstand, times1, x1, y1, z1,h,delta_h,tmax, delta_t):
    #slope only for Helga #fit is pushed across the curve
    slope=np.zeros(x1.shape[0])
    intercept=np.zeros(x1.shape[0])
    r_value=np.zeros(x1.shape[0])
    p_value=np.zeros(x1.shape[0])
    std_err=np.zeros(x1.shape[0])
    Ntimes=x1.shape[0]

    N=int(Ntimes/10) #number of averaged points
    #first point in time is skipped
    for n in range(N+1,x1.shape[0],1):
        slope[n], intercept[n], r_value[n], p_value[n], std_err[n] = stats.linregress(times1[n-N:n], ln_abstand[n-N:n,1])###

    lyapunov_manu=slope[N:x1.shape[0]] #the first entries stay empty
    times1=times1[N:x1.shape[0]] #slope is assigned to the last time of the fitting
    times=times1
    return lyapunov_manu/(2.*np.pi), times, N, slope, intercept

def delta_h_variieren(option,h,delta_h,tmax, Ntimes, delta_t,steps):
    #initial distance is changed, the resulting lyapunov exponent should stay nearly the same
    #steps: how often delta_h should be changed
    #since the saturation threshhold and the initial area threshhold stays constant, less time steps can be used for the fit if the initial distance gets larger
    lyapunov_manu=np.zeros(steps)
    delta_h_n=np.zeros(steps)
    delta_h=delta_h/10
    for n in range(0,steps,1):
        delta_h=delta_h*3
        abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1 = Abstand_ln(option,h,delta_h,tmax, Ntimes, delta_t)
        #only the end lyapunov value is necessary
        lyapunov_manu[n] = lyapunov_manuell_gesamt_Endwert(ln_abstand, times1)
        delta_h_n[n] = delta_h
    #plot
    fig, ax1 = plt.subplots(1,1)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set(ylabel ='lyapunov exponent', xlabel = 'delta_h')
    ax1.plot(np.transpose(delta_h_n),lyapunov_manu,'o-')
    plt.title('Variation distance Helga and shadow Helga')
    fig.savefig('lyapunov_delta_h.png')

def h_variieren(n,start,step, option ,h, delta_h, tmax, Ntimes, delta_t):
    #to find the resonant orbit sizes, it is necessary to compare different orbits and its lyapunov exponents
    #n: how many orbits should be calculated
    #start: starting semimajoraxis orbit ratio
    #step: deifference of two considered orbit ratios
    lyapunov_manu=np.zeros(n)
    H=np.zeros(n)
    for i,h in enumerate(np.arange(start,start+n*step,step)):
        abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1 = Abstand_ln(option ,h, delta_h, tmax, Ntimes, delta_t)
        lyapunov_manu[i]=lyapunov_manuell_gesamt_Endwert(ln_abstand, times1)
        H[i]=h
    #plot
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='lyapunov exponent', xlabel = '$a$/$a_{Jupiter}$')
    ax1.plot(np.transpose(H),lyapunov_manu[:],'o-')
    ax1.grid()
    plt.title('lyapunov exponent for different semimajoraxes')
    fig.savefig('lyapunov_a_manuell.png')

def lyapunov_ueber_aktuellem_a(a,lyapunov,N):
    #lyapunov over current semimajoraxis
    A=np.zeros(a.shape[0]-N)
    #Semimajoraxis is averaged over the number of points where the slope is determined
    for i,t in enumerate(a[1:a.shape[0]-N+1,1]):
        A[i]=mean(a[i:i+N,1])
    #plot
    fig, ax1 = plt.subplots(1,1)
    ax1.set(ylabel ='Lyapunov-Exponent', xlabel = 'a aktuell')
    ax1.plot(A,lyapunov,"o")
    plt.title('Lyapunov über gemitteltem a')
    fig.savefig('lyapunov_manuell_ueber_aktuellem_a.png')

def plot_ln_abstand(ln_abstand,times,slope,intercept,h,tmax,Ntimes,indexstart,indexend):
        #plot of ln_abstand and its fit
        #indexstart: index for the start of the fit or end of the skipped initial region
        #indexend: index for the end of the fit or the start of the saturation
        fig, ax1 = plt.subplots(1,1)
        ax1.set(ylabel ='ln(d)', xlabel = 'time t in Jovian years')
        times=times/11.863 #time in Jovian years
        ax1.grid()

        if np.any(slope!=0) or np.any(intercept!=0):
            m=1
            if type(slope) is tuple:
                m=slope.shape[0]
            for i in range(2*m):
                fit=np.zeros(len(times))
                fit=slope[i]*times*11.863+intercept[i]
                ax1.plot(times[indexstart[i]:indexend[i]],fit[indexstart[i]:indexend[i]],'b-')
                ax1.plot(times,ln_abstand[:,i],'-')
                ax1.plot(times[indexstart[i]],ln_abstand[indexstart[i],i],'ro')
                ax1.plot(times[indexend[i]],ln_abstand[indexend[i],i],'ro')
        plt.title('Development of ln(distance)')
        #ax1.legend(['nonresonant orbit','resonant orbit'],loc='lower right')
        fig.savefig(os.path.join('ZAbstand_ln_' +'_'+ str(h) +'_'+ str(tmax) +'_' + str(Ntimes)+ '.png'))


h=0.6967
tmax=0.8e6
Ntimes=400
#Example:
abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1 = Abstand_ln('Helga', h, 1e-10,tmax, Ntimes, 1.)
lyapunov_manu,m1,o1, slope1, intercept1 = lyapunov_manuell_gesamt(ln_abstand, times1, x1, y1, z1, h, 1e-10, tmax, 1.)
#or
h_variieren(31,0.696,0.0001,'Helga',0.697, 1e-10, tmax, Ntimes, 1.)

# abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1 = Abstand_ln('Helga', 0.696, 1e-10,tmax, Ntimes, 1.)
# lyapunov_manu,m1,o1, slope1, intercept1 = lyapunov_manuell_gesamt(ln_abstand, times1, x1, y1, z1, 0.696, 1e-10, tmax, 1.)
# LN=np.zeros((len(times1),2))
# LN[:,0]=ln_abstand[:,1]
# abstand, ln_abstand, sim1, sim2, times1, x1, y1, z1 = Abstand_ln('Helga', 0.6967, 1e-10,tmax, Ntimes, 1.)
# lyapunov_manu,m2,o2, slope2, intercept2 = lyapunov_manuell_gesamt(ln_abstand, times1, x1, y1, z1, 0.6967, 1e-10, tmax, 1.)
# LN[:,1]=ln_abstand[:,1]
#
# plot_ln_abstand(LN,times1, [slope1[len(intercept1)-1],slope2[len(intercept2)-1]], [intercept1[len(intercept1)-1],intercept2[len(intercept2)-1]],h,tmax,Ntimes,[m1,m2],[o1,o2])
