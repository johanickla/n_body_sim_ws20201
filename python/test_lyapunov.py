import rebound
import numpy as np
import warnings
import matplotlib.pyplot as plt
import visualize_orbit

h = 0.696

def test_1():
    sim1 = visualize_orbit.setup('Helga', h)
    print('Number of particles before "init_megno" :', sim1.N)
    sim1.init_megno()
    print('Number of particles after "init_megno" :', sim1.N)
    sim1.status()
    sim1.integrate(10)
    l1 = sim1.calculate_lyapunov()
    # sim1.status()

    sim2 = visualize_orbit.setup('Helga', h)
    sim2.init_megno()
    sim2.integrate(5)
    # sim2.status()
    sim2.integrate(5)
    l2 = sim2.calculate_lyapunov()
    # sim2.status()
    print('L_exp 1:  ', round(l1,3), '   L_exp 2:  ', round(l2,3))

def test_2():
    # a = [16,64,256,1024,4098]
    a = [2,4,8,16,32,64]
    Lyapunov=[]
    S=[]
    T = []
    for i in range(len(a)):
        times = np.linspace(0, 100, a[i])
        l = np.zeros(len(times))
        sim = visualize_orbit.setup('Helga', h)
        sim.integrator = "whfast"
        sim.dt = 0.01
        sim.init_megno()
        for k in range(a[i]) :
            time = times[k]
            sim.integrate(time)
            exp = round(sim.calculate_lyapunov(),4)
            l[k] = exp
            # sim.status()
        Lyapunov.append(l)
        T.append(times)
        s = round(sum(Lyapunov[i]),3)
        S.append(s)
    # print('Lyapunov:' ,Lyapunov)
    # print('rowwise sum over Lyapunov array:', S)
    # print('Times: ', T)
    # sim.status()

    fig, ax1 = plt.subplots(1,1,figsize=(10,6))
    ax1.grid()
    for i in range(len(a)):
        ax1.plot(T[i],Lyapunov[i],'o-', label = '%d' %a[i])
        # ax1.set_xscale('log')
        # ax1.set_yscale('log')
        ax1.legend()
    fig.savefig('test_lyapunov.png')

def test_3():
    sim = visualize_orbit.setup('Helga', h)
    sim.init_megno()
    for time in np.linspace(0,100,10):
        sim.integrate(time)
        exp = round(sim.calculate_lyapunov(),4)
        print(exp)
        sim.status()


if __name__ == "__main__":
    # test_1()
    test_2()
    # test_3()
