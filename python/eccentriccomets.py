import rebound
import numpy as np
sim = rebound.Simulation()
np.random.seed(42)

#integrator options
sim.integrator = "mercurius"
sim.dt = 1
sim.testparticle_type = 1

#collision and boundary options
sim.collision = "direct"
sim.collision_resolve = "merge"
sim.collision_resolve_keep_sorted = 1
sim.boundary = "open"
boxsize = 200.
sim.configure_box(boxsize)
sim.track_energy_offset = 1

#simulation time
tmax = 1e4

#massive bodies
sim.add(m=1., r=0.005)                     # Sun
a_neptune = 30.05
sim.add(m=5e-5,r=2e-4,a=a_neptune,e=0.01)  # Neptune

# semi-active bodies
n_comets = 100
a = np.random.random(n_comets)*10 + a_neptune
e = np.random.random(n_comets)*0.009 + 0.99
inc = np.random.random(n_comets)*np.pi/2.
m = 1e-10
r = 1e-7

for i in xrange(0,n_comets):
    rand = np.random.random()*2*np.pi
    sim.add(m=m, r=r, a=a[i], e=e[i], inc=inc[i], Omega=0, omega=rand, f=rand)

sim.move_to_com()
E0 = sim.calculate_energy()
#%matplotlib inline
fig = rebound.OrbitPlot(sim,Narc=300)

sim.getWidget(size=(500,300),scale=1.8*a_neptune)
sim.integrate(tmax)
dE = abs((sim.calculate_energy() - E0)/E0)
print(dE)
