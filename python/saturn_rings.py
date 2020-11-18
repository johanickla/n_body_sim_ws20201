import rebound
import numpy as np
sim = rebound.Simulation()
OMEGA = 0.00013143527     # [1/s]
sim.ri_sei.OMEGA = OMEGA

# Finally, let us define the surface density
# of the ring and the particle density.

surface_density = 400.    # kg/m^2
particle_density = 400.   # kg/m^3

# The gravitational constant in SI units is

sim.G = 6.67428e-11       # N m^2 / kg^2

# We choose a timestep of 1/1000th of the orbital period.

sim.dt = 1e-3*2.*np.pi/OMEGA

# We enable gravitational softening to smear out any
# potential numerical artefacts at very small scales.

sim.softening = 0.2       # [m]

# Next up, we configure the simulation box. By default REBOUND used
# no boundary conditions, but here we have shear periodic
# boundaries and a finite simulation domain,
# so we need to let REBOUND know about the simulation boxsize
# (note that it is significantly smaller than a,
# so our local approximation is very good. In this example we’ll work in SI units.

boxsize = 200.            # [m]
sim.configure_box(boxsize)

# Because we have shear-periodic boundary conditions, we use ghost boxes
# to simulate the gravity of neighbouring ring patches. The more ghostboxes we use, the smoother the gravitational force accross the boundary. Here, two layers of ghost boxes in the x and y direction are enough (this is a total of 24 ghost boxes). We don’t need ghost boxes in the z direction because a rings is a two dimensional system.

sim.configure_ghostboxes(2,2,0)

# We can now setup which REBOUND modules we want to use for our simulation.
# Besides the SEI integrator and the shear-periodic boundary conditions
# mentioned above, we select the tree modules for both gravity and collisions.
# This speeds up the code from O(N2) to O(Nlog(N)) for
# large numbers of particles N


sim.integrator = "sei"
sim.boundary   = "shear"
sim.gravity    = "tree"
sim.collision  = "tree"
sim.collision_resolve = "hardsphere"

# When two ring particles collide, they loose energy during
# their the bounce.
# We here use a velocity dependent Bridges et. al. coefficient of restitution.
#  It is implemented as a python function
# (a C implementation would be faster!). We let REBOUND know which function
# we want to use by setting the coefficient_of_restitution
# function pointer in the simulation instance.

def cor_bridges(r, v):
        eps = 0.32*pow(abs(v)*100.,-0.234)
        if eps>1.:
            eps=1.
        if eps<0.:
            eps=0.
        return eps
sim.coefficient_of_restitution = cor_bridges

# To initialize the particles, we will draw random
# numbers from a power law distribution.

def powerlaw(slope, min_v, max_v):
    y = np.random.uniform()
    pow_max = pow(max_v, slope+1.)
    pow_min = pow(min_v, slope+1.)
    return pow((pow_max-pow_min)*y + pow_min, 1./(slope+1.))

# Now we can finally add particles to REBOUND.
# Note that we initialize particles so that they have
# initially no velocity relative to the mean shear flow.

total_mass = 0.
while total_mass < surface_density*(boxsize**2):
    radius = powerlaw(slope=-3, min_v=1, max_v=4)  # [m]
    mass = particle_density*4./3.*np.pi*(radius**3)
    x = np.random.uniform(low=-boxsize/2., high=boxsize/2.)
    sim.add(
        m=mass,
        r=radius,
        x=x,
        y=np.random.uniform(low=-boxsize/2., high=boxsize/2.),
        z=np.random.normal(),
        vx = 0.,
        vy = -3./2.*x*OMEGA,
        vz = 0.)
    total_mass += mass

# To see what is going on in our simulation,
# we create a function to plot the current positions of
# particles and call it once to visualise the initial conditions.

#%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib.patches as patches
def plotParticles(sim):
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(111,aspect='equal')
    ax.set_ylabel("radial coordinate [m]")
    ax.set_xlabel("azimuthal coordinate [m]")
    ax.set_ylim(-boxsize/2.,boxsize/2.)
    ax.set_xlim(-boxsize/2.,boxsize/2.)
    for i, p in enumerate(sim.particles):
        circ = patches.Circle((p.y, p.x), p.r, facecolor='darkgray', edgecolor='black')
        ax.add_patch(circ)
    plt.show()

#-------------

plotParticles(sim)
plt.show()
plt.savefig('saturnringe.png')
