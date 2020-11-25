# Import the rebound module
import rebound

# Create Simulation object
sim = rebound.Simulation()
sim.ri_sei.OMEGA =1
# Add particle to rebound

sim.add( m=1.) #sun
sim.add( m=1e-3, a=1., e=0.1 ) # Planet 1 (Jupiter)
sim.add( a=1.4, e=0.1 )       # Massless test particle

o = sim.calculate_orbits()

# Output orbits in Jacobi coordinates
for o in sim.calculate_orbits():
    print(o)

# Output orbits in Heliocentric coordinates
for o in sim.calculate_orbits(primary=sim.particles[0]):
    print(o)

# Output cartesian coordinates
for p in sim.particles:
    print(p)
