import rebound
import numpy as np
def setupSimulation():
    sim = rebound.Simulation()
    sim.add(m=1., hash="Sun")
    sim.add(x=0.4,vx=5., hash="Mercury")
    sim.add(a=0.7, hash="Venus")
    sim.add(a=1., hash="Earth")
    sim.move_to_com()
    return sim

sim = setupSimulation()
sim.status()
