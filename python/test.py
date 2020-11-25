import numpy as np
import rebound
sim = rebound.Simulation()
sim.add(m=1)
sim.add(a=1)
fig, ax_main = rebound.OrbitPlot(sim)
fig.savefig("image.png") # save figure to file
fig.show() # show figure on screen
