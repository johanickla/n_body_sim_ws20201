import rebound
import numpy as np
import warnings
import matplotlib
import matplotlib.pyplot as plt
import visualize_orbit

h = 0.696
sim = visualize_orbit.setup('Helga', h)

l = sim.calculate_lyapunov()
print(l)
