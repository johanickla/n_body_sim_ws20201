import sys
sys.path.append('C:/Users/hannes/Documents/GitHub/n_body_sim_ws20201/python')
import numpy as np
RW = np.zeros((1,8))
n = 8
molp= [1.00000000e-05, 1.29154967e-05,  2.78255940e-05, 3.59381366e-05, 4.64158883e-05, 5.99484250e-05, 7.74263683e-05, 1.00000000e-04]
RW[0,:]=[2.84774637e-05, 4.75368311e-05, 3.26286186e-05, 4.26223922e-05, 4.40024877e-05, 6.71254385e-05, 3.64044614e-05, 3.45124728e-05]
RWA = np.mean(RW,axis=0)
print(RWA)

from resonance_width import plot_res_widths,multi_resonance_width_fit
fig = plot_res_widths(RWA,molp)
fig.savefig('resonance_widths_3_4.png')
