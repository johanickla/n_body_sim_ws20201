# finding_resonances_test.py
import numpy as np
import matplotlib.pyplot as plt
#---------------------------------------------------------
a_resonants = np.array([0.697956,0.697910,0.697700,0.697498,0.697155,0.696579,0.6959,0.6955,0.6952])
Mu = np.array([0.0001,0.00015849,0.00025119,.000398,.000631,.001,0.00178,0.00316,0.00562])

# def a_resonant(mu):
#     b = 0.6981
#     a = -0.0012/0.001
#     a0 = 5.204
#     c0 = 0.02
#     # a_res = b*(1-c0*mu**(1/3))
#     c1 = 0.15
#     # a_res = b*(1-c1*mu**(2/3))
#     c2 = 0.06
#     # a_res = b*(1-c2*mu**(1/2))
#     c4 = 0.01
#     # a_res = b*(1-c4*mu**(1/4))
#     a_res = b+a*mu
#     return  a_res
# B = a_resonant(Mu)
#
# fig, ax1 = plt.subplots()
# ax1.set(xscale = 'linear',
#         ylabel = 'Resonanzstelle',
#         xlabel = 'Massenverhältnis $\mu = m_{Planet}/m_{Star}$'
#         )
# ax1.plot(Mu,a_resonants,'-o')
# ax1.plot(Mu,B,'-')
# plt.show()
#-----------------------------------------------------
def semimajors_surrounding(mu,a_steps):
    a_res = a_resonant(mu)
    C = np.linspace(a_res, a_res+0.001, int(a_steps/2)+1, endpoint = True)
    D = np.linspace(a_res-0.001, a_res, int(a_steps/2)-1, endpoint = False)
    A = np.concatenate((C,D), axis=0)
    return A

# print(semimajors_surrounding(0.696, 20))
#-------------------------------------------------
# def pres_res_widths(m):
#     return 0.16*(m**(2/3))
# M = np.array([0.0001,0.00056,0.001,0.00178,0.00316])
# C_goal = [4*1e-4,8*1e-4,15*1e-4,20*1e-4,30*1e-4]
# C_calc =  pres_res_widths(M)
#
# fig2, ax1 = plt.subplots()
# ax1.set(xscale = 'linear',
#         ylabel = 'notwendige Breite',
#         xlabel = 'Massenverhältnis $\mu = m_{Planet}/m_{Star}$'
#         )
# ax1.plot(M,C_goal,'-o')
# ax1.plot(M,C_calc,'-')
# plt.show()
#-----------------------------------------------------------------------
k = 0.0034
print(k)
print('%e' %(k))
print('%.2e' %(k))
print(('der', r'$\alpha > \beta$')
