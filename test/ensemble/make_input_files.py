import os
import numpy as np
from scipy.stats import qmc
#
with open('input_base.py','r',encoding='utf-8') as file:
    data = file.readlines()


N = 1900
sampler = qmc.LatinHypercube(d=5)
sample = sampler.random(n=N) # (N,5)
# bounds: [visc_melt_vent_dense, vesic, phi_inf, log10(max_cryst_rate), log10(yield_strength_crust)]
l_bounds = np.array([70,  0.1, 0.4,  np.log10(0.05e-4), np.log10(1e3)])
u_bounds = np.array([250, 0.4, 0.59, np.log10(0.5e-4),  np.log10(1e5)])
sample = l_bounds[None,:] + (u_bounds - l_bounds)[None,:]*sample
sample[:,-1] = 10**sample[:,-1]
sample[:,-2] = 10**sample[:,-2]

for i in range(N):
    visc_melt_vent_dense, vesic, phi_inf, max_cryst_rate, yield_strength_crust = sample[i]
    data[2] = 'visc_melt_vent_dense = {:f}\n'.format(visc_melt_vent_dense)
    data[3] = 'vesic = {:f}\n'.format(vesic)
    data[5] = 'phi_inf = {:f}\n'.format(phi_inf)
    data[6] = 'max_cryst_rate = {:f}\n'.format(max_cryst_rate)
    data[7] = 'yield_strength_crust = {:f}\n'.format(yield_strength_crust)
    data[8] = 'sim_num = {}\n'.format(i)
    #
    with open('input_sim_{}.py'.format(i),'w',encoding='utf-8') as file:
        file.writelines(data)
