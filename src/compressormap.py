import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import numpy as np
import itertools

from turbomachine import gridmap
from compressor import Compressor
from mfp2mach import mach2mfp

def plot_compressor_map(ax, P0_min=1, P0_max=4, samples=25, c=Compressor()):

    MFP_choke = mach2mfp(1,c.gam)
    sol_cr = c.general_explicit_map({'MFP1': MFP_choke, 'MFP2': MFP_choke})
    sol_P1 = c.general_explicit_map({'MFP2': MFP_choke, 'P0_ratio': P0_min})

    def MFP_max(P0_ratio):
        MFP = np.interp(P0_ratio, 
                       [sol_P1.params['P0_ratio'], sol_cr.params['P0_ratio']],
                       [sol_P1.params['MFP1'], sol_cr.params['MFP1']])
        return MFP
  
    P0_ratio = np.linspace(P0_min, P0_max, samples)
    P0_grid = np.empty((samples, samples))
    MFP_grid = np.empty_like(P0_grid)
    for i, p in enumerate(P0_ratio):
        MFP = np.linspace(1e-6,MFP_max(p), samples)
        MFP_grid[:, i] = MFP
        P0_grid[:, i] = p

    params = gridmap(c, MFP_grid, P0_grid, 'MFP1')
    params['eff'] = (c.gam-1)/c.gam*np.log(P0_grid)/np.log(params['T0_ratio'])

    ax.plot(MFP_grid, P0_grid, 'k,')
    CS2 = ax.contour(MFP_grid, P0_grid, params['Mb'], levels=np.arange(0,2,0.1), 
                        colors='k', linewidths=1.5)
    ax.clabel(CS2, CS2.levels, fmt='%.1f')
    CS2 = ax.contour(MFP_grid, P0_grid, params['eff'], colors='k', linewidths=0.5,
                      levels=[0.5,0.8,0.9,0.95])
    ax.clabel(CS2, CS2.levels, fmt='%.2f')
    ax.plot([MFP_choke, MFP_choke, sol_P1.params['MFP1']],
             [P0_max, sol_cr.params['P0_ratio'], P0_min], 'k--', label='choke limit')
    ax.set_xlim((0,0.6))
    ax.set_ylim((1,4))
    ax.set_xlabel("MFP")
    ax.set_ylabel(r"$\frac{P_{03}}{P_{02}}$")

plt.figure(figsize=(6,8),dpi=150)
plot_compressor_map(plt.gca())
plt.savefig('map.pdf', dpi=300)

plt.show()
