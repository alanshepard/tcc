import matplotlib.pyplot as plt
import numpy as np

from compressor import Compressor
from turbine import Turbine, TurbineExtendedMap
from engine import Engine

fig, ax = plt.subplots(1,2)

Compressor().plot_map(ax[0])
TurbineExtendedMap().plot_map(ax[1])

e = Engine()
e.initial_guess = e.general_explicit_map({'M_flight': 0, 'MFP': 0.3}).params

wline = e.working_line(np.linspace(600,1033,10))

ax[0].plot([p['MFP'] for p in wline], [p['P0_ratio_c'] for p in wline], '-ko')
ax[1].plot([p['MFP4']*p['Mb_t'] for p in wline], [1/p['P0_ratio_t'] for p in wline], '-ko')

for xy, label in [((p['MFP'], p['P0_ratio_c']), p['T04']) for p in wline]:
    ax[0].annotate(label, xy)

plt.show()