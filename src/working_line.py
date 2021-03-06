import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import itertools

from compressor import Compressor
from turbine import Turbine, TurbineExtendedMap
from engine import Engine
import tccsty

samples = 144

fig, ax = plt.subplots(1,2, figsize=((297-50)/25.4, (210-50-30)/25.4), sharey=True)

Compressor().plot_map(ax[0], P0_max=5, samples=samples)
TurbineExtendedMap().plot_map(ax[1], extent=(0,1.2,1,3),samples=samples)

# plot work line
e = Engine()
e.initial_guess = e.general_explicit_map({'M_flight': 0, 'MFP': 0.3}).params

T04 = np.arange(600,1001,50)
wline = e.working_line(T04)

ax[0].plot([p['MFP'] for p in wline], [p['P0_ratio_c'] for p in wline], 'k.-', linewidth=tccsty.verythick, label='35mm nozzle')
ax[1].plot([p['MFP4']*p['Mb_t'] for p in wline], [1/p['P0_ratio_t'] for p in wline], 'k.-', linewidth=tccsty.verythick)

for xy, label in [((p['MFP'], p['P0_ratio_c']), p['T04']) for p in wline]:
    ax[0].annotate("{:.0f}".format(label), xy, horizontalalignment='right',
                                               verticalalignment='bottom')
# for xy, label in [((p['MFP4']*p['Mb_t'], 1/p['P0_ratio_t']), p['T04']) for p in wline]:
#     ax[1].annotate("{:.0f}K".format(label), xy, horizontalalignment='left',
#             verticalalignment='top')

#45mm nozzle 
wline45 = pd.read_csv('wline45mm.csv')
T04 = np.arange(200,1001,10)
T04 = T04[(T04<=wline45['T04'].max()) & (T04>=wline45['T04'].min())]

P0_ratio_c = interp1d(wline45['T04'], wline45['P0_ratio_c'])(T04)
MFP_c = interp1d(wline45['T04'], wline45['MFP'])(T04)
P0_ratio_t = interp1d(wline45['T04'], wline45['P0_ratio_t'])(T04)
MbMFP_t = interp1d(wline45['T04'], wline45['MFP4']*wline45['Mb_t'])(T04)

ax[0].plot(MFP_c, P0_ratio_c, 'k-^', linewidth=tccsty.verythick, label='45mm nozzle')
ax[1].plot(MbMFP_t, 1/P0_ratio_t, 'k-^', linewidth=tccsty.verythick)

for xy, label in zip(zip(MFP_c+0.01, P0_ratio_c), T04):
    ax[0].annotate("{:.0f}".format(label), xy, horizontalalignment='left',
                                               verticalalignment='top')
#for xy, label in zip(zip(MbMFP_t, 1/P0_ratio_t), T04):
#    ax[1].annotate("{:.0f}".format(label), xy, horizontalalignment='left',
#                                               verticalalignment='top')


ax[0].legend(loc='upper left')
axins = zoomed_inset_axes(ax[1], 5.7, loc=1)  # zoom = 6

# sub region of the original image
x1, x2, y1, y2 = 0.05, 0.145,1.11,1.31
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
TurbineExtendedMap().plot_map(axins)
axins.plot([p['MFP4']*p['Mb_t'] for p in wline], [1/p['P0_ratio_t'] for p in wline], 'k.-', linewidth=tccsty.verythick)
for (x, y, label), ha in zip(
        [(p['MFP4']*p['Mb_t'], 1/p['P0_ratio_t'], p['T04']) for p in wline],
        itertools.cycle(['left','right'])
        ):
    axins.text(x+ (2e-3 if ha is 'left' else -2e-3),y,"{:.0f}".format(label),horizontalalignment=ha,
            verticalalignment='center')

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax[1], axins, loc1=2, loc2=4, facecolor="none", edgecolor="0.5", linewidth=plt.rcParams['axes.linewidth'])
axins.set_xlabel('')
axins.set_ylabel('')

ax[0].set_ylim(1,3)
ax[0].set_xlim(0,0.62)
ax[0].set_title('Compressor')
ax[1].set_ylim(1,3)
ax[1].set_title('Turbine')

plt.savefig('wline.pdf')
