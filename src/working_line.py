import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
plt.rcParams['text.usetex'] = True
import numpy as np
import itertools

from compressor import Compressor
from turbine import Turbine, TurbineExtendedMap
from engine import Engine

fig, ax = plt.subplots(1,2, figsize=(7.3*np.sqrt(2), 7.3), sharey=True)

Compressor().plot_map(ax[0], P0_max=3, samples=200)
TurbineExtendedMap().plot_map(ax[1], extent=(0,1.2,1,3),samples=200)

e = Engine()
e.initial_guess = e.general_explicit_map({'M_flight': 0, 'MFP': 0.3}).params

wline = e.working_line(np.linspace(600,1033,10))

ax[0].plot([p['MFP'] for p in wline], [p['P0_ratio_c'] for p in wline], 'k.-', linewidth=1.2)
ax[1].plot([p['MFP4']*p['Mb_t'] for p in wline], [1/p['P0_ratio_t'] for p in wline], 'k.-', linewidth=1.2)

for xy, label in [((p['MFP'], p['P0_ratio_c']), p['T04']) for p in wline]:
    ax[0].annotate("{:.0f}".format(label), xy, horizontalalignment='right',
                                               verticalalignment='bottom')
# for xy, label in [((p['MFP4']*p['Mb_t'], 1/p['P0_ratio_t']), p['T04']) for p in wline]:
#     ax[1].annotate("{:.0f}K".format(label), xy, horizontalalignment='left',
#             verticalalignment='top')

axins = zoomed_inset_axes(ax[1], 5.5, loc=1)  # zoom = 6

# sub region of the original image
x1, x2, y1, y2 = 0.05, 0.14,1.11,1.31
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
TurbineExtendedMap().plot_map(axins)
axins.plot([p['MFP4']*p['Mb_t'] for p in wline], [1/p['P0_ratio_t'] for p in wline], 'k.-', linewidth=1.2)
for (xy, label), ha in zip(
        [((p['MFP4']*p['Mb_t'], 1/p['P0_ratio_t']), p['T04']) for p in wline],
        itertools.cycle(['left','right'])
        ):
    axins.annotate("{:.0f}".format(label), xy,horizontalalignment=ha,
            verticalalignment='center')

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax[1], axins, loc1=2, loc2=4, fc="none", ec="0.5")
axins.set_xlabel('')
axins.set_ylabel('')

ax[0].set_ylim(1,3)
ax[0].set_title('Compressor')
ax[1].set_ylim(1,3)
ax[1].set_title('Turbine')

plt.savefig('wline.pdf', dpi=600)

plt.show()
