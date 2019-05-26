import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import pchip

from engine import Engine
from mfp2mach import mfp2mach
import tccsty

e = Engine()

sol = e.general_explicit_map({'M_flight': 0, 'T04': 600})
for T04 in np.linspace(600,1000,10):
    sol = e.general_explicit_map({'M_flight': 0, 'T04': T04}, initial_guesses=sol.params)
sol.params

p = e.dimensionalize(**sol.params)

# Data for plot
x = np.linspace(0,1,5)
mach = np.array(
        [mfp2mach(sol.params['MFP'], e.compressor.gam),
         mfp2mach(sol.params['MFP3'], e.compressor.gam),
         mfp2mach(sol.params['MFP4'], e.compressor.gam),
         mfp2mach(sol.params['MFP5'], e.turbine.gam),
         mfp2mach(sol.params['MFP5']*e.turbine.geom['A2']/e.A8, e.turbine.gam)
        ])
pr = np.array((p.P02, p.P03, p.P04, p.P05, p.P8))/e.P01
tr = np.array((p.T02, p.T03, p.T04, p.T05, p.T8))/e.T01

#interpolate data
x_interp = np.linspace(0,1,100)
mach_interp = pchip(x, mach)(x_interp)
pr_interp = pchip(x, pr)(x_interp)
tr_interp = pchip(x, tr)(x_interp)

#disable spines
plt.rcParams.update({
        'axes.spines.left'   : False,   # display axis spines
        'axes.spines.bottom' : False,
        'axes.spines.top'    : False,
        'axes.spines.right'  : False,
        'xtick.bottom'       : False,
        'axes.facecolor'     :'none'})

textwidth = (21-5)/2.54
fig, ax = plt.subplots(1, 1, figsize=(0.7*textwidth, textwidth/2))

# Keep pressure and temperature ratios constants. 
# Adapted from https://matplotlib.org/examples/subplots_axes_and_figures/fahrenheit_celsius_scales.html

def pr2p(pr):
    """
    Converts pressure ratio to absplute pressure
    """
    return pr*e.P01*1e-5

def tr2t(tr):
    """
    Converts pressure ratio to absplute pressure
    """
    return tr*e.T01


def convert_yaxis(ax_from, ax_to, conversion_function):
    """
    Update second axis according with first axis.
    """
    y1, y2 = ax_from.get_ylim()
    ax_to.set_ylim(conversion_function(y1), conversion_function(y2))
    ax_to.figure.canvas.draw()


#set up axes
ax_ratio = ax
ax_P = ax_ratio.twinx()
ax_P.yaxis.set_label_position('left')
ax_P.yaxis.set_ticks_position('left')
ax_T = ax_ratio.twinx()
ax_ratio.callbacks.connect("ylim_changed", lambda ax_from: convert_yaxis(ax_from, ax_T, tr2t))
ax_ratio.callbacks.connect("ylim_changed", lambda ax_from: convert_yaxis(ax_from, ax_P, pr2p))

# plot
ax_ratio.plot(x, pr, 'k.')
ax_ratio.plot(x_interp, pr_interp, 'k-')
ax_ratio.plot(x, tr, 'k.')
ax_ratio.plot(x_interp, tr_interp, 'k-')
ax_ratio.set_yticks([]) #Remove y axis from axis_ratio

ax_ratio.set_xticks(x)
ax_ratio.set_xticklabels(('02','03','04','05','8'))
ax_ratio.set_xlabel('Station')

ax_ratio.set_ylim(None, 5)

ax_P.set_yticks(np.arange(1,3,0.5))
ax_P.spines['left'].set_visible(True)
ax_P.spines['left'].set_bounds(1,2.5)
ax_P.set_ylabel('$P$ [bar]')
ax_P.yaxis.set_label_coords(-0.1, 1.75, transform=ax_P.yaxis.get_ticklabels()[0].get_transform())

ax_T.set_yticks(np.arange(200,1001,100))
ax_T.spines['right'].set_visible(True)
ax_T.spines['right'].set_bounds(200,1000)
ax_T.set_ylabel('$T$ [K]')
ax_T.yaxis.set_label_coords(1.1,600, transform=ax_T.yaxis.get_ticklabels()[0].get_transform())

ax_mach_bottom = 0.6
mach_rect = list(ax_ratio.get_position().bounds)
mach_rect[1] = ax_mach_bottom
mach_rect[3] = 1-ax_mach_bottom
ax_mach = fig.add_axes(mach_rect , sharex=ax_ratio)
ax_mach.plot(x, mach, 'k.')
ax_mach.plot(x_interp, mach_interp, 'k-')
ax_mach.set_yticks(np.arange(0,1.1,0.2))
ax_mach.spines['top'].set_visible(True)
ax_mach.spines['left'].set_visible(True)
ax_mach.spines['left'].set_bounds(0,1)
ax_mach.set_ylabel('$M$')
ax_mach.tick_params(labelbottom=False)

plt.savefig('stations{}K.pdf'.format(int(sol.params['T04'])))

