import numpy as np
import matplotlib.pyplot as plt

from mfp2mach import mach2mfp
import tccsty

gam=1.3
mach = np.linspace(0,4,100)
mfp = mach2mfp(mach,gam)

plt.rcParams.update({
        'axes.spines.left'   : True,   # display axis spines
        'axes.spines.bottom' : True,
        'axes.spines.top'    : False,
        'axes.spines.right'  : False,
        'axes.labelpad'      : 2,
        })


plt.figure(figsize=(2,1.5))
plt.tick_params(bottom=True, top=False, left=True, right=False)
plt.plot(mach,mfp,'k', lw=tccsty.thick)
plt.xlim(0,4)
plt.ylim(0,0.6)
plt.xlabel('M')
plt.ylabel('MFP')
plt.xticks(range(5))

plt.savefig('mfp_fig.pdf')

plt.plot(2*[1],[0,mach2mfp(1,gam)], 'r--')
plt.plot([0,1],2*[mach2mfp(1,gam)], 'r--')
plt.annotate(r'$\text{{MFP}}={:.3f}, \gamma={:.1f}$'.format(mach2mfp(1,gam),gam), 
            (1.01,mach2mfp(1,gam)+1e-2),
            color='r'
            )

#plt.savefig('mfp_fig2.pdf')

