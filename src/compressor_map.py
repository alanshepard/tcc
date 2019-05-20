import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import numpy as np
import itertools

from turbomachine import gridmap
from compressor import Compressor
from mfp2mach import mach2mfp


# In[2]:

c = Compressor()
c.implicit_map(MFP1=-0.1, MFP2=0.1, Mb=0.1, P0_ratio=1., T0_ratio=1)


# In[3]:

c.general_explicit_map(params={'MFP2': 0.3, 'MFP1': 0.3})#, initial_guesses={'MFP2':0.2, 'Mb':1, 'T0_ratio':1})


# In[4]:

# Define map corners
P0_max = 4
P0_min = 1
MFP_choke = mach2mfp(1,c.gam)
sol_cr = c.general_explicit_map({'MFP1': MFP_choke, 'MFP2': MFP_choke})
sol_P1 = c.general_explicit_map({'MFP2': MFP_choke, 'P0_ratio': P0_min})
plt.plot([0, sol_P1.params['MFP1'], MFP_choke], [1, 1, sol_cr.params['P0_ratio']])


# In[5]:

# Define grid
samples= 25
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

plt.plot(MFP_grid, P0_grid)
plt.plot(MFP_grid[0,0], P0_grid[0,0], 'b*')
plt.plot(MFP_grid[-1,-1], P0_grid[-1,-1], 'r*')
# Calculate coarse grid
# Interpolate to fine grid
# Calculate fine grid


# In[6]:

params = gridmap(c, MFP_grid, P0_grid, 'MFP1')


# In[7]:

eff = (c.gam-1)/c.gam*np.log(P0_grid)/np.log(params['T0_ratio'])


# In[8]:

plt.imshow(success, origin='lower')


# In[9]:

plt.imshow(nfev, origin='lower')


# In[10]:

plt.imshow(MFP2, origin='lower')


# In[ ]:

plt.imshow(Mb, origin='lower')


# In[ ]:

plt.imshow(T0_ratio, origin='lower')


# In[ ]:

plt.imshow(eff, origin='lower')


# In[ ]:




# In[ ]:

plt.figure(figsize=(6,8),dpi=150)
# plt.plot(MFP_grid, P0_grid, 'k,')
CS2 = plt.contour(MFP_grid, P0_grid, Mb, levels=np.arange(0,2,0.1), 
                    colors='k', linewidths=1.5)
plt.clabel(CS2, CS2.levels, fmt='%.1f')
CS2 = plt.contour(MFP_grid, P0_grid, eff, colors='k', linewidths=0.5,
                  levels=[0.5,0.8,0.9,0.95])
plt.clabel(CS2, CS2.levels, fmt='%.2f')
plt.plot([MFP_choke, MFP_choke, sol_P1.params['MFP1']],
         [P0_max, sol_cr.params['P0_ratio'], P0_min], 'k--', label='choke limit')
plt.xlim((0,0.6))
plt.ylim((1,4))
plt.xlabel("MFP")
plt.ylabel(r"$\frac{P_{03}}{P_{02}}$")
# plt.legend(("Mb", "eta_p"))
plt.legend()
plt.savefig('map.pdf', dpi=300)

ax = plt.gca()


# In[ ]:

# np.savetxt('MFP_grid.gz', MFP_grid)
# np.savetxt('P0_grid.gz', P0_grid)
# np.savetxt('nfev_grid.gz', nfev)
# np.savetxt('success_grid.gz', success)
# np.savetxt('MFP2_grid.gz', MFP2)
# np.savetxt('eff_grid.gz', eff)
# np.savetxt('T0_grid.gz', T0_ratio)


# In[ ]:

from turbine import Turbine


# In[ ]:

t = Turbine()


# In[ ]:

t.implicit_map(0.2,0.2,0.4,1,1)


# In[ ]:

t.general_explicit_map(params={'P0_ratio': 1, 'MFP1': 0.5})


# In[ ]:

from engine import Engine


# In[ ]:

e = Engine()


# In[ ]:

e.implicit_map(**e.__class__.DEFAULT_PARAMS)


# In[ ]:

sol=e.general_explicit_map(params={'M_flight':0, 'T04': 1500})
print(sol.success)
print(sol.nfev)
print(sol.fun)
sol.params


# In[ ]:

plt.figure(figsize=(6,8),dpi=150)
sol = []
for T04 in np.linspace(273.15, 1033.15, 10):
   sol.append(e.general_explicit_map(params={'M_flight':0, 'T04': T04}))
   plt.annotate("{:.0f}K".format(T04), (sol[-1].params['MFP'], sol[-1].params['P0_ratio_c']))
    
MFP = [s.params['MFP'] for s in sol]
P0_ratio = [s.params['P0_ratio_c'] for s in sol]

# plt.plot(MFP_grid, P0_grid, 'k,')
CS2 = plt.contour(MFP_grid, P0_grid, Mb, levels=np.arange(0,2,0.1), 
                    colors='k', linewidths=0.7*1.625)
plt.clabel(CS2, CS2.levels, fmt='%.1f')
CS2 = plt.contour(MFP_grid, P0_grid, eff, colors='k', linewidths=0.7,
                  levels=[0.5,0.8,0.9,0.94,0.95])
plt.clabel(CS2, CS2.levels, fmt='%.2f')
plt.plot([MFP_choke, MFP_choke, sol_P1.params['MFP1']],
         [P0_max, sol_cr.params['P0_ratio'], P0_min], 'k--', label='choke limit')

plt.plot(MFP, P0_ratio, 'k.-', linewidth=0.7*1.62**2)
plt.xlim((0,0.6))
plt.ylim((1,2))
plt.xlabel("MFP")
plt.ylabel(r"$\frac{P_{03}}{P_{02}}$")
# plt.legend(("Mb", "eta_p"))
plt.legend()
plt.savefig('map.pdf', dpi=300)


# In[ ]:




# In[ ]:



