# coding: utf-8
from turbine import *
t=TurbineExtendedMap()
from mfp2mach import mach2mfp
import matplotlib.pyplot as plt
MFP_choke = mach2mfp(1, t.gam)
sol_cr = t.general_explicit_map({'MFP1':MFP_choke, 'MFP2':MFP_choke})
sol_almost_cr = t.general_explicit_map({'MFP1':0.4, 'MFP2':0.4})
print(sol_almost_cr)
t.initial_guess = sol_almost_cr.params

fig, (ax1, ax2) = plt.subplots(1,2)

ax1.plot(sol_cr.params['MFP1'], sol_cr.params['P0_ratio'], 's')

ax1.plot(sol_almost_cr.params['MFP1'], sol_almost_cr.params['P0_ratio'], 's')
assert sol_almost_cr.success
print(sol_cr)


sol=sol_almost_cr
for MbMFP in reversed(np.linspace(0, sol_cr.params['MbMFP'])):
    sol = t.general_explicit_map({'MbMFP': MbMFP, 'MFP1': MFP_choke}, initial_guesses=sol.params)
    ax1.plot(sol.params['MbMFP'], 1/sol.params['P0_ratio'], '+')
    ax2.plot(sol.params['MFP1'], sol.params['MFP2'], '+')
    if not sol.success:
        break

sol=sol_cr
for MbMFP in np.linspace(sol_cr.params['MFP1'], 0, 100):
    last_sol=sol
    sol = t.general_explicit_map({'MFP1': MbMFP, 'MFP2': MFP_choke}, initial_guesses=sol.params)
    ax1.plot(sol.params['MbMFP'], 1/sol.params['P0_ratio'], 'v')
    ax2.plot(sol.params['MFP1'], sol.params['MFP2'], 's')
    if not sol.success:        
        sol=last_sol
        break


for MbMFP in np.arange(sol.params['MbMFP'], 5, 0.1):
    last_sol=sol
    sol = t.general_explicit_map({'MbMFP': MbMFP, 'MFP2': MFP_choke}, initial_guesses=sol.params)
    ax1.plot(sol.params['MbMFP'], 1/sol.params['P0_ratio'], '^')
    ax2.plot(sol.params['MFP1'], sol.params['MFP2'], '^')
    if not sol.success:
        sol=last_sol
        break

sol = t.general_explicit_map({'MbMFP':sol.params['MbMFP'], 'MFP1':sol.params['MFP1']*0.7})
ax1.plot(sol.params['MbMFP'], sol.params['P0_ratio'], 'g*')

for MbMFP in np.linspace(sol_cr.params['MFP2'], 0, 100):
    sol = t.general_explicit_map({'MFP2': MbMFP, 'MFP1': MFP_choke}, initial_guesses=sol.params)
    ax1.plot(sol.params['MbMFP'], 1/sol.params['P0_ratio'], 'o')
    ax2.plot(sol.params['MFP1'], sol.params['MFP2'], 'o')
    if not sol.success:        break

plt.show()
