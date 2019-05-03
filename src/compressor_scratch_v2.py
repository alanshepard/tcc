from math import tan, pi, sin, exp
import numpy as np

import compressor_losses as losses
from mfp2mach import mach2mfp, mfp2mach

# TODO
# slip_factor
# other losses
# choke line

def compressor(state, geom, tol=1e-6):

    MFP1 = state['MFP1']
    Mb = state['Mb']
    gam = state['gam']
    beta2 = geom['beta2']
    A_ratio = geom['A_ratio']
    Z = geom['Z']
    R_ratio = geom['R_ratio']

    slip_factor = slip_wiesner(beta2, Z, R_ratio)

    MFP_choke = mach2mfp(1,gam)
    M1 = mfp2mach(MFP1, gam)

    # Initial_guess on MFP2. Assuming MFP2 as choked is safer because it allows us to approach the true value from above,
    #  without the risk of MFP2 appearing to be larger than the choke value.
    MFP2 = 0.1#MFP_choke

    # Assuming compressor is isentropic and choked, we can get a good guess T0_ratio
    T0_ratio = (A_ratio*MFP2/MFP1)**(-2*(gam-1)/(gam+1)) 

    # Now we can begin the iterative part 
    success = np.full_like(MFP1, False, dtype=bool)
    choked = np.full_like(MFP1, False, dtype=bool)
    for i in range(1000):
        M2 = mfp2mach(MFP2, gam)
        phi1 = MFP1/Mb*(1+(gam-1)/2*M1**2)**(1/(gam-1))
        phi2 = MFP2/Mb*(1+(gam-1)/2*M2**2)**(1/(gam-1))*T0_ratio**0.5
        psi_euler = slip_factor*(1-phi2*tan(beta2))
        loss_internal = losses.internal(slip_factor=slip_factor, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **geom)
        loss_parasitic = losses.parasitic(slip_factor=slip_factor, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **geom)
        psi_isen = psi_euler-sum(loss_internal.values())
        psi_actual = psi_euler+sum(loss_parasitic.values())
        T0_ratio = (gam-1)*psi_actual*Mb**2 + 1
        P0_ratio = ((gam-1)*psi_isen*Mb**2 + 1)**(gam/(gam-1))

        MFP2_new = MFP1/(A_ratio*P0_ratio)*T0_ratio**0.5

        err_MFP = np.abs(MFP2_new - MFP2)

        MFP2 = MFP2_new

        success[err_MFP<tol] = True
        choked[MFP2>=MFP_choke] = True

        if np.all(success[~choked]):
            break

    if np.any(~choked & ~success):
        raise RuntimeError("Problem converging and not choked!!! Check if Mb>0 (Mb={})".format(Mb))


    return {
            'success': success,
            'choked': choked,
            'slip_factor': slip_factor,
            'MFP2': MFP2,
            'phi1': phi1,
            'phi2': phi2,
            'psi_euler': psi_euler,
            'psi_isen': psi_isen,
            'psi_actual': psi_actual,
            'T0_ratio': T0_ratio,
            'P0_ratio': P0_ratio,
            'eff': np.log(P0_ratio**((gam-1)/gam))/np.log(T0_ratio),
            'iter': i,
            'internal_losses': loss_internal,
            'parasitic_losses': loss_parasitic,
            }

def slip_wiesner(beta2, Z, R_ratio):
    slip_factor = 1 - beta2**0.5/Z**0.7

    e = 1/R_ratio
    e_limit = 1/exp(8.16*sin(beta2/Z))

    if e>e_limit:
        slip_factor *= 1 - ((e-e_limit)/(1-e_limit))**3

    return slip_factor

def choke_line(Mb, gam, geom, tol=1e-6):
    """
    This function returns the part of the map where the outlet is choked (i.e. M2==1)
    """

    beta2 = geom['beta2']
    A_ratio = geom['A_ratio']
    Z = geom['Z']
    R_ratio = geom['R_ratio']

    slip_factor = slip_wiesner(beta2, Z, R_ratio)

    MFP_choke = mach2mfp(1,gam)

    # Initial_guess on MFP1.
    MFP1 = MFP_choke

    # Assuming compressor is isentropic and choked from inlet to outlet, we can get a good guess T0_ratio
    T0_ratio = (A_ratio)**(-2*(gam-1)/(gam+1)) 

    # Now we can begin the iterative part 
    success = np.full_like(Mb, False, dtype=bool)
    inlet_choked = np.full_like(Mb, False, dtype=bool)
    for i in range(10000):
        M1 = mfp2mach(MFP1, gam, tol/100)
        phi1 = MFP1/Mb*(1+(gam-1)/2*M1**2)**(1/(gam-1))
        phi2 = MFP_choke/Mb*((gam+1)/2)**(1/(gam-1))*T0_ratio**0.5
        psi_euler = slip_factor*(1-phi2*tan(beta2))
        loss_internal = losses.internal(slip_factor=slip_factor, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **geom)
        loss_parasitic = losses.parasitic(slip_factor=slip_factor, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **geom)
        psi_isen = psi_euler-sum(loss_internal.values())
        psi_actual = psi_euler+sum(loss_parasitic.values())
        T0_ratio = (gam-1)*psi_actual*Mb**2 + 1
        P0_ratio = ((gam-1)*psi_isen*Mb**2 + 1)**(gam/(gam-1))

        MFP1_new = MFP_choke*(A_ratio*P0_ratio)/T0_ratio**0.5

        err_MFP = np.abs(MFP1_new - MFP1)

        MFP1 = MFP1_new

        success[err_MFP<tol] = True
        inlet_choked[MFP1>MFP_choke] = True

        if np.all(success[~inlet_choked]):
            break

    if np.any(~inlet_choked & ~success):
        pass
        raise RuntimeError("Problem converging and not choked!!! Check if Mb>0 (Mb={})".format(Mb))


    print(i)
    return {
            'success': success,
            'inlet_choked': inlet_choked,
            'slip_factor': slip_factor,
            'MFP1': MFP1,
            'phi1': phi1,
            'phi2': phi2,
            'psi_euler': psi_euler,
            'psi_isen': psi_isen,
            'psi_actual': psi_actual,
            'T0_ratio': T0_ratio,
            'P0_ratio': P0_ratio,
            'eff': np.log(P0_ratio**((gam-1)/gam))/np.log(T0_ratio),
            'iter': i,
            'internal_losses': loss_internal,
            'parasitic_losses': loss_parasitic,
            }

def critical_choke_point(gam, geom, tol=1e-10):
    """Find Mb for which MFP1=MFP2=MFP_choke"""
    beta2 = geom['beta2']
    A_ratio = geom['A_ratio']
    Z = geom['Z']
    R_ratio = geom['R_ratio']

    slip_factor = slip_wiesner(beta2, Z, R_ratio)

    MFP_choke = mach2mfp(1,gam)

    #Guess Mb
    # Assuming compressor is isentropic and choked from inlet to outlet, we can take a good guess
    # at T0_ratio and Mb
    eff=1
    T0_ratio = 1.1# A_ratio**(-(2*gam-1)/((2*eff-1)*gam-1))
    Mb = 1.1 #((T0_ratio-1)/(gam-1))**0.5

    # Now we can begin the iterative part 
    success = False
    for i in range(1000):
        phi1 = MFP_choke/Mb*((gam+1)/2)**(1/(gam-1))
        phi2 = MFP_choke/Mb*((gam+1)/2)**(1/(gam-1))*T0_ratio**0.5
        psi_euler = slip_factor*(1-phi2*tan(beta2))
        loss_internal = losses.internal(slip_factor=slip_factor, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **geom)
        loss_parasitic = losses.parasitic(slip_factor=slip_factor, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **geom)
        psi_isen = psi_euler-sum(loss_internal.values())
        psi_actual = psi_euler+sum(loss_parasitic.values())
        T0_ratio = (gam-1)*psi_actual*Mb**2 + 1
        P0_ratio = ((gam-1)*psi_isen*Mb**2 + 1)**(gam/(gam-1))
        eff = np.log(P0_ratio**((gam-1)/gam))/np.log(T0_ratio)
        if eff<0.6:
            eff=0.6

        Mb_new = ((A_ratio**(-(2*(gam-1)/((2*eff-1)*gam-1))) -1)/((gam-1)*psi_actual))**0.5

        err = np.abs(Mb_new - Mb)
  
        Mb = Mb_new

        if err<tol:
            success = True
        #     break

    return {
            'success': success,
            'slip_factor': slip_factor,
            'phi1': phi1,
            'phi2': phi2,
            'psi_euler': psi_euler,
            'psi_isen': psi_isen,
            'psi_actual': psi_actual,
            'T0_ratio': T0_ratio,
            'P0_ratio': P0_ratio,
            'eff': eff,
            'iter': i,
            'internal_losses': loss_internal,
            'parasitic_losses': loss_parasitic,
            }



state = {
        'gam':  1.4,
        'Mb':   0.8,
        'MFP1': np.linspace(0, mach2mfp(1,1.4)) 
        }

geom = {
        #Measured parameters
        'b': 5e-3,
        'D2': 64e-3,
        'D1t': 42.8e-3,
        'D1h': 17e-3,
        'beta2': 10*pi/180,
        'beta1': 40*pi/180,
        'alfa3': 69*pi/180,
        'Z': 10,
        }

# Add derived geometry to geom
geom['A1'] = pi/4*(geom['D1t']**2-geom['D1h']**2)
geom['A2'] = pi*geom['D2']*geom['b']
geom['A_ratio'] = geom['A2']/geom['A1']
geom['R_ratio'] = geom['D2']/((geom['D1h']**2 + geom['D1t']**2)/2)**0.5
geom['D1t_D2'] = geom['D1t']/geom['D2']
geom['D1h_D2'] = geom['D1h']/geom['D2']

import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt

sol = compressor(state, geom)
fig, ax = plt.subplots(1,1)
for k, v in sol['internal_losses'].items():
    ax.plot(state['MFP1'],v, label=k)
ax.legend()


fig, (ax1, ax2) = plt.subplots(2,1)
ax1.set_ylim((0,1))
fig, ax3 = plt.subplots(1,1)
sol = critical_choke_point(1.4, geom)
ax2.plot(mach2mfp(1,1.4), sol['P0_ratio'], '*')
ax1.plot(mach2mfp(1,1.4), sol['eff'], '*')
ax2.plot(2*[mach2mfp(1, 1.4)], [1,3], 'k--')

choked = choke_line(np.linspace(0.1,1.5,1000), 1.4, geom)
MFP1 = choked['MFP1']
success = choked['success']
P0_ratio = choked['P0_ratio']
eff = choked['eff']
ax2.plot(MFP1[success], P0_ratio[success], 'k--')
ax1.plot(MFP1[success&(0<=eff)&(eff<=1)], eff[success&(0<=eff)&(eff<=1)],'k--')
print(mach2mfp(1,1.4))
print(max(MFP1[success]))


Mb = np.linspace(0.1, 1, 10)
MFP1 = np.linspace(0, mach2mfp(1,1.4))
Mb_grid, MFP1_grid = np.meshgrid(Mb, MFP1)
state = {
        'gam': 1.4,
        'Mb': Mb_grid,
        'MFP1': MFP1_grid,
}

# sol = compressor(state, geom)
# ax1.contour(MFP1_grid, sol['P0_ratio'], Mb_grid)
# ax1.contour(MFP1_grid, sol['MFP2'], Mb_grid)

# plt.show()

for Mb in np.linspace(0.1,1,10):
    MFP1 = np.linspace(0, mach2mfp(1,1.4),100)
    state = {
            'gam':  1.4,
            'Mb':   Mb,
            'MFP1': MFP1, 
    }
        
    sol = compressor(state, geom)
    print(sol['iter'])
    success = sol['success']
    MFP2 = sol['MFP2']
    P0_ratio = sol['P0_ratio']
    phi2=sol['phi2']
    eff = sol['eff']
    ax2.plot(MFP1[success], P0_ratio[success], '-+')
    ax1.plot(MFP1[success & (P0_ratio>1)], eff[success & (P0_ratio>1)])
    ax3.plot(MFP1[success], MFP2[success])

np.savetxt('comp.txt', np.c_[MFP1[success], P0_ratio[success]])

plt.show()

