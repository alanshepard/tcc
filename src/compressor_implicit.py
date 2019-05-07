from math import tan, pi

import compressor_losses as losses
from mfp2mach import mach2mfp, mfp2mach

def compressor(x, state, geom):
    """
    x: vector of variables which depend on iteration.
    """

    # Unpack vector of variables which depend on iteration
    M2, T0_ratio = x

    ## Explicit equations
    # NOTE FOR MFP2: 
    # MFP2 could be calculated from MFP1, but this often leads to results larger than the choke limit,
    # and that is non-physical. Thus it was preferred to leave M2 as a guess and just report an error on it.

    MFP2 = mach2mfp(M2, state['gam'])
    phi1 = state['MFP1']/state['Mb']*(1+(state['gam']-1)/2*state['M1']**2)
    phi2 = MFP2/state['Mb']*(1+(state['gam']-1)/2*M2**2)*T0_ratio**0.5
    psi_euler = 1-phi2*tan(geom['beta2'])
    loss_internal = sum(losses.internal(
                                        MFP2=MFP2, 
                                        phi1=phi1, 
                                        phi2=phi2, 
                                        psi_euler=psi_euler, 
                                        slip_factor=0.9,
                                        **geom, 
                                        **state
                                        ).values()
                        )
    loss_parasitic = sum(losses.parasitic(
                                        MFP2=MFP2, 
                                        phi1=phi1, 
                                        phi2=phi2,
                                        psi_euler=psi_euler,
                                        slip_factor=0.9,
                                        **geom, 
                                        **state
                                        ).values()
                        )
    psi_isen = psi_euler - loss_internal
    psi_actual = psi_euler + loss_parasitic
    P0_ratio = ((state['gam']-1)*psi_isen*state['Mb']**2 + 1)**(state['gam']/(state['gam']-1))

    ## Implicit equations
    err_MFP = state['MFP1']/MFP2 - geom['A_ratio']*P0_ratio/T0_ratio
    err_T0_ratio = T0_ratio - ((state['gam']-1)*psi_actual*state['Mb']**2 + 1)

    return err_MFP, err_T0_ratio

state = {
        'gam':  1.4,
        'Mb':   0.3,
        'MFP1': 0.5 
        }

# Add M1 to state to avoid iterating within the non linear solver
state['M1'] = mfp2mach(state['MFP1'], state['gam'])

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

x0 = (0.5, 1)
print(compressor(x0, state, geom))

import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(2,1)
for Mb in np.linspace(0,1,10):
    T0_ratio_list = []
    M2_list = []
    M1_list = []
    for M1 in np.linspace(0,1,100):
        state = {
                'gam':  1.4,
                'Mb':   Mb,
                'MFP1': mach2mfp(M1, state['gam']), 
                'M1': M1
        }
        
        sol = root(compressor, x0, (state, geom))
        print(sol)
        if sol.success and sol.x[0]<1:
            M1_list.append(M1)
            M2_list.append(sol.x[0])
            T0_ratio_list.append(sol.x[1])

    ax1.plot(M1_list, T0_ratio_list, '-+')
    ax2.plot(M1_list, M2_list, '-+')

plt.show()

