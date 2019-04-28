from math import tan, pi
import numpy as np

import compressor_losses as losses
from mfp2mach import mach2mfp, mfp2mach

def compressor(state, geom, ):

    MFP1 = state['MFP1']
    Mb = state['Mb']
    gam = state['gam']
    beta2 = geom['beta2']
    A_ratio = geom['A_ratio']

    MFP_choke = mach2mfp(1,gam)
    M1 = mfp2mach(MFP1, gam)
    
    # Initial_guess on MFP2. Assuming MFP2 as choked is safer because it allows us to approach the true value from above,
    #  without the risk of MFP2 appearing to be larger than the choke value.
    MFP2 = MFP_choke
    
    # Assuming compressor is isentropic and choked, we can get a good guess T0_ratio
    T0_ratio = (A_ratio*MFP2/MFP1)**(1-gam) 
    
    # Now we can begin the iterative part 
    success = np.full_like(MFP1, False, dtype=bool)
    choked = np.full_like(MFP1, False, dtype=bool)
    for iter in range(50):
        print(MFP2)
        M2 = mfp2mach(MFP2, gam)
        phi1 = MFP1/Mb*(1+(gam-1)/2*M1**2)**(1/(gam-1))
        phi2 = MFP2/Mb*(1+(gam-1)/2*M2**2)**(1/(gam-1))*T0_ratio**0.5
        psi_euler = 1-phi2*tan(beta2)
        loss_internal = losses.internal(slip_factor=0.9, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **geom)
        loss_parasitic = losses.parasitic(slip_factor=0.9, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **geom)
        psi_isen = psi_euler-sum(loss_internal.values())
        psi_actual = psi_euler+sum(loss_parasitic.values())
        T0_ratio = (gam-1)*psi_actual*Mb**2 + 1
        P0_ratio = ((gam-1)*psi_isen*Mb**2 + 1)**(gam/(gam-1))

        MFP2_new = MFP1/(A_ratio*P0_ratio)*T0_ratio

        err_MFP = MFP2_new - MFP2

        MFP2 = MFP2_new

        success[err_MFP<1e-8] = True
        choked[MFP2>MFP_choke] = True

        if np.all(np.abs(err_MFP[~choked]) < 1e-8):
            break

    if np.any(~choked & ~success):
        raise RuntimeError("Problem converging and not choked!!! Check if Mb>0 (Mb={})".format(Mb))


    return {
            'success': success,
            'choked': choked,
            'MFP2': MFP2,
            'phi1': phi1,
            'phi2': phi2,
            'psi_euler': psi_euler,
            'loss_internal': loss_internal,
            'loss_parasitic': loss_parasitic,
            'psi_isen': psi_isen,
            'psi_actual': psi_actual,
            'T0_ratio': T0_ratio,
            'P0_ratio': P0_ratio,
            }

def choke_line():
        """
        This function returns the part of the map where the outlet is choked (i.e. M2==1)
        """

        MFP_choke = mach2mfp(1, gam)

        phi1 = MFP1/Mb*(1+(gam-1)/2*M1**2)**(1/(gam-1))
        phi2 = MFP_choke/Mb*(1+(gam-1)/2*M2**2)**(1/(gam-1))*T0_ratio**0.5
        psi_euler = 1-phi2*tan(beta2)
        loss_internal = losses.internal(slip_factor=0.9, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **geom)
        loss_parasitic = losses.parasitic(slip_factor=0.9, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **geom)
        psi_isen = psi_euler-sum(loss_internal.values())
        psi_actual = psi_euler+sum(loss_parasitic.values())
        T0_ratio = (gam-1)*psi_actual*Mb**2 + 1
        P0_ratio = ((gam-1)*psi_isen*Mb**2 + 1)**(gam/(gam-1))

        MFP2_new = MFP1/(A_ratio*P0_ratio)*T0_ratio



state = {
        'gam':  1.4,
        'Mb':   1.3,
        'MFP1': 0.1 
        }

geom = {
        #Measured parameters
        'b': 5e-3,
        'D2': 64e-3,
        'D1t': 42.8e-3,
        'D1h': 17e-3,
        'beta2': 10*pi/180,
        'beta1': 40*pi/180,
        'Z': 10,
}

# Add derived geometry to geom
geom['A1'] = pi/4*(geom['D1t']**2-geom['D1h']**2)
geom['A2'] = pi*geom['D2']*geom['b']
geom['A_ratio'] = geom['A2']/geom['A1']
geom['R_ratio'] = geom['D2']/((geom['D1h']**2 + geom['D1t']**2)/2)**0.5
geom['D1t_D2'] = geom['D1t']/geom['D2']
geom['D1h_D2'] = geom['D1h']/geom['D2']

print(compressor(state, geom))

import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(2,1)

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
    P0_ratio_list = []
    MFP2_list = []
    MFP1_list = []
    for MFP1 in np.linspace(0,mach2mfp(1,1.4),100):
        state = {
                'gam':  1.4,
                'Mb':   Mb,
                'MFP1': MFP1, 
        }
        
        sol = compressor(state, geom)
        print(sol)
        if sol['success']:
            MFP1_list.append(MFP1)
            MFP2_list.append(sol['MFP2'])
            P0_ratio_list.append(sol['P0_ratio'])

    ax1.plot(MFP1_list, P0_ratio_list, '-+')
    ax2.plot(MFP1_list, MFP2_list, '-+')

plt.show()
