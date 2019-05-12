from math import tan, pi
import numpy as np

from mfp2mach import mach2mfp, mfp2mach
import compressor_losses as losses
from turbomachine import Turbomachine


class Compressor(Turbomachine):
    def __init__(self):
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
        
        self.geom=geom
        
        self.gam=1.4
    
    def implicit_map(self, MFP1, MFP2, Mb, T0_ratio, P0_ratio, tol=1e-13):
        #Unpack constants
        gam = self.gam
        slip_factor = 0.9
        beta2 = self.geom['beta2']
        A_ratio = self.geom['A_ratio']        
    
        #Non linear model
        M1 = mfp2mach(MFP1, gam, tol)
        M2 = mfp2mach(MFP2, gam, tol)
        
        phi1 = MFP1/Mb*(1+(gam-1)/2*M1**2)**(1/(gam-1))
        phi2 = MFP2/Mb*(1+(gam-1)/2*M2**2)**(1/(gam-1))*T0_ratio**0.5
        
        psi_euler = slip_factor*(1-phi2*tan(beta2))
        loss_internal = losses.internal(slip_factor=slip_factor, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **self.geom)
        loss_parasitic = losses.parasitic(slip_factor=slip_factor, phi1=phi1, phi2=phi2, psi_euler=psi_euler, **self.geom)
        psi_isen = psi_euler-sum(loss_internal.values())
        psi_actual = psi_euler+sum(loss_parasitic.values())
        
        res_T0_ratio = T0_ratio - ((gam-1)*psi_actual*Mb**2 + 1)
        res_P0_ratio = P0_ratio - ((gam-1)*psi_isen*Mb**2 + 1)**(gam/(gam-1))
        res_MFP = MFP2 - MFP1*T0_ratio**0.5/(A_ratio*P0_ratio)
        
        return res_MFP, res_T0_ratio, res_P0_ratio
    