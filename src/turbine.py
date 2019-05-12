from math import tan, pi

from mfp2mach import mach2mfp, mfp2mach
from turbomachine import Turbomachine


class Turbine(Turbomachine):
    def __init__(self):
        self.gam = 1.4
        self.geom = {'beta5': -56*pi/180,
                     'alfa4':  60*pi/180,
                     'D1t': 65.3e-3,
                     'D1h': 45e-3,
                     'D2t': 67e-3,
                     'D2h': 46e-3,
                    }
        self.geom['A1'] = pi/4*(self.geom['D1t']**2 - self.geom['D1h']**2)
        self.geom['A2'] = pi/4*(self.geom['D2t']**2 - self.geom['D2h']**2)
        self.geom['A_ratio'] = self.geom['A2']/self.geom['A1']

    DEFAULT_PARAMS = {'MFP1':0.3, 'MFP2':0.3, 'Mb':0.5, 'T0_ratio':1, 'P0_ratio':1}
    N_FREE_PARAMS = 2
    
    def implicit_map(self, MFP1, MFP2, Mb, T0_ratio, P0_ratio, tol=1e-12):

        # Expose instance variables
        gam     = self.gam
        beta5   = self.geom['beta5']
        alfa4   = self.geom['alfa4']
        A_ratio = self.geom['A_ratio']

        # Calculate auxiliary variables
        M4 = mfp2mach(MFP1, gam, tol)
        M5 = mfp2mach(MFP2, gam, tol)

        phi4 = MFP1/Mb*(1+(gam-1)/2*M4**2)**(1/(gam-1))
        phi5 = MFP2/Mb*(1+(gam-1)/2*M5**2)**(1/(gam-1))*T0_ratio**0.5

        psi = 1 + phi5*tan(beta5) - phi4*tan(alfa4)

        # Non-linear system
        res_MFP = MFP2 - MFP1*T0_ratio**0.5/(A_ratio*P0_ratio)
        res_P0_ratio = P0_ratio - (1+(gam-1)*psi*Mb**2)**(gam/(gam-1))
        res_T0_ratio = T0_ratio - (1+(gam-1)*psi*Mb**2)

        return res_MFP, res_T0_ratio, res_P0_ratio

