from math import tan, pi
import numpy as np

from mfp2mach import mach2mfp, mfp2mach
from turbomachine import Turbomachine, gridmap


class Turbine(Turbomachine):
    def __init__(self):
        super().__init__()
        self.gam = 1.3
        self.geom = {'beta5': -56*pi/180,
                     'alfa4': 60*pi/180,
                     'D1t': 65.3e-3,
                     'D1h': 45e-3,
                     'D2t': 67e-3,
                     'D2h': 46e-3,
                    }
        self.geom['A1'] = pi/4*(self.geom['D1t']**2 - self.geom['D1h']**2)
        self.geom['A2'] = pi/4*(self.geom['D2t']**2 - self.geom['D2h']**2)
        self.geom['A_ratio'] = self.geom['A2']/self.geom['A1']
        self.geom['Dt_ratio'] = self.geom['D2t']/self.geom['D1t']

    DEFAULT_PARAMS = {'MFP1':0.3, 'MFP2':0.3, 'Mb':0.5, 'T0_ratio':1, 'P0_ratio':1}
    N_FREE_PARAMS = 2
    
    def implicit_map(self, MFP1, MFP2, Mb, T0_ratio, P0_ratio, tol=1e-12):

        # Expose instance variables
        gam     = self.gam
        beta5   = self.geom['beta5']
        alfa4   = self.geom['alfa4']
        A_ratio = self.geom['A_ratio']
        Dt_ratio = self.geom['Dt_ratio']

        # Calculate auxiliary variables
        M4 = mfp2mach(MFP1, gam, tol)
        M5 = mfp2mach(MFP2, gam, tol)

        phi4 = MFP1/Mb*(1+(gam-1)/2*M4**2)**(1/(gam-1))
        phi5 = MFP2/Mb*(1+(gam-1)/2*M5**2)**(1/(gam-1))*T0_ratio**0.5

        psi = 1 + phi5*tan(beta5) - 1/(Dt_ratio)**2*phi4*tan(alfa4)

        # Non-linear system
        res_MFP = MFP2 - MFP1*T0_ratio**0.5/(A_ratio*P0_ratio)
        res_P0_ratio = P0_ratio - (1+(gam-1)*psi*Mb**2)**(gam/(gam-1))
        res_T0_ratio = T0_ratio - (1+(gam-1)*psi*Mb**2)

        return res_MFP, res_T0_ratio, res_P0_ratio

class TurbineExtendedMap(Turbine):
    """
    This class adds the parameter MbMFP=Mb*MFP1 to the turbine map
    """
    DEFAULT_PARAMS = Turbine.DEFAULT_PARAMS.copy()
    DEFAULT_PARAMS['MbMFP'] = Turbine.DEFAULT_PARAMS['Mb']*Turbine.DEFAULT_PARAMS['MFP1']

    def implicit_map(self, MFP1, MFP2, Mb, T0_ratio, P0_ratio, MbMFP, tol=1e-12):
        res_turbine = super().implicit_map(MFP1, MFP2, Mb, T0_ratio, P0_ratio, tol)

        return (*res_turbine, MbMFP-Mb*MFP1)

    def plot_map(self, ax, extent=None, samples=25, grid=False):

        if extent is None:
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
        else:
            xmin, xmax, ymin, ymax = extent

        #Choke line
        MFP_choke = mach2mfp(1, self.gam)
        sol = self.general_explicit_map({'MFP1': 0.9*MFP_choke, 'MFP2': 0.9*MFP_choke})
        MbMFP = np.linspace(sol.params['MbMFP'], xmax, samples*5)
        pr_choke = np.empty_like(MbMFP)
        for i, MbMFP_ in enumerate(MbMFP):
            if sol.success:
                old_sol = sol
            sol = self.general_explicit_map({'MbMFP': MbMFP_, 'MFP2': MFP_choke}, 
                                            initial_guesses=old_sol.params)
            pr_choke[i] = sol.params['P0_ratio'] if sol.success else np.nan

        ax.plot(MbMFP, 1/pr_choke, 'k--', linewidth=0.8)


        #Map

        MbMFP, P0_ratio = np.meshgrid(np.linspace(xmin, xmax, samples), np.linspace(1/ymax, 1/ymin, samples))
        
        params = gridmap(self, MbMFP, P0_ratio, 'MbMFP')
        
        params['eff'] = (self.gam-1)/self.gam*np.log(P0_ratio)/np.log(params['T0_ratio'])
    
        if grid:
            ax.plot(MbMFP, 1/P0_ratio, 'k,')

        CS = ax.contour(MbMFP, 1/P0_ratio, params['Mb'], levels=np.arange(0,2,0.1), 
                            colors='k', linewidths=0.8)
        ax.clabel(CS, CS.levels, fmt='%.1f', rightside_up=False)
        CS.collections[0].set_label('$M_{bt}$')

        #CS2 = ax.contour(MbMFP, 1/P0_ratio, params['eff'], colors='k', linewidths=0.4,
        #                  levels=[0.5,0.8,0.9,0.95])
        #CS2.collections[0].set_label(r'$\eta_p$')
        #ax.clabel(CS2, CS2.levels, fmt='%.2f')

        ax.set_xlabel(r"$M_{bt}\text{MFP}_4$")
        ax.set_ylabel(r"$\frac{P_{04}}{P_{05}}$", rotation=0)
