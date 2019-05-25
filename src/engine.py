from math import pi
from collections import namedtuple

from turbomachine import Turbomachine
from compressor import Compressor
from turbine import Turbine
from mfp2mach import mach2mfp

DimensionalParameters = namedtuple('DimensionalEngineParameters',
             'T02 T03 T04 T05 P02 P03 P04 P05 omega_c omega_t mdot2 mdot3 mdot4 mdot5')

class Engine(Turbomachine):
    def __init__(self):
        self.compressor = Compressor()
        self.turbine = Turbine()
        self.initial_guess = self.__class__.DEFAULT_PARAMS.copy()
        self.cpc = 1004 # J/(kg K)
        self.cpt = 1225 # J/(kg K)
        self.T01 = 273.15 # K
        self.P01 = 101e3  # Pa
        self.R_c = (1-1/self.compressor.gam)*self.cpc # from gam and cp
        self.R_t = (1-1/self.turbine.gam)*self.cpt
        self.A8 = pi/4*35e-3**2 # Measured: D8 = 45mm
        self.FHV = 43e6 #J/kg Fuel lower calorific value
        self.eta_combustion = 0.5


    DEFAULT_PARAMS = {'P0_ratio_c': 1.2347241010782356,
             'Mb_t': 0.27072466251896482,
              'MFP4': 0.15310678698124691,
               'MFP3': 0.14924987675288623,
                'M_flight': 0,
                 'mdotf': 0.0019637316313999187,
                  'P0_ratio_t': 0.92680225375718472,
                   'MFP5': 0.15451547834023707,
                    'Mb_c': 0.44527731422329964,
                     'MFP': 0.14801196525452109,
                      'T0_ratio_t': 0.98261102832775471,
                       'T0_ratio_c': 1.0669702796449636,
                        'T04': 840.63888888888891}

    N_FREE_PARAMS = 2

    def implicit_map(self,
                     M_flight,                                # flight mach number
                     MFP, MFP3, Mb_c, T0_ratio_c, P0_ratio_c, # compressor
                     mdotf, T04,                             # burner
                     MFP4, MFP5, Mb_t, T0_ratio_t, P0_ratio_t # turbine 
                     ):
        """
        This function implements burner, nozzle and spool dynamics and integrates it with
        turbine and compressor maps

        Variable balance: 
            12 variables
          - 11 equations
        --------
             2 free choices (e.g. M_flight and mdotf)
             
        """

        # Expose constants
        A8 = self.A8
        FHV = self.FHV
        eta_combustion = self.eta_combustion
        A5 = self.turbine.geom['A2']
        gam_c = self.compressor.gam
        gam_t = self.turbine.gam
        cpc = self.cpc
        cpt = self.cpt
        T01 = self.T01 
        P01 = self.P01  
        R_c = self.R_c
        R_t = self.R_t

        ### Dimensionalize everything ###
        dim_params = self.dimensionalize(
                     M_flight,                                # flight mach number
                     MFP, MFP3, Mb_c, T0_ratio_c, P0_ratio_c, # compressor
                     mdotf, T04,                             # burner
                     MFP4, MFP5, Mb_t, T0_ratio_t, P0_ratio_t # turbine 
                     )
        T02 = dim_params.T02
        T03 = dim_params.T03
        T05 = dim_params.T05 
        P02 = dim_params.P02 
        P03 = dim_params.P03 
        P04 = dim_params.P04 
        P05 = dim_params.P05 
        omega_c = dim_params.omega_c
        omega_t = dim_params.omega_t
        mdot2 = dim_params.mdot2
        mdot3 = dim_params.mdot3
        mdot4 = dim_params.mdot4
        mdot5 = dim_params.mdot5

        ### CONSTRAINTS ###
        # * Energy addition in ther burner
        res_T04 = T03 + mdotf*FHV*eta_combustion/(mdot3*cpc) - T04
        # * Spool has constant speed
        res_omega = omega_c-omega_t
        
        # * Conservation of energy in the spool
        res_energy = mdot2*cpc*(T03-T02) - mdot4*cpt*(T04-T05)

        # * Consevation of mass in the combustor
        res_mdot = mdot3 + mdotf - mdot4

        # * Nozzle exit is either choked or at ambient pressure
        Pa = P01/(1+(gam_c-1)/2*M_flight**2)**(gam_c/(gam_c-1))
        MFP8 = MFP5*A5/A8
        res_MFP_nozzle = (P05/Pa-1)*MFP8 - (P05/Pa-1)*mach2mfp(min(1, (2/(gam_t-1)*((P05/Pa)**((gam_t-1)/gam_t)-1))**0.5) if P05/Pa>1 else 0 ,gam_t)

        ### remove dimensions of residuals ###
        # References for removing dimensions of residuals
        
        mdotref = (T01*R_c)**0.5/(P01*gam_c**0.5*self.compressor.geom['A1'])
        h_ref = cpc*T01
        omega_ref = (gam_c*R_c*T01)**0.5/self.compressor.geom['D2']
        
        #remove dimensions
        res_omega /= omega_ref
        res_energy /= h_ref*mdotref
        res_mdot /= mdotref
        res_T04 /= T04

        ### COMPRESSOR MAP ###
        res_MFP_c, res_T0_ratio_c, res_P0_ratio_c = self.compressor.implicit_map(MFP, MFP3, Mb_c, T0_ratio_c, P0_ratio_c, tol=1e-13)
        ### TURBINE MAP ###
        res_MFP_t, res_T0_ratio_t, res_P0_ratio_t = self.turbine.implicit_map(MFP4, MFP5, Mb_t, T0_ratio_t, P0_ratio_t, tol=1e-13)

        return (res_omega, res_energy, res_mdot, res_MFP_nozzle,
                res_T04,
                res_MFP_c, res_T0_ratio_c, res_P0_ratio_c, 
                res_MFP_t, res_T0_ratio_t, res_P0_ratio_t)

    def dimensionalize(self,
                     M_flight,                                # flight mach number
                     MFP, MFP3, Mb_c, T0_ratio_c, P0_ratio_c, # compressor
                     mdotf, T04,                             # burner
                     MFP4, MFP5, Mb_t, T0_ratio_t, P0_ratio_t # turbine 
                     ):
        gam_c = self.compressor.gam
        gam_t = self.turbine.gam
        T01 = self.T01
        P01 = self.P01
        R_c = self.R_c
        R_t = self.R_t

        ### Dimensionalize everything ###
        a01 = (gam_c*R_c*T01)**0.5
        a04 = (gam_t*R_t*T04)**0.5
        T02 = T01
        T03 = T02*T0_ratio_c
        T04 = T04
        T05 = T04*T0_ratio_t
        P02 = P01
        P03 = P02*P0_ratio_c
        P04 = P03
        P05 = P03*P0_ratio_t
        omega_c = Mb_c * a01 / (self.compressor.geom['D2']/2)
        omega_t = Mb_t * a04 / (self.turbine.geom['D1t']/2)
        mdot2 = MFP * self.compressor.geom['A1']*P01*gam_c**0.5/(T01*R_c)**0.5
        mdot3 = MFP3 * self.compressor.geom['A2']*P03*gam_c**0.5/(T03*R_c)**0.5
        mdot4 = MFP4 * self.turbine.geom['A1']*P04*gam_t**0.5/(T04*R_t)**0.5
        mdot5 = MFP5 * self.turbine.geom['A2']*P05*gam_t**0.5/(T05*R_t)**0.5

        dim_params = DimensionalParameters(
                                           T02     = T02    ,
                                           T03     = T03    ,
                                           T04     = T04    ,
                                           T05     = T05    ,
                                           P02     = P02    ,
                                           P03     = P03    ,
                                           P04     = P04    ,
                                           P05     = P05    ,
                                           omega_c = omega_c,
                                           omega_t = omega_t,
                                           mdot2   = mdot2  ,
                                           mdot3   = mdot3  ,
                                           mdot4   = mdot4  ,
                                           mdot5   = mdot5  
                                           )

        return dim_params

    def working_line(self, T04_grid):
        
        wline=[]
        params = self.initial_guess
        for T04 in T04_grid:
            sol = self.general_explicit_map({'M_flight': 0, 'T04': T04}, params)
            params = sol.params
            wline.append(params)

        return wline
            


    def ode_fun(self,t,y):
        P4, mdot, omega = y

        self.compressor.explicit_map(mdot, omega)
        T04 = self.throtle(t, mdot, T03)
        self.turbine.explicit_map(P4, omega)

        P4_dot = a01**2*(mdot-mdott)
        mdotdot = A1/Lc*(P3-P4)
        omega_dot = 1/J*torque_t-torque_c

        return P4_dot, mdotdot, omega_dot


if __name__ == '__main__':
    e=Engine()
    e.implicit_map(**e.__class__.DEFAULT_PARAMS)
    sol = e.general_explicit_map(params={'M_flight': 0, 'T04': 1000})#, 'MFP': 0.4})
    print(sol)
