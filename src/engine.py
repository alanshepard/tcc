from math import pi

from turbomachine import Turbomachine
from compressor import Compressor
from turbine import Turbine
from mfp2mach import mach2mfp

# Outputs: N, F

class Engine(Turbomachine):
    def __init__(self):
        self.compressor = Compressor()
        self.turbine = Turbine()


    DEFAULT_PARAMS = {'P0_ratio_c': 1.2347241010792356,
             'Mb_t': 0.27072466251896482,
              'MFP4': 0.15310678698124691,
               'MFP3': 0.14924987675288623,
                'M_flight': 0,
                 'mdot_f': 0.0019637316313999187,
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
                     mdot_f, T04,                             # burner
                     MFP4, MFP5, Mb_t, T0_ratio_t, P0_ratio_t # turbine 
                     ):
        """
        This function implements burner, nozzle and spool dynamics and integrates it with
        turbine and compressor maps

        Variable balance: 
            12 variables
          - 11 equations
        --------
             2 free choices (e.g. M_flight and mdot_f)
             
        """

        # Expose constants
        A5 = self.turbine.geom['A2']
        A9 = pi/4*30e-3**2 #TODO: check
        gam_c = self.compressor.gam
        gam_t = self.turbine.gam
        cpc = 1004 # J/(kg K)
        cpt = 1225 # J/(kg K)
        T01 = 273.15 # K
        P01 = 101e3  # Pa
        FHV = 43e6 #J/kg Fuel lower calorific value
        eta_combustion = 0.5
        R_c = (1-1/gam_c)*cpc # from gam and cp
        R_t = (1-1/gam_t)*cpt

        ### Dimensionalize everything ###
        T02 = T01
        T03 = T02*T0_ratio_c
        T05 = T04*T0_ratio_t
        P02 = P01
        P03 = P02*P0_ratio_c
        P04 = P03
        P05 = P03*P0_ratio_t
        a01 = (gam_c*R_c*T01)**0.5
        omega_c = Mb_c * a01 / (self.compressor.geom['D2']/2)
        a04 = (gam_t*R_t*T04)**0.5
        omega_t = Mb_t * a04 / (self.turbine.geom['D1t']/2)
        mdot_2 = MFP * self.compressor.geom['A1']*P01*gam_c**0.5/(T01*R_c)**0.5
        mdot_3 = MFP3 * self.compressor.geom['A2']*P03*gam_c**0.5/(T03*R_c)**0.5
        mdot_4 = MFP4 * self.turbine.geom['A1']*P04*gam_t**0.5/(T04*R_t)**0.5
        mdot_5 = MFP5 * self.turbine.geom['A2']*P05*gam_t**0.5/(T05*R_t)**0.5

        ### CONSTRAINTS ###
        # * Energy addition in ther burner
        res_T04 = T03 + mdot_f*FHV*eta_combustion/(mdot_3*cpc) - T04
        # * Spool has constant speed
        res_omega = omega_c-omega_t
        
        # * Conservation of energy in the spool
        res_energy = mdot_2*cpc*(T03-T02) - mdot_4*cpt*(T04-T05)

        # * Consevation of mass in the combustor
        res_mdot = mdot_3 + mdot_f - mdot_4

        # * Nozzle exit is either choked or at ambient pressure
        Pa = P01/(1+(gam_c-1)/2*M_flight**2)**(gam_c/(gam_c-1))
        MFP9 = MFP5*A5/A9
        res_MFP_nozzle = (P05/Pa-1)*MFP9 - (P05/Pa-1)*mach2mfp(min(1, (2/(gam_t-1)*((P05/Pa)**((gam_t-1)/gam_t)-1))**0.5) if P05/Pa>1 else 0 ,gam_t)

        ### remove dimensions of residuals ###
        # References for removing dimensions of residuals
        
        mdot_ref = (T01*R_c)**0.5/(P01*gam_c**0.5*self.compressor.geom['A1'])
        h_ref = cpc*T01
        omega_ref = (gam_c*R_c*T01)**0.5/self.compressor.geom['D2']
        
        #remove dimensions
        res_omega /= omega_ref
        res_energy /= h_ref*mdot_ref
        res_mdot /= mdot_ref
        res_T04 /= T04

        ### COMPRESSOR MAP ###
        res_MFP_c, res_T0_ratio_c, res_P0_ratio_c = self.compressor.implicit_map(MFP, MFP3, Mb_c, T0_ratio_c, P0_ratio_c, tol=1e-13)
        ### TURBINE MAP ###
        res_MFP_t, res_T0_ratio_t, res_P0_ratio_t = self.turbine.implicit_map(MFP4, MFP5, Mb_t, T0_ratio_t, P0_ratio_t, tol=1e-13)

        return (res_omega, res_energy, res_mdot, res_MFP_nozzle,
                res_T04,
                res_MFP_c, res_T0_ratio_c, res_P0_ratio_c, 
                res_MFP_t, res_T0_ratio_t, res_P0_ratio_t)


if __name__ == '__main__':
    e=Engine()
    e.implicit_map(**e.__class__.DEFAULT_PARAMS)
    sol = e.general_explicit_map(params={'M_flight': 0, 'T04': 1000})#, 'MFP': 0.4})
    print(sol)
