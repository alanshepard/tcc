from compressor import Compressor
from turbine import Turbine

# Outputs: N, F

class Engine:
    def __init__(self):
        self.compressor = Compressor()
        self.turbine = Turbine()

    def implicit_map(self,
                     M_flight                                 # flight mach number
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
        gam_c = self.compressor.gam
        gam_t = self.turbine.gam
        cpc = 1004e3 # J/(kg K)
        cpt = 1225e3 # J/(kg K)
        T01 = 273.15 # K
        P01 = 101e3  # Pa
        FHV = #Fuel lower calorific value
        eta_combustion = 
        R_c = # from gam and cp
        R_t = 

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
        omega_t = Mb_t * a04 / (self.turbine.geom['D1']/2)
        mdot_c = MFP * 
        mdot_t = MFP4 *

        ### CONSTRAINTS ###
        # * Energy addition in ther burner
        res_T04 = T03 + mdot_f*FHV*eta_combustion/(mdot_c*cpc) - T04
        # * Spool has constant speed
        res_omega = omega_c-omega_t
        
        # * Conservation of energy in the spool
        res_energy = mdot_c*cpc*(T03-T02) - mdot_t*cpt(T04-T05)

        # * Consevation of mass in the combustor
        res_mdot = mdot_3 + mdot_f - mdot_4

        # * Nozzle exit is either choked or at ambient pressure
        A_ratio_star = MFP5*((gam_t+1)/2)**((gam_t+1)/(2*(gam-1)))
        if A9/A5 > A_ratio_star:
            # Flow is subsonic
            P09 = P05 #isentropic nozzle
            P9 = P01/(1+(gam_c-1)/2*M_flight**2)**(gam_c/(gam_c-1))
            assert P9<P09
            M9 = (((P09/P9)**((gam_t-1)/gam_t) - 1)*2/(gam_t-1))**0.5
            assert M9<1
        else:
            # Flow is sonic
            M9 = 1
        
        MFP9 = mach2mfp(M9, gam_t)
        # isentropic nozzle
        res_MFP_nozzle = MFP9-MFP5*A5/A9

        ### remove dimensions of residuals ###
        # References for removing dimensions of residuals
        
        mdot_ref = (T01*R_c)**0.5/(P01*gam_c**0.5*self.compressor.geom['A1'])
        h_ref = cpc*T01
        omega_ref = (gam_c*R_c*T01)**0.5/self.compressor.geom['D1']
        
        #remove dimensions
        res_omega /= omega_ref
        res_energy /= h_ref*mdot_ref
        res_mdot /= mdot_ref
        res_T04 /= T04

        ### COMPRESSOR MAP ###
        res_MFP_c, res_T0_ratio_c, res_P0_ratio_c = self.compressor.implicit_map(MFP, MFP3, Mb_c, T0_ratio_c, P0_ratio_c, tol=1e-13)
        ### TURBINE MAP ###
        res_MFP_t, res_T0_ratio_t, res_P0_ratio_t = self.turbine.implicit_map(MFP4, MFP5, Mb_t, T0_ratio_t, P0_ratio_t, tol=1e-13)

        return (res_omega, res_energy, res_mdot, res_MFP_nozzle, res_T04
                res_MFP_c, res_T0_ratio_c, res_P0_ratio_c, 
                res_MFP_t, res_T0_ratio_t, res_P0_ratio_t)
