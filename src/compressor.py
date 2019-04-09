from math import pi, tan, cos
from numpy import sqrt

def norm(x, y, z)
    return sqrt(x**2+y**2+z**2)

class Compressor:
    def __init__():
        # Geometry
        beta1 # blade leading edge angle with respect to axial direction
        beta2 # blade trailing edge angle with respect to radial direction. 
              # Positive in direction oposite to rotation
        D1t   # impeler blade leading edge diameter at hub
        D1h   # impeler blade leading edge diameter at hub
        D2    # impeler blade trailing edge diameter
        A1    # circular entrance area
        A2    # ring shaped exit area
        z_fb  # full length blades
        z_sb  # splitter blades
        L     # blade mean streamline length
        L_sb  # splitter blade mean streamline length

    @property
    def z():
        """ 
        Effective number of blades (AUNGIER, 1995)
        """

        return self.z_fb + self.z_sb*self.L_sb/self.L

    @property
    def slip_factor():
        """
        AUNGIER, 1995
        """
        s = 1 - sqrt(beta2)/z**0.7

        #TODO: adicionar limite

        return s

    def impeller_distortion_factor(phi2):
        """
        AUNGIER, 1995
        """

        def loss_hub_to_shroud():
            """
            Assumptions:
            * Streamlines aligned with axis on entrance and 90 deg on exit
            """
            pass

        #To avoid iterating, calculate W2/U2 withot correction factors
        W2_U2 = abs(phi2)/cos(self.beta2)

        B2 = 
        l = 1/(1-B2)
        return l

    def load_param_euler(flow_param_exit):
        pass

    def flow_param_exit(mdot, omega):
        pass

    def eff_isen(flow_param_exit):
        losses = [self.incidence_loss()]
        eff = 1.-sum(losses)/load_param_euler(flow_param_exit)
        return eff


    ### Decomposition of W2 ###
    def  W2x_U2(phi2):
        return 0.

    def W2r_U2(phi2):
        """
        W2r = mdot/(rho*A2) = phi2*U2
        """
        return phi2 # BY DEFINITION!

    def W2u_U2(phi2):
        return self.V2u_U2(phi2)-1

    def W2_U2(phi2):
        return norm(self.W2r_U2(phi2),
                    self.W2x_U2(phi2),
                    self.W2u_U2(phi2)
                    )

    ### Decomposition of V2 ###
    
    def V2x_U2():
        return 0.

    def V2u_U2(phi2):    
        """
        V2u = sigma*(1-lambda*phi2*tan(beta2))*U2
        where sigma is the slip factor 
        and lambda is the impeller distortion factor
        """

        return self.slip_factor*(1-self.impeller_distortion_factor*phi2*tan(beta2))
    
    def V2r_U2(phi2):
        return self.W2r_U2(phi2)

    def V2_U2(phi2):
        return norm(self.V2x_U2(),
                    self.V2u_U2(phi2),
                    self.V2r_U2(phi2)
                    )

    ##

    ### Decomposition of W1t ###
    def W1tr_U2(): return 0.

    def W1tu_U2():
        """
        W1tu = U1t = D1t/D2*U2
        """
        return self.D1t/self.D2

    def W1tx_U2(phi2):
        """
        W1tx = V1x
        """
        return self.V1x_U2(phi2)

    def W1t_U2(phi2):
        return norm(self.W1tr_U2(),
                    self.W1tu_U2(),
                    self.W1tx_U2(phi2)
                    )

    ##

    ### Decomposition of W1h
    def W1hr_U2():
        return 0.

    def W1hu_U2():
        """
        W1hu = omega*D1h = D1h/D2*U2
        """
        return self.D1h/self.D2

    def W1hx_U2(phi2):
        """
        W1hx = V1x
        """
        return self.V1x_U2(phi2)

    def W1rht_U2(phi2):
        return norm(self.W1hr_U2(),
                    self.W1hx_U2(phi2),
                    self.W1hu_U2()
                    )

    ##
    
    ### Decomposition of V1 ###
    
    def V1r_U2(): 
        return 0.

    def V1u_U2(): 
        return 0.

    def V1x_U2(phi2):
        """
        V1x=mdot/(rho*A1)=A1/A2*phi2*U2
        """
        return self.A1/self.A2*phi2
    
    def V1_U2(phi2)
        return norm(self.V1x_U2(phi2),
                    self.V1u_U2(phi2),
                    self.V1r_U2(phi2)
                    )
        
    ##

    ### Losses ###
    #All losses written as a decrement on the load parameter
    def incidence_loss():
        f_inc = 0.7 # 0.5 to 0.7
        
        Wui_U2 = self.W1u_U2() self.D1t/self.D2

        delta_load_factor = f_inc *(W1u_U2)**2/2
        
        return delta_load_factor

    def blade_loading_loss(load_param_euler):

        D1t_D2 = self.D1t/self.D2

        W2_W1t = (phi2/cos(beta2))/(sqrt((self.D1t/self.D2)**2 + (phi2*(self.A1/self.A2))**2)
        
        Df = 1 - W2/W1t*(1 - (0.75*load_param_euler)
                             /(self.n_blades/pi*(1-D1t_D2)+2*(D1t_D2)))

        delta_load_factor = 0.05*Df**2

        return delta_load_factor


    def skin_friction_loss():
        
        delta_load_factor = 2*Cf*Lb/Dhyd*W_U2**2

        return delta_load_factor
        

    def clearence_loss():
        pass
    def mixing_loss():
        pass
    def vaneless_difuser_loss():
        pass
    def disc_friction_loss():
        pass
    def recirculation_loss():
        pass
    def leakage_loss():
        pass

    ##
