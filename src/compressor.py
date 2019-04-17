from math import pi, tan, cos, sin
from numpy import sqrt
import numpy as np

def norm(x, y, z):
    return sqrt(x**2+y**2+z**2)

class Compressor:
    def __init__(self):
        #Fluid (air)
        self.gam = 1.4
        self.Cp = 1004.
        self.T1 = 275.

        # Geometry for VTR 160L
        
        #Impeller
        self.beta1 = 40*pi/180 # blade leading edge angle with respect to axial direction
        self.beta2 = 0.         # blade trailing edge angle with respect to radial direction
                                # Positive in direction opposite to rotation
        self.D1t = 0.106        # blade leading edge diameter at hub
        self.D1h = 0.054        # blade leading edge diameter at hub
        self.D2  = 0.180        # blade trailing edge diameter
        self.b2  = 0.007        # exit blade width
        self.z_fb  = 20.        # full length blades
        self.z_sb  = 0.         # splitter blades
        self.L     = 1.         # blade mean streamline length
        self.L_sb  = 1.         # splitter blade mean streamline length

        # Diffuser
        self.alfa2 = -62*pi/180 # vane leading edge angle with respect to radial direction
                                # Positive in direction opposite to rotation
        self.alfa3 = -50*pi/180 # vane leading edge angle with respect to radial direction
                                # Positive in direction opposite to rotation
        self.D3 = 0.215         # Inlet diameter
        self.D4 = 0.258         # Outlet diameter

    def __call__(self, flow_rate, omega):

        # Run adimensional model
        phi = self.flow_param_exit(flow_rate, omega)
        psi = self.load_param_euler(phi)
        eff_isen = self.eff_isen(phi)

        # Dimensionalize
        U2  = omega * self.D2/2
        delta_h = psi*U2**2

        # Get pressure ratio
        T0_ratio = 1. + delta_h/(self.Cp*self.T1)
        P0_ratio = (eff_isen*T0_ratio)**(self.gam/(self.gam-1))

        return P0_ratio, T0_ratio, eff_isen

    @property
    def D1(self):
        """
        Area averaged impeller diameter
        """
        return (self.D1t**2+self.D1h**2)**0.5/2

    @property
    def A1(self):
        """
        Anular entrance area
        """
        return pi*((self.D1t/2)**2 - (self.D1h/2)**2)

    @property
    def A2(self):
        """
        ring shaped exit area
        """
        return 2*pi*(self.D2/2)**2*self.b2

    @property
    def z(self):
        """ 
        Effective number of blades (AUNGIER, 1995)
        """

        if self.z_sb == 0:
            return self.z_fb

        return self.z_fb + self.z_sb*self.L_sb/self.L

    @property
    def slip_factor(self):
        """
        AUNGIER, 1995
        """
        s = 1 - sqrt(cos(self.beta2))/self.z**0.7

        #TODO: adicionar limite

        return s

    def load_param_euler(self, flow_param_exit):
        
        psi = self.slip_factor * (1-flow_param_exit*tan(self.beta2))

        return psi

    def flow_param_exit(self, flow_rate, omega):
        
        W2r = flow_rate/self.A2
        U2 = omega*self.D2/2
        phi = W2r/U2
        
        return phi

    def eff_isen(self, flow_param_exit):
        losses = [self.impeller_incidence_loss(flow_param_exit)]#, self.diffuser_incidence_loss(flow_param_exit)]
        eff = self.load_param_euler(flow_param_exit)/(sum(losses)+self.load_param_euler(flow_param_exit))
        return eff


    def impeller_distortion_factor(self, phi2):
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

        B2 = 0 #TODO!
        l = 1/(1-B2)
        return l


    ### Decomposition of W2 ###
    def  W2x_U2(self, phi2):
        return 0.

    def W2r_U2(self, phi2):
        """
        W2r = mdot/(rho*A2) = phi2*U2
        """
        return phi2 # BY DEFINITION!

    def W2u_U2(self, phi2):
        return self.V2u_U2(phi2)-1

    def W2_U2(self, phi2):
        return norm(self.W2r_U2(phi2),
                    self.W2x_U2(phi2),
                    self.W2u_U2(phi2)
                    )

    ### Decomposition of V2 ###
    
    def V2x_U2(self):
        return 0.

    def V2u_U2(self, phi2):    
        """
        V2u = sigma*(1-lambda*phi2*tan(beta2))*U2
        where sigma is the slip factor 
        and lambda is the impeller distortion factor
        """

        return self.slip_factor*(1-self.impeller_distortion_factor(phi2)*phi2*tan(self.beta2))
    
    def V2r_U2(self, phi2):
        return self.W2r_U2(phi2)

    def V2_U2(self, phi2):
        return norm(self.V2x_U2(),
                    self.V2u_U2(phi2),
                    self.V2r_U2(phi2)
                    )

    ##

    ### Decomposition of W1 ###
    def W1r_U2(self):
        return 0.

    def W1u_U2(self):
        """
        W1tu = U1t = D1/D2*U2
        """
        return self.D1/self.D2

    def W1x_U2(self, phi2):
        """
        W1x = V1x
        """
        return self.V1x_U2(phi2)

    def W1_U2(self, phi2):
        return norm(self.W1r_U2(),
                    self.W1u_U2(),
                    self.W1x_U2(phi2)
                    )

    ### Decomposition of W1t ###
    def W1tr_U2(self):
        return 0.

    def W1tu_U2(self):
        """
        W1tu = U1t = D1t/D2*U2
        """
        return self.D1t/self.D2

    def W1tx_U2(self, phi2):
        """
        W1tx = V1x
        """
        return self.V1x_U2(phi2)

    def W1t_U2(self, phi2):
        return norm(self.W1tr_U2(),
                    self.W1tu_U2(),
                    self.W1tx_U2(phi2)
                    )

    ##

    ### Decomposition of W1h
    def W1hr_U2(self):
        return 0.

    def W1hu_U2(self):
        """
        W1hu = omega*D1h = D1h/D2*U2
        """
        return self.D1h/self.D2

    def W1hx_U2(self, phi2):
        """
        W1hx = V1x
        """
        return self.V1x_U2(phi2)

    def W1rht_U2(self, phi2):
        return norm(self.W1hr_U2(),
                    self.W1hx_U2(phi2),
                    self.W1hu_U2()
                    )

    ##
    
    ### Decomposition of V1 ###
    
    def V1r_U2(self): 
        return 0.

    def V1u_U2(self): 
        return 0.

    def V1x_U2(self, phi2):
        """
        V1x=mdot/(rho*A1)=A1/A2*phi2*U2
        """
        return self.A1/self.A2*phi2
    
    def V1_U2(self, phi2):
        return norm(self.V1x_U2(phi2),
                    self.V1u_U2(),
                    self.V1r_U2()
                    )
        
    ##

    ### Losses ###
    #All losses written as a decrement on the load parameter
    def impeller_incidence_loss(self, phi2):
        """
        GRAVDAHL 1999, GALVAS 1973
        NASA shock loss theory: kinetic energy normal to blade is destroyed
        """
    
        #TODO: FIX PROBLEM HERE!!!
        beta1_flow = np.arctan2(self.W1u_U2(), self.W1x_U2(phi2))
        W1normal_U2 = self.W1_U2(phi2)*sin(beta1_flow-self.beta1)
        assert abs(W1normal_U2)<= self.W1_U2(phi2)
        
        delta_load_factor = (W1normal_U2)**2/2
        
        return delta_load_factor

    def diffuser_incidence_loss(self, phi2):
        """
        RAVDAHL 1999, GALVAS 1973
        NASA shock loss theory: kinetic energy normal to blade is destroyed
        """
        Unormal_U2 = self.V2r_U2(phi2)*sin(self.alfa2)+self.V2u_U2(phi2)*cos(self.beta2)
        
        delta_load_factor = Unormal_U2**2/2

        return delta_load_factor

    def blade_loading_loss(self, load_param_euler):

        D1t_D2 = self.D1t/self.D2

        W2_W1u = (phi2/cos(beta2))/(sqrt((self.D1t/self.D2)**2 + (phi2*(self.A1/self.A2))**2))
        
        Df = 1 - W2/W1u*(1 - (0.75*load_param_euler)
                             /(self.n_blades/pi*(1-D1t_D2)+2*(D1t_D2)))

        delta_load_factor = 0.05*Df**2

        return delta_load_factor


    def skin_friction_loss(self):
        
        delta_load_factor = 2*Cf*Lb/Dhyd*W_U2**2

        return delta_load_factor
        

    def clearence_loss(self):
        pass
    def mixing_loss(self):
        pass
    def vaneless_diffuser_loss(self):
        pass
    def disc_friction_loss(self):
        pass
    def recirculation_loss(self):
        pass
    def leakage_loss(self):
        pass

    ##


import matplotlib.pyplot as plt
c = Compressor()
R = 287
P01 = 101e3
x=[]
y=[]
for q in np.linspace(0, 0.005):
    P_ratio, T_ratio, eff_isen = c(q, 35e3*2*pi/60)
    T2 = c.T1*T_ratio
    P2 = P01*P_ratio
    rho = P2/(R*T2)
    mdot = q*rho
    print(mdot, P_ratio)
    x.append(mdot)
    y.append(P_ratio)

plt.plot(x, y)
plt.show()