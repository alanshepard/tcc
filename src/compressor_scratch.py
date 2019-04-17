from math import tan, pi
from warnings import warn
import numpy as np
import matplotlib.pyplot as plt

from mfp2mach import mfp2mach, mach2mfp

def phi1(MFP, Mrotor, gam):
    
    M1 = mfp2mach(MFP, gam)
    phi1_ = MFP/Mrotor*(1+(gam-1)/2*M1**2)**(1/(gam-1))
    
    return phi1_

def MFP2(MFP, gam, A_ratio, P0_ratio, T0_ratio):
    MFP_ratio = T0_ratio**0.5/(P0_ratio*A_ratio)
    MFP2_ = MFP_ratio*MFP
    # MFP2_ = np.minimum(MFP2_, mach2mfp(1,gam))
    return MFP2_

def phi2(MFP2_, Mrotor, gam, T0_ratio):
    
    M2rel = mfp2mach(MFP2_, gam)
    phi2_ = MFP2_/Mrotor * np.sqrt(T0_ratio) * (1+(gam-1)/2*M2rel**2)**(1/(gam-1))

    return phi2_


def psi2T0_ratio(psi_actual, Mrotor, gam):
    return 1+ psi_actual*Mrotor*2*(gam-1)


def psi2P0_ratio(psi_isen, Mrotor, gam):
    return psi2T0_ratio(psi_isen, Mrotor, gam)**(gam/(gam-1))

def compressor_dimensionless(MFP, Mrotor, gam, beta1, beta2, A_ratio, R_ratio):
    slip_factor = 0.9

    # guess P0_ratio and T0_ratio
    P0_ratio = 1
    T0_ratio = 1

    tol = 1e-8
    iterations = 0
    considering_losses = False
    while True:
        MFP2_ = MFP2(MFP, gam, A_ratio, P0_ratio, T0_ratio)
        phi2_ = phi2(MFP2_, Mrotor, gam, T0_ratio)
        phi1_ = phi1(MFP, Mrotor, gam)

        psi_euler = slip_factor*(1-phi2_*tan(beta2))

        internal_losses = {
            "incidence": incidence_loss(phi1_, beta1, R_ratio),
            "skin friction": skin_friction_loss(0.004, 20, slip_factor, beta2, phi1_, phi2_,1/R_ratio , 1/R_ratio)
        }

        parasitic_losses ={}

        if considering_losses:
            psi_isen = psi_euler - sum(internal_losses.values())
            psi_actual = psi_euler + sum(parasitic_losses.values())
        else:
            psi_isen=psi_euler
            psi_actual = psi_euler
        
        P0_ratio_new = psi2P0_ratio(psi_isen, Mrotor, gam)
        T0_ratio_new = psi2T0_ratio(psi_actual, Mrotor, gam)

        if    (np.all(np.abs((P0_ratio_new - P0_ratio)/P0_ratio)<=tol)
           and np.all(np.abs((T0_ratio_new - T0_ratio)/T0_ratio)<=tol)):
            if considering_losses:
                break
            else:
                considering_losses=True
                print("Converged without losses, iter={}".format(iterations))
        
        if iterations>=50:
            warn("Exceeded number of iterations")
            break

        iterations = iterations+1
        P0_ratio = P0_ratio_new
        T0_ratio = T0_ratio_new
    
    MFP_choke = mach2mfp(1,gam)
    choked = (MFP>MFP_choke) | (MFP2_ >MFP_choke)

    print(iterations)
    return P0_ratio, T0_ratio, choked, phi2_, internal_losses, parasitic_losses

def plot_contour_map(MFP, M0rotor, gam, beta1, beta2, A_ratio, R_ratio):
    MFP_grid, M0rotor_grid = np.meshgrid(MFP, M0rotor)
    P0_ratio, T0_ratio, choked, phi2_, internal_losses, parasitic_losses = compressor_dimensionless(MFP_grid, M0rotor_grid, gam, beta1, beta2, A_ratio, R_ratio)
    T0_ratio_isen = P0_ratio**((gam-1)/gam)
    eff = np.log(T0_ratio_isen)/np.log(T0_ratio)

    #P0_ratio[choked] = 0
    #eff[choked] = 0

    cs_P0 = plt.contour(MFP_grid, P0_ratio, M0rotor_grid)
    plt.clabel(cs_P0)
    cs_eff = plt.contour(MFP_grid, P0_ratio, eff, 30)
    plt.clabel(cs_eff)
    #cs_choked = plt.contour(MFP_grid, P0_ratio, choked, 1)
    #print(choked)

    plt.title("Compressor map")
    plt.xlabel("MFP")
    plt.ylabel("P02/P01")
    plt.legend(("M0rotor", "eta_p"))
    
    return P0_ratio, eff


def plot_map(MFP, M0rotor, gam, beta1, beta2, A_ratio, R_ratio):

    fig, ax= plt.subplots(2,1)
    fig2, ax2 = plt.subplots(1, 1)

    for M in M0rotor:
        P0_ratio, T0_ratio, choked, phi2_, internal_losses, parasitic_losses = compressor_dimensionless(MFP, M, gam, beta1, beta2, A_ratio, R_ratio)
        T0_ratio_isen = P0_ratio ** ((gam - 1) / gam)
        eff = np.log(T0_ratio_isen) / np.log(T0_ratio)
        ax[0].plot(MFP, eff)
        ax[1].plot(MFP, P0_ratio)
        ax2.plot(MFP, MFP2(MFP, gam, A_ratio, P0_ratio, T0_ratio))


def incidence_loss(phi1_, beta1, R_ratio):
    
    beta1_opt = np.arctan(1/(R_ratio*phi1_))
    delta_psi = ((1/R_ratio)**2+phi1_**2)*np.sin(beta1-beta1_opt)**2

    return delta_psi

def skin_friction_loss(Cf, Lb_Dhyd, slip_factor, beta2, phi1_, phi2_, D1t_D2, D1h_D2):
    #every speed is adimentionalized by U2
    V2 = np.sqrt(phi2_**2 + (slip_factor*(1-phi2_*tan(beta2)))**2)
    W1t = np.sqrt(D1t_D2**2 + phi1_**2)
    W1h = np.sqrt(D1h_D2**2 + phi1_**2)
    W2 = phi2_*np.sqrt(1+tan(beta2)**2)

    Wbar = (phi1_ + V2+ W1t + 2*W1h + 3*W2)/8

    return 2*Cf*Lb_Dhyd*Wbar**2


slip_factor=0.9
gam = 1.4
MFP_choke = mach2mfp(1,gam)
MFP = np.linspace(0,MFP_choke,1000)
M = mfp2mach(MFP, gam)
Mrotor = 0.3
beta2 = 30*pi/180
beta1 = 40*pi/180
A_ratio = 2
R_ratio = 2
P0_ratio, T0_ratio, choked, phi2_, internal_losses, parasitic_losses = compressor_dimensionless(MFP, Mrotor, gam, beta1, beta2, A_ratio, R_ratio)
MFP2_ = MFP2(MFP, gam, A_ratio, P0_ratio, T0_ratio)
phi2_ = phi2(MFP2_, Mrotor, gam, T0_ratio)
phi1_ = phi1(MFP, Mrotor, gam)

for loss_name, delta_psi in internal_losses.items():
    plt.plot(phi2_, delta_psi, label=loss_name)
    plt.legend()

plt.figure()
plt.plot(MFP, P0_ratio)
plt.plot(MFP, T0_ratio)
plt.figure()
# plt.plot(MFP, phi1_)
plt.plot(MFP, phi2_)
plt.figure()
plt.plot(MFP, MFP2_)
plt.xlabel("MFP1")
plt.xlabel("MFP2")
# plt.plot(MFP, M2(MFP, Mrotor, gam, P0_ratio, T0_ratio, A_ratio, slip_factor, beta2))
plt.hlines(MFP_choke, 0, MFP_choke)
plt.vlines(MFP_choke, 0, MFP_choke)
plt.figure()
P0_ratio, eff = plot_contour_map(MFP, np.linspace(0.3,1.5,100), gam, beta1, beta2, A_ratio, R_ratio)
plot_map(MFP, np.linspace(0.3,1,10), gam, beta1, beta2, A_ratio, R_ratio)

plt.show()



class Compressor:
    def __init__(self):
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
    




#TODO: plot individal losses, return losses in a map, return relevant params (e.g. phi2)
#Improve visualization of what is happening

