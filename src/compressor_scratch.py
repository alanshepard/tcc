from math import tan, pi, sqrt
from warnings import warn
import numpy as np
import matplotlib.pyplot as plt

from mfp2mach import mfp2mach, mach2mfp
import compressor_losses

def fill_nan(a):
    
    ind = np.where(~np.isnan(a))[0]
    first, last = ind[0], ind[-1]
    a[:first] = a[first]
    a[last + 1:] = a[last]
    
    return a

def phi1(MFP, Mrotor, gam):
    
    M1 = mfp2mach(MFP, gam)
    phi1_ = MFP/Mrotor*(1+(gam-1)/2*M1**2)**(1/(gam-1))
    
    return phi1_

def MFP2(MFP, gam, A_ratio, P0_ratio, T0_ratio):
    MFP_ratio = T0_ratio**0.5/(P0_ratio*A_ratio)
    MFP2_ = MFP_ratio*MFP

    # Choke handling
    #MFP_choke = mach2mfp(1,gam)
    #MFP[MFP2_>=MFP_choke] = np.nan #MFP_choke/MFP_ratio[MFP2_>=MFP_choke]
    # MFP2_ = np.minimum(MFP2_, MFP_choke)
    
    return MFP2_

def phi2(MFP2_, Mrotor, gam, T0_ratio):
    
    M2rel = mfp2mach(MFP2_, gam)
    phi2_ = MFP2_/Mrotor * np.sqrt(T0_ratio) * (1+(gam-1)/2*M2rel**2)**(1/(gam-1))

    return phi2_


def psi2T0_ratio(psi_actual, Mrotor, gam):
    return 1+ psi_actual*Mrotor**2*(gam-1)


def psi2P0_ratio(psi_isen, Mrotor, gam):
    return psi2T0_ratio(psi_isen, Mrotor, gam)**(gam/(gam-1))

def compressor_dimensionless(MFP, Mrotor, gam, beta1, beta2, Z, A_ratio, R_ratio):
    slip_factor = 0.9

    # guess P0_ratio and T0_ratio
    MFP2_ = np.ones_like(MFP)
    P0_ratio = np.ones_like(MFP)

    T0_ratio = MFP2_/MFP*A_ratio*P0_ratio

    tol = 1e-8
    iterations = 0
    while True:

        fill_nan(P0_ratio)
        fill_nan(T0_ratio)

        MFP2_ = MFP2(MFP, gam, A_ratio, P0_ratio, T0_ratio)
        phi2_ = phi2(MFP2_, Mrotor, gam, T0_ratio)
        phi1_ = phi1(MFP, Mrotor, gam)

        psi_euler = slip_factor*(1-phi2_*tan(beta2))

        internal_losses = compressor_losses.internal(psi_euler=psi_euler, Z=Z, slip_factor=slip_factor, beta1=beta1, beta2=beta2, phi1=phi1_, phi2=phi2_, R_ratio=R_ratio)

        parasitic_losses = compressor_losses.parasitic()

        psi_isen = psi_euler - sum(internal_losses.values())
        psi_actual = psi_euler + sum(parasitic_losses.values())
        
        P0_ratio_new = psi2P0_ratio(psi_isen, Mrotor, gam)
        T0_ratio_new = psi2T0_ratio(psi_actual, Mrotor, gam)

        if    (np.all(np.abs((P0_ratio_new - P0_ratio)/P0_ratio)<=tol)
           and np.all(np.abs((T0_ratio_new - T0_ratio)/T0_ratio)<=tol)):
            break
        
        if iterations>=50:
            warn("Exceeded number of iterations")
            break

        iterations = iterations+1
        P0_ratio = P0_ratio_new
        T0_ratio = T0_ratio_new
    
    print(iterations)
    return P0_ratio_new, T0_ratio_new, phi2_, psi_euler, MFP2_, internal_losses, parasitic_losses

def plot_contour_map(MFP, M0rotor, gam, beta1, beta2, A_ratio, R_ratio):
    MFP_grid, M0rotor_grid = np.meshgrid(MFP, M0rotor)
    P0_ratio, T0_ratio, phi2_, psi_euler, MFP2_, internal_losses, parasitic_losses = compressor_dimensionless(MFP_grid, M0rotor_grid, gam, beta1, beta2, A_ratio, R_ratio)
    T0_ratio_isen = P0_ratio**((gam-1)/gam)
    eff = np.log(T0_ratio_isen)/np.log(T0_ratio)

    cs_P0 = plt.contour(MFP_grid, P0_ratio, M0rotor_grid)
    plt.clabel(cs_P0)
    cs_eff = plt.contour(MFP_grid, P0_ratio, psi_euler, 30)
    plt.clabel(cs_eff)

    plt.title("Compressor map")
    plt.xlabel("MFP")
    plt.ylabel("P02/P01")
    plt.legend(("M0rotor", "eta_p"))
    
    return P0_ratio, eff


def plot_map(MFP, M0rotor, gam, beta1, beta2, Z, A_ratio, R_ratio):

    fig, ax= plt.subplots(2,1)
    fig2, ax2 = plt.subplots(1, 1)
    fig3, ax3 = plt.subplots(1, 1)

    # Loop through isospeed lines
    for M in M0rotor:

        # Calculate performance
        P0_ratio, T0_ratio, phi2_, psi, MFP2_, internal_losses, parasitic_losses = compressor_dimensionless(MFP, M, gam, beta1, beta2, Z, A_ratio, R_ratio)
        T0_ratio_isen = P0_ratio ** ((gam - 1) / gam)
        eff = np.log(T0_ratio_isen) / np.log(T0_ratio)

        # Isospeeds MFP x P0ratio
        i = np.nanargmax(P0_ratio)
        j = np.nanargmax(eff) 
        ax[1].plot(MFP[i], P0_ratio[i], 'ro')
        ax[1].plot(MFP[j], P0_ratio[j], 'go')
        ax[1].plot(MFP[i:], P0_ratio[i:])
        ax[1].grid(True)
        ax[1].set_ylabel("P02/P01")
        ax[1].set_xlabel("MFP")

        # Plot MFP x polytropic efficiency
        ax[0].plot(MFP[i:], eff[i:])
        ax[0].set_ylim([0,1])
        ax[0].set_ylabel("Polytropic efficiency")


        # Isospeeds MFP x MFP2
        ax2.plot(MFP, MFP2_)
        ax2.set_ylabel("MFP2")
        ax2.set_xlabel("MFP")
        # Choke box
        MFP_choke = mach2mfp(1,gam)
        ax2.hlines(MFP_choke, 0, MFP_choke)
        ax2.vlines(MFP_choke, 0, MFP_choke)

        # flow vs work coefficients
        ax3.plot(MFP, psi)

    return ax, ax2




b = 5e-3
D2 = 64e-3
D1t = 42.8e-3
D1h = 17e-3

A1 = pi/4*(D1t**2-D1h**2)
A2 = pi*D2*b

slip_factor=0.9
gam = 1.4
MFP_choke = mach2mfp(1,gam)
MFP = np.linspace(0,MFP_choke,300)
M = mfp2mach(MFP, gam)
beta2 = 10*pi/180
beta1 = 40*pi/180
A_ratio = A2/A1
Z=10
print("A_ratio", A_ratio)
R_ratio = D2/sqrt((D1h**2 + D1t**2)/2)
print("R_ratio", R_ratio)

# for loss_name, delta_psi in internal_losses.items():
#     plt.plot(phi2_, delta_psi, label=loss_name)
#     plt.legend()


# plt.figure()
# plot_contour_map(MFP, np.linspace(0.9,1.5,100), gam, beta1, beta2, A_ratio=2, R_ratio=2)

M0rotor = np.linspace(0.2, 1.3, 10)
# MFP_choke, P0_ratio, T0_ratio = choke_line(M0rotor, gam, beta1, beta2, A_ratio, R_ratio)
ax, ax2 = plot_map(MFP, M0rotor, gam, beta1, beta2, Z, A_ratio, R_ratio)

# ax[1].plot(MFP_choke, P0_ratio, "o")
# ax2.plot(MFP_choke, MFP2(MFP_choke, gam, A_ratio, P0_ratio, T0_ratio), "o")

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

