from math import tan, pi
from warnings import warn
import numpy as np


from mfp2mach import mfp2mach, mach2mfp

def phi1(MFP, Mrotor, gam):
    
    M1 = mfp2mach(MFP, gam)
    phi1_ = MFP/Mrotor*(1+(gam-1)/2*M1**2)**(-1/gam-1)
    
    return phi1_

def MFP2(MFP, gam, A_ratio, P0_ratio, T0_ratio):
    MFP_ratio = T0_ratio**0.5/(P0_ratio*A_ratio)
    MFP2_ = MFP_ratio*MFP
    MFP2_ = np.minimum(MFP2_, mach2mfp(1,gam))
    return MFP2_

def phi2(MFP, Mrotor, gam, A_ratio, P0_ratio, T0_ratio):
    
    MFP2_ = MFP2(MFP, gam, A_ratio, P0_ratio, T0_ratio)
    M2rel = mfp2mach(MFP2_, gam)
    phi2_ = MFP/Mrotor * T0_ratio/(P0_ratio*A_ratio) * (1+(gam-1)/2*M2rel**2)**(-1/(gam-1))

    return phi2_

#TODO: fix or remove this
def M2(MFP, Mrotor, gam, P0_ratio, T0_ratio, A_ratio, slip_factor, beta2):
    tol=1e-8

    MFP2_ = MFP2(MFP, gam, P0_ratio, T0_ratio, A_ratio)
    M2rel = mfp2mach(MFP2_, gam)
    phi2_ = phi2(MFP, Mrotor, gam, P0_ratio, T0_ratio, A_ratio)
 
    M2_ = M2rel #initial_guess
    while True:
        M2new = Mrotor/np.sqrt(T0_ratio*(1+(gam-1)/2*M2_**2))*np.sqrt(phi2_**2 + (slip_factor*(1-phi2_*tan(beta2)))**2)

        if np.max(np.abs(M2new-M2_))<=tol:
            break

        M2_=M2new

    return M2new



def psi_euler(phi2, beta2, slip_factor):
    
    psi = slip_factor*(1-phi2*tan(beta2))
    
    return psi


def psi2T0_ratio(psi_actual, Mrotor, gam):
    return 1+ psi_actual*Mrotor*2*(gam-1)


def psi2P0_ratio(psi_isen, Mrotor, gam):
    return psi2T0_ratio(psi_isen, Mrotor, gam)**(gam/(gam-1))

def compressor_dimensionless(MFP, Mrotor, gam, beta2, A_ratio):
    slip_factor = 0.9

    phi1_ = phi1(MFP, Mrotor, gam)

    # guess P0_ratio and T0_ratio
    P0_ratio = 1
    T0_ratio = 1

    tol = 1e-8
    iterations = 0
    while True:
        phi2_ = phi2(MFP, Mrotor, gam, A_ratio, P0_ratio, T0_ratio)

        psi = psi_euler(phi2_, beta2, slip_factor)

        P0_ratio_new = psi2P0_ratio(psi, Mrotor, gam)
        T0_ratio_new = psi2T0_ratio(psi, Mrotor, gam)

        if    (np.all(np.abs((P0_ratio_new - P0_ratio)/P0_ratio)<=tol)
           and np.all(np.abs((T0_ratio_new - T0_ratio)/T0_ratio)<=tol)):
            break
        
        if iterations>=50:
            warn("Exceeded number of iterations")
            break

        iterations = iterations+1
        P0_ratio = P0_ratio_new
        T0_ratio = T0_ratio_new
    
    return P0_ratio, T0_ratio



import matplotlib.pyplot as plt

print(compressor_dimensionless(0, 0.6, 1.4, 30*pi/180, 1))

slip_factor=0.9
gam = 1.4
MFP_choke = mach2mfp(1,gam)
MFP = np.linspace(0,MFP_choke,100)
M = mfp2mach(MFP, gam)
Mrotor = 0.12
beta2 = 15*pi/180
A_ratio = 1
P0_ratio, T0_ratio = compressor_dimensionless(MFP, Mrotor, gam, beta2, A_ratio)
phi2_ = phi2(MFP, Mrotor, gam, A_ratio, P0_ratio, T0_ratio)
phi1_ = phi1(MFP, Mrotor, gam)
plt.plot(MFP, P0_ratio)
plt.plot(MFP, T0_ratio)
plt.figure()
plt.plot(MFP, phi1_)
plt.plot(MFP, phi2_)
plt.figure()
plt.plot(MFP, MFP2(MFP, gam, A_ratio, P0_ratio, T0_ratio))
plt.plot(MFP, M2(MFP, Mrotor, gam, P0_ratio, T0_ratio, A_ratio, slip_factor, beta2))
plt.hlines(MFP_choke, 0, MFP_choke)
plt.vlines(MFP_choke, 0, MFP_choke)

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
    





