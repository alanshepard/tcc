import numpy as np
from math import tan, pi
from scipy.optimize import minimize, root

from mfp2mach import mach2mfp

#Inputs: geometry, P0_ratio = P04/P05, M0rotor, gamma
gam = 1.4
beta5 = -56*pi/180
alfa4 = 60*pi/180
A_ratio = (65.3**2-45**2)/(67.**2-46.**2)
M0rotor = 0.8

def f(x, P0_ratio):
    M4, M5 = x
    MFP5 = mach2mfp(M5, gam)
    MFP4 = mach2mfp(M4, gam) 

    sqrt_T0_ratio = MFP5/MFP4 * A_ratio * P0_ratio 

    phi4 = MFP4/M0rotor*(1+(gam-1)/2*M4**2)**(1/(gam-1))
    phi5 = MFP5/M0rotor*(1+(gam-1)/2*M5**2)**(1/(gam-1))*sqrt_T0_ratio

    psi = 1 + phi5*tan(beta5) - phi4*tan(alfa4)

    P0_ratio_new = (1+(gam-1)*psi*M0rotor**2)**(gam/(gam-1))

    err_P0 = P0_ratio_new - P0_ratio
    err_MFP = MFP4/MFP5 - A_ratio*P0_ratio**((2*gam-1)/(gam))

    return P0_ratio_new, err_MFP, psi


import matplotlib.pyplot as plt

M = np.linspace(5e-2,1,10)
grid = np.meshgrid(M, M)
print(np.shape(grid))
P0_ratio=0.9
P0, err_MFP, psi = f(grid, P0_ratio)

ex = P0

c = plt.imshow(P0, cmap='RdBu', vmin=0, vmax=1.2,
              extent=[5e-2,1,5e-2,1],
              interpolation='nearest', origin='lower')
plt.colorbar(c)
plt.contour(M,M,psi, levels=[0], color='red')
plt.contour(M,M,P0, levels=[P0_ratio], color='red')
plt.contour(M,M,err_MFP, levels=[0], color='red')
plt.figure(2)
d = plt.imshow(psi, cmap='RdBu', vmin=np.min(psi), vmax=np.max(psi),
              extent=[5e-2,1,5e-2,1],
              interpolation='nearest', origin='lower')
plt.colorbar(d)
# plt.contour(M,M,errP0, levels=[0], color='blue')
plt.contour(M,M,err_MFP, levels=[0], color='red')
plt.contour(M,M,P0) 

def g(x, P):
    P0_ratio_new, err_MFP, psi = f(x, P)
    return P-P0_ratio_new, err_MFP

fig, ax = plt.subplots(2,2)
for M0rotor in np.linspace(0.1, 0.5, 10):
    MFP1=[]
    MFP2=[]
    P=[]
    x0=(0.5,0.5)
    for p in np.linspace(1,1/3,100):
        sol = root(g,x0,p)
        print(sol)
        if not sol.success or np.any(sol.x>=1):
            break
        x0 = sol.x
        P.append(p)
        MFP1.append(mach2mfp(sol.x[0], gam))
        MFP2.append(mach2mfp(sol.x[1], gam))
    # P[-1]=0.5
    ax[0][0].plot(MFP1,MFP2)
    ax[1][0].plot(np.array(list(reversed(MFP1))),1/np.array(list(reversed(P))))
    ax[0][1].plot(1/np.array(list(reversed(P))),np.array(list(reversed(MFP2))))
    

# Choke box
MFP_choke = mach2mfp(1,gam)
ax[0][0].hlines(MFP_choke, 0, MFP_choke)
ax[0][0].vlines(MFP_choke, 0, MFP_choke)
ax[0][0].set_ylabel("MFP_5")
ax[1][0].vlines(MFP_choke, 0, 3)
ax[1][0].set_xlabel("MFP_4")
ax[1][0].set_ylabel("P04/P05")
ax[0][1].set_xlabel("P04/P05")

plt.delaxes(ax[1][1])

# print(minimize(f, [1,1], bounds=[(0.1,1), (0.1,1)]))
# print(root(f,[1,1]))

plt.show()




