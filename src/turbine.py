import numpy as np
from math import tan, pi
from scipy.optimize import minimize

from mfp2mach import mach2mfp

#Inputs: geometry, P0_ratio = P04/P05, M0rotor, gamma
gam = 1.4
beta5 = -20*pi/180
alfa4 =  20*pi/180
A_ratio = 1.2
M0rotor = 0.3
P0_ratio = 0.8

def f(x):
    M4, M5 = x
    MFP5 = mach2mfp(M5, gam)
    MFP4 = mach2mfp(M4, gam) 

    sqrt_T0_ratio = MFP5/MFP4 * A_ratio * P0_ratio 

    phi4 = MFP4/M0rotor*(1+(gam-1)/2*M4**2)**(1/(gam-1))
    phi5 = MFP5/M0rotor*(1+(gam-1)/2*M5**2)**(1/(gam-1))*sqrt_T0_ratio

    psi = 1 - phi4*tan(alfa4) - phi5*tan(beta5)

    P0_ratio_new = (1-(gam-1)*psi*M0rotor**2)**(gam/(gam-1))

    err_P0 = P0_ratio_new - P0_ratio
    err_MFP = MFP4/MFP5 - A_ratio*P0_ratio**((2*gam-1)/gam)

    return err_P0, err_MFP


import matplotlib.pyplot as plt
M = np.linspace(0.3,1,100)
grid = np.meshgrid(M, M)
print(np.shape(grid))

err_P0, err_MFP = f(grid)

ex = 1/(err_P0+P0_ratio)

c = plt.imshow(ex, cmap='RdBu', vmin=np.nanmin(ex), vmax=np.nanmax(ex),
              extent=[0.3,1,0.3,1],
              interpolation='nearest', origin='lower')
plt.colorbar(c)
plt.contour(M,M,err_P0, levels=[0], color='b')
plt.contour(M,M,err_MFP, levels=[0], color='r')

# print(minimize(f, [1,1], bounds=[(0.1,1), (0.1,1)]))

plt.show()




