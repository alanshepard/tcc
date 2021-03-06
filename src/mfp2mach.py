import numpy as np
from warnings import warn

def mach2mfp(M, gam):
    """
    Mass flow parameter

    MFP = mdot*sqrt(R*T0)/(P0*A*sqrt(gam)) = mdot/(rho0*A*a0)
    """
    return M*(1+(gam-1)/2*M**2)**(-(gam+1)/(2*(gam-1)))

def mfp2mach(MFP, gam, tol=1e-13):
    """ 
    Only works for subsonic
    
    Reference:DER 1974
    """
    
    MFP_choked = mach2mfp(1,gam)
    #Anything with MFP>MFP_choked is choked, so
    # MFP[MFP>MFP_choked] = MFP_choked

    def z(MFP):
        z_= -np.sqrt(1-MFP/MFP_choked)
        return z_
    
    def dM_dz(M):
        num = 2*z(mach2mfp(M,gam))*((2+(gam-1)*M**2)/(gam+1))**((gam+1)/(2*(gam-1)))
        den = ((gam+1)*M**2)/(2+(gam-1)*M**2)-1

        return np.where(M<1, np.divide(num, den, where=M!=1), 2/(gam+1))
    
    #Newton method
    iterations = 0
    M = 1-z(MFP)**2 #initial guess
    while(np.max(np.abs(z(mach2mfp(M,gam))-z(MFP))) > tol):
        
        if iterations >= 100:
            warn("Exceeded maximum iterations")
            break

        iterations = iterations+1
        
        M = M-(z(mach2mfp(M,gam))-z(MFP))*dM_dz(M)

    return M

if __name__ == "__main__":
    M=np.linspace(0,1,1000)
    gam=1.4
    M_err = M - mfp2mach(mach2mfp(M,gam),gam,tol=1e-13)
    print(np.c_[M,M_err])
    assert np.max(np.abs(M_err))<=1e-13
    print("TEST OK")
    
    import matplotlib.pyplot as plt
    M=np.linspace(0,2,100)
    plt.plot(M, mach2mfp(M, gam))
    plt.show()
