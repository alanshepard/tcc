import numpy as np
from math import tan, pi
import inspect
from warnings import warn

def internal(**kwargs):

    #Add kwargs
    kwargs["Cf"] = 0.005
    kwargs["Lb_Dhyd"] = 5
    kwargs["D1t_D2"] = 1/kwargs["R_ratio"]
    kwargs["D1h_D2"] = 1/kwargs["R_ratio"]


    internal_loss_func_dict = {
        "inducer_incidence_loss": inducer_incidence,
        "diffuser_incidence_loss": diffuser_incidence,
        "skin_friction_loss": skin_friction,
        "blade_loading_loss": blade_loading,
    }

    internal_loss_dict = eval_loss_func_dict(internal_loss_func_dict, **kwargs)
    # assert not np.any(np.isnan(list(internal_loss_dict.values())))
    
    return internal_loss_dict

def parasitic(**kwargs):
    parasitic_loss_func_dict = {
    }

    parasitic_loss_dict = eval_loss_func_dict(parasitic_loss_func_dict, **kwargs)

    return parasitic_loss_dict
    
def eval_loss_func_dict(loss_func_dict, **available_kwargs):
    loss_dict = {}
    for name, func in loss_func_dict.items():
        expected_kwargs = inspect.getargspec(func).args
        kwargs = {k:available_kwargs[k] for k in expected_kwargs if k in available_kwargs}
        loss_dict[name] = func(**kwargs)
    return loss_dict


def inducer_incidence(phi1, beta1, R_ratio):
    
    beta1_opt = np.arctan(1/(R_ratio*phi1))
    delta_psi = ((1/R_ratio)**2+phi1**2)*np.sin(beta1-beta1_opt)**2

    return delta_psi

def diffuser_incidence(psi_euler, phi2, alfa3):
    alfa3_opt = np.arctan(psi_euler/phi2)
    delta_psi = (psi_euler**2+phi2**2)*np.sin(alfa3-alfa3_opt)**2

    return delta_psi

def skin_friction(Cf, Lb_Dhyd, slip_factor, beta2, phi1, phi2, D1t_D2, D1h_D2):
    #every speed is non-dimensionalized by U2
    V2 = np.sqrt(phi2**2 + (slip_factor*(1-phi2*tan(beta2)))**2)
    W1t = np.sqrt(D1t_D2**2 + phi1**2)
    W1h = np.sqrt(D1h_D2**2 + phi1**2)
    W2 = phi2*np.sqrt(1+tan(beta2)**2)

    Wbar = (phi1 + V2+ W1t + 2*W1h + 3*W2)/8

    delta_psi = 2*Cf*Lb_Dhyd*Wbar**2

    return delta_psi

def blade_loading(psi_euler, phi1, phi2, beta2, Z, D1t_D2, D1h_D2):
    W1t = np.sqrt(D1t_D2**2 + phi1**2)
    W1h = np.sqrt(D1h_D2**2 + phi1**2)
    W2 = phi2*np.sqrt(1+tan(beta2)**2)

    Df = 1 - W2/W1t + (0.75*psi_euler)/((W1t/W2)*(Z/pi)*(1-D1t_D2)+2*D1t_D2)
    
    delta_psi = 0.05*Df**2

    return delta_psi

def clearence():
    pass