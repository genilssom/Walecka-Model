import numpy as np
from scipy.optimize import fsolve
from .constants import *
from .arrays import *
from .constant_couplings import *

def fmn_star(mn_star):
    EF = np.sqrt(kF**2 + mn_star**2)
    ns = g*mn_star/(2*pi**2) * (kF*EF - mn_star**2*np.log((kF+EF)/mn_star))
    funct = mn_star - (m_N - gs2*ns)
    return funct
#
for i in range(num_nB):
    kF = ((3*(np.pi**2)*nB[i])/g)**(1/3)
    if (i==0):
        MN_star = fsolve(fmn_star, m_N)
    else:
        MN_star = fsolve(fmn_star, MN_star)
    m_Ns[i] = MN_star/ m_N
    m_Nsapp[i] = 1. / (1 + gs2*(kF**2/pi**2))
    EF = np.sqrt(kF**2 + MN_star**2)
    ns = g*MN_star/(2*pi**2) * (kF*EF - MN_star**2*np.log((kF+EF)/MN_star))
    epsilon[i] = 0.5*gw2*nB[i]**2 + 0.5*gs2*ns**2 + (g/(8*pi**2)) \
                * ((2*kF**3 + MN_star**2 * kF)*EF - MN_star**4*np.log((kF+EF)/MN_star))
    
    E0[i] = ((epsilon[i]/nB[i]) - m_N)*hbc
    iemin = np.argmin(E0)
    n0[i] = nB[i]/nb_sat
    Emin = E0[iemin]
    nBmin = nB[iemin]