import numpy as np
from scipy.optimize import fsolve
from .constants import *
from .arrays import *

def fmns(mns,*dat_gS2):
    kF  = dat_gS2[0]
    gS2 = dat_gS2[1]
    EF = np.sqrt(kF**2 + mns**2)
    ns = (g*mns/(2*pi**2)) * ( kF*EF - mns**2*np.log((kF+EF)/mns) )
    func = mns - (m_N - gS2*ns)
    return func
#
def mstar(kF,m_guess,gS2):
    dat_gS2 = (kF,gS2)
    mstar = fsolve(fmns, m_guess, args = dat_gS2)
    mnstar = mstar[0]
    return mnstar
#
def energy(kF,gS2,gW2):
    nb = g*kF**3/3/pi**2
    dat_gS2 = (kF,gS2)
    if kF < 1:
        m_guess = m_N
    else:
        m_guess = 0.1*m_N
    mns = fsolve(fmns, m_guess, args = dat_gS2)
    mnstar = mns[0]
    EF = np.sqrt(kF**2 + mnstar**2)
    ns = (g*mnstar/(2*pi**2)) * ( kF*EF - mnstar**2*np.log((kF+EF)/mnstar) )
    ener = 0.5*gW2*nb**2 + 0.5*gS2*ns**2 \
         + (g/(8*pi**2)) * (2*kF*EF**3 - mnstar**2 * kF*EF \
         - mnstar**4*np.log((kF+EF)/mnstar))    
    return ener
#
def pressure(kF,gS2,gW2):
    nb = g*kF**3/3/pi**2
    dat_gS2 = (kF,gS2)
    if kF < 1:
        m_guess = m_N
    else:
        m_guess = 0.1*m_N
    mns = fsolve(fmns, m_guess, args = dat_gS2)
    mnstar = mns[0]
    EF = np.sqrt(kF**2 + mnstar**2)
    ns = (g*mnstar/(2*pi**2)) * ( kF*EF - mnstar**2*np.log((kF+EF)/mnstar) )
    pres = 0.5*gW2*nb**2 - 0.5*gS2*ns**2 \
         + (1/(g*8*pi**2)) * (((2/3)*kF**3 - mnstar**2 * kF)*EF \
         + mnstar**4*np.log((kF+EF)/mnstar))
    return pres
#
nb = nb_sat
def coupls(gm2,*data):
    kF = data[0]
    gS2 = gm2[0]
    gW2 = gm2[1]
    delkF = kF/500.
    delen = energy(kF+delkF,gS2,gW2) - energy(kF-delkF,gS2,gW2)
    der_en = delen/2./delkF - g*kF**2/pi**2 * energy(kF,gS2,gW2)/nb
    FeP = np.empty((2))
    FeP[0] = der_en
    FeP[1] = energy(kF,gS2,gW2)/nb - m_N + 16.3/hbc
    return FeP
#
target_energy = -16.3
tolerance = 1e-2 
kF = ((3*pi**2*nb)/g)**(1/3)
mnstar = m_N

gS2 = 2 
gW2 = 2
while True:
    data = (kF)
    couplGuess = np.array([gS2, gW2])
    coupl = fsolve(coupls, couplGuess, args=data)
    current_energy = (energy(kF, coupl[0], coupl[1]) / nb - m_N) * hbc

    if abs(current_energy - target_energy) < tolerance:
        break

    gS2 = gS2 + 0.1 * (gS2 - coupl[0])
    gW2 = gW2 + 0.1 * (gW2 - coupl[1])

gs2 = coupl[0]
gw2 = coupl[1]



