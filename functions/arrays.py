import numpy as np
from .constants import *
#
nB       = np.logspace(np.log10(nB_min), np.log10(nB_max), num_nB)
epsilon  = np.array([0. for n in range(0, num_nB)])
E0       = np.array([0. for n in range(0, num_nB)])
E0_coupl = np.array([0. for n in range(0, num_nB)])
m_Ns     = np.array([0. for n in range(0, num_nB)])
m_Nsapp  = np.array([0. for n in range(0, num_nB)])
m_NsI    = np.array([0. for n in range(0, num_nB)])
n0       = np.array([0. for n in range(0, num_nB)])
gm2      = np.array([0. for n in range(0,1)])