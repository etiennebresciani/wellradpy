# -*- coding: utf-8 -*-
"""
Created on 2020/03/20

@author: Etienne Bresciani
"""

from wellradpy import recovery as re
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Parameters
###############################################################################

T = 10. # Transmissivity
S = 1.e-4 # Storativity
Q = 100. # Pumping rate
tp = 1. # Pumping duration
sc = 0.01 # Apparent resolution

###############################################################################
# Radius of investigation and associated quantities
###############################################################################

# Maximum radius of investigation
rinvmax = re.rinvmax(T, S, Q, tp, sc)

# Time at which the radius of investigation is maximum
tmax = re.tmax(T, Q, tp, sc)

# Termination time of the recovery test
tend = re.tend(T, Q, tp, sc)

# Radius of investigation as a function of time (counted from the beginning of
# pumping)
Nt = 1000
t = np.linspace(1.0001*tp, tp+0.99999*(tend-tp), Nt)
rinv = np.full(t.shape, np.nan)
for i, ti in enumerate(t):
    rinv[i] = re.rinv(ti, T, S, Q, tp, sc)

# Plot
fig, ax = plt.subplots(1,1)
ax.plot(t, rinv)
ax.set_xlabel('t')
ax.set_ylabel(r'$r_{inv}$')
ax.set_xlim(tp, np.max(t))
ax.set_ylim(0., 1.1*rinvmax)
