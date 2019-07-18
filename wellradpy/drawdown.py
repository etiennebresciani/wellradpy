# -*- coding: utf-8 -*-
"""
Created on 2019/06/26

@author: Etienne Bresciani
"""

import numpy as np # version 1.16.2
from .tools import E1, E1inv, Xinv, Yinv, Zinv

###############################################################################
# Radius of influence functions
###############################################################################

def rinfl_absdrawd(t, T, S, Q, sc):
    """
    Calculate radius of influence based on absolute drawdown criterion.
    """
    C = 2 * np.sqrt(E1inv(4*np.pi*T*sc/Q))
    return C * np.sqrt(T*t/S)

def rinfl_reldrawd(t, T, S, alpha, rw):
    """
    Calculate radius of influence based on relative drawdown criterion.
    """
    C = 2 * np.sqrt(E1inv(alpha*E1(S*rw**2/(4*T*t))))
    return C * np.sqrt(T*t/S)

def rinfl_relflow(t, T, S, alpha):
    """
    Calculate radius of influence based on relative flow rate criterion.
    """
    C = 2 * np.sqrt(-np.log(alpha))
    return C * np.sqrt(T*t/S)

def rinfl_relvol(t, T, S, alpha):
    """
    Calculate radius of influence based on relative volume criterion.
    """
    C = 2 * np.sqrt(Xinv(alpha))
    return C * np.sqrt(T*t/S)

def rinfl_quasisteady(t, T, S):
    """
    Calculate radius of influence based on quasi-steady state model.
    """
    C = 2
    return C * np.sqrt(T*t/S)

def rinfl_jones(t, T, S):
    """
    Calculate radius of influence based on Jones'formula.
    """
    C = 4
    return C * np.sqrt(T*t/S)

def rinfl_closedres(t, T, S):
    """
    Calculate radius of influence based on extension of closed reservoir regime.
    """
    C = 2.83
    return C * np.sqrt(T*t/S)

def rinfl_impulse(t, T, S):
    """
    Calculate radius of influence based on impulse response peak.
    """
    C = 2
    return C * np.sqrt(T*t/S)

def rinfl_log(t, T, S):
    """
    Calculate radius of influence based on extension of logarithmic regime.
    """
    C = 1.5
    return C * np.sqrt(T*t/S)

def rinfl_dim(t, T, S):
    """
    Calculate radius of influence based on dimensional analysis.
    """
    C = 1
    return C * np.sqrt(T*t/S)

###############################################################################
# Radius of investigation functions
###############################################################################

def rinv_absdrawddiff_normal(t, T, S, Q, sc):
    """
    Calculate radius of investigation based on absolute drawdown difference criterion (analysis in normal scale).
    """
    C = np.sqrt(E1inv(4*np.pi*T*sc/Q))
    return C * np.sqrt(T*t/S)

def rinv_absdrawddiff_log(t, T, S, Q, sc, rw):
    """
    Calculate radius of investigation based on absolute drawdown difference criterion (analysis in log scale).
    """
    uw = S * rw**2 / (4*T*t)
    C = np.sqrt(Yinv(4*np.pi*T*sc/Q, uw))
    return C * np.sqrt(T*t/S)

def rinv_absdrawdderivdiff_normal(t, T, S, Q, sc, delta):
    """
    Calculate radius of investigation based on absolute drawdown derivative difference criterion (analysis in normal scale).
    """
    C = np.sqrt(-np.log(4*np.pi*T*np.sqrt(2)*sc/(Q*delta)))
    return C * np.sqrt(T*t/S)

def rinv_absdrawdderivdiff_log(t, T, S, Q, sc, delta, rw):
    """
    Calculate radius of investigation based on absolute drawdown derivative difference criterion (analysis in log scale).
    """
    uw = S * rw**2 / (4*T*t)
    C = np.sqrt(Zinv(4*np.pi*T*np.sqrt(2)*sc/(Q*delta), uw))
    return C * np.sqrt(T*t/S)

def rinv_reldrawdderivdiff_log(t, T, S, alpha):
    """
    Calculate radius of investigation based on relative drawdown derivative difference criterion (analysis in log scale).
    """
    C = np.sqrt(-np.log(np.power(2, alpha)-1))
    return C * np.sqrt(T*t/S)

def rinv_consthead(t, T, S):
    """
    Calculate radius of investigation based on semi-empirical start of constant-head boundary effect.
    """
    C = 2.64
    return C * np.sqrt(T*t/S)

def rinv_closedres(t, T, S):
    """
    Calculate radius of investigation based on intersection of unbounded and closed boundary regimes.
    """
    C = 2
    return C * np.sqrt(T*t/S)

def rinv_linearbarr(t, T, S):
    """
    Calculate radius of investigation based on intersection of unbounded and linear barrier regimes.
    """
    C = 0.75
    return C * np.sqrt(T*t/S)

def rinv_impulse(t, T, S):
    """
    Calculate radius of investigation based on impulse response difference peak.
    """
    C = 1
    return C * np.sqrt(T*t/S)
