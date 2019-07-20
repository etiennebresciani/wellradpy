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

def rinfl_absdrawd(t, T, S, Q, sc=0.05):
    """
    Calculate radius of influence during pumping based on an absolute drawdown criterion.

    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    Q: float
        Pumping rate.
    sc: float, optional
        Absolute drawdown threshold.
    
    Returns
    -------
    rinfl: float
        Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2 * np.sqrt(E1inv(4*np.pi*T*sc/Q))
    return C * np.sqrt(T*t/S)

def rinfl_reldrawd(t, T, S, rw, alpha=0.01):
    """
    Calculate radius of influence during pumping based on an relative drawdown criterion.

    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    rw: float
        Well radius.
    alpha: float, optional
        Drawdown threshold relative to drawdown at the well.
    
    Returns
    -------
    rinfl: float
        Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2 * np.sqrt(E1inv(alpha*E1(S*rw**2/(4*T*t))))
    return C * np.sqrt(T*t/S)

def rinfl_relflow(t, T, S, alpha=0.01):
    """
    Calculate radius of influence during pumping based on an relative flow rate criterion.

    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    alpha: float, optional
        Flow rate threshold relative to pumping rate.
    
    Returns
    -------
    rinfl: float
        Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2 * np.sqrt(-np.log(alpha))
    return C * np.sqrt(T*t/S)

def rinfl_relvol(t, T, S, alpha=0.01):
    """
    Calculate radius of influence during pumping based on an relative volume criterion.

    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    alpha: float, optional
        Volume threshold relative to volume of cone of depression.
    
    Returns
    -------
    rinfl: float
        Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2 * np.sqrt(Xinv(alpha))
    return C * np.sqrt(T*t/S)

def rinfl_quasisteady(t, T, S):
    """
    Calculate radius of influence based on quasi-steady state model.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    
    Returns
    -------
    rinfl: float
        Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2
    return C * np.sqrt(T*t/S)

def rinfl_jones(t, T, S):
    """
    Calculate radius of influence based on Jones'formula.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    
    Returns
    -------
    rinfl: float
        Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 4
    return C * np.sqrt(T*t/S)

def rinfl_closedres(t, T, S):
    """
    Calculate radius of influence based on extension of closed reservoir regime.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    
    Returns
    -------
    rinfl: float
        Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2.83
    return C * np.sqrt(T*t/S)

def rinfl_impulse(t, T, S):
    """
    Calculate radius of influence based on impulse response peak.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    
    Returns
    -------
    rinfl: float
        Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2
    return C * np.sqrt(T*t/S)

def rinfl_log(t, T, S):
    """
    Calculate radius of influence based on extension of logarithmic regime.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    
    Returns
    -------
    rinfl: float
        Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 1.5
    return C * np.sqrt(T*t/S)

def rinfl_dim(t, T, S):
    """
    Calculate radius of influence based on dimensional analysis.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    
    Returns
    -------
    rinfl: float
        Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 1
    return C * np.sqrt(T*t/S)

###############################################################################
# Radius of investigation functions
###############################################################################

def rinv_absdrawddiff_lin(t, T, S, Q, sc=0.05):
    """
    Calculate radius of investigation based on absolute drawdown difference criterion when data are analyzed in lin scale.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    Q: float
        Pumping rate.
    sc: float, optional
        Absolute drawdown difference threshold.
    
    Returns
    -------
    rinfl: float
        Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = np.sqrt(E1inv(4*np.pi*T*sc/Q))
    return C * np.sqrt(T*t/S)

def rinv_absdrawddiff_log(t, T, S, Q, rw, sc=0.05):
    """
    Calculate radius of investigation based on absolute drawdown difference criterion when data are analyzed in log scale.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    Q: float
        Pumping rate.
    sc: float, optional
        Absolute drawdown difference threshold.
    
    Returns
    -------
    rinfl: float
        Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    uw = S * rw**2 / (4*T*t)
    C = np.sqrt(Yinv(4*np.pi*T*sc/Q, uw))
    return C * np.sqrt(T*t/S)

def rinv_absdrawdderivdiff_lin(t, T, S, Q, delta, sc=0.05):
    """
    Calculate radius of investigation based on absolute drawdown derivative difference criterion when data are analyzed in lin scale.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    Q: float
        Pumping rate.
    delta: float
        Window size in derivative calculation (note: restrictions apply; see manuscript for details).
    sc: float, optional
        Absolute drawdown difference threshold.
    
    Returns
    -------
    rinfl: float
        Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = np.sqrt(-np.log(4*np.pi*T*np.sqrt(2)*sc/(Q*delta)))
    return C * np.sqrt(T*t/S)

def rinv_absdrawdderivdiff_log(t, T, S, Q, rw, delta, sc=0.05):
    """
    Calculate radius of investigation based on absolute drawdown derivative difference criterion when data are analyzed in log scale.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    Q: float
        Pumping rate.
    delta: float
        Window size in derivative calculation (note: restrictions apply; see manuscript for details).
    sc: float, optional
        Absolute drawdown difference threshold.
    
    Returns
    -------
    rinfl: float
        Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    uw = S * rw**2 / (4*T*t)
    C = np.sqrt(Zinv(4*np.pi*T*np.sqrt(2)*sc/(Q*delta), uw))
    return C * np.sqrt(T*t/S)

def rinv_reldrawdderivdiff_log(t, T, S, alpha=0.5):
    """
    Calculate radius of investigation based on relative drawdown derivative difference criterion when data are analyzed in log scale.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    alpha: float, optional
        Confidence level at which the presence of a linear barrier would be detected.
    
    Returns
    -------
    rinfl: float
        Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = np.sqrt(-np.log(np.power(2, alpha)-1))
    return C * np.sqrt(T*t/S)

def rinv_consthead(t, T, S):
    """
    Calculate radius of investigation based on semi-empirical start of constant-head boundary effect.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    
    Returns
    -------
    rinfl: float
        Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2.64
    return C * np.sqrt(T*t/S)

def rinv_closedres(t, T, S):
    """
    Calculate radius of investigation based on intersection of unbounded and closed boundary regimes.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    
    Returns
    -------
    rinfl: float
        Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2
    return C * np.sqrt(T*t/S)

def rinv_linearbarr(t, T, S):
    """
    Calculate radius of investigation based on intersection of unbounded and linear barrier regimes.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    
    Returns
    -------
    rinfl: float
        Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 0.75
    return C * np.sqrt(T*t/S)

def rinv_impulse(t, T, S):
    """
    Calculate radius of investigation based on impulse response difference peak.
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    
    Returns
    -------
    rinfl: float
        Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 1
    return C * np.sqrt(T*t/S)
