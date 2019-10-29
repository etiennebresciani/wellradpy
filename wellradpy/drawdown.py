# -*- coding: utf-8 -*-
"""
Created on 2019/06/26

@author: Etienne Bresciani
"""

import numpy as np # version 1.16.2
from .tools import E1, E1inv, Finv, Ginv, Hinv

###############################################################################
# Radius of influence functions
###############################################################################

def rinfl_absdraw(t, T, S, Q, sc=0.05):
    """
    Calculate radius of influence during drawdown based on an absolute drawdown criterion.

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
    Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    sc_star = 4*np.pi*T*sc/Q
    C = 2 * np.sqrt(E1inv(sc_star))
    return C * np.sqrt(T*t/S)

def rinfl_reldraw(t, T, S, rw, alpha=0.01):
    """
    Calculate radius of influence during drawdown based on a relative drawdown criterion.

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
    Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    uw = S*rw**2/(4*T*t)
    C = 2 * np.sqrt(E1inv(alpha*E1(uw)))
    return C * np.sqrt(T*t/S)

def rinfl_relflow(t, T, S, alpha=0.01):
    """
    Calculate radius of influence during drawdown based on a relative flow rate criterion.

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
    Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2 * np.sqrt(-np.log(alpha))
    return C * np.sqrt(T*t/S)

def rinfl_relvol(t, T, S, alpha=0.01):
    """
    Calculate radius of influence during drawdown based on a relative volume criterion.

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
    Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2 * np.sqrt(Finv(alpha))
    return C * np.sqrt(T*t/S)

def rinfl_quasisteady(t, T, S):
    """
    Calculate radius of influence during drawdown based on quasi-steady state model.
    
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
    Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2
    return C * np.sqrt(T*t/S)

def rinfl_jones(t, T, S):
    """
    Calculate radius of influence during drawdown based on Jones'formula.
    
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
    Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 4
    return C * np.sqrt(T*t/S)

def rinfl_closedres(t, T, S):
    """
    Calculate radius of influence during drawdown based on extension of closed reservoir regime.
    
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
    Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2.83
    return C * np.sqrt(T*t/S)

def rinfl_impulse(t, T, S):
    """
    Calculate radius of influence during drawdown based on impulse response peak.
    
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
    Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2
    return C * np.sqrt(T*t/S)

def rinfl_log(t, T, S):
    """
    Calculate radius of influence during drawdown based on extension of logarithmic regime.
    
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
    Radius of influence.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 1.5
    return C * np.sqrt(T*t/S)

###############################################################################
# Radius of investigation functions
###############################################################################

def rinv_absdrawdiff(t, T, S, Q, sc=0.05):
    """
    Calculate radius of investigation during drawdown based on an absolute drawdown difference criterion.
    
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
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    sc_star = 4*np.pi*T*sc/Q
    C = np.sqrt(E1inv(sc_star))
    return C * np.sqrt(T*t/S)

def rinv_absdrawderivdiff(t, T, S, Q, delta, sc=0.05):
    """
    Calculate radius of investigation during drawdown based on an absolute drawdown derivative difference criterion.
    
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
        Window size used to calculate derivative (note: restrictions apply; see manuscript for details).
    sc: float, optional
        Absolute drawdown derivative difference threshold.
    
    Returns
    -------
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    sc_star = 4*np.pi*T*sc/Q
    C = np.sqrt(-np.log(np.sqrt(2)*sc_star/(delta)))
    return C * np.sqrt(T*t/S)

def rinv_reldrawdiff(t, T, S, rw, alpha=0.01):
    """
    Calculate radius of investigation during drawdown based on a relative drawdown difference criterion.
    
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
        Relative drawdown difference threshold.
    
    Returns
    -------
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    uw = S*rw**2/(4*T*t)
    C = np.sqrt(E1inv(alpha*E1(uw)))
    return C * np.sqrt(T*t/S)

def rinv_reldrawderivdiff(t, T, S, rw, alpha=0.01):
    """
    Calculate radius of investigation during drawdown based on a relative drawdown derivative difference criterion.
    
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
        Relative drawdown derivative difference threshold.
    
    Returns
    -------
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = np.sqrt(S*rw**2/(4*T*t) - np.log(alpha))
    return C * np.sqrt(T*t/S)

def rinv_reldrawave(t, T, S, rw, alpha=0.01):
    """
    Calculate radius of investigation during drawdown based on a relative drawdown averaging criterion.
    
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
        Relative drawdown averaging threshold.
    
    Returns
    -------
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    uw = S*rw**2/(4*T*t)
    C = 2 * np.sqrt(Ginv(alpha, uw))
    return C * np.sqrt(T*t/S)

def rinv_reldrawderivave(t, T, S, rw, alpha=0.01):
    """
    Calculate radius of investigation during drawdown based on a relative drawdown derivative averaging criterion.
    
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
        Relative drawdown derivative averaging threshold.
    
    Returns
    -------
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    uw = S*rw**2/(4*T*t)
    C = 2 * np.sqrt(Hinv(alpha, uw))
    return C * np.sqrt(T*t/S)

def rinv_propbarrierregime_lin(t, T, S, alpha=0.5):
    """
    Calculate radius of investigation during drawdown based on a proportion of linear barrier regime (linear scale analysis).
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    alpha: float, optional
        Confidence level at which the presence of a linear barrier would be detected using drawdown derivative.
    
    Returns
    -------
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = np.sqrt(-np.log(alpha))
    return C * np.sqrt(T*t/S)

def rinv_propbarrierregime_log(t, T, S, alpha=0.5):
    """
    Calculate radius of investigation during drawdown based on a proportion of linear barrier regime (logarithmic scale analysis).
    
    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    alpha: float, optional
        Confidence level at which the presence of a linear barrier would be detected using drawdown derivative.
    
    Returns
    -------
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = np.sqrt(-np.log(np.power(2, alpha)-1))
    return C * np.sqrt(T*t/S)

def rinv_consthead(t, T, S):
    """
    Calculate radius of investigation during drawdown based on semi-empirical start of constant-head boundary effect.
    
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
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2.64
    return C * np.sqrt(T*t/S)

def rinv_closedres(t, T, S):
    """
    Calculate radius of investigation during drawdown based on intersection of unbounded and closed boundary regimes.
    
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
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 2
    return C * np.sqrt(T*t/S)

def rinv_linearbarr(t, T, S):
    """
    Calculate radius of investigation during drawdown based on intersection of unbounded and linear barrier regimes.
    
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
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 0.75
    return C * np.sqrt(T*t/S)

def rinv_impulse(t, T, S):
    """
    Calculate radius of investigation during drawdown based on impulse response difference peak.
    
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
    Radius of investigation.
    
    Notes
    -----
    Units as you wish, but must be consistent for all the parameters.
    
    """
    
    C = 1
    return C * np.sqrt(T*t/S)
