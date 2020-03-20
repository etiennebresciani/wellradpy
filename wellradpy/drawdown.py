# -*- coding: utf-8 -*-
"""
Created on 2019/06/26

@author: Etienne Bresciani
"""

import numpy as np # version 1.16.2
from .utils import E1, E1inv, whittaker
import scipy.optimize as opt # version 1.2.1
import scipy.integrate as integrate # version 1.2.1

###############################################################################
# Radius of influence functions
###############################################################################

def rinfl_absdraw(t, T, S, Q, sc=0.05):
    """
    Calculate radius of influence during drawdown based on an absolute drawdown
    criterion.

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
    Calculate radius of influence during drawdown based on a relative drawdown
    criterion.

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
    Calculate radius of influence during drawdown based on a relative flow rate
    criterion.

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

def _F(u):
    """
    F function defined by F(u) = exp(-u) - u*E1(u).

    Parameters
    ----------
    u: float
        Any positive real number.

    Returns
    -------
    F(u).

    """
    return np.exp(-u) - u*E1(u)

def _func_root_F(u, x):
    return _F(u) - x

def _Finv(x):
    """
    Inverse F function.

    Parameters
    ----------
    x: float
        Any positive real number.

    Returns
    -------
    Finv(x).

    """
    # Note: Another method (e.g. Newton) could be more efficient, but bisection
    # is simple and robust, and we expect that efficiency will not be an issue
    # in practice
    sol = opt.root_scalar(_func_root_F, args=(x), method='bisect',
                          bracket=(1e-12, 1e2), rtol=1e-5)
    return sol.root

def rinfl_relvol(t, T, S, alpha=0.01):
    """
    Calculate radius of influence during drawdown based on a relative volume
    criterion.

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
    C = 2 * np.sqrt(_Finv(alpha))
    return C * np.sqrt(T*t/S)

def rinfl_quasisteady(t, T, S):
    """
    Calculate radius of influence during drawdown based on quasi-steady state
    model.

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
    Calculate radius of influence during drawdown based on extension of closed
    reservoir regime.

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
    Calculate radius of influence during drawdown based on impulse response
    peak.

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
    Calculate radius of influence during drawdown based on extension of
    logarithmic regime.

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
    Calculate radius of investigation during drawdown based on an absolute
    drawdown difference criterion.

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
    Calculate radius of investigation during drawdown based on an absolute
    drawdown derivative difference criterion.

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
        Window size used to calculate derivative (note: restrictions apply;
        see manuscript for details).
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
    Calculate radius of investigation during drawdown based on a relative
    drawdown difference criterion.

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
    Calculate radius of investigation during drawdown based on a relative
    drawdown derivative difference criterion.

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

def _w_aux(u):
    res = np.exp(-2*u) * whittaker(4*u)
    return float(res)

def _w(u):
    """
    Weighting function that measures the contribution of transmissivity
    variations to drawdown.

    Parameters
    ----------
    u: float
        Any positive real number.

    Returns
    -------
    w(u).

    """
    [aux_integral, aux_integral_err] = integrate.quad(_w_aux, u, np.inf,
                                                      epsrel=1e-5)
    return np.sqrt(np.pi) / u * aux_integral

def _G(u, uw):
    """
    G function defined by G(u, uw) = int_u^inf(w)/int_uw^inf(w).

    Parameters
    ----------
    u: float
        Any positive real number.

    Returns
    -------
    G(u, uw).

    """
    # use this value instead of np.inf otherwise the integration results appear
    # to be wrong; this value is safe for a realistic range of values of u
    # (typically u < 5, while w decreases extremely fast beyond that)
    effective_inf = 10
    [cumul_tail, err_tail] = integrate.quad(_w, u, effective_inf, epsrel=1e-5)
    if uw==1e-2:
        cumul_all = 1.520 # pre-calculated to speed up this particular case
    elif uw==1e-16:
        cumul_all = 17.632 # pre-calculated to speed up this particular case
    else:
        [cumul_all, err_all] = integrate.quad(_w, uw, effective_inf,
                                              epsrel=1e-5)
    return cumul_tail/cumul_all

def _func_root_G(u, x, uw):
    return _G(u, uw) - x

def _Ginv(x, uw):
    """
    Inverse G function for a fixed uw.

    Parameters
    ----------
    x: float
        Any positive real number.

    Returns
    -------
    Ginv(x).

    """
    # Set up bounds
    if uw==0:
        b1 = 1e-16
    else:
        b1 = 1.00001 * uw
    b2 = 1e1
    # Note: Another method (e.g. Newton) could be more efficient, but bisection
    # is simple and robust, and we expect that efficiency will not be an issue
    # in practice
    sol = opt.root_scalar(_func_root_G, args=(x, uw), method='bisect',
                          bracket=(b1, b2), rtol=1e-5)
    return sol.root

def rinv_reldrawave(t, T, S, rw, alpha=0.01):
    """
    Calculate radius of investigation during drawdown based on a relative
    drawdown averaging criterion.

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
    C = 2 * np.sqrt(_Ginv(alpha, uw))
    return C * np.sqrt(T*t/S)

def _wprime(u):
    """
    Weighting function that measures the contribution of transmissivity
    variations to drawdown derivative.

    Parameters
    ----------
    u: float
        Any positive real number.

    Returns
    -------
    wprime(u).

    """
    res = np.sqrt(np.pi) * np.exp(-2*u) * whittaker(4*u)
    return float(res)

def _H(u, uw):
    """
    H function defined by H(u, uw) = int_u^inf(wprime)/int_uw^inf(wprime).

    Parameters
    ----------
    u: float
        Any positive real number.

    Returns
    -------
    H(u, uw).

    """
    [cumul_tail, err_tail] = integrate.quad(_wprime, u, np.inf, epsrel=1e-5)
    [cumul_all, err_all] = integrate.quad(_wprime, uw, np.inf, epsrel=1e-5)
    return cumul_tail/cumul_all

def _func_root_H(u, x, uw):
    return _H(u, uw) - x

def _Hinv(x, uw):
    """
    Inverse H function for a fixed uw.

    Parameters
    ----------
    x: float
        Any positive real number.

    Returns
    -------
    Hinv(x).

    """
    # Set up bounds
    if uw==0:
        b1 = 1e-16
    else:
        b1 = 1.00001 * uw
    b2 = 1e1
    # Note: Another method (e.g. Newton) could be more efficient, but bisection
    # is simple and robust, and we expect that efficiency will not be an issue
    # in practice
    sol = opt.root_scalar(_func_root_H, args=(x, uw), method='bisect',
                          bracket=(b1, b2), rtol=1e-5)
    return sol.root

def rinv_reldrawderivave(t, T, S, rw, alpha=0.01):
    """
    Calculate radius of investigation during drawdown based on a relative
    drawdown derivative averaging criterion.

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
    C = 2 * np.sqrt(_Hinv(alpha, uw))
    return C * np.sqrt(T*t/S)

def rinv_propbarrierregime_lin(t, T, S, alpha=0.5):
    """
    Calculate radius of investigation during drawdown based on a proportion of
    linear barrier regime (linear scale analysis).

    Parameters
    ----------
    t: float
        Time from beginning of pumping.
    T: float
        Transmissivity.
    S: float
        Storativity.
    alpha: float, optional
        Confidence level at which the presence of a linear barrier would be
        detected using drawdown derivative.

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
    Calculate radius of investigation during drawdown based on a proportion of
    linear barrier regime (logarithmic scale analysis).

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
    Calculate radius of investigation during drawdown based on semi-empirical
    start of constant-head boundary effect.

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
    Calculate radius of investigation during drawdown based on intersection of
    unbounded and closed boundary regimes.

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
    Calculate radius of investigation during drawdown based on intersection of
    unbounded and linear barrier regimes.

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
    Calculate radius of investigation during drawdown based on impulse response
    difference peak.

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
