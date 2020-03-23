# -*- coding: utf-8 -*-
"""
Created on 2020/03/20

@author: Etienne Bresciani
"""

import numpy as np # version 1.16.2
from .utils import E1, E1inv
import scipy.optimize as opt # version 1.2.1

def _barrier_effect_star(rinv_star, t_star):
    return E1(rinv_star**2/t_star) - E1(rinv_star**2/(t_star-1))

def _func_root_rinv_star(rinv_star_unknown, sc_star_target, t_star):
    return _barrier_effect_star(rinv_star_unknown, t_star) - sc_star_target

def _rinv_star(sc_star, t_star):
    """
    Calculate the dimensionless radius of investigation during recovery.

    Parameters
    ----------
    sc_star: float
        Dimensionless apparent resolution.
    t_star: float
        Dimensionless time from beginning of pumping.

    Returns
    -------
    rinv_star

    """
    # Note: Another method (e.g. Newton) could be more efficient, but bisection
    # is simple and robust, and we expect that efficiency will not be an issue
    # in practice
    sol = opt.root_scalar(_func_root_rinv_star, args=(sc_star,t_star),
                          method='bisect', bracket=(1e-12, 1e3), rtol=1e-5)
    return sol.root

def rinv(t, T, S, Q, tp, sc=0.05):
    """
    Calculate the radius of investigation during recovery.

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
    tp: float
        Pumping duration.
    sc: float, optional
        Apparent resolution.

    Returns
    -------
    rinv

    """
    sc_star = 4*np.pi*T*sc/Q
    t_star = t/tp
    rinv_star = _rinv_star(sc_star, t_star)
    rp = np.sqrt(T*tp/S)
    rinv = rinv_star*rp
    return rinv

def _barrier_effect_at_tmax_star(tmax_star):
    return E1((tmax_star-1)*np.log(tmax_star/(tmax_star-1))) - \
           E1(tmax_star*np.log(tmax_star/(tmax_star-1)))

def _func_root_tmax_star(tmax_star_unknown, sc_star_target):
    return _barrier_effect_at_tmax_star(tmax_star_unknown) - sc_star_target

def _tmax_star(sc_star):
    """
    Calculate the dimensionless time at which radius of investigation is
    maximum during recovery.

    Parameters
    ----------
    sc_star: float
        Dimensionless apparent resolution.

    Returns
    -------
    tmax_star

    """
    # Note: Another method (e.g. Newton) could be more efficient, but bisection
    # is simple and robust, and we expect that efficiency will not be an issue
    # in practice
    sol = opt.root_scalar(_func_root_tmax_star, args=(sc_star),
                          method='bisect', bracket=(1.0000000001,1e5), rtol=1e-5)
    return sol.root

def tmax(T, Q, tp, sc=0.05):
    """
    Calculate the time at which the radius of investigation is maximum during
    recovery.

    Parameters
    ----------
    T: float
        Transmissivity.
    Q: float
        Pumping rate.
    tp: float
        Pumping duration.
    sc: float, optional
        Apparent resolution.

    Returns
    -------
    tmax

    """
    sc_star = 4*np.pi*T*sc/Q
    tmax_star = _tmax_star(sc_star)
    tmax = tmax_star*tp
    return tmax

def _rinvmax_star(sc_star):
    """
    Calculate the dimensionless maximum radius of investigation during
    recovery.

    Parameters
    ----------
    sc: float, optional
        Dimensionless apparent resolution.

    Returns
    -------
    tmax

    """
    tmax_star = _tmax_star(sc_star)
    rinvmax_star = np.sqrt(tmax_star * (tmax_star-1) * np.log(tmax_star / \
                           (tmax_star-1)))
    return rinvmax_star

def rinvmax(T, S, Q, tp, sc=0.05):
    """
    Calculate the maximum radius of investigation during recovery.

    Parameters
    ----------
    T: float
        Transmissivity.
    S: float
        Storativity.
    Q: float
        Pumping rate.
    tp: float
        Pumping duration.
    sc: float, optional
        Apparent resolution.

    Returns
    -------
    tmax

    """
    sc_star = 4*np.pi*T*sc/Q
    rinvmax_star = _rinvmax_star(sc_star)
    rp = np.sqrt(T*tp/S)
    rinvmax = rinvmax_star*rp
    return rinvmax

def _barrier_effect_at_tend_star(tend_star):
    return E1(1.e-10/tend_star) - E1(1.e-10/(tend_star-1))

def _func_root_tend_star(tend_star_unknown, sc_star_target):
    return _barrier_effect_at_tend_star(tend_star_unknown) - sc_star_target

def _tend_star(sc_star):
    """
    Calculate the dimensionless time at which the radius of investigation
    becomes zero during recovery (i.e., when the recovery test may be
    considered terminated).

    Parameters
    ----------
    sc_star: float
        Dimensionless apparent resolution.

    Returns
    -------
    tend_star

    """
    # Note: Another method (e.g. Newton) could be more efficient, but bisection
    # is simple and robust, and we expect that efficiency will not be an issue
    # in practice
    sol = opt.root_scalar(_func_root_tend_star, args=(sc_star),
                          method='bisect', bracket=(1.00000001,1e5), rtol=1e-5)
    return sol.root

def tend(T, Q, tp, sc=0.05):
    """
    Calculate the time at which the radius of investigation becomes zero during
    recovery (i.e., when the recovery test may be considered terminated).

    Parameters
    ----------
    T: float
        Transmissivity.
    Q: float
        Pumping rate.
    tp: float
        Pumping duration.
    sc: float, optional
        Apparent resolution.

    Returns
    -------
    tend

    """
    sc_star = 4*np.pi*T*sc/Q
    tend_star = _tend_star(sc_star)
    tend = tend_star*tp
    return tend
