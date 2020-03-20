# -*- coding: utf-8 -*-
"""
Created on 2019/06/26

@author: Etienne Bresciani
"""

import numpy as np # version 1.16.2
import scipy.special as spe # version 1.2.1
import scipy.optimize as opt # version 1.2.1

def E1(u):
    """
    Exponential integral function.

    Parameters
    ----------
    u: float
        Any positive real number.

    Returns
    -------
    Exponential integral of u.

    Notes
    -----
    This is simply a wrapper to improve code readability.

    """
    return spe.expn(1, u)

def _func_root_E1(u, x):
    return E1(u) - x

def E1inv(x):
    """
    Inverse exponential integral function.

    Parameters
    ----------
    x: float
        Any positive real number.

    Returns
    -------
    Inverse exponential integral of x.

    """
    # Note: Another method (e.g. Newton) could be more efficient, but bisection
    # is simple and robust, and we expect that efficiency will not be an issue
    # in practice
    sol = opt.root_scalar(_func_root_E1, args=(x), method='bisect',
                          bracket=(1e-12, 1e2), rtol=1e-5)
    return sol.root

def E1inv_appr(x):
    """
    Analytical approximation of inverse exponential integral function.

    Parameters
    ----------
    x: float
        Any positive real number.

    Returns
    -------
    Inverse exponential integral of x.

    """
    return 3.656*np.power(x, -0.1295) - 3.445

def whittaker(z):
    """
    Whittaker function for kapa=1/2 and mu=1/2.

    Parameters
    ----------
    z: float
        Any positive real number.
    """
    return np.exp(-0.5*z) * z * spe.hyperu(0.5, 2, z)
