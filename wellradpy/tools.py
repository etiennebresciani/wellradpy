# -*- coding: utf-8 -*-
"""
Created on 2019/06/26

@author: Etienne Bresciani
"""

import numpy as np # version 1.16.2
import scipy.special as spe # version 1.2.1
import scipy.optimize as opt # version 1.2.1
import scipy.integrate as integrate # version 1.2.1
import mpmath as mp

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

def func_root_E1(u, x):
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
    
    # Note: Another method (e.g. Newton) could be more efficient, but bisection is simple and robust, and we expect that efficiency will not be an issue in practice
    sol = opt.root_scalar(func_root_E1, args=x, method='bisect', bracket=(1e-12, 1e2), rtol=1e-5)
    return sol.root

def F(u):
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

def func_root_F(u, x):
    return F(u) - x

def Finv(x):
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
    
    # Note: Another method (e.g. Newton) could be more efficient, but bisection is simple and robust, and we expect that efficiency will not be an issue in practice
    sol = opt.root_scalar(func_root_F, args=x, method='bisect', bracket=(1e-12, 1e2), rtol=1e-5)
    return sol.root

def w_aux(u):
    res = np.exp(-2*u) * mp.whitw(0.5, 0.5, 4*u)
    return float(res)

def w(u):
    """
    Weighting function that measures the contribution of transmissivity variations to drawdown.

    Parameters
    ----------
    u: float
        Any positive real number.
    
    Returns
    -------
    w(u).
    
    """
    
    [aux_integral, aux_integral_err] = integrate.quad(w_aux, u, np.inf, epsrel=1e-5)
    return np.sqrt(np.pi) / u * aux_integral

def G(u, uw):
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
    
    effective_inf = 10 # use this value instead of np.inf otherwise the integration results appear to be wrong; this value is safe for a realistic range of values of u (typically u < 5, while w decreases extremely fast beyond that)
    [cumul_tail, err_tail] = integrate.quad(w, u, effective_inf, epsrel=1e-5)
    if uw==1e-2:
        cumul_all = 1.520 # pre-calculated to speed up this particular case
    elif uw==1e-16:
        cumul_all = 17.632 # pre-calculated to speed up this particular case
    else:
        [cumul_all, err_all] = integrate.quad(w, uw, effective_inf, epsrel=1e-5)
    return cumul_tail/cumul_all

def func_root_G(u, x, uw):
    return G(u, uw) - x

def Ginv(x, uw):
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
    # Note: Another method (e.g. Newton) could be more efficient, but bisection is simple and robust, and we expect that efficiency will not be an issue in practice
    sol = opt.root_scalar(func_root_G, args=(x, uw), method='bisect', bracket=(b1, b2), rtol=1e-5)
    return sol.root

def wprime(u):
    """
    Weighting function that measures the contribution of transmissivity variations to drawdown derivative.

    Parameters
    ----------
    u: float
        Any positive real number.
    
    Returns
    -------
    wprime(u).
    
    """
    
    res = np.sqrt(np.pi) * np.exp(-2*u) * mp.whitw(0.5,0.5,4*u)
    return float(res)

def H(u, uw):
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
    [cumul_tail, err_tail] = integrate.quad(wprime, u, np.inf, epsrel=1e-5)
    [cumul_all, err_all] = integrate.quad(wprime, uw, np.inf, epsrel=1e-5)
    return cumul_tail/cumul_all

def func_root_H(u, x, uw):
    return H(u, uw) - x

def Hinv(x, uw):
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
    # Note: Another method (e.g. Newton) could be more efficient, but bisection is simple and robust, and we expect that efficiency will not be an issue in practice
    sol = opt.root_scalar(func_root_H, args=(x, uw), method='bisect', bracket=(b1, b2), rtol=1e-5)
    return sol.root
