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
    Wrapper of exponential integral function.

    Parameters
    ----------
    u: float
        Any positive real number.
    
    Returns
    -------
    Exponential integral of u.
    
    Notes
    -----
    The purpose of this function is only to improve code readability.
    
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

def X(u):
    """
    X function defined by X(u) = exp(-u) - u*E1(u).

    Parameters
    ----------
    u: float
        Any positive real number.
    
    Returns
    -------
    X(u).
    
    """
    
    return np.exp(-u) - u*E1(u)

def func_root_X(u, x):
    return X(u) - x

def Xinv(x):
    """
    Inverse X function.

    Parameters
    ----------
    x: float
        Any positive real number.
    
    Returns
    -------
    Xinv(x).
    
    """
    
    # Note: Another method (e.g. Newton) could be more efficient, but bisection is simple and robust, and we expect that efficiency will not be an issue in practice
    sol = opt.root_scalar(func_root_X, args=x, method='bisect', bracket=(1e-12, 1e2), rtol=1e-5)
    return sol.root

def Y(u, uw):
    """
    Y function defined by Y(u,uw) = (E1(uw)+E1(u)) * log((E1(uw)+E1(u))/E1(uw)).

    Parameters
    ----------
    u: float
        Any positive real number.
    uw: float
        Any positive real number.
    
    Returns
    -------
    Y(u,uw).
    
    """
    
    return (E1(uw)+E1(u)) * np.log((E1(uw)+E1(u))/E1(uw))

def func_root_Y(u, x, uw):
    return Y(u, uw) - x

def Yinv(x, uw):
    """
    Inverse Y function.

    Parameters
    ----------
    x: float
        Any positive real number.
    
    Returns
    -------
    Yinv(x).
    
    """
    
    # Note: Another method (e.g. Newton) could be more efficient, but bisection is simple and robust, and we expect that efficiency will not be an issue in practice
    sol = opt.root_scalar(func_root_Y, args=(x, uw), method='bisect', bracket=(1e-12, 1e2), rtol=1e-5)
    return sol.root

def Z(u, uw):
    """
    Z function defined by Z(u,uw) = (exp(-uw)+exp(-u)) * np.log((exp(-uw)+exp(-u))/exp(-uw)).

    Parameters
    ----------
    u: float
        Any positive real number.
    uw: float
        Any positive real number.
    
    Returns
    -------
    Z(u,uw).
    
    """
    
    return (np.exp(-uw)+np.exp(-u)) * np.log((np.exp(-uw)+np.exp(-u))/np.exp(-uw))

def func_root_Z(u, x, uw):
    return Z(u,uw) - x

def Zinv(x, uw):
    """
    Inverse Z function.

    Parameters
    ----------
    x: float
        Any positive real number.
    
    Returns
    -------
    Zinv(x).
    
    """
    
    # Note: Another method (e.g. Newton) could be more efficient, but bisection is simple and robust, and we expect that efficiency will not be an issue in practice
    sol = opt.root_scalar(func_root_Z, args=(x, uw), method='bisect', bracket=(1e-12, 1e2), rtol=1e-5)
    return sol.root
