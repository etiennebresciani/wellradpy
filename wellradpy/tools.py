# -*- coding: utf-8 -*-
"""
Created on 2019/06/26

@author: Etienne Bresciani
"""

import numpy as np # version 1.16.2
import scipy.special as spe # version 1.2.1
import scipy.optimize as opt # version 1.2.1

def E1(u):
    return spe.expn(1, u)

def func_root_E1(u, x):
    return E1(u) - x

def E1inv(x):
    # Note: Another method (e.g. Newton) could be more efficient, but bisection is simple and robust, and we expect that efficiency will not be an issue in practice
    sol = opt.root_scalar(func_root_E1, args=x, method='bisect', bracket=(1e-12, 1e2), rtol=1e-5)
    return sol.root

def X(u):
    return np.exp(-u) - u*E1(u)

def func_root_X(u, x):
    return X(u) - x

def Xinv(x):
    # Note: Another method (e.g. Newton) could be more efficient, but bisection is simple and robust, and we expect that efficiency will not be an issue in practice
    sol = opt.root_scalar(func_root_X, args=x, method='bisect', bracket=(1e-12, 1e2), rtol=1e-5)
    return sol.root

def Y(u, u_w):
    return (E1(u_w)+E1(u)) * np.log((E1(u_w)+E1(u))/E1(u_w))

def func_root_Y(u, x, u_w):
    return Y(u, u_w) - x

def Yinv(x, u_w):
    # Note: Another method (e.g. Newton) could be more efficient, but bisection is simple and robust, and we expect that efficiency will not be an issue in practice
    sol = opt.root_scalar(func_root_Y, args=(x, u_w), method='bisect', bracket=(1e-12, 1e2), rtol=1e-5)
    return sol.root

def Z(u, u_w):
    return (np.exp(-u_w)+np.exp(-u)) * np.log((np.exp(-u_w)+np.exp(-u))/np.exp(-u_w))

def func_root_Z(u, x, u_w):
    return Z(u,u_w) - x

def Zinv(x, u_w):
    # Note: Another method (e.g. Newton) could be more efficient, but bisection is simple and robust, and we expect that efficiency will not be an issue in practice
    sol = opt.root_scalar(func_root_Z, args=(x, u_w), method='bisect', bracket=(1e-12, 1e2), rtol=1e-5)
    return sol.root
