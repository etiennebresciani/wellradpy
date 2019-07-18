# -*- coding: utf-8 -*-
"""
Created on 2019/06/26

@author: Etienne Bresciani
"""

from wellradpy import drawdown as dr

# Parameters (not all needed, depending on the method; units as you wish, but consistent)
t = 1
T = 10
S = 1e-4
Q = 30
rw = 0.15
sc = 0.05
alpha = 0.01
alpha_confidence_level = 0.5
delta = 0.4 # Restrictions apply (see manuscript)

###############################################################################
# Radius of influence
###############################################################################

# Method based on absolute drawdown criterion
rinfl_absdrawd = dr.rinfl_absdrawd(t, T, S, Q, sc)

# Method based on relative drawdown criterion
rinfl_reldrawd = dr.rinfl_reldrawd(t, T, S, alpha, rw)

# Method based on relative flow rate criterion
rinfl_relflow = dr.rinfl_relflow(t, T, S, alpha)

# Method based on relative volume criterion
rinfl_relvol = dr.rinfl_relvol(t, T, S, alpha)

# Quasi-steady state model
rinfl_quasisteady = dr.rinfl_quasisteady(t, T, S)

# Jones'formula
rinfl_jones = dr.rinfl_jones(t, T, S)

# Extension of closed reservoir regime
rinfl_closedres = dr.rinfl_closedres(t, T, S)

# Impulse response peak
rinfl_impulse = dr.rinfl_impulse(t, T, S)

# Extension of logarithmic regime
rinfl_log = dr.rinfl_log(t, T, S)

# Dimensional analysis
rinfl_dim = dr.rinfl_dim(t, T, S)

###############################################################################
# Radius of investigation
###############################################################################

# Method based on absolute drawdown difference criterion (analysis in normal scale)
rinv_absdrawddiff_normal = dr.rinv_absdrawddiff_normal(t, T, S, Q, sc)

# Method based on absolute drawdown difference criterion (analysis in log scale)
rinv_absdrawddiff_log = dr.rinv_absdrawddiff_log(t, T, S, Q, sc, rw)

# Method based on absolute drawdown derivative difference criterion (analysis in normal scale)
rinv_absdrawdderivdiff_normal = dr.rinv_absdrawdderivdiff_normal(t, T, S, Q, sc, delta)

# Method based on absolute drawdown derivative difference criterion (analysis in log scale)
rinv_absdrawdderivdiff_log = dr.rinv_absdrawdderivdiff_log(t, T, S, Q, sc, delta, rw)

# Method based on relative drawdown derivative difference criterion (analysis in log scale)
rinv_reldrawdderivdiff_log = dr.rinv_reldrawdderivdiff_log(t, T, S, alpha_confidence_level)

# Semi-empirical start of constant-head boundary effect
rinv_consthead = dr.rinv_consthead(t, T, S)

# Intersection of unbounded and closed boundary regimes
rinv_closedres = dr.rinv_closedres(t, T, S)

# Intersection of unbounded and linear barrier regimes
rinv_linearbarr = dr.rinv_linearbarr(t, T, S)

# Impulse response difference peak
rinv_impulse = dr.rinv_impulse(t, T, S)
