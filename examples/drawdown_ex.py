# -*- coding: utf-8 -*-
"""
Created on 2019/06/26

@author: Etienne Bresciani
"""

from wellradpy import drawdown as dr

###############################################################################
# Parameters (not all needed, depending on the method)
###############################################################################

t = 1 # Time from beginning of pumping
T = 10 # Transmissivity
S = 1e-4 # Storativity
Q = 30 # Pumping rate
rw = 0.15 # Well radius
sc = 0.05 # Absolute drawdown threshold criterion
alpha = 0.01 # Relative threshold criterion
alpha_conf = 0.5 # Confidence level criterion
delta = 0.4 # Window size in derivative calculation

###############################################################################
# Radius of influence
###############################################################################

# Based on absolute drawdown criterion
rinfl_absdrawd = dr.rinfl_absdrawd(t, T, S, Q, sc)

# Based on relative drawdown criterion
rinfl_reldrawd = dr.rinfl_reldrawd(t, T, S, rw, alpha)

# Based on relative flow rate criterion
rinfl_relflow = dr.rinfl_relflow(t, T, S, alpha)

# Based on relative volume criterion
rinfl_relvol = dr.rinfl_relvol(t, T, S, alpha)

# Based on quasi-steady state model
rinfl_quasisteady = dr.rinfl_quasisteady(t, T, S)

# Based on Jones'formula
rinfl_jones = dr.rinfl_jones(t, T, S)

# Based on extension of closed reservoir regime
rinfl_closedres = dr.rinfl_closedres(t, T, S)

# Based on impulse response peak
rinfl_impulse = dr.rinfl_impulse(t, T, S)

# Based on extension of logarithmic regime
rinfl_log = dr.rinfl_log(t, T, S)

# Based on dimensional analysis
rinfl_dim = dr.rinfl_dim(t, T, S)

###############################################################################
# Radius of investigation
###############################################################################

# Based on absolute drawdown difference criterion (analysis in lin scale)
rinv_absdrawddiff_lin = dr.rinv_absdrawddiff_lin(t, T, S, Q, sc)

# Based on absolute drawdown difference criterion (analysis in log scale)
rinv_absdrawddiff_log = dr.rinv_absdrawddiff_log(t, T, S, Q, rw, sc)

# Based on absolute drawdown derivative difference criterion (analysis in lin scale)
rinv_absdrawdderivdiff_lin = dr.rinv_absdrawdderivdiff_lin(t, T, S, Q, delta, sc)

# Based on absolute drawdown derivative difference criterion (analysis in log scale)
rinv_absdrawdderivdiff_log = dr.rinv_absdrawdderivdiff_log(t, T, S, Q, rw, delta, sc)

# Based on relative drawdown derivative difference criterion (analysis in log scale)
rinv_reldrawdderivdiff_log = dr.rinv_reldrawdderivdiff_log(t, T, S, alpha_conf)

# Based on semi-empirical start of constant-head boundary effect
rinv_consthead = dr.rinv_consthead(t, T, S)

# Based on intersection of unbounded and closed boundary regimes
rinv_closedres = dr.rinv_closedres(t, T, S)

# Based on intersection of unbounded and linear barrier regimes
rinv_linearbarr = dr.rinv_linearbarr(t, T, S)

# Based on impulse response difference peak
rinv_impulse = dr.rinv_impulse(t, T, S)
