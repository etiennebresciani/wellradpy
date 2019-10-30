# -*- coding: utf-8 -*-
"""
Created on 2019/06/26

@author: Etienne Bresciani
"""

from wellradpy import drawdown as dr

###############################################################################
# Parameters (not all always needed, depending on the definition)
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

# Based on an absolute drawdown criterion
rinfl_absdraw = dr.rinfl_absdraw(t, T, S, Q, sc)

# Based on a relative drawdown criterion
rinfl_reldraw = dr.rinfl_reldraw(t, T, S, rw, alpha)

# Based on a relative flow rate criterion
rinfl_relflow = dr.rinfl_relflow(t, T, S, alpha)

# Based on a relative volume criterion
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

###############################################################################
# Radius of investigation
###############################################################################

# Based on an absolute drawdown difference criterion
rinv_absdrawdiff = dr.rinv_absdrawdiff(t, T, S, Q, sc)

# Based on an absolute drawdown derivative difference criterion
rinv_absdrawderivdiff = dr.rinv_absdrawderivdiff(t, T, S, Q, delta, sc)

# Based on a relative drawdown difference criterion
rinv_reldrawdiff = dr.rinv_reldrawdiff(t, T, S, rw, alpha)

# Based on a relative drawdown derivative difference criterion
rinv_reldrawderivdiff = dr.rinv_reldrawderivdiff(t, T, S, rw, alpha)

# Based on a relative drawdown averaging criterion (RUNTIME ~HOURS)
#rinv_reldrawave = dr.rinv_reldrawave(t, T, S, rw, alpha)

# Based on a relative drawdown derivative averaging criterion (RUNTIME ~SECONDS)
rinv_reldrawderivave = dr.rinv_reldrawderivave(t, T, S, rw, alpha)

# Based on a proportion of linear barrier regime (linear scale analysis)
rinv_propbarrierregime_lin = dr.rinv_propbarrierregime_lin(t, T, S, alpha_conf)

# Based on a proportion of linear barrier regime (logarithmic scale analysis)
rinv_propbarrierregime_log = dr.rinv_propbarrierregime_log(t, T, S, alpha_conf)

# Based on semi-empirical start of constant-head boundary effect
rinv_consthead = dr.rinv_consthead(t, T, S)

# Based on intersection of unbounded and closed boundary regimes
rinv_closedres = dr.rinv_closedres(t, T, S)

# Based on intersection of unbounded and linear barrier regimes
rinv_linearbarr = dr.rinv_linearbarr(t, T, S)

# Based on impulse response difference peak
rinv_impulse = dr.rinv_impulse(t, T, S)
