# -*- coding: utf-8 -*-
"""
Created on 2019/10/29

@author: Etienne Bresciani
"""

import numpy as np # version 1.16.2
import sys
sys.path.insert(0, '..\examples')
from drawdown_ex import * # Execute the script and retrieve all the variables

np.savez('drawdown_ex_results',
         rinfl_absdraw=rinfl_absdraw,
         rinfl_reldraw=rinfl_reldraw,
         rinfl_relflow=rinfl_relflow,
         rinfl_relvol=rinfl_relvol,
         rinfl_quasisteady=rinfl_quasisteady,
         rinfl_jones=rinfl_jones,
         rinfl_closedres=rinfl_closedres,
         rinfl_impulse=rinfl_impulse,
         rinfl_log=rinfl_log,
         rinv_absdrawdiff=rinv_absdrawdiff,
         rinv_absdrawderivdiff=rinv_absdrawderivdiff,
         rinv_reldrawdiff=rinv_reldrawdiff,
         rinv_reldrawderivdiff=rinv_reldrawderivdiff,
         rinv_reldrawave=rinv_reldrawave,
         rinv_reldrawderivave=rinv_reldrawderivave,
         rinv_propbarrierregime_lin=rinv_propbarrierregime_lin,
         rinv_propbarrierregime_log=rinv_propbarrierregime_log,
         rinv_consthead=rinv_consthead,
         rinv_closedres=rinv_closedres,
         rinv_linearbarr=rinv_linearbarr,
         rinv_impulse=rinv_impulse)
