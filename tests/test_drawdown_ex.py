# -*- coding: utf-8 -*-
"""
Created on 2019/10/29

@author: Etienne Bresciani
"""

import numpy as np # version 1.16.2
import sys
sys.path.insert(0, '..\examples')
from drawdown_ex import * # Execute the script and retrieve all the variables

# Load references results
reference_res = np.load('drawdown_ex_results.npz')

# Calculate the difference
diff = np.empty(0)
diff = np.append(diff, rinfl_absdraw - reference_res['rinfl_absdraw'])
diff = np.append(diff, rinfl_reldraw - reference_res['rinfl_reldraw'])
diff = np.append(diff, rinfl_relflow - reference_res['rinfl_relflow'])
diff = np.append(diff, rinfl_relvol - reference_res['rinfl_relvol'])
diff = np.append(diff, rinfl_quasisteady - reference_res['rinfl_quasisteady'])
diff = np.append(diff, rinfl_jones - reference_res['rinfl_jones'])
diff = np.append(diff, rinfl_closedres - reference_res['rinfl_closedres'])
diff = np.append(diff, rinfl_impulse - reference_res['rinfl_impulse'])
diff = np.append(diff, rinfl_log - reference_res['rinfl_log'])
diff = np.append(diff, rinv_absdrawdiff - reference_res['rinv_absdrawdiff'])
diff = np.append(diff, rinv_absdrawderivdiff - reference_res['rinv_absdrawderivdiff'])
diff = np.append(diff, rinv_reldrawdiff - reference_res['rinv_reldrawdiff'])
diff = np.append(diff, rinv_reldrawderivdiff - reference_res['rinv_reldrawderivdiff'])
#diff = np.append(diff, rinv_reldrawave - reference_res['rinv_reldrawave'])
#diff = np.append(diff, rinv_reldrawderivave - reference_res['rinv_reldrawderivave'])
diff = np.append(diff, rinv_propbarrierregime_lin - reference_res['rinv_propbarrierregime_lin'])
diff = np.append(diff, rinv_propbarrierregime_log - reference_res['rinv_propbarrierregime_log'])
diff = np.append(diff, rinv_consthead - reference_res['rinv_consthead'])
diff = np.append(diff, rinv_closedres - reference_res['rinv_closedres'])
diff = np.append(diff, rinv_linearbarr - reference_res['rinv_linearbarr'])
diff = np.append(diff, rinv_impulse - reference_res['rinv_impulse'])

error_indices = np.nonzero(diff)[0]
if error_indices.size==0:
    print('Test passed successfully')
else:
    print('!!!!!!!!!!!!! TEST FAILED !!!!!!!!!!!!!')
