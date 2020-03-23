# -*- coding: utf-8 -*-
"""
Created on 2019/10/29

@author: Etienne Bresciani
"""

import numpy as np # version 1.16.2
import sys
sys.path.insert(0, '..\examples')
from recovery_ex import * # Execute the script and retrieve all the variables

# Load references results
reference_res = np.load('recovery_ex_results.npz')

# Calculate the difference
diff = np.empty(0)
diff = np.append(diff, rinvmax - reference_res['rinvmax'])
diff = np.append(diff, tmax - reference_res['tmax'])
diff = np.append(diff, tend - reference_res['tend'])
diff = np.append(diff, rinv[0] - reference_res['rinvfirst'])

error_indices = np.nonzero(diff)[0]
if error_indices.size==0:
    print('Test passed successfully')
else:
    print('!!!!!!!!!!!!! TEST FAILED !!!!!!!!!!!!!')
