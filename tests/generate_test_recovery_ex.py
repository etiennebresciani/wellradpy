# -*- coding: utf-8 -*-
"""
Created on 2019/10/29

@author: Etienne Bresciani
"""

import numpy as np # version 1.16.2
import sys
sys.path.insert(0, '..\examples')
from recovery_ex import * # Execute the script and retrieve all the variables

np.savez('recovery_ex_results',
         rinvmax=rinvmax,
         tmax=tmax,
         tend=tend,
         rinvfirst=rinv[0])
