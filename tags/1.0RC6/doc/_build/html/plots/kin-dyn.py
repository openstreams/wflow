# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 10:09:44 2012

@author: schelle
"""

import matplotlib.pyplot as plt
import numpy as np
from pcrut import *
  
runkin, d = readtss('run.tss')
rundyn, d = readtss('rundyn.tss')
 
plt.figure()
ax =plt.plot(rundyn[:,6],label="Dynamic")
ax =plt.plot(runkin[:,6],label="Kinematic")
leg = plt.legend() 
plt.title("Difference between kinematic and dynamic wave")
