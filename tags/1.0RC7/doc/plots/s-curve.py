# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 10:09:44 2012

@author: schelle
"""

import matplotlib.pyplot as plt
import numpy as np
   
def sCurve(X,a=0.0,b=1.0,c=1.0):

    s = 1.0/(b + np.exp(-c * (X-a)))
    return s

dem50 = 300
dem90 = 313
dem10 = 275
perc = 0.9

C = -np.log(1.0/(perc) - 1)/(dem90 - dem50)
CC = -np.log(1.0/(1-perc) - 1)/(dem10 - dem50)
CCC = (C + CC ) * 0.5
zz = np.array(np.arange(250,380))
S = sCurve(zz,a=dem50,c=C)
SS = sCurve(zz,a=dem50,c=CC)
SSS = sCurve(zz,a=dem50,c=CCC)
plt.plot(zz,S,label='fitted to 90%')
plt.plot(zz,SS,label='fitted to 10%')
plt.plot(zz,SSS,label='fitted to average')
plt.plot(dem90,0.9,'p')
plt.plot(dem10,0.1,'p')
plt.legend()
plt.show()
