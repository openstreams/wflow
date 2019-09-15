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

C=8
CC=4
CCC=3
CCCC=1

swredu = 0.038
b=1/(1-swredu)
zz = np.array(np.arange(-3,3,0.1))
S = sCurve(zz,a=0,b=b,c=C) + swredu
SS = sCurve(zz,a=0,b=b,c=CC) + swredu
SSS = sCurve(zz,a=0,b=b,c=CCC) + swredu
SSSS = sCurve(zz,a=0,b=b,c=CCCC) + swredu
plt.plot(zz,S,label='c=' + str(C))
plt.plot(zz,SS,label='c=' + str(CC))
plt.plot(zz,SSS,label='c=' + str(CCC))
plt.plot(zz,SSSS,label='c=' + str(CCCC))
plt.title('Infiltration reduction for frozen soil')
plt.ylabel('Reduction factor (cf_soil)')
plt.xlabel('Temperature $^\circ$C')
plt.legend()
plt.show()
