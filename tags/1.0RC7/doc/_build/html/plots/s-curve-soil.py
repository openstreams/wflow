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

C=-8000
CC=-1
CCC=-0.5
CCCC=-0.3
zz = np.array(np.arange(250,300))
S = sCurve(zz,a=275,c=C)
SS = sCurve(zz,a=275,c=CC)
SSS = sCurve(zz,a=275,c=CCC)
SSSS = sCurve(zz,a=275,c=CCCC)
plt.plot(zz,S,label='c=' + str(C))
plt.plot(zz,SS,label='c=' + str(CC))
plt.plot(zz,SSS,label='c=' + str(CCC))
plt.plot(zz,SSSS,label='c=' + str(CCCC))
plt.title('Wet roots fraction for a rooting depth of 275 mm')
plt.ylabel('Fraction of wet roots')
plt.xlabel('Water table depth below surface (zi in mm)')
plt.legend()
plt.show()
