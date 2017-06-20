# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 16:10:57 2013

@author: ir. Frederiek Sperna Weiland
"""

import  shutil

src ="prcp_giss_aom_a1b_1952-02-28_03min.map"

startYear = src[18:22]

print startYear

for i in range (0,52,4):

    dst = src[0:18] + str(int(startYear)+i) + src[22:26] + str(29) + src[28:38]

    print dst

    shutil.copy(src,dst)