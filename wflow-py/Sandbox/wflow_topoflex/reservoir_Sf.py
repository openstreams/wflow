# -*- coding: utf-8 -*-
"""
Created on Fri Apr 04 17:15:38 2014

@author: TEuser

List all function versions
"""

import numpy
import pdb
from copy import copy as copylist

try:
    from  wflow.wf_DynamicFramework import *
except ImportError:
    from  wf_DynamicFramework import *
import scipy

def selectSfR(i):
    """
    not all functions are still in this file, the older functions can be found
    (with the same numbering) in h:\My Documents\memo's\python scripts\wflow\
    """
    if i == 1:
        name = 'fastRunoff_lag'
    if i == 2:
        name = 'fastRunoff_lag2'
    return name

def fastRunoff_no_reservoir(self, k):
    """
    This function is used when no unsaturated zone reservoir is used and only
    passes fluxes from the upper reservoirs to the lower
    Qf = 0.
    Storage in fast reservoir = 0.
    """
    self.Qfin_[k] = self.Qu * (1 - self.D)
    self.Qf_[k] = self.Qfin_[k]
    self.Sf_[k] = 0.
    self.wbSf_[k] = self.Qfin_[k] - self.Qf_[k] - self.Sf[k] + self.Sf_t[k]

def fastRunoff_lag2(self, k):
    """
    - Lag is applied before inflow into the fast reservoir 
    - Lag formula is derived from Fenicia (2011)
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato
    - not a semi analytical solution for Sf anymore
    - Code for ini-file: 2
    """

    if self.FR_L:
        self.Qu = areatotal(self.Qu_[k] * self.percentArea ,nominal(self.TopoId)) 
    else:
        self.Qu = self.Qu_[k]
    
    self.Qfin = (1 - self.D[k]) * self.Qu
      
    if self.D[k] < 1.00: 
        if self.convQu[k]:
            self.QfinLag = self.convQu[k][-1]
            self.Qf = self.Sf[k] * self.Kf[k]                
            self.Sf[k] = self.Sf[k] + self.QfinLag - self.Qf
                
            self.convQu[k].insert(0, 0 * scalar(self.TopoId))            #convolution Qu for following time steps
            self.Tfmap = self.Tf[k] * scalar(self.TopoId)
            del self.convQu[k][-1]
            temp = [self.convQu[k][i] + (2/self.Tfmap-2/(self.Tfmap*(self.Tfmap+1))*(self.Tfmap-i))*self.Qfin for i in range(len(self.convQu[k]))]
            self.convQu[k] = temp
            
#            self.Qfin_[k] = self.Qfin
#            self.Qfinput_[k] = self.QfinLag
                
        else:
            self.Qf = self.Sf[k] * self.Kf[k]                
            self.Sf[k] = self.Sf[k] + self.Qfin - self.Qf
            
#            self.Qfin_[k] = self.Qfin
#            self.Qfinput_[k] = self.Qfin
    else:
        self.Qf = self.ZeroMap
        self.Qfinput_[k] = self.ZeroMap
        self.Qfin_[k] = self.ZeroMap

    self.wbSf_[k] = self.Qfin - self.Qf - self.Sf[k] + self.Sf_t[k] - sum(self.convQu[k]) + sum(self.convQu_t[k])    
    
    self.Qf_[k] = self.Qf
#    self.QuA_[k] = self.Qu

    
def routingQf_combined(self):       
    """
    - Routing of fluxes from fast reservoir 
    - Qf is devided over the reservoir numbers for the timesteps matching with the average
    travel time for a calcultation cell
    """
    if sum(pcr2numpy(self.Transit,NaN)) > 0:    
        self.Qflag = self.trackQ[0]    
        self.trackQ.append(0 * scalar(self.TopoId))    
        del self.trackQ[0]
        temp = [self.trackQ[i] + ifthenelse(rounddown(self.Transit) == i*scalar(self.TopoId), (self.Transit - i*scalar(self.TopoId)) * self.Qfcub, ifthenelse(roundup(self.Transit) == i*scalar(self.TopoId), (i*scalar(self.TopoId) - self.Transit) * self.Qfcub, 0)) for i in range(len(self.trackQ))]
        self.trackQ = temp    
    
    else:
        self.Qflag = self.Qfcub
    
    self.wbSfrout = self.Qfcub - self.Qflag - sum(self.trackQ) + sum (self.trackQ_t)
    
    self.Qflag_ = self.Qflag