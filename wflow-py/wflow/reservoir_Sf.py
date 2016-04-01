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

def fastRunoff_no_reservoir(self, k):
    """
    This function is used when no unsaturated zone reservoir is used and only
    passes fluxes from the upper reservoirs to the lower
    Qf = Qf_in.
    Storage in fast reservoir = 0.
    """
    self.Qfin_[k] = self.Qu * (1 - self.D)
    self.Qf_[k] = self.Qfin_[k]
    self.Sf_[k] = 0.
    self.wbSf_[k] = self.Qfin_[k] - self.Qf_[k] - self.Sf[k] + self.Sf_t[k]
    
def fastAgriRunoff_no_reservoir(self, k):
    """
    This function is used when no unsaturated zone reservoir is used and only
    passes fluxes from the upper reservoirs to the lower
    Qfa = Qfa_in.
    Storage in fast reservoir = 0.
    """

    self.Qfa_[k] = self.Qa_[k]
    self.Sfa[k] = 0.
    self.wbSfa_[k] = self.Qfa_[k] - self.Qfa_[k] - self.Sfa[k] + self.Sfa_t[k] - sum(self.convQa[k]) + sum(self.convQa_t[k])

def fastRunoff_lag2(self, k):
    """
    - Lag is applied before inflow into the fast reservoir 
    - Lag formula is derived from Fenicia (2011)
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato
    - not a semi analytical solution for Sf anymore
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
            
        else:
            self.Qf = self.Sf[k] * self.Kf[k]                
            self.Sf[k] = self.Sf[k] + self.Qfin - self.Qf
            
    else:
        self.Qf = self.ZeroMap
        self.Qfinput_[k] = self.ZeroMap
        self.Qfin_[k] = self.ZeroMap

    self.wbSf_[k] = self.Qfin - self.Qf - self.Sf[k] + self.Sf_t[k] - sum(self.convQu[k]) + sum(self.convQu_t[k])    
    
    self.Qf_[k] = self.Qf


def fastRunoff_lag_forAgri_combined(self,k):
    """
    - Lag is applied before inflow into the fast reservoir 
    - Lag formula is derived from Fenicia (2011)
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato
    - not a semi analytical solution for Sf anymore
    - When separate ditch/road fast reservoir is taken into account, this option should not be used
    """

    if self.FR_L:
        self.Qu = areatotal(self.Qu_[k] * self.percentArea ,nominal(self.TopoId)) 
        self.Qa = areatotal(self.Qa_[k] * self.percentArea ,nominal(self.TopoId)) 
    else:
        self.Qu = self.Qu_[k]
        self.Qa = self.Qa_[k]
    
    self.Qfin = (1 - self.D[k]) * self.Qu + self.Qa

     
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
            
                
        else:
            self.Qf = self.Sf[k] * self.Kf[k]                
            self.Sf[k] = self.Sf[k] + self.Qfin - self.Qf
            

    else:                              
        self.Qf = self.ZeroMap
        self.Qfinput_[k] = self.ZeroMap
        self.Qfin_[k] = self.ZeroMap

    self.wbSf_[k] = self.Qfin - self.Qf - self.Sf[k] + self.Sf_t[k] - sum(self.convQu[k]) + sum(self.convQu_t[k])    
    
    self.Qf_[k] = self.Qf
    
def fastRunoff_lag_agriDitch(self,k):
    """
    - Lag is applied before inflow into the fast reservoir 
    - Lag formula is derived from Fenicia (2011)
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato
    - not a semi analytical solution for Sf anymore
    - very fast responding reservoir to represent fast drainage via roads and ditches
    """

    if self.FR_L:
        self.Qa = areatotal(self.Qa_[k] * self.percentArea ,nominal(self.TopoId)) 
    else:
        self.Qa = self.Qa_[k]
    
    self.Qfain = self.Qa

    if self.convQa[k]:
        self.QfainLag = self.convQa[k][-1]
        self.Qfa = self.Sfa[k] * self.Kfa[k]                
        self.Sfa[k] = self.Sfa[k] + self.QfainLag - self.Qfa
            
        self.convQa[k].insert(0, 0 * scalar(self.TopoId))            #convolution Qu for following time steps
        self.Tfmap = self.Tfa[k] * scalar(self.TopoId)
        del self.convQa[k][-1]
        temp = [self.convQa[k][i] + (2/self.Tfmap-2/(self.Tfmap*(self.Tfmap+1))*(self.Tfmap-i))*self.Qfain for i in range(len(self.convQa[k]))]
        self.convQa[k] = temp

    else:
        self.Qfa = self.Sfa[k] * self.Kfa[k]                
        self.Sfa[k] = self.Sfa[k] + self.Qfain - self.Qfa

    self.wbSfa_[k] = self.Qfain - self.Qfa - self.Sfa[k] + self.Sfa_t[k] - sum(self.convQa[k]) + sum(self.convQa_t[k])    
    
    self.Qfa_[k] = self.Qfa
    
def routingQf_combined(self):       
    """
    - Routing of fluxes from fast reservoir 
    - Qf is devided over the reservoir numbers for the timesteps matching with the average
    travel time for a calcultation cell
    """
           
    if nansum(pcr2numpy(self.Transit,NaN)) > 0:
        self.Qflag = self.trackQ[0]  # first bucket is transferred to outlet
        self.trackQ.append(0 * scalar(self.TopoId))  # add new bucket for present time step
        del self.trackQ[0]  # remove first bucket (transferred to outlet)
        temp = [self.trackQ[i] + ifthenelse(rounddown(self.Transit) == i*scalar(self.TopoId),
                                            (self.Transit - i*scalar(self.TopoId)) * self.Qftotal / 1000 * self.surfaceArea,
                                            ifthenelse(roundup(self.Transit) == i*scalar(self.TopoId),
                                                       (i*scalar(self.TopoId) - self.Transit) * self.Qftotal / 1000 * self.surfaceArea, 0)) for i in range(len(self.trackQ))]
        self.trackQ = temp    
    
    else:
        self.Qflag = self.Qftotal
    
    self.wbSfrout = self.Qftotal - self.Qflag - sum(self.trackQ) + sum (self.trackQ_t)
    
    self.Qflag_ = self.Qflag
    
    self.Qtlag = self.Qflag_ / self.timestepsecs + self.Qs_ / 1000 * self.surfaceArea / self.timestepsecs
    self.QLagTot = areatotal(self.Qtlag, nominal(self.TopoId))  # catchment total runoff with looptijd

   
def routingQf_Qs_grid(self):
    """
    - Routing of both Qf and Qs
    - based on a velocity map
    """
    self.Qtot = self.Qftotal + self.Qs_  # total local discharge in mm/hour
    self.Qtotal = self.Qtot / 1000 * self.surfaceArea / self.timestepsecs  # total local discharge in m3/s
    self.Qstate_t = self.Qstate
    self.Qrout = accutraveltimeflux(self.TopoLdd, self.Qstate + self.Qtotal, self.velocity)
    self.Qstate = accutraveltimestate(self.TopoLdd, self.Qstate + self.Qtotal, self.velocity)
    
    self.Qtlag = self.Qrout
    self.QLagTot = self.Qrout
    
    # water balance of flux routing
    self.dSdt = self.Qstate-self.Qstate_t
    self.WB_rout = (accuflux(self.TopoLdd, self.Qtotal - self.dSdt)-self.Qrout)/accuflux(self.TopoLdd, self.Qtotal)