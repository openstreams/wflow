# -*- coding: utf-8 -*-
"""
Created on Wed Feb 04 14:52:30 2015

@author: teuser
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 03 16:31:35 2014

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
import JarvisCoefficients

def selectSaR(i):
    """
    not all functions are still in this file, the older functions can be found
    (with the same numbering) in h:\My Documents\memo's\python scripts\wflow\
    """
    if i == 1:
        name = 'agriZone_Jarvis'
    elif i == 2:
        name = 'agriZone_Ep'
    elif i == 3:
        name = 'agriZone_Ep_Sa'
    elif i == 4:
        name = 'agriZone_Ep_Sa_cropG'
    return name

    
      
def agriZone_Jarvis(self,k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on Jarvis stress functions
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa
    - Code for ini-file: 1
    """
    self.Qa = max(self.Pe - (self.samax[k] - self.Sa_t[k]),0)    
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa)
    self.SaN = self.Sa[k] / self.samax[k]
    self.SuN = self.Su[k] / self.sumax[k]
    
    JarvisCoefficients.calcEu(self,k,1)           #calculation of Ea based on Jarvis stress functions
    self.Ea1 = self.Eu
    
#    if self.teller == 45:
#        pdb.set_trace()
    self.Fa1 = self.Fmin[k] + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * self.SuN)
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Fa1 - self.Ea1

    self.Sa_diff = ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Fa = self.Fa1 + (self.Fa1/ifthenelse(self.Fa1 + self.Ea1 > 0 , self.Fa1 + self.Ea1 , 1)) * self.Sa_diff
    self.Ea = self.Ea1 + (self.Ea1/ifthenelse(self.Fa1 + self.Ea1 > 0 , self.Fa1 + self.Ea1 , 1)) * self.Sa_diff
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Ea - self.Fa   
    self.Sa[k] = ifthenelse(self.Sa[k] < 0, 0 , self.Sa[k])    
    self.Sa_diff2 = ifthen(self.Sa[k] < 0, self.Sa[k]) 

    self.wbSa_[k] = self.Pe - self.Ea - self.Qa - self.Fa - self.Sa[k] + self.Sa_t[k]
    
    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa
    self.Fa_[k] = self.Fa
    

def agriZone_Ep(self,k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa
    - Code for ini-file: 2
    """
    JarvisCoefficients.calcEp(self,k)
    self.PotEvaporation = self.EpHour    
    
    self.Qa = max(self.Pe - (self.samax[k] - self.Sa_t[k]),0)    
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa)
    self.SaN = self.Sa[k] / self.samax[k]
    self.SuN = self.Su[k] / self.sumax[k]
    
    self.Ea1 = max((self.PotEvaporation - self.Ei),0) * min(self.Sa[k] / (self.samax[k] * self.LP[k]),1)        
    
    self.Fa1 = self.Fmin[k] + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * self.SuN)
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Fa1 - self.Ea1

    self.Sa_diff = ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Fa = self.Fa1 + (self.Fa1/ifthenelse(self.Fa1 + self.Ea1 > 0 , self.Fa1 + self.Ea1 , 1)) * self.Sa_diff
    self.Ea = self.Ea1 + (self.Ea1/ifthenelse(self.Fa1 + self.Ea1 > 0 , self.Fa1 + self.Ea1 , 1)) * self.Sa_diff
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Ea - self.Fa   
    self.Sa[k] = ifthenelse(self.Sa[k] < 0, 0 , self.Sa[k])    
    self.Sa_diff2 = ifthen(self.Sa[k] < 0, self.Sa[k]) 

    self.wbSa_[k] = self.Pe - self.Ea - self.Qa - self.Fa - self.Sa[k] + self.Sa_t[k]
    
    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa
    self.Fa_[k] = self.Fa

def agriZone_Ep_Sa(self,k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa
    - Fa is based on storage in Sa
    - Code for ini-file: 3
    """
    JarvisCoefficients.calcEp(self,k)
    self.PotEvaporation = self.EpHour    
    
    self.Qa = max(self.Pe - (self.samax[k] - self.Sa_t[k]),0)    
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa)
    self.SaN = self.Sa[k] / self.samax[k]
    self.SuN = self.Su[k] / self.sumax[k]
    
    self.Ea1 = max((self.PotEvaporation - self.Ei),0) * min(self.Sa[k] / (self.samax[k] * self.LP[k]),1)        
    
    self.Fa1 = ifthenelse(self.SaN > 0,self.Fmin[k] + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),0)
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Fa1 - self.Ea1

    self.Sa_diff = ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Fa = self.Fa1 + (self.Fa1/ifthenelse(self.Fa1 + self.Ea1 > 0 , self.Fa1 + self.Ea1 , 1)) * self.Sa_diff
    self.Ea = self.Ea1 + (self.Ea1/ifthenelse(self.Fa1 + self.Ea1 > 0 , self.Fa1 + self.Ea1 , 1)) * self.Sa_diff
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Ea - self.Fa   
    self.Sa[k] = ifthenelse(self.Sa[k] < 0, 0 , self.Sa[k])    
    self.Sa_diff2 = ifthen(self.Sa[k] < 0, self.Sa[k]) 

    self.wbSa_[k] = self.Pe - self.Ea - self.Qa - self.Fa - self.Sa[k] + self.Sa_t[k]
    
    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa
    self.Fa_[k] = self.Fa
    
    
def agriZone_Ep_Sa_cropG(self,k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa
    - Fa is based on storage in Sa
    - Code for ini-file: 4
    """
    JarvisCoefficients.calcEp(self,k)
    self.PotEvaporation = self.EpHour    
    
    self.samax2 = self.samax[k] * self.cropG   
    self.Qaadd = max(self.Sa_t[k] - self.samax2,0)
    
#    if self.teller == 66:    
#        pdb.set_trace()    
    
    self.Qa = max(self.Pe - (self.samax2 - self.Sa_t[k]),0) + self.Qaadd
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa)
    self.SaN = self.Sa[k] / self.samax2
    self.SuN = self.Su[k] / self.sumax[k]
    
    self.Ea1 = max((self.PotEvaporation - self.Ei),0) * min(self.Sa[k] / (self.samax2 * self.LP[k]),1)        
    
    self.Fa1 = ifthenelse(self.SaN > 0,self.Fmin[k] + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),0)
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Fa1 - self.Ea1

    self.Sa_diff = ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Fa = self.Fa1 + (self.Fa1/ifthenelse(self.Fa1 + self.Ea1 > 0 , self.Fa1 + self.Ea1 , 1)) * self.Sa_diff
    self.Ea = self.Ea1 + (self.Ea1/ifthenelse(self.Fa1 + self.Ea1 > 0 , self.Fa1 + self.Ea1 , 1)) * self.Sa_diff
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Ea - self.Fa   
    self.Sa[k] = ifthenelse(self.Sa[k] < 0, 0 , self.Sa[k])    
    self.Sa_diff2 = ifthen(self.Sa[k] < 0, self.Sa[k]) 

    self.wbSa_[k] = self.Pe - self.Ea - self.Qa - self.Fa - self.Sa[k] + self.Sa_t[k]
    
    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa
    self.Fa_[k] = self.Fa
