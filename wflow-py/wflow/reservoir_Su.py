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

def unsatZone_no_reservoir(self, k):
    """
    This function is used when no unsaturated zone reservoir is used and only
    passes fluxes from the upper reservoirs to the lower
    Qu = Pe
    Eu = 0.
    Perc = 0.
    Cap = 0.
    Storage in unsaturated zone = 0.
    """
    self.Qu_[k] = max(self.Pe_[k], 0)
    self.Eu_[k] = 0.
    self.Perc_[k] = 0.
    self.Su[k] = 0.
    self.Cap_[k] = 0.
    self.wbSu_[k] = self.Pe - self.Eu - self.Qu - self.Perc + self.Cap - self.Su[k] + self.Su_t[k]

def unsatZone_LP_beta(self,k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation linear until LP, from than with potential rate
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato
    """
    
    self.Su[k] = ifthenelse(self.Su_t[k] + self.Pe > self.sumax[k], self.sumax[k], self.Su_t[k] + self.Pe) 
    self.Quadd = ifthenelse(self.Su_t[k] + self.Pe > self.sumax[k], self.Su_t[k] + self.Pe - self.sumax[k], 0)
    self.SuN = self.Su[k] / self.sumax[k]
    self.SiN = self.Si[k] / self.imax[k]    
    
    self.Eu1 = max((self.PotEvaporation - self.Ei),0) * min(self.Su[k] / (self.sumax[k] * self.LP[k]),1)
        
    self.Qu1 = (self.Pe - self.Quadd) * (1 - (1 - self.SuN) ** self.beta[k])
    self.Perc1 = self.perc[k] * self.SuN
    self.Su[k] = self.Su_t[k] + (self.Pe - self.Quadd) - self.Qu1 - self.Eu1 - self.Perc1

    self.Su_diff = ifthenelse(self.Su[k] < 0, self.Su[k], 0)
    self.Eu = self.Eu1 + (self.Eu1 / ifthenelse(self.Qu1 + self.Eu1 +self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Qu = self.Qu1 + (self.Qu1/ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Perc = ifthenelse (self.Perc1 > 0, self.Perc1 + (self.Perc1/ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff, self.Perc1)
    self.Su[k] = self.Su_t[k] + (self.Pe - self.Quadd) - self.Eu - self.Qu - self.Perc   
    self.Su[k] = ifthenelse(self.Su[k] < 0, 0 , self.Su[k])    
    self.Su_diff2 = ifthen(self.Su[k] < 0, self.Su[k]) 

    self.Cap = min(self.cap[k] * (1 - self.Su[k] / self.sumax[k]), self.Ss)    
    self.Su[k] = self.Su[k] + self.Cap    
    
    self.wbSu_[k] = self.Pe - self.Eu - self.Qu - self.Quadd - self.Perc + self. Cap - self.Su[k] + self.Su_t[k]
    
    self.Eu_[k] = self.Eu
    self.Qu_[k] = self.Qu + self.Quadd
    self.Cap_[k] = self.Cap
    self.Perc_[k] = self.Perc
   
      
def unsatZone_LP_beta_Jarvis(self,k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on Jarvis stress functions
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qu is determined with a beta function (same as in HBV?)
    """
    self.Su[k] = ifthenelse(self.Su_t[k] + self.Pe > self.sumax[k], self.sumax[k], self.Su_t[k] + self.Pe) 
    self.Quadd = ifthenelse(self.Su_t[k] + self.Pe > self.sumax[k], self.Su_t[k] + self.Pe - self.sumax[k], 0)
    self.SuN = self.Su[k] / self.sumax[k]
    self.SiN = self.Si[k] / self.imax[k]    
    
    JarvisCoefficients.calcEu(self,k,1)           #calculation of Eu based on Jarvis stress functions
            
    self.Qu1 = (self.Pe - self.Quadd) * (1 - (1 - self.SuN) ** self.beta[k])
    self.Perc1 = self.perc[k] * self.SuN
    self.Su[k] = self.Su_t[k] + (self.Pe - self.Quadd) - self.Qu1 - self.Eu - self.Perc1

    self.Su_diff = ifthenelse(self.Su[k] < 0, self.Su[k], 0)
    self.Qu = self.Qu1 + (self.Qu1/ifthenelse(self.Qu1 + self.Perc1 > 0 , self.Qu1 + self.Perc1 , 1)) * self.Su_diff
    self.Perc = ifthenelse (self.Perc1 > 0, self.Perc1 + (self.Perc1/ifthenelse(self.Qu1 + self.Perc1 > 0 , self.Qu1 + self.Perc1 , 1)) * self.Su_diff, self.Perc1)
    self.Su[k] = self.Su_t[k] + (self.Pe - self.Quadd) - self.Eu - self.Qu - self.Perc   
    self.Su[k] = ifthenelse(self.Su[k] < 0, 0 , self.Su[k])    
    self.Su_diff2 = ifthen(self.Su[k] < 0, self.Su[k]) 

    self.Cap = min(self.cap[k] * (1 - self.Su[k] / self.sumax[k]), self.Ss)    
    self.Su[k] = self.Su[k] + self.Cap    
    
    self.wbSu_[k] = self.Pe - self.Eu - self.Qu - self.Quadd - self.Perc + self. Cap - self.Su[k] + self.Su_t[k]
    
    self.Eu_[k] = self.Eu
    self.Qu_[k] = self.Qu + self.Quadd
    self.Cap_[k] = self.Cap
    self.Perc_[k] = self.Perc
    
    
def unsatZone_LP_beta_Ep(self,k):
    """
    - Potential evaporation is calculated with formula in 'JarvisCoefficients', but without
    using the Jarvis stress functions
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation linear until LP, from than with potential rate
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato
    - Qu is determined with a beta function (same as in HBV?)
    """

    JarvisCoefficients.calcEp(self,k)
    self.PotEvaporation = cover(ifthenelse(self.EpHour >= 0, self.EpHour, 0),0)  
    
    self.Su[k] = ifthenelse(self.Su_t[k] + self.Pe > self.sumax, self.sumax, self.Su_t[k] + self.Pe) 
    self.Quadd = ifthenelse(self.Su_t[k] + self.Pe > self.sumax, self.Su_t[k] + self.Pe - self.sumax, 0)
    self.SuN = self.Su[k] / self.sumax
    self.SiN = self.Si[k] / self.imax[k]    
    
    self.Eu1 = max((self.PotEvaporation - self.Ei),0) * min(self.Su[k] / (self.sumax * self.LP[k]),1)
    
    self.Qu1 = (self.Pe - self.Quadd) * (1 - (1 - self.SuN) ** self.beta[k])
    self.Perc1 = self.perc[k] * self.SuN
    self.Su[k] = self.Su_t[k] + (self.Pe - self.Quadd) - self.Qu1 - self.Eu1 - self.Perc1

    self.Su_diff = ifthenelse(self.Su[k] < 0, self.Su[k], 0)
    self.Eu = self.Eu1 + (self.Eu1 / ifthenelse(self.Qu1 + self.Eu1 +self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Qu = self.Qu1 + (self.Qu1/ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Perc = ifthenelse (self.Perc1 > 0, self.Perc1 + (self.Perc1/ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff, self.Perc1)
    self.Su[k] = self.Su_t[k] + (self.Pe - self.Quadd) - self.Eu - self.Qu - self.Perc   
    self.Su[k] = ifthenelse(self.Su[k] < 0, 0 , self.Su[k])    
    self.Su_diff2 = ifthen(self.Su[k] < 0, self.Su[k]) 

    self.Cap = min(self.cap[k] * (1 - self.Su[k] / self.sumax), self.Ss)    
    self.Su[k] = self.Su[k] + self.Cap    
    
    self.wbSu_[k] = self.Pe - self.Eu - self.Qu - self.Quadd - self.Perc + self. Cap - self.Su[k] + self.Su_t[k]
    
    self.Eu_[k] = self.Eu
    self.Qu_[k] = self.Qu + self.Quadd
    self.Cap_[k] = self.Cap
    self.Perc_[k] = self.Perc

def unsatZone_LP_beta_Ep_cropG(self,k):
    """
    - Potential evaporation is calculated with formula in 'JarvisCoefficients', but without
    using the Jarvis stress functions
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation linear until LP, from than with potential rate
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato
    - Qu is determined with a beta function (same as in HBV?)
    - root zone storage for crop land is decreased in autumn and winter
    """
    
    JarvisCoefficients.calcEp(self,k)
    self.PotEvaporation = cover(ifthenelse(self.EpHour >= 0, self.EpHour, 0),0)  

    self.cropG_scal = pcr2numpy(self.cropG,NaN)
    if any(self.cropG_scal == 1):
        self.sumax2 = self.sumax[k]
    elif any(self.cropG_scal > 0):      
        self.sumax2 = self.sumax[k] * (1 - numpy.max(self.cropG_scal[self.cropG_scal >= 0]) * (1-self.redsu[k]))
    else:
        self.sumax2 = self.sumax[k] * self.redsu[k]
        
    self.Su[k] = ifthenelse(self.Su_t[k] + self.Pe > self.sumax2, self.sumax2, self.Su_t[k] + self.Pe) 
    self.Quadd = ifthenelse(self.Su_t[k] + self.Pe > self.sumax2, self.Su_t[k] + self.Pe - self.sumax2, 0)
    self.SuN = self.Su[k] / self.sumax2
    self.SiN = self.Si[k] / self.imax[k]    
    
    self.Eu1 = max((self.PotEvaporation - self.Ei),0) * min(self.Su[k] / (self.sumax2 * self.LP[k]),1)
        
    self.Qu1 = (self.Pe - self.Quadd) * (1 - (1 - self.SuN) ** self.beta[k])
    self.Perc1 = self.perc[k] * self.SuN
    self.Su[k] = self.Su_t[k] + (self.Pe - self.Quadd) - self.Qu1 - self.Eu1 - self.Perc1

    self.Su_diff = ifthenelse(self.Su[k] < 0, self.Su[k], 0)
    self.Eu = self.Eu1 + (self.Eu1 / ifthenelse(self.Qu1 + self.Eu1 +self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Qu = self.Qu1 + (self.Qu1/ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Perc = ifthenelse (self.Perc1 > 0, self.Perc1 + (self.Perc1/ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff, self.Perc1)
    self.Su[k] = self.Su_t[k] + (self.Pe - self.Quadd) - self.Eu - self.Qu - self.Perc   
    self.Su[k] = ifthenelse(self.Su[k] < 0, 0 , self.Su[k])    
    self.Su_diff2 = ifthen(self.Su[k] < 0, self.Su[k]) 

    self.Cap = min(self.cap[k] * (1 - self.Su[k] / self.sumax2), self.Ss)    
    self.Su[k] = self.Su[k] + self.Cap    
    
    self.wbSu_[k] = self.Pe - self.Eu - self.Qu - self.Quadd - self.Perc + self. Cap - self.Su[k] + self.Su_t[k]
    
    self.Eu_[k] = self.Eu
    self.Qu_[k] = self.Qu + self.Quadd
    self.Cap_[k] = self.Cap
    self.Perc_[k] = self.Perc

def unsatZone_forAgri_Jarvis(self,k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on Jarvis stress functions
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qu is determined with a beta function (same as in HBV?)
    - inflow is infiltration from agriculture reservoir
    """
    self.Su[k] = ifthenelse(self.Su_t[k] + self.Fa > self.sumax[k], self.sumax[k], self.Su_t[k] + self.Fa) 
    self.Quadd = ifthenelse(self.Su_t[k] + self.Fa > self.sumax[k], self.Su_t[k] + self.Fa - self.sumax[k], 0)
    self.SuN = self.Su[k] / self.sumax[k]
    self.SiN = self.Si[k] / self.imax[k]    
    
    JarvisCoefficients.calcEu(self,k,2)           #calculation of Eu based on Jarvis stress functions
    self.Eu1 = self.Eu
            
    self.Qu1 = (self.Fa - self.Quadd) * (1 - (1 - self.SuN) ** self.beta[k])
    self.Perc1 = self.perc[k] * self.SuN
    self.Su[k] = self.Su_t[k] + (self.Fa - self.Quadd) - self.Qu1 - self.Eu - self.Perc1

    self.Su_diff = ifthenelse(self.Su[k] < 0, self.Su[k], 0)
    self.Eu = self.Eu1 + (self.Eu1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Qu = self.Qu1 + (self.Qu1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Perc = ifthenelse (self.Perc1 > 0, self.Perc1 + (self.Perc1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff, self.Perc1)
    self.Su[k] = self.Su_t[k] + (self.Fa - self.Quadd) - self.Eu - self.Qu - self.Perc   
    self.Su[k] = ifthenelse(self.Su[k] < 0, 0 , self.Su[k])    
    self.Su_diff2 = ifthen(self.Su[k] < 0, self.Su[k]) 

    self.Cap = min(self.cap[k] * (1 - self.Su[k] / self.sumax[k]), self.Ss)    
    self.Su[k] = self.Su[k] + self.Cap    
    
    self.wbSu_[k] = self.Fa - self.Eu - self.Qu - self.Quadd - self.Perc + self. Cap - self.Su[k] + self.Su_t[k]
    
    self.Eu_[k] = self.Eu
    self.Qu_[k] = self.Qu + self.Quadd
    self.Cap_[k] = self.Cap
    self.Perc_[k] = self.Perc

def unsatZone_forAgri_Ep(self,k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on beta/LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qu is determined with a beta function (same as in HBV?)
    - inflow is infiltration from agriculture reservoir
    """
    
    JarvisCoefficients.calcEp(self,k)
    self.PotEvaporation = cover(ifthenelse(self.EpHour >= 0, self.EpHour, 0),0)  

    self.Su[k] = ifthenelse(self.Su_t[k] + self.Fa > self.sumax[k], self.sumax[k], self.Su_t[k] + self.Fa) 
    self.Quadd = ifthenelse(self.Su_t[k] + self.Fa > self.sumax[k], self.Su_t[k] + self.Fa - self.sumax[k], 0)
    self.SuN = self.Su[k] / self.sumax[k]
    self.SiN = self.Si[k] / self.imax[k]    
    
    self.Eu1 = ifthenelse(self.Ft_[k] == 1, max((self.PotEvaporation - self.Ei - self.Ea),0) * min(self.Su[k] / (self.sumax * self.LP[k]),1), 0)         # no transpiration in case of frozen soil. Added on 22 feb 2016        
    
    self.Qu1 = (self.Fa - self.Quadd) * (1 - (1 - self.SuN) ** self.beta[k])
    self.Perc1 = self.perc[k] * self.SuN
    self.Su[k] = self.Su_t[k] + (self.Fa - self.Quadd) - self.Qu1 - self.Eu - self.Perc1

    self.Su_diff = ifthenelse(self.Su[k] < 0, self.Su[k], 0)
    self.Eu = self.Eu1 + (self.Eu1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Qu = self.Qu1 + (self.Qu1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Perc = ifthenelse (self.Perc1 > 0, self.Perc1 + (self.Perc1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff, self.Perc1)
    self.Su[k] = self.Su_t[k] + (self.Fa - self.Quadd) - self.Eu - self.Qu - self.Perc   
    self.Su[k] = ifthenelse(self.Su[k] < 0, 0 , self.Su[k])    
    self.Su_diff2 = ifthen(self.Su[k] < 0, self.Su[k]) 

    self.Cap = min(self.cap[k] * (1 - self.Su[k] / self.sumax[k]), self.Ss)    
    self.Su[k] = self.Su[k] + self.Cap    
    
    self.wbSu_[k] = self.Fa - self.Eu - self.Qu - self.Quadd - self.Perc + self. Cap - self.Su[k] + self.Su_t[k]
    
    self.Eu_[k] = self.Eu
    self.Qu_[k] = self.Qu + self.Quadd
    self.Cap_[k] = self.Cap
    self.Perc_[k] = self.Perc

def unsatZone_forAgri_hourlyEp(self,k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on beta/LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qu is determined with a beta function (same as in HBV?)
    - inflow is infiltration from agriculture reservoir
    """
    
    self.Su[k] = ifthenelse(self.Su_t[k] + self.Fa > self.sumax[k], self.sumax[k], self.Su_t[k] + self.Fa) 
    self.Quadd = ifthenelse(self.Su_t[k] + self.Fa > self.sumax[k], self.Su_t[k] + self.Fa - self.sumax[k], 0)
    self.SuN = self.Su[k] / self.sumax[k]
    self.SiN = self.Si[k] / self.imax[k]    
    
    self.Eu1 = ifthenelse(self.Ft_[k] == 1, max((self.PotEvaporation - self.Ei - self.Ea),0) * min(self.Su[k] / (self.sumax[k] * self.LP[k]),1), 0)         # no transpiration in case of frozen soil. Added on 31 mrt 2016        
    
    self.Qu1 = (self.Fa - self.Quadd) * (1 - (1 - self.SuN) ** self.beta[k])
    self.Perc1 = self.perc[k] * self.SuN
    self.Su[k] = self.Su_t[k] + (self.Fa - self.Quadd) - self.Qu1 - self.Eu - self.Perc1

    self.Su_diff = ifthenelse(self.Su[k] < 0, self.Su[k], 0)
    self.Eu = self.Eu1 + (self.Eu1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Qu = self.Qu1 + (self.Qu1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Perc = ifthenelse (self.Perc1 > 0, self.Perc1 + (self.Perc1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff, self.Perc1)
    self.Su[k] = self.Su_t[k] + (self.Fa - self.Quadd) - self.Eu - self.Qu - self.Perc   
    self.Su[k] = ifthenelse(self.Su[k] < 0, 0 , self.Su[k])    
    self.Su_diff2 = ifthen(self.Su[k] < 0, self.Su[k]) 

    self.Cap = min(self.cap[k] * (1 - self.Su[k] / self.sumax[k]), self.Ss)    
    self.Su[k] = self.Su[k] + self.Cap    
    
    self.wbSu_[k] = self.Fa - self.Eu - self.Qu - self.Quadd - self.Perc + self. Cap - self.Su[k] + self.Su_t[k]
    
    self.Eu_[k] = self.Eu
    self.Qu_[k] = self.Qu + self.Quadd
    self.Cap_[k] = self.Cap
    self.Perc_[k] = self.Perc


def unsatZone_forAgri_Jarvis_cropG(self,k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on Jarvis stress functions
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qu is determined with a beta function (same as in HBV?)
    - inflow is infiltration from agriculture reservoir
    """
    self.cropG_scal = pcr2numpy(self.cropG,NaN)
    if any(self.cropG_scal == 1):
        self.sumax2 = self.sumax[k]
    else:
        self.sumax2 = self.sumax[k] * self.redsu[k]
    
    self.Su[k] = ifthenelse(self.Su_t[k] + self.Fa > self.sumax2, self.sumax2, self.Su_t[k] + self.Fa) 
    self.Quadd = ifthenelse(self.Su_t[k] + self.Fa > self.sumax2, self.Su_t[k] + self.Fa - self.sumax2, 0)
    self.SuN = self.Su[k] / self.sumax2
    self.SiN = self.Si[k] / self.imax[k]    
    
    JarvisCoefficients.calcEu(self,k,2)           #calculation of Eu based on Jarvis stress functions
    self.Eu1 = self.Eu
            
    self.Qu1 = (self.Fa - self.Quadd) * (1 - (1 - self.SuN) ** self.beta[k])
    self.Perc1 = self.perc[k] * self.SuN
    self.Su[k] = self.Su_t[k] + (self.Fa - self.Quadd) - self.Qu1 - self.Eu - self.Perc1

    self.Su_diff = ifthenelse(self.Su[k] < 0, self.Su[k], 0)
    self.Eu = self.Eu1 + (self.Eu1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Qu = self.Qu1 + (self.Qu1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Perc = ifthenelse (self.Perc1 > 0, self.Perc1 + (self.Perc1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff, self.Perc1)
    self.Su[k] = self.Su_t[k] + (self.Fa - self.Quadd) - self.Eu - self.Qu - self.Perc   
    self.Su[k] = ifthenelse(self.Su[k] < 0, 0 , self.Su[k])    
    self.Su_diff2 = ifthen(self.Su[k] < 0, self.Su[k]) 

    self.Cap = min(self.cap[k] * (1 - self.Su[k] / self.sumax2), self.Ss)    
    self.Su[k] = self.Su[k] + self.Cap    
    
    self.wbSu_[k] = self.Fa - self.Eu - self.Qu - self.Quadd - self.Perc + self. Cap - self.Su[k] + self.Su_t[k]
    
    self.Eu_[k] = self.Eu
    self.Qu_[k] = self.Qu + self.Quadd
    self.Cap_[k] = self.Cap
    self.Perc_[k] = self.Perc

def unsatZone_forAgri_Ep_cropG(self,k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on beta/LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qu is determined with a beta function (same as in HBV?)
    - inflow is infiltration from agriculture reservoir
    """
    
    JarvisCoefficients.calcEp(self,k)
    self.PotEvaporation = cover(ifthenelse(self.EpHour >= 0, self.EpHour, 0),0)  
    
    self.cropG_scal = pcr2numpy(self.cropG,NaN)
    if any(self.cropG_scal == 1):
        self.sumax2 = self.sumax[k]
    else:
        self.sumax2 = self.sumax[k] * self.redsu[k]
        
    self.Su[k] = ifthenelse(self.Su_t[k] + self.Fa > self.sumax2, self.sumax2, self.Su_t[k] + self.Fa) 
    self.Quadd = ifthenelse(self.Su_t[k] + self.Fa > self.sumax2, self.Su_t[k] + self.Fa - self.sumax2, 0)
    self.SuN = self.Su[k] / self.sumax2
    self.SiN = self.Si[k] / self.imax[k]    
    
    self.Eu1 = max((self.PotEvaporation - self.Ei),0) * min(self.Su[k] / (self.sumax2 * self.LP[k]),1)        
    
    self.Qu1 = (self.Fa - self.Quadd) * (1 - (1 - self.SuN) ** self.beta[k])
    self.Perc1 = self.perc[k] * self.SuN
    self.Su[k] = self.Su_t[k] + (self.Fa - self.Quadd) - self.Qu1 - self.Eu - self.Perc1

    self.Su_diff = ifthenelse(self.Su[k] < 0, self.Su[k], 0)
    self.Eu = self.Eu1 + (self.Eu1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Qu = self.Qu1 + (self.Qu1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff
    self.Perc = ifthenelse (self.Perc1 > 0, self.Perc1 + (self.Perc1 / ifthenelse(self.Qu1 + self.Eu1 + self.Perc1 > 0 , self.Qu1 + self.Eu1 + self.Perc1 , 1)) * self.Su_diff, self.Perc1)
    self.Su[k] = self.Su_t[k] + (self.Fa - self.Quadd) - self.Eu - self.Qu - self.Perc   
    self.Su[k] = ifthenelse(self.Su[k] < 0, 0 , self.Su[k])    
    self.Su_diff2 = ifthen(self.Su[k] < 0, self.Su[k]) 

    self.Cap = min(self.cap[k] * (1 - self.Su[k] / self.sumax2), self.Ss)    
    self.Su[k] = self.Su[k] + self.Cap    
    
    self.wbSu_[k] = self.Fa - self.Eu - self.Qu - self.Quadd - self.Perc + self. Cap - self.Su[k] + self.Su_t[k]
    
    self.Eu_[k] = self.Eu
    self.Qu_[k] = self.Qu + self.Quadd
    self.Cap_[k] = self.Cap
    self.Perc_[k] = self.Perc

   
