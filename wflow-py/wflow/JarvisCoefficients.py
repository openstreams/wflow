# -*- coding: utf-8 -*-
"""
Created on Wed Oct 01 14:15:26 2014

@author: TEuser

TO DO:
compile new input time series
calculate Jarvis equations based on daily data 
calculate Jarvis coefficient for soil moisture based on hourly data
calculate daily reference open water 
calculate daily potential evaporation
donwscale dialy potential evaporation to hourly values
use hourly evaporation as model input

Calculations based on Zhou et al 2006 and Stewart 1998
"""

import numpy
import pdb
from copy import copy as copylist

try:
    from  wflow.wf_DynamicFramework import *
except ImportError:
    from  wf_DynamicFramework import *
import scipy


def calcEp(self,k):
    """
    REQUIRED INPUT:
    """
#    resistenceAeroD(self)
#    potential_evaporation(self,k)
    downscale_evaporation(self,k)
    
def calcEpSnow(self,k):
    """
    REQUIRED INPUT:
	daily input data
    """
#    resistenceAeroD(self)
#    potential_evaporation(self,k)
    self.EpDaySnow = self.EpDay
    downscale_evaporation_snow(self,k)

def calcEpSnowHour(self,k):
    """
    REQUIRED INPUT:
	hourly input data
    """
#    resistenceAeroD(self)
#    potential_evaporation(self,k)
    self.EpHour = self.PotEvaporation / self.lamda * self.lamdaS
    
def calcEu_laiFixed(self,k):
    """
    REQUIRED INPUT:
    """
    JC_temperature(self,k)
    JC_vapourDeficit(self,k)
    JC_LAIeffective(self,k)
    JC_solarRadiation(self,k)
    JC_soilMoisture(self,k)
    resistenceAeroD(self)
    resistenceTotal(self,k)
#    potential_evaporation(self,k)
    downscale_evaporation(self,k)
    self.EupHour = self.EpHour - self.Ei_[k]
    self.Eu = self.EupHour * self.k
    
def calcEu(self,k,n):
    """
    REQUIRED INPUT:
    n = number of evaporation fluxes above the reservoir to be calculated (to decrease Ep)
    """
    JC_temperature(self,k)
    JC_vapourDeficit(self,k)
    JC_solarRadiation(self,k)
    JC_soilMoisture(self,k)
    resistenceAeroD(self)
    resistenceTotal_laiHRU(self,k)
    downscale_evaporation(self,k)
    if n == 1:
        self.EupHour = self.EpHour - self.Ei_[k]
    elif n == 2:
        self.EupHour = self.EpHour - self.Ei_[k] - self.Ea_[k]
    self.Eu = self.EupHour * self.k


def JC_temperature(self,k):
    """
    REQUIRED INPUT:
        - mean daily temperature (K)
    PARAMETERS:
        - Topt (optimum temperature for transpiration, based on elevation and latitude (Cui et al, 2012))
    """
    self.JC_ToptC = self.JC_Topt - 273.15   # changed on 9/4/2015, before Topt in K 
    JC_temp1 = ifthenelse(self.Tmean < 273.15, 0, ifthenelse(pcrand(self.Tmean >= self.JC_Topt - 1, self.Tmean <= self.JC_Topt + 1), 1, 1 - self.JC_Topt ** -2 * (self.Tmean - self.JC_Topt) ** 2))
    self.JC_temp = ifthenelse(JC_temp1 < 0, 0, JC_temp1)
    self.JC_temp_[k] = self.JC_temp

def JC_vapourDeficit(self,k):
    """
    REQUIRED INPUT:
        - vapour pressure deficit (kPa)
    PARAMETERS:
        - D05 (vapourpressure deficit halfway between 1 and Cd2, fixed at 1.5 kPa (Matsumoto et al. 2008))
        - Cd1 (first vapour pressure parameter, fixed at 3 (Matsumoto et al. 2008))
        - Cd2 (second vapour pressure parameter, fixed at 0.1 (Matsumoto et al. 2008))
        """
    
    denom = 1 + (self.vpd / self.JC_D05[k]) ** self.JC_cd1[k]
    JC_vpd1 = (1 / denom) * (1 - self.JC_cd2[k]) + self.JC_cd2[k]    
    self.JC_vpd = ifthenelse(JC_vpd1 < 0, 0, JC_vpd1)   
    self.JC_vpd_[k] = self.JC_vpd
   
def JC_LAIeffective(self,k):
    """
    REQUIRED INPUT:
        - LAI (-)
    PARAMETERS:
        - none (Allen et al., 2006 & Zhou et al., 2006)   
    """
    
    self.JC_laiEff = self.LAI / (0.2 * self.LAI + 1)        

    
def JC_solarRadiation(self,k):
    """
    REQUIRED INPUT:
        - incoming short wave radiation (W/m2)
    PARAMETERS:
        - Cr (radiation stress parameter, fixed at 100 (Zhou et al. 2006))
    """    
    
    rad_si_Wm2 = self.rad_si / 86400
    JC_rad1 = rad_si_Wm2 * (1 + self.JC_cr[k]/1000) * (1 / (self.JC_cr[k] + rad_si_Wm2))
    self.JC_rad = ifthenelse(JC_rad1 < 0, 0, JC_rad1) 
    self.JC_rad_[k] = self.JC_rad
    
def JC_soilMoisture(self,k):
    """
    REQUIRED INPUT:
        - storage unsaturated zone
    PARAMETERS:
        - Cuz (soil moisture stress parameter, fixed at 0.07 (Matsumoto et al., 2008))
        - SuFC (level of saturation in soil at field capacity)
        - SuWP (level of saturation in soil at wilting point)
    """
    SuN = self.Su[k] / self.sumax[k]   
    JC_sm1 = ifthenelse(SuN < self.SuWP[k], 0, ifthenelse(SuN > self.SuFC[k], 1,\
        (SuN - self.SuWP[k]) * (self.SuFC[k] - self.SuWP[k] + self.JC_cuz[k]) / ((self.SuFC[k] - self.SuWP[k]) * (SuN - self.SuWP[k] + self.JC_cuz[k]))))
    self.JC_sm = ifthenelse(JC_sm1 < 0, 0, JC_sm1)
    self.JC_sm_[k] = self.JC_sm

def resistenceAeroD(self):
    """
    REQUIRED INPUT:
        - all Jarvis stress functions
        - windspeed at 2m height (m/s)
        - sgamma (slope of saturation vapour function)
    PARAMETERS:
        - rstmin (minimal stomatal resistence, depending on land use, s/m)
        - gamma (psychrometric constant)
    """
    self.aeroRes = 245 / (86400 * (0.54 * self.wind2m + 0.5))

def resistenceTotal(self,k):
    """
    REQUIRED INPUT:
        - all Jarvis stress functions
        - windspeed at 2m height (m/s)
        - sgamma (slope of saturation vapour function)
    PARAMETERS:
        - rstmin (minimal stomatal resistence, depending on land use, s/m)
        - gamma (psychrometric constant)
    """
    JC_all = self.JC_laiEff * self.JC_rad * self.JC_vpd * self.JC_temp * self.JC_sm
    stomRes1 = ifthenelse(JC_all == 0, 50000, self.JC_rstmin[k] / JC_all)
    stomRes2 = ifthenelse(stomRes1 > 50000, 50000, stomRes1)
    self.stomRes = stomRes2 / (3600*24)
    
    self.k = 1 / (1 + self.stomRes * self.gamma / (self.aeroRes * (self.sgamma + self.gamma)))    
    self.JC_k_[k] = self.k 
    
def resistenceTotal_laiHRU(self,k):
    """
    REQUIRED INPUT:
        - all Jarvis stress functions
        - windspeed at 2m height (m/s)
        - sgamma (slope of saturation vapour function)
    PARAMETERS:
        - rstmin (minimal stomatal resistence, depending on land use, s/m)
        - gamma (psychrometric constant)
    """
    JC_all = self.JC_rad * self.JC_vpd * self.JC_temp * self.JC_sm
    stomRes1 = ifthenelse(JC_all == 0, 50000, self.rst_lai[k] / JC_all)
    stomRes2 = ifthenelse(stomRes1 > 50000, 50000, stomRes1)
    self.stomRes = stomRes2 / (3600*24)
    
    self.k = 1 / (1 + self.stomRes * self.gamma / (self.aeroRes * (self.sgamma + self.gamma)))    
    self.JC_k_[k] = self.k 

def potential_evaporation(self,k):
    """
    REQUIRED INPUT:
        - net radiation
        - vapour pressure deficit
        - sgamma (slope of staturation vapour function)
        - aeroRes (aerodynamic resistance)
    PARAMETERS:
        - gamma (psychrometric constant)
        - Cp (specific heat of moist air, 1.01 MJkg-1K-1)
        - rhoA (density of air)
        - rhoW (density of water)
        - lamda (latent heat of vaporisation)
    """
    nom = self.sgamma * self.Rn + self.rhoA * self.Cp * self.vpd / self.aeroRes 
    denom = self.rhoW * self.lamda * (self.sgamma + self.gamma)
    self.EpDay = (nom / denom) * 1000
    self.EpD_[k] = self.EpDay
    

def downscale_evaporation(self,k):
    """
    REQUIRED INPUT:
        daily evaporation (EpDay)
        - hour of the day (x; derived from self.thestep)
        - start of the day (derived from global radiation)
        - end of the day (derived from global radiation)
    PARAMETERS:
    - 
    """
    
    teller = numpy.max(pcr2numpy(self.thestep,nan))
    x = teller - floor(teller/24) * 24 * scalar(self.TopoId)
    DL = self.DE - self.DS + 1
    P = 2 * pi / (2 * DL)                              # period
    SH = DL - 12                                        #horizontal shift of new signal
    aN = -1 * self.EpDay * P                       # nominator of the amplitude
    aDN = sin((P * (self.DE + SH)) * 180 / pi) - sin((P * (self.DS + SH)) * 180 / pi)     # denominator of the amplitude
    ampl = aN / aDN                                     # amplitude of new signal
    
    self.EpHour = max(ifthenelse(pcrand(x >= self.DS, x <= self.DE), -1 * ampl * cos((P * (x + SH)) * 180 / pi), 0), 0)
    

def downscale_evaporation_snow(self,k):
    """
    REQUIRED INPUT:
        daily evaporation (EpDay)
        - hour of the day (x; derived from self.teller)
        - start of the day (derived from global radiation)
        - end of the day (derived from global radiation)
    PARAMETERS:
    - 
    """
    
    teller = numpy.max(pcr2numpy(self.thestep,nan))
    x = teller - floor(teller/24) * 24 * scalar(self.TopoId)
    DL = self.DE - self.DS + 1
    P = 2 * pi / (2 * DL)                              # period
    SH = DL - 12                                         #horizontal shift of new signal  
    aN = -1 * self.EpDaySnow * P                       # nominator of the amplitude
    aDN = sin((P * (self.DE + SH)) * 180 / pi) - sin((P * (self.DS + SH)) * 180 / pi)     # denominator of the amplitude
    ampl = aN / aDN                                     # amplitude of new signal
    
    self.EpHour = max(ifthenelse(pcrand(x >= self.DS, x <= self.DE), -1 * ampl * cos((P * (x + SH)) * 180 / pi), 0), 0)    
    
      
    
    
    
    