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

import pcraster as pcr
try:
    from wflow.wf_DynamicFramework import *
except ImportError:
    from .wf_DynamicFramework import *
from . import JarvisCoefficients


def selectSaR(i):
    """
    not all functions are still in this file, the older functions can be found
    (with the same numbering) in h:\My Documents\memo's\python scripts\wflow\
    """
    if i == 1:
        name = "agriZone_Jarvis"
    elif i == 2:
        name = "agriZone_Ep"
    elif i == 3:
        name = "agriZone_Ep_Sa"
    elif i == 4:
        name = "agriZone_Ep_Sa_cropG"
    elif i == 5:
        name = "agriZone_Ep_Sa_cropG_beta"
    elif i == 6:
        name = "agriZone_Ep_Sa_beta"
    elif i == 7:
        name = "agriZone_Ep_Sa_beta_frost"
    elif i == 8:
        name = "agriZone_Ep_Sa_beta_Fvar"
    elif i == 9:
        name = "agriZone_hourlyEp_Sa_beta_Fvar"
    elif i == 10:
        name = "agriZone_hourlyEp_Sa_beta_frost"
    elif i == 11:
        name = "agriZone_hourlyEp_Sa_beta_frostSamax"
    elif i == 12:
        name = "agriZone_Ep_Sa_beta_frostSamax"
    elif i == 13:
        name = "agriZone_Ep_Sa_beta_frostSamax_surfTemp"
    return name


def agriZone_no_reservoir(self, k):
    """
    This function is used when no unsaturated zone reservoir is used and only
    passes fluxes from the upper reservoirs to the lower
    self.Qa_[k] = 0.
    self.Ea_[k] = 0.
    self.Sa[k] = 0.
    self.Fa_[k] = Pe
    Storage in unsaturated zone = 0.
    """
    self.Qa_[k] = 0.0
    self.Ea_[k] = 0.0
    self.Sa[k] = 0.0
    self.Fa_[k] = pcr.max(self.Pe_[k], 0)
    self.wbSa_[k] = (
        self.Pe_[k]
        - self.Ea_[k]
        - self.Qa_[k]
        - self.Fa_[k]
        - self.Sa[k]
        + self.Sa_t[k]
    )


def agriZone_Jarvis(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on Jarvis stress functions
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa
    - Code for ini-file: 1
    """
    self.Qa = pcr.max(self.Pe - (self.samax[k] - self.Sa_t[k]), 0)
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    JarvisCoefficients.calcEu(
        self, k, 1
    )  # calculation of Ea based on Jarvis stress functions
    self.Ea1 = self.Eu

    self.Fa1 = self.Fmin[k] + (self.Fmax[k] - self.Fmin[k]) * e ** (
        -self.decF[k] * self.SuN
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Fa = (
        self.Fa1
        + (self.Fa1 / pcr.ifthenelse(self.Fa1 + self.Ea1 > 0, self.Fa1 + self.Ea1, 1))
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (self.Ea1 / pcr.ifthenelse(self.Fa1 + self.Ea1 > 0, self.Fa1 + self.Ea1, 1))
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Ea - self.Fa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = self.Pe - self.Ea - self.Qa - self.Fa - self.Sa[k] + self.Sa_t[k]

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa
    self.Fa_[k] = self.Fa


def agriZone_Ep(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa
    - Code for ini-file: 2
    """
    JarvisCoefficients.calcEp(self, k)
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour >= 0, self.EpHour, 0), 0)

    self.Qa = pcr.max(self.Pe - (self.samax[k] - self.Sa_t[k]), 0)
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax[k] * self.LP[k]), 1
    )

    self.Fa1 = self.Fmin[k] + (self.Fmax[k] - self.Fmin[k]) * e ** (
        -self.decF[k] * self.SuN
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Fa = (
        self.Fa1
        + (self.Fa1 / pcr.ifthenelse(self.Fa1 + self.Ea1 > 0, self.Fa1 + self.Ea1, 1))
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (self.Ea1 / pcr.ifthenelse(self.Fa1 + self.Ea1 > 0, self.Fa1 + self.Ea1, 1))
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Ea - self.Fa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = self.Pe - self.Ea - self.Qa - self.Fa - self.Sa[k] + self.Sa_t[k]

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa
    self.Fa_[k] = self.Fa


def agriZone_Ep_Sa(self, k):
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
    JarvisCoefficients.calcEp(self, k)
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour >= 0, self.EpHour, 0), 0)

    self.Qa = pcr.max(self.Pe - (self.samax[k] - self.Sa_t[k]), 0)
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax[k] * self.LP[k]), 1
    )

    self.Fa1 = pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),
        0,
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Fa = (
        self.Fa1
        + (self.Fa1 / pcr.ifthenelse(self.Fa1 + self.Ea1 > 0, self.Fa1 + self.Ea1, 1))
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (self.Ea1 / pcr.ifthenelse(self.Fa1 + self.Ea1 > 0, self.Fa1 + self.Ea1, 1))
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Ea - self.Fa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = self.Pe - self.Ea - self.Qa - self.Fa - self.Sa[k] + self.Sa_t[k]

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa
    self.Fa_[k] = self.Fa


def agriZone_Ep_Sa_cropG(self, k):
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
    JarvisCoefficients.calcEp(self, k)
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour >= 0, self.EpHour, 0), 0)

    self.samax2 = self.samax[k] * self.cropG
    self.Qaadd = pcr.max(self.Sa_t[k] - self.samax2, 0)

    self.Qa = pcr.max(self.Pe - (self.samax2 - self.Sa_t[k]), 0) + self.Qaadd
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax2 * self.LP[k]), 1
    )

    self.Fa1 = pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),
        0,
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Fa = (
        self.Fa1
        + (self.Fa1 / pcr.ifthenelse(self.Fa1 + self.Ea1 > 0, self.Fa1 + self.Ea1, 1))
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (self.Ea1 / pcr.ifthenelse(self.Fa1 + self.Ea1 > 0, self.Fa1 + self.Ea1, 1))
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qa) - self.Ea - self.Fa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = self.Pe - self.Ea - self.Qa - self.Fa - self.Sa[k] + self.Sa_t[k]

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa
    self.Fa_[k] = self.Fa


def agriZone_Ep_Sa_cropG_beta(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa --> incorporation of beta function
    - Fa is based on storage in Sa
    - Code for ini-file: 5
    """

    JarvisCoefficients.calcEp(self, k)
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour >= 0, self.EpHour, 0), 0)

    self.samax2 = self.samax[k] * self.cropG
    self.Qaadd = pcr.max(self.Sa_t[k] + self.Pe - self.samax2, 0)

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax2 * self.LP[k]), 1
    )
    self.Qa1 = (self.Pe - self.Qaadd) * (1 - (1 - self.SaN) ** self.beta[k])
    self.Fa1 = pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),
        0,
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Qa1 - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Qa = (
        self.Qa1
        + (
            self.Qa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Fa = (
        self.Fa1
        + (
            self.Fa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (
            self.Ea1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Ea - self.Fa - self.Qa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = (
        self.Pe - self.Ea - self.Qa - self.Qaadd - self.Fa - self.Sa[k] + self.Sa_t[k]
    )

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa + self.Qaadd
    self.Fa_[k] = self.Fa


def agriZone_Ep_Sa_beta(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa --> incorporation of beta function
    - Fa is based on storage in Sa
    - Code for ini-file: 6
    """

    JarvisCoefficients.calcEp(self, k)
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour >= 0, self.EpHour, 0), 0)

    self.samax2 = self.samax[k] * pcr.scalar(self.catchArea)
    self.Qaadd = pcr.max(self.Sa_t[k] + self.Pe - self.samax2, 0)

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd)
    self.SaN = pcr.min(pcr.max(self.Sa[k] / self.samax2, 0), 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax2 * self.LP[k]), 1
    )
    self.Qa1 = (self.Pe - self.Qaadd) * (1 - (1 - self.SaN) ** self.beta[k])
    self.Fa1 = pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),
        0,
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Qa1 - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Qa = (
        self.Qa1
        + (
            self.Qa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Fa = (
        self.Fa1
        + (
            self.Fa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (
            self.Ea1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Ea - self.Fa - self.Qa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = (
        self.Pe - self.Ea - self.Qa - self.Qaadd - self.Fa - self.Sa[k] + self.Sa_t[k]
    )

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa + self.Qaadd
    self.Fa_[k] = self.Fa


def agriZone_hourlyEp_Sa_beta(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa --> incorporation of beta function
    - Fa is based on storage in Sa
    - Code for ini-file:
    """

    # JarvisCoefficients.calcEp(self,k)
    # self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour >= 0, self.EpHour, 0),0)

    self.samax2 = self.samax[k] * pcr.scalar(self.catchArea)
    self.Qaadd = pcr.max(self.Sa_t[k] + self.Pe - self.samax2, 0)

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd)
    self.SaN = pcr.min(pcr.max(self.Sa[k] / self.samax2, 0), 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax2 * self.LP[k]), 1
    )
    self.Qa1 = (self.Pe - self.Qaadd) * (1 - (1 - self.SaN) ** self.beta[k])
    self.Fa1 = pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),
        0,
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Qa1 - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Qa = (
        self.Qa1
        + (
            self.Qa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Fa = (
        self.Fa1
        + (
            self.Fa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (
            self.Ea1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Ea - self.Fa - self.Qa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = (
        self.Pe - self.Ea - self.Qa - self.Qaadd - self.Fa - self.Sa[k] + self.Sa_t[k]
    )

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa + self.Qaadd
    self.Fa_[k] = self.Fa


def agriZone_Ep_Sa_beta_frost(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa --> incorporation of beta function
    - Fa is based on storage in Sa
    - Fa is decreased in case of frozen soil
    - Code for ini-file: 7
    """

    JarvisCoefficients.calcEp(self, k)
    self.PotEvaporation = self.EpHour

    self.samax2 = self.samax[k] * pcr.scalar(self.catchArea)
    self.Qaadd = pcr.max(self.Sa_t[k] + self.Pe - self.samax2, 0)
    self.FrDur[k] = pcr.min(
        self.FrDur[k]
        + (self.Tmean - 273.15) / 86400 * self.timestepsecs * self.dayDeg[k],
        0,
    )

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax2 * self.LP[k]), 1
    )
    self.Qa1 = (self.Pe - self.Qaadd) * (1 - (1 - self.SaN) ** self.beta[k])
    self.Ft = pcr.min(
        pcr.max(
            self.FrDur[k] / (self.FrDur1[k] - self.FrDur0[k])
            - self.FrDur0[k] / (self.FrDur1[k] - self.FrDur0[k]),
            0,
        ),
        1,
    )
    self.Fa1 = self.Ft * pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),
        0,
    )

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Qa1 - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Qa = (
        self.Qa1
        + (
            self.Qa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Fa = (
        self.Fa1
        + (
            self.Fa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (
            self.Ea1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Ea - self.Fa - self.Qa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = (
        self.Pe - self.Ea - self.Qa - self.Qaadd - self.Fa - self.Sa[k] + self.Sa_t[k]
    )

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa + self.Qaadd
    self.Fa_[k] = self.Fa
    self.Ft_[k] = self.Ft


def agriZone_hourlyEp_Sa_beta_frost(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa --> incorporation of beta function
    - Fa is based on storage in Sa
    - Fa is decreased in case of frozen soil
    - Code for ini-file: 10
    """

    # JarvisCoefficients.calcEp(self,k)
    # self.PotEvaporation = self.EpHour

    self.samax2 = self.samax[k] * pcr.scalar(self.catchArea)
    self.Qaadd = pcr.max(self.Sa_t[k] + self.Pe - self.samax2, 0)
    self.FrDur[k] = pcr.min(
        self.FrDur[k] + (self.Temperature) / 86400 * self.timestepsecs * self.dayDeg[k],
        0,
    )

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax2 * self.LP[k]), 1
    )
    self.Qa1 = (self.Pe - self.Qaadd) * (1 - (1 - self.SaN) ** self.beta[k])
    self.Ft = pcr.min(
        pcr.max(
            self.FrDur[k] / (self.FrDur1[k] - self.FrDur0[k])
            - self.FrDur0[k] / (self.FrDur1[k] - self.FrDur0[k]),
            0,
        ),
        1,
    )
    self.Fa1 = self.Ft * pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),
        0,
    )

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Qa1 - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Qa = (
        self.Qa1
        + (
            self.Qa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Fa = (
        self.Fa1
        + (
            self.Fa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (
            self.Ea1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Ea - self.Fa - self.Qa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = (
        self.Pe - self.Ea - self.Qa - self.Qaadd - self.Fa - self.Sa[k] + self.Sa_t[k]
    )

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa + self.Qaadd
    self.Fa_[k] = self.Fa
    self.Ft_[k] = self.Ft


def agriZone_hourlyEp_Sa_beta_frostSamax(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa --> incorporation of beta function
    - Fa is based on storage in Sa
    - Fa is decreased in case of frozen soil
    - Code for ini-file: 11
    """

    # JarvisCoefficients.calcEp(self,k)
    # self.PotEvaporation = self.EpHour

    self.FrDur[k] = pcr.min(self.FrDur[k] + (self.Temperature) * self.dayDeg[k], 0)
    self.Ft = pcr.min(
        pcr.max(
            self.FrDur[k] / (self.FrDur1[k] - self.FrDur0[k])
            - self.FrDur0[k] / (self.FrDur1[k] - self.FrDur0[k]),
            0.1,
        ),
        1,
    )

    self.samax2 = self.samax[k] * pcr.scalar(self.catchArea) * self.Ft
    self.Qaadd = pcr.max(self.Sa_t[k] + self.Pe - self.samax2, 0)

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax2 * self.LP[k]), 1
    )
    self.Qa1 = (self.Pe - self.Qaadd) * (1 - (1 - self.SaN) ** self.beta[k])
    self.Fa1 = pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * pcr.exp((-self.decF[k] * (1 - self.SaN))),
        0,
    )

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Qa1 - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Qa = (
        self.Qa1
        + (
            self.Qa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Fa = (
        self.Fa1
        + (
            self.Fa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (
            self.Ea1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Ea - self.Fa - self.Qa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = (
        self.Pe - self.Ea - self.Qa - self.Qaadd - self.Fa - self.Sa[k] + self.Sa_t[k]
    )

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa + self.Qaadd
    self.Fa_[k] = self.Fa
    self.Ft_[k] = self.Ft


def agriZone_Ep_Sa_beta_frostSamax(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa --> incorporation of beta function
    - Fa is based on storage in Sa
    - Fa is decreased in case of frozen soil
    - Code for ini-file: 12
    """

    JarvisCoefficients.calcEp(self, k)
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour >= 0, self.EpHour, 0), 0)

    self.FrDur[k] = pcr.min(
        self.FrDur[k]
        + pcr.ifthenelse(
            self.Temperature > 0, self.ratFT[k] * self.Temperature, self.Temperature
        )
        * self.dayDeg[k],
        0,
    )
    self.Ft = pcr.min(
        pcr.max(
            self.FrDur[k] / (self.FrDur1[k] - self.FrDur0[k])
            - self.FrDur0[k] / (self.FrDur1[k] - self.FrDur0[k]),
            self.samin[k],
        ),
        1,
    )

    self.samax2 = self.samax[k] * pcr.scalar(self.catchArea) * self.Ft
    self.Qaadd = pcr.max(self.Sa_t[k] + self.Pe - self.samax2, 0)

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax2 * self.LP[k]), 1
    )
    self.Qa1 = (self.Pe - self.Qaadd) * (1 - (1 - self.SaN) ** self.beta[k])

    self.Fa1 = pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),
        0,
    )

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Qa1 - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Qa = (
        self.Qa1
        + (
            self.Qa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Fa = (
        self.Fa1
        + (
            self.Fa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (
            self.Ea1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Ea - self.Fa - self.Qa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = (
        self.Pe - self.Ea - self.Qa - self.Qaadd - self.Fa - self.Sa[k] + self.Sa_t[k]
    )

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa + self.Qaadd
    self.Fa_[k] = self.Fa
    self.Ft_[k] = self.Ft


def agriZone_Ep_Sa_beta_frostSamax_surfTemp(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa --> incorporation of beta function
    - Fa is based on storage in Sa
    - Fa is decreased in case of frozen soil
    - Code for ini-file: 13
    """

    JarvisCoefficients.calcEp(self, k)
    self.PotEvaporation = self.EpHour

    self.FrDur[k] = pcr.min(
        self.FrDur[k]
        + pcr.ifthenelse(
            self.TempSurf > 0, self.ratFT[k] * self.TempSurf, self.TempSurf
        )
        * self.dayDeg[k],
        0,
    )
    self.Ft = pcr.min(
        pcr.max(
            self.FrDur[k] / (self.FrDur1[k] - self.FrDur0[k])
            - self.FrDur0[k] / (self.FrDur1[k] - self.FrDur0[k]),
            self.samin[k],
        ),
        1,
    )

    self.samax2 = self.samax[k] * pcr.scalar(self.catchArea) * self.Ft
    self.Qaadd = pcr.max(self.Sa_t[k] + self.Pe - self.samax2, 0)

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax2 * self.LP[k]), 1
    )
    self.Qa1 = (self.Pe - self.Qaadd) * (1 - (1 - self.SaN) ** self.beta[k])

    self.Fa1 = pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),
        0,
    )

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Qa1 - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Qa = (
        self.Qa1
        + (
            self.Qa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Fa = (
        self.Fa1
        + (
            self.Fa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (
            self.Ea1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Ea - self.Fa - self.Qa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = (
        self.Pe - self.Ea - self.Qa - self.Qaadd - self.Fa - self.Sa[k] + self.Sa_t[k]
    )

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa + self.Qaadd
    self.Fa_[k] = self.Fa
    self.Ft_[k] = self.Ft


def agriZone_Ep_Sa_beta_Fvar(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa --> incorporation of beta function
    - Fa is based on storage in Sa
    - Code for ini-file: 8
    """

    JarvisCoefficients.calcEp(self, k)
    self.PotEvaporation = self.EpHour

    self.samax2 = self.samax[k] * pcr.scalar(self.catchArea)
    self.Qaadd = pcr.max(self.Sa_t[k] + self.Pe - self.samax2, 0)

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax2 * self.LP[k]), 1
    )
    self.Qa1 = (self.Pe - self.Qaadd) * (1 - (1 - self.SaN) ** self.beta[k])
    self.Fa1 = self.cropG * pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),
        0,
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Qa1 - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Qa = (
        self.Qa1
        + (
            self.Qa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Fa = (
        self.Fa1
        + (
            self.Fa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (
            self.Ea1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Ea - self.Fa - self.Qa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = (
        self.Pe - self.Ea - self.Qa - self.Qaadd - self.Fa - self.Sa[k] + self.Sa_t[k]
    )

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa + self.Qaadd
    self.Fa_[k] = self.Fa


def agriZone_hourlyEp_Sa_beta_Fvar(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Formula for evaporation based on LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato --> Eu is 
    no longer taken into account for this correction
    - Qa u is determined from overflow from Sa --> incorporation of beta function
    - Fa is based on storage in Sa
    - Code for ini-file: 9
    """

    #    JarvisCoefficients.calcEp(self,k)
    #    self.PotEvaporation = self.EpHour

    self.samax2 = self.samax[k] * pcr.scalar(self.catchArea)
    self.Qaadd = pcr.max(self.Sa_t[k] + self.Pe - self.samax2, 0)

    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd)
    self.SaN = pcr.min(self.Sa[k] / self.samax2, 1)
    self.SuN = self.Su[k] / self.sumax[k]

    self.Ea1 = pcr.max((self.PotEvaporation - self.Ei), 0) * pcr.min(
        self.Sa[k] / (self.samax2 * self.LP[k]), 1
    )
    self.Qa1 = (self.Pe - self.Qaadd) * (1 - (1 - self.SaN) ** self.beta[k])
    self.Fa1 = self.cropG * pcr.ifthenelse(
        self.SaN > 0,
        self.Fmin[k]
        + (self.Fmax[k] - self.Fmin[k]) * e ** (-self.decF[k] * (1 - self.SaN)),
        0,
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Qa1 - self.Fa1 - self.Ea1

    self.Sa_diff = pcr.ifthenelse(self.Sa[k] < 0, self.Sa[k], 0)
    self.Qa = (
        self.Qa1
        + (
            self.Qa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Fa = (
        self.Fa1
        + (
            self.Fa1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Ea = (
        self.Ea1
        + (
            self.Ea1
            / pcr.ifthenelse(
                self.Fa1 + self.Ea1 + self.Qa1 > 0, self.Fa1 + self.Ea1 + self.Qa1, 1
            )
        )
        * self.Sa_diff
    )
    self.Sa[k] = self.Sa_t[k] + (self.Pe - self.Qaadd) - self.Ea - self.Fa - self.Qa
    self.Sa[k] = pcr.ifthenelse(self.Sa[k] < 0, 0, self.Sa[k])
    self.Sa_diff2 = pcr.ifthen(self.Sa[k] < 0, self.Sa[k])

    self.wbSa_[k] = (
        self.Pe - self.Ea - self.Qa - self.Qaadd - self.Fa - self.Sa[k] + self.Sa_t[k]
    )

    self.Ea_[k] = self.Ea
    self.Qa_[k] = self.Qa + self.Qaadd
    self.Fa_[k] = self.Fa
