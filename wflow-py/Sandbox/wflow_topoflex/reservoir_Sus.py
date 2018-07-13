# -*- coding: utf-8 -*-
"""
Created on Wed Jul 02 15:05:31 2014

@author: teuser

List all function versions

- Combination of unsaturated and saturated zone
"""


import numpy
import pdb
from copy import copy as copylist

try:
    from wflow.wf_DynamicFramework import *
except ImportError:
    from wf_DynamicFramework import *
import scipy


def selectSusR(i):
    if i == 1:
        name = "unsatSatZone_GWout"
    elif i == 2:
        name = "unsatSatZone_noGWout"
    elif i == 3:
        name = "unsatSatZone_noGWout_VSA"
    return name


def unsatSatZone_noGWout_VSA(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Evaporation rate constrained by moisture stress via LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato
    - the saturated part of the reservoir has only drainage outflow
    - the overflow in case of full saturation increases with contributing area (contributing 
    area is estimated, linear relation is derived from HAND distribution in wetland)
    - Code for ini-file: 3
    """

    #    pdb.set_trace()
    A = min(
        self.Co[k]
        * (
            1
            + ((self.Sus_t[k] - self.susmax2[k]) * (1 - self.Co[k]))
            / ((self.susmax3[k] - self.susmax2[k]) * self.Co[k])
        ),
        1,
    )
    self.Qo = ifthenelse(self.Sus_t[k] >= self.susmax2[k], self.Pe * A, 0)
    self.Sus[k] = self.Sus_t[k] + (self.Pe - self.Qo)
    self.Eu1 = max((self.PotEvaporation - self.Ei), 0) * min(
        self.Sus[k] / (self.susmax2[k] * self.LP[k]), 1
    )
    self.Qd1 = self.Kd[k] * max((self.Sus[k] - self.susmax1[k]), 0)

    self.Sus[k] = self.Sus_t[k] + (self.Pe - self.Qo) - self.Eu1 - self.Qd1

    self.Su_diff = ifthenelse(self.Sus[k] < 0, self.Sus[k], 0)
    self.Eu = (
        self.Eu1
        + (self.Eu1 / ifthenelse(self.Qd1 + self.Eu1 > 0, self.Qd1 + self.Eu1, 1))
        * self.Su_diff
    )
    self.Qd = (
        self.Qd1
        + (self.Qd1 / ifthenelse(self.Qd1 + self.Eu1 > 0, self.Qd1 + self.Eu1, 1))
        * self.Su_diff
    )
    self.Sus[k] = self.Sus_t[k] + self.Pe - self.Qo - self.Eu - self.Qd
    self.Sus[k] = ifthenelse(self.Sus[k] < 0, 0, self.Sus[k])
    self.Su_diff2 = ifthen(self.Sus[k] < 0, self.Sus[k])

    self.wbSus_[k] = self.Pe - self.Qo - self.Qd - self.Eu - self.Sus[k] + self.Sus_t[k]

    self.Eu_[k] = self.Eu
    self.Qo_[k] = self.Qo
    self.Qd_[k] = self.Qd
    self.Su_diff_[k] = self.Su_diff


def unsatSatZone_noGWout(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Evaporation rate constrained by moisture stress via LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato
    - the saturated part of the reservoir has only drainage outflow
    - Code for ini-file: 2
    """

    self.Qo = max(self.Pe - (self.susmax2[k] - self.Sus_t[k]), 0)
    self.Sus[k] = self.Sus_t[k] + self.Pe - self.Qo
    self.Eu1 = max((self.PotEvaporation - self.Ei), 0) * min(
        self.Sus[k] / (self.susmax2[k] * self.LP[k]), 1
    )
    self.Qd1 = self.Kd[k] * max((self.Sus[k] - self.susmax1[k]), 0)

    self.Sus[k] = self.Sus_t[k] + self.Pe - self.Qo - self.Eu1 - self.Qd1

    self.Su_diff = ifthenelse(self.Sus[k] < 0, self.Sus[k], 0)
    self.Eu = (
        self.Eu1
        + (self.Eu1 / ifthenelse(self.Qd1 + self.Eu1 > 0, self.Qd1 + self.Eu1, 1))
        * self.Su_diff
    )
    self.Qd = (
        self.Qd1
        + (self.Qd1 / ifthenelse(self.Qd1 + self.Eu1 > 0, self.Qd1 + self.Eu1, 1))
        * self.Su_diff
    )
    self.Sus[k] = self.Sus_t[k] + self.Pe - self.Qo - self.Eu - self.Qd
    self.Sus[k] = ifthenelse(self.Sus[k] < 0, 0, self.Sus[k])
    self.Su_diff2 = ifthen(self.Sus[k] < 0, self.Sus[k])

    self.wbSus_[k] = self.Pe - self.Qo - self.Qd - self.Eu - self.Sus[k] + self.Sus_t[k]

    self.Eu_[k] = self.Eu
    self.Qo_[k] = self.Qo
    self.Qd_[k] = self.Qd
    self.Su_diff_[k] = self.Su_diff


def unsatSatZone_GWout(self, k):
    """
    - Potential evaporation is decreased by energy used for interception evaporation    
    - Evaporation rate constrained by moisture stress via LP
    - Outgoing fluxes are determined based on (value in previous timestep + inflow) 
    and if this leads to negative storage, the outgoing fluxes are corrected to rato
    - the saturated part of the reservoir has gw outflow, in addition to the drainage outflow
    - NIET UITGEBREID GETEST VOOR SLUITEN WB, NU OOK CONFLICT MET Qs
    - Code for ini-file: 1
    """

    self.Qo = max(self.Pe - (self.susmax2[k] - self.Sus_t[k]), 0)
    self.Sus[k] = self.Sus_t[k] + self.Pe - self.Qo
    self.Eu1 = max((self.PotEvaporation - self.Ei), 0) * min(
        self.Sus[k] / (self.susmax2[k] * self.LP[k]), 1
    )
    self.Qd1 = self.Kd[k] * max((self.Sus[k] - self.susmax1[k]), 0)

    self.Sus[k] = self.Sus_t[k] + self.Pe - self.Qo - self.Eu1 - self.Qd1

    self.Su_diff = ifthenelse(self.Sus[k] < 0, self.Sus[k], 0)
    self.Eu = (
        self.Eu1
        + (self.Eu1 / ifthenelse(self.Qd1 + self.Eu1 > 0, self.Qd1 + self.Eu1, 1))
        * self.Su_diff
    )
    self.Qd = (
        self.Qd1
        + (self.Qd1 / ifthenelse(self.Qd1 + self.Eu1 > 0, self.Qd1 + self.Eu1, 1))
        * self.Su_diff
    )
    self.Sus[k] = self.Sus_t[k] + self.Pe - self.Qo - self.Eu - self.Qd
    self.Sus[k] = ifthenelse(self.Sus[k] < 0, 0, self.Sus[k])
    self.Su_diff2 = ifthen(self.Sus[k] < 0, self.Sus[k])

    self.Qus = self.Ks[k] * self.Sus[k]
    self.Sus[k] = self.Sus[k] - self.Qus

    self.wbSus_[k] = (
        self.Pe - self.Qo - self.Qd - self.Eu - self.Qus - self.Sus + self.Sus_t
    )

    self.Eu_[k] = self.Eu
    self.Qo_[k] = self.Qo
    self.Qd_[k] = self.Qd
    self.Qus_[k] = self.Qus
    self.Su_diff_[k] = self.Su_diff
