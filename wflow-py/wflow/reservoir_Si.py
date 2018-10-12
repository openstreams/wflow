# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 18:01:39 2014

@author: teuser

List all function versions
"""
import numpy

try:
    from wflow.wf_DynamicFramework import *
except ImportError:
    from .wf_DynamicFramework import *
from . import JarvisCoefficients


def selectSiR(i):
    """
    not all functions are still in this file, the older functions can be found
    (with the same numbering) in h:\My Documents\memo's\python scripts\wflow\
    """
    if i == 1:
        name = "interception_overflow"
    if i == 2:
        name = "interception_overflow2"
    if i == 3:
        name = "interception_overflow_Ep"
    return name


def interception_no_reservoir(self, k):
    """
    Effective rainfall = rainfall
    Interception evaporation = 0.
    Storage in interception = 0.
    """
    self.Pe_[k] = max(self.Precipitation, 0)
    self.Ei_[k] = 0.0
    self.Si_[k] = 0.0
    self.wbSi_[k] = (
        self.Precipitation - self.Ei_[k] - self.Pe_[k] - self.Si[k] + self.Si_t[k]
    )


def interception_overflow2(self, k):
    """
    - Effective rainfall is all that does not fit into the interception reservoir
    - Outgoing fluxes are determined separately
    - Code for ini-file: 2
    """

    self.Pe = max(self.Precipitation - (self.imax[k] - self.Si_t[k]), 0)
    self.Si[k] = self.Si_t[k] + (self.Precipitation - self.Pe)
    self.Ei = ifthenelse(
        self.Sw[k] == 0, min(self.PotEvaporation, self.Si[k]), 0
    )  # ifstatement added on 3-11-2015 for snow module
    self.Si[k] = self.Si[k] - self.Ei

    self.wbSi_[k] = self.Precipitation - self.Ei - self.Pe - self.Si[k] + self.Si_t[k]

    self.Pe = self.Pe + self.Qw_[k]  # added on 3-11-2015 for snow module
    self.Ei = self.Ei + (
        self.Ew_[k] / self.lamda * self.lamdaS
    )  # lambda added on 31-3-2016

    self.Ei_[k] = self.Ei
    self.Pe_[k] = self.Pe
    self.PotEvaporation_ = self.PotEvaporation

    if self.URFR_L:
        self.Ei = areatotal(self.Ei * self.percentArea, nominal(self.TopoId))
        self.Pe = areatotal(self.Pe * self.percentArea, nominal(self.TopoId))
        self.PotEvaporation = areatotal(
            self.PotEvaporation * self.percentArea, nominal(self.TopoId)
        )
        self.Si[k] = areatotal(self.Si[k] * self.percentArea, nominal(self.TopoId))


def interception_overflow_Ep(self, k):
    """
    - Effective rainfall is all that does not fit into the interception reservoir
    - Outgoing fluxes are determined separately
    - this version cannot be run with Su averaged (for the current code)
    - Code for ini-file: 3
    """

    JarvisCoefficients.calcEp(
        self, k
    )  # this line indicates that hourly profiles are made out of daily pot evap data based on start day (DS.tss) and day end (DE)
    self.PotEvaporation = cover(ifthenelse(self.EpHour >= 0, self.EpHour, 0), 0)

    self.Pe = max(self.Precipitation - (self.imax[k] - self.Si_t[k]), 0)
    self.Si[k] = self.Si_t[k] + (self.Precipitation - self.Pe)
    self.Ei = ifthenelse(
        self.Sw[k] == 0,
        min(
            (self.PotEvaporation - (self.Ew_[k] / self.lamda * self.lamdaS)), self.Si[k]
        ),
        0,
    )  # ifstatement added on 3-11-2015 for snow module, '-self.Ew_[k]' added on 17-2-2016
    self.Si[k] = self.Si[k] - self.Ei

    self.wbSi_[k] = self.Precipitation - self.Ei - self.Pe - self.Si[k] + self.Si_t[k]

    self.Pe = self.Pe + self.Qw_[k]  # added on 3-11-2015 for snow module
    self.Ei = self.Ei + (
        self.Ew_[k] / self.lamda * self.lamdaS
    )  # lambda added on 17-2-2016

    self.Ei_[k] = self.Ei
    self.Pe_[k] = self.Pe
    self.Ep_[k] = self.EpHour

    if self.URFR_L:
        self.Ei = areatotal(self.Ei * self.percentArea, nominal(self.TopoId))
        self.Pe = areatotal(self.Pe * self.percentArea, nominal(self.TopoId))
        self.PotEvaporation = areatotal(
            self.PotEvaporation * self.percentArea, nominal(self.TopoId)
        )
        self.Si[k] = areatotal(self.Si[k] * self.percentArea, nominal(self.TopoId))
