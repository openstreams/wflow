# -*- coding: utf-8 -*-
"""
Created on Thu Apr 03 17:06:26 2014

@author: TEuser

List all function versions

not all functions are still in this file, the older functions can be found
in h:\My Documents\memo's\python scripts\wflow\
"""


try:
    from wflow.wf_DynamicFramework import *
except ImportError:
    from .wf_DynamicFramework import *


def groundWater_no_reservoir(self):
    """
    Recharge = Qu * D (per class)
    QsinClass = sum of recharge, percolation and capillary rise per class
    Qsin = QsinClass, summed over areal fraction and all classes

    """

    self.Recharge = [x * y for x, y in zip(self.Qu_, self.D)]
    self.QsinClass = [
        a + b - c for a, b, c in zip(self.Recharge, self.Perc_, self.Cap_)
    ]  # all fluxes are summed per class
    self.QsinTemp = sum(
        [a * b for a, b in zip(self.QsinClass, self.percent)]
    )  # fluxes are summed per cell according to percentage of class
    self.Qsin = pcr.areatotal(
        self.QsinTemp * self.percentArea, pcr.nominal(self.TopoId)
    )  # areatotal is taken, according to area percentage of cell

    self.Qs = self.Qsin
    self.Ss = 0.0
    self.wbSs = self.Qsin - self.Qs - self.Ss + self.Ss_t


def groundWaterCombined3(self):
    """
    - When using this option, the declaration of Ss and Qs should not be divided into separate classes
    - Stocks are lumped for the entire catchment
    - Differences with 'groundWaterCombined' are the order of calculations to try to decrease the numerical error --> turns out there is no difference
    - cannot be used combined with Sus for wetland areas
    """

    self.Recharge = [x * y for x, y in zip(self.Qu_, self.D)]
    self.QsinClass = [
        a + b - c for a, b, c in zip(self.Recharge, self.Perc_, self.Cap_)
    ]  # all fluxes are summed per class
    self.QsinTemp = sum(
        [a * b for a, b in zip(self.QsinClass, self.percent)]
    )  # fluxes are summed per cell according to percentage of class
    self.Qsin = pcr.areatotal(
        self.QsinTemp * self.percentArea, pcr.nominal(self.TopoId)
    )  # areatotal is taken, according to area percentage of cell

    self.Qs = self.Ss * self.Ks[0]
    # self.Ss = self.Ss * pcr.exp(-self.Ks[0]) + self.Qsin

    # add a gain/lossterm based on constant value
    self.Gain = pcr.ifthenelse(
        self.Closure < 0.0, -pcr.min(self.Ss, -self.Closure), self.Closure
    )  # if negative, limited by storage in Ss
    self.Ss = self.Ss * pcr.exp(-self.Ks[0]) + self.Qsin + self.Gain

    self.wbSs = self.Qsin - self.Qs - self.Ss + self.Ss_t + self.Gain

    self.Qs_ = self.Qs
    self.Recharge_ = self.Recharge
    self.QsinClass_ = self.QsinClass
    self.QsinTemp_ = self.QsinTemp
    self.Qsin_ = self.Qsin
    self.Gain_ = self.Gain
