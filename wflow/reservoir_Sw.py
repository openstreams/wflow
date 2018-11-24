# -*- coding: utf-8 -*-
"""
Created on Thu Apr 03 16:31:35 2014

@author: Teuser

List all function versions
"""


try:
    from wflow.wf_DynamicFramework import *
except ImportError:
    from .wf_DynamicFramework import *
from . import JarvisCoefficients


def selectSwR(i):

    if i == 1:
        name = "snow"
    if i == 2:
        name = "snowHour"

    return name


def snow_no_reservoir(self, k):
    """
    This function is used when no snow zone reservoir is used and only
    passes fluxes from the upper reservoirs to the lower
    Qw = Psnow
    Ew = 0.
    Storage in snow zone = 0.

    !!!still needs a final check!!!

    k is the class indication
    self contains all the variables of the model

    """
    try:
        JarvisCoefficients.calcEpSnow(self, k)
    except:
        JarvisCoefficients.calcEpSnowHour(self, k)
    self.PotEvaporation = self.EpHour
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour > 0, self.EpHour, 0), 0)

    self.Qw_[k] = pcr.max(self.PrecipitationSnow, 0)
    self.Ew_[k] = 0.0
    self.Sw[k] = 0.0
    self.wbSw_[k] = (
        self.Precipitation - self.Ew_[k] - self.Qw_[k] - self.Sw[k] + self.Sw_t[k]
    )


def snow(self, k):
    """
    - snow melt based on degree day factor and
    -
    - Code for ini-file: 1
    """
    JarvisCoefficients.calcEpSnow(self, k)
    self.PotEvaporation = self.EpHour
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour > 0, self.EpHour, 0), 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow

    self.Ew1 = pcr.max(pcr.min(self.PotEvaporation, self.Sw[k]), 0)
    self.Qw1 = pcr.max(self.Fm[k] * (self.Temperature - self.Tm[k]), 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew1 - self.Qw1

    self.Sw_diff = pcr.ifthenelse(self.Sw[k] < 0, self.Sw[k], 0)
    self.Ew = (
        self.Ew1
        + (self.Ew1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Qw = (
        self.Qw1
        + (self.Qw1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew - self.Qw
    self.Sw[k] = pcr.ifthenelse(self.Sw[k] < 0, 0, self.Sw[k])
    self.Sw_diff2 = pcr.ifthen(self.Sw[k] < 0, self.Sw[k])

    self.wbSw_[k] = (
        self.PrecipitationSnow - self.Ew - self.Qw - self.Sw[k] + self.Sw_t[k]
    )

    self.Ew_[k] = self.Ew
    self.Qw_[k] = self.Qw


def snow_rain(self, k):
    """
    - snow melt based on degree day factor and minimum surface temperature
    - meltfactor increases with temperature
    -
    - Code for ini-file: 6
    """

    JarvisCoefficients.calcEpSnow(self, k)
    # self.PotEvaporation = self.EpHour
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour > 0, self.EpHour, 0), 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow

    self.Fm2 = pcr.max(self.Fm[k] * self.Precipitation, self.Fm[k])
    self.Ew1 = pcr.max(pcr.min(self.PotEvaporation, self.Sw[k]), 0)
    self.Qw1 = pcr.max(
        pcr.min(self.Fm2 * (self.Temperature - self.Tm[k]), self.Sw[k]), 0
    )

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew1 - self.Qw1

    self.Sw_diff = pcr.ifthenelse(self.Sw[k] < 0, self.Sw[k], 0)
    self.Ew = (
        self.Ew1
        + (self.Ew1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Qw = (
        self.Qw1
        + (self.Qw1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew - self.Qw
    self.Sw[k] = pcr.ifthenelse(self.Sw[k] < 0, 0, self.Sw[k])
    self.Sw_diff2 = pcr.ifthen(self.Sw[k] < 0, self.Sw[k])

    self.wbSw_[k] = (
        self.PrecipitationSnow - self.Ew - self.Qw - self.Sw[k] + self.Sw_t[k]
    )

    self.Ew_[k] = self.Ew
    self.Qw_[k] = self.Qw


def snow_rain_hourlyEp(self, k):
    """
    - snow melt based on degree day factor and minimum surface temperature
    - meltfactor increases with temperature
    -
    - Code for ini-file: 6
    """

    JarvisCoefficients.calcEpSnowHour(self, k)
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour > 0, self.EpHour, 0), 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow

    self.Fm2 = pcr.max(self.Fm[k] * self.Precipitation, self.Fm[k])
    self.Ew1 = pcr.max(pcr.min(self.PotEvaporation, self.Sw[k]), 0)
    self.Qw1 = pcr.max(
        pcr.min(self.Fm2 * (self.Temperature - self.Tm[k]), self.Sw[k]), 0
    )

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew1 - self.Qw1

    self.Sw_diff = pcr.ifthenelse(self.Sw[k] < 0, self.Sw[k], 0)
    self.Ew = (
        self.Ew1
        + (self.Ew1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Qw = (
        self.Qw1
        + (self.Qw1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew - self.Qw
    self.Sw[k] = pcr.ifthenelse(self.Sw[k] < 0, 0, self.Sw[k])
    self.Sw_diff2 = pcr.ifthen(self.Sw[k] < 0, self.Sw[k])

    self.wbSw_[k] = (
        self.PrecipitationSnow - self.Ew - self.Qw - self.Sw[k] + self.Sw_t[k]
    )

    self.Ew_[k] = self.Ew
    self.Qw_[k] = self.Qw


def snow_rain_Tsurf(self, k):
    """
    - snow melt based on degree day factor and minimum surface temperature
    - meltfactor increases with temperature
    -
    - Code for ini-file: 3
    """

    JarvisCoefficients.calcEpSnow(self, k)
    # self.PotEvaporation = self.EpHour
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour > 0, self.EpHour, 0), 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow

    self.Fm2 = pcr.max(self.Fm[k] * self.Precipitation, self.Fm[k])
    self.Ew1 = pcr.max(pcr.min(self.PotEvaporation, self.Sw[k]), 0)
    self.Qw1 = pcr.max(pcr.min(self.Fm2 * (self.TempSurf - self.Tm[k]), self.Sw[k]), 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew1 - self.Qw1

    self.Sw_diff = pcr.ifthenelse(self.Sw[k] < 0, self.Sw[k], 0)
    self.Ew = (
        self.Ew1
        + (self.Ew1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Qw = (
        self.Qw1
        + (self.Qw1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew - self.Qw
    self.Sw[k] = pcr.ifthenelse(self.Sw[k] < 0, 0, self.Sw[k])
    self.Sw_diff2 = pcr.ifthen(self.Sw[k] < 0, self.Sw[k])

    self.wbSw_[k] = (
        self.PrecipitationSnow - self.Ew - self.Qw - self.Sw[k] + self.Sw_t[k]
    )

    self.Ew_[k] = self.Ew
    self.Qw_[k] = self.Qw


def snow_rain_TsurfAir(self, k):
    """
    - snow melt based on degree day factor and minimum surface temperature
    - meltfactor increases with temperature
    -
    - Code for ini-file: 4
    """
    JarvisCoefficients.calcEpSnow(self, k)
    # self.PotEvaporation = self.EpHour
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour > 0, self.EpHour, 0), 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow
    self.Temp = (self.TempSurf + self.Temperature) / 2

    self.Fm2 = pcr.max(self.Fm[k] * self.Precipitation, self.Fm[k])
    self.Ew1 = pcr.max(pcr.min(self.PotEvaporation, self.Sw[k]), 0)
    self.Qw1 = pcr.max(pcr.min(self.Fm2 * (self.Temp - self.Tm[k]), self.Sw[k]), 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew1 - self.Qw1

    self.Sw_diff = pcr.ifthenelse(self.Sw[k] < 0, self.Sw[k], 0)
    self.Ew = (
        self.Ew1
        + (self.Ew1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Qw = (
        self.Qw1
        + (self.Qw1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew - self.Qw
    self.Sw[k] = pcr.ifthenelse(self.Sw[k] < 0, 0, self.Sw[k])
    self.Sw_diff2 = pcr.ifthen(self.Sw[k] < 0, self.Sw[k])

    self.wbSw_[k] = (
        self.PrecipitationSnow - self.Ew - self.Qw - self.Sw[k] + self.Sw_t[k]
    )

    self.Ew_[k] = self.Ew
    self.Qw_[k] = self.Qw


def snow_rain_Tsurf_noEw(self, k):
    """
    - snow melt based on degree day factor and minimum surface temperature
    - meltfactor increases with temperature
    -
    - Code for ini-file: 5
    """
    JarvisCoefficients.calcEpSnow(self, k)
    # self.PotEvaporation = self.EpHour
    self.PotEvaporation = pcr.cover(pcr.ifthenelse(self.EpHour > 0, self.EpHour, 0), 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow

    self.Fm2 = pcr.max(self.Fm[k] * self.Precipitation, self.Fm[k])
    self.Ew1 = 0
    self.Qw1 = pcr.max(pcr.min(self.Fm2 * (self.TempSurf - self.Tm[k]), self.Sw[k]), 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew1 - self.Qw1

    self.Sw_diff = pcr.ifthenelse(self.Sw[k] < 0, self.Sw[k], 0)
    self.Ew = (
        self.Ew1
        + (self.Ew1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Qw = (
        self.Qw1
        + (self.Qw1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew - self.Qw
    self.Sw[k] = pcr.ifthenelse(self.Sw[k] < 0, 0, self.Sw[k])
    self.Sw_diff2 = pcr.ifthen(self.Sw[k] < 0, self.Sw[k])

    self.wbSw_[k] = (
        self.PrecipitationSnow - self.Ew - self.Qw - self.Sw[k] + self.Sw_t[k]
    )

    self.Ew_[k] = self.Ew
    self.Qw_[k] = self.Qw


def snowHour(self, k):
    """
    - snow melt based on degree day factor and
    - for hourly input  data
    - Code for ini-file: 2
    """
    #    JarvisCoefficients.calcEpSnowHour(self,k)
    #    self.PotEvaporation = self.EpHour
    #    self.PotEvaporation = pcr.ifthenelse(self.EpHour > 0, self.EpHour, 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow

    self.Ew1 = pcr.max(pcr.min(self.PotEvaporation, self.Sw[k]), 0)
    #    self.Ew1 = 0
    self.Qw1 = pcr.max(self.Fm[k] * (self.Temperature - self.Tm[k]), 0)

    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew1 - self.Qw1

    self.Sw_diff = pcr.ifthenelse(self.Sw[k] < 0, self.Sw[k], 0)
    self.Ew = (
        self.Ew1
        + (self.Ew1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Qw = (
        self.Qw1
        + (self.Qw1 / pcr.ifthenelse(self.Ew1 + self.Qw1 > 0, self.Ew1 + self.Qw1, 1))
        * self.Sw_diff
    )
    self.Sw[k] = self.Sw_t[k] + self.PrecipitationSnow - self.Ew - self.Qw
    self.Sw[k] = pcr.ifthenelse(self.Sw[k] < 0, 0, self.Sw[k])
    self.Sw_diff2 = pcr.ifthen(self.Sw[k] < 0, self.Sw[k])

    #    if any(pcr.pcr2numpy(self.Sw[k],np.nan) > 0):
    #        pdb.set_trace()
    self.wbSw_[k] = (
        self.PrecipitationSnow - self.Ew - self.Qw - self.Sw[k] + self.Sw_t[k]
    )

    self.Ew_[k] = self.Ew
    self.Qw_[k] = self.Qw


#    if any(pcr.pcr2numpy(self.Qw,np.nan) > 0):
#        pdb.set_trace()
