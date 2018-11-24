# Copyright (c) J. Schellekens 2005-2011
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
wflow_funcs -  hydrological functions library
---------------------------------------------

In addition this library contain a number of hydrological functions
that may be used within the wflow models

"""

import pcraster as pcr


def rainfall_interception_hbv(Rainfall, PotEvaporation, Cmax, InterceptionStorage):
    """
    Returns:
    TF, Interception, IntEvap,InterceptionStorage
    """
    Interception = pcr.min(
        Rainfall, Cmax - InterceptionStorage
    )  #: Interception in mm/timestep

    InterceptionStorage = (
        InterceptionStorage + Interception
    )  #: Current interception storage
    TF = Rainfall - Interception
    IntEvap = pcr.min(
        InterceptionStorage, PotEvaporation
    )  #: Evaporation from interception storage
    InterceptionStorage = InterceptionStorage - IntEvap

    return TF, Interception, IntEvap, InterceptionStorage


def rainfall_interception_gash(
    Cmax, EoverR, CanopyGapFraction, Precipitation, CanopyStorage, maxevap=9999
):
    """
    Interception according to the Gash model (For daily timesteps). 
    """
    # TODO:  add other rainfall interception method (lui)
    # TODO: Include subdaily Gash model
    # TODO: add LAI variation in year
    # Hack for stemflow

    pt = 0.1 * CanopyGapFraction

    P_sat = pcr.max(
        pcr.scalar(0.0),
        pcr.cover(
            (-Cmax / EoverR) * pcr.ln(1.0 - (EoverR / (1.0 - CanopyGapFraction - pt))),
            pcr.scalar(0.0),
        ),
    )

    # large storms P > P_sat
    largestorms = Precipitation > P_sat

    Iwet = pcr.ifthenelse(
        largestorms,
        ((1 - CanopyGapFraction - pt) * P_sat) - Cmax,
        Precipitation * (1 - CanopyGapFraction - pt),
    )
    Isat = pcr.ifthenelse(largestorms, (EoverR) * (Precipitation - P_sat), 0.0)
    Idry = pcr.ifthenelse(largestorms, Cmax, 0.0)
    Itrunc = 0

    StemFlow = pt * Precipitation

    ThroughFall = Precipitation - Iwet - Idry - Isat - Itrunc - StemFlow
    Interception = Iwet + Idry + Isat + Itrunc

    # Non corect for area without any Interception (say open water Cmax -- zero)
    CmaxZero = Cmax <= 0.0
    ThroughFall = pcr.ifthenelse(CmaxZero, Precipitation, ThroughFall)
    Interception = pcr.ifthenelse(CmaxZero, pcr.scalar(0.0), Interception)
    StemFlow = pcr.ifthenelse(CmaxZero, pcr.scalar(0.0), StemFlow)

    # Now corect for maximum potential evap
    OverEstimate = pcr.ifthenelse(
        Interception > maxevap, Interception - maxevap, pcr.scalar(0.0)
    )
    Interception = pcr.min(Interception, maxevap)
    # Add surpluss to the thoughdfall
    ThroughFall = ThroughFall + OverEstimate

    return ThroughFall, Interception, StemFlow, CanopyStorage


def rainfall_interception_modrut(
    Precipitation, PotEvap, CanopyStorage, CanopyGapFraction, Cmax
):
    """
    Interception according to a modified Rutter model. The model is solved
    explicitly and there is no drainage below Cmax.
    
    Returns:
        - NetInterception: P - TF - SF (may be different from the actual wet canopy evaporation)
        - ThroughFall:
        - StemFlow:
        - LeftOver: Amount of potential eveporation not used
        - Interception: Actual wet canopy evaporation in this thimestep
        - CanopyStorage: Canopy storage at the end of the timestep
    
    """

    ##########################################################################
    # Interception according to a modified Rutter model with hourly timesteps#
    ##########################################################################

    p = CanopyGapFraction
    pt = 0.1 * p

    # Amount of P that falls on the canopy
    Pfrac = pcr.max((1 - p - pt), 0) * Precipitation

    # S cannot be larger than Cmax, no gravity drainage below that
    DD = pcr.ifthenelse(CanopyStorage > Cmax, CanopyStorage - Cmax, 0.0)
    CanopyStorage = CanopyStorage - DD

    # Add the precipitation that falls on the canopy to the store
    CanopyStorage = CanopyStorage + Pfrac

    # Now do the Evap, make sure the store does not get negative
    dC = -1 * pcr.min(CanopyStorage, PotEvap)
    CanopyStorage = CanopyStorage + dC

    LeftOver = PotEvap + dC
    # Amount of evap not used

    # Now drain the canopy storage again if needed...
    D = pcr.ifthenelse(CanopyStorage > Cmax, CanopyStorage - Cmax, 0.0)
    CanopyStorage = CanopyStorage - D

    # Calculate throughfall
    ThroughFall = DD + D + p * Precipitation
    StemFlow = Precipitation * pt

    # Calculate interception, this is NET Interception
    NetInterception = Precipitation - ThroughFall - StemFlow
    Interception = -dC

    return NetInterception, ThroughFall, StemFlow, LeftOver, Interception, CanopyStorage


# baseflow seperation methods
# see http://mssanz.org.au/MODSIM97/Vol%201/Chapman.pdf


def bf_oneparam(discharge, k):
    bf = list(range(0, len(discharge)))
    for i in range(1, len(discharge)):
        bf[i] = (k * bf[i - 1] / (2.0 - k)) + ((1.0 - k) * discharge[i] / (2.0 - k))
        if bf[i] > discharge[i]:
            bf[i] = discharge[i]

    return bf


def bf_twoparam(discharge, k, C):
    bf = list(range(0, len(discharge)))
    for i in range(1, len(discharge)):
        bf[i] = (k * bf[i - 1] / (1.0 + C)) + ((C) * discharge[i] / (1.0 + C))
        if bf[i] > discharge[i]:
            bf[i] = discharge[i]

    return bf


def bf_threeparam(discharge, k, C, a):
    bf = list(range(0, len(discharge)))
    for i in range(1, len(discharge)):
        bf[i] = (k * bf[i - 1] / (1.0 + C)) + (
            (C) * discharge[i] + a * discharge[i - 1] / (1.0 + C)
        )
        if bf[i] > discharge[i]:
            bf[i] = discharge[i]

    return bf
