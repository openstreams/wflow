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

It contains both the kinematic wave, interception, snow/glaciers and
reservoirs/lakes modules.

"""

from numba import jit
import math
import numpy as np
import pcraster as pcr
from numba.errors import NumbaPendingDeprecationWarning
import warnings

warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)


# ldd definitie
_ldd = np.array([[7, 8, 9], [4, 5, 6], [1, 2, 3]])
_ldd_us = np.fliplr(np.flipud(_ldd))
_pits = 5

@jit(nopython=True)
def _up_nb(ldd_f, idx0, shape, _ldd_us=_ldd_us):
    """returns a numpy array with 1d indices of upstream neighbors on a ldd
    """
    nrow, ncol = shape
    r = idx0 // ncol
    c = idx0 % ncol
    wdw_idx = list()
    for dr in range(-1, 2):
        row = r + dr
        for dc in range(-1, 2):
            col = c + dc
            if dr == 0 and dc == 0: # skip pit -> return empty array
                continue
            elif row < 0 or row >= nrow or col < 0 or col >= ncol: # out of bounds
                pass
            else:
                idx = idx0 + dc + dr*ncol
                if ldd_f[idx] == _ldd_us[dr+1, dc+1]:
                    wdw_idx.append(idx)
    return np.array(wdw_idx, dtype=np.int64)

@jit(nopython=True)
def set_dd(ldd, _ldd_us=_ldd_us, pit_value=_pits):
    """set drainage direction network from downstream to upstream
    """
    shape = ldd.shape
    ldd = ldd.ravel()
    nodes = list()
    nodes_up = list()

    # most downstream indices
    idx_ds = np.where(ldd==np.array(pit_value).astype(ldd.dtype))[0].astype(np.int64)
           
    # move upstream
    while True:
        nodes.append(idx_ds)
        idx_next = list()
        nbs_up = np.ones((idx_ds.size, 8), dtype=np.int64)*-1
        for i, idx in enumerate(idx_ds):
            idx_up = _up_nb(ldd, idx, shape, _ldd_us)
            if np.any(idx_up):
                idx_next.extend(idx_up)
                nbs_up[i, :idx_up.size] = idx_up
        nodes_up.append(nbs_up)
        if len(idx_next) == 0:
            break
        idx_ds = np.array(idx_next, dtype=np.int64)
    return nodes[::-1], nodes_up[::-1]


###############################################################################
# Kinematic wave modules                                                      #
###############################################################################

def estimate_iterations_kin_wave(Q, Beta, alpha, timestepsecs, dx, mv):
    
    celerity = pcr.ifthen(Q > 0.0, 1.0 / (alpha * Beta * Q**(Beta-1)))
    courant = (timestepsecs / dx) * celerity
    np_courant = pcr.pcr2numpy(courant, mv)
    np_courant[np_courant==mv] = np.nan
    try:
        it_kin = int(np.ceil(1.25*(np.nanpercentile(np_courant,95))))
    except:
        it_kin = 1
    
    return it_kin

@jit(nopython=True)
def kinematic_wave(Qin,Qold,q,alpha,beta,deltaT,deltaX):
    
    epsilon = 1e-12
    MAX_ITERS = 3000

    if ((Qin+Qold+q) == 0.):
        return 0.
    else:
        #common terms
        ab_pQ = alpha*beta*pow(((Qold+Qin)/2.),beta-1.)
        deltaTX = deltaT/deltaX
        C = deltaTX*Qin + alpha*pow(Qold,beta) + deltaT*q
        
        Qkx   = (deltaTX * Qin + Qold * ab_pQ + deltaT * q) / (deltaTX + ab_pQ)
        
        if math.isnan(Qkx):
            Qkx = 0.
        
        Qkx   = max(Qkx, 1e-30)
        fQkx  = deltaTX * Qkx + alpha * pow(Qkx, beta) - C
        dfQkx = deltaTX + alpha * beta * pow(Qkx, beta - 1.)
        Qkx   = Qkx - fQkx / dfQkx
        Qkx   = max(Qkx, 1e-30)
        count = 0
        
        while abs(fQkx) > epsilon and count < MAX_ITERS:
            fQkx  = deltaTX * Qkx + alpha * pow(Qkx, beta) - C
            dfQkx = deltaTX + alpha * beta * pow(Qkx, beta - 1.)
            Qkx  =  Qkx - fQkx / dfQkx
            Qkx   = max(Qkx, 1e-30)
            count = count + 1
          
        return Qkx


@jit(nopython=True)
def kin_wave(rnodes, rnodes_up, Qold, q, Alpha, Beta, DCL, River, Bw, AlpTermR, AlpPow, deltaT, it=1):
    
    acc_flow = np.zeros(Qold.size, dtype=np.float64)
    acc_flow = np.concatenate((acc_flow, np.array([0], dtype=np.float64)))


    for v in range(0,it):
        shape = Qold.shape
        # flat new state
        Qnew = np.zeros(Qold.size, dtype=np.float64)
        # append zero to end to deal with nodata (-1) in indices
        Qnew = np.concatenate((Qnew, np.array([0], dtype=np.float64)))

        for i in range(len(rnodes)):
            for j in range(len(rnodes[i])):
                idx = rnodes[i][j]
                nbs = rnodes_up[i][j]
    
                Qin = np.sum(Qnew[nbs])
                Qnew[idx] = kinematic_wave(Qin, Qold[idx], q[idx], Alpha[idx], Beta[idx], deltaT/it, DCL[idx])
        
                acc_flow[idx] = acc_flow[idx] + Qnew[idx] * (deltaT/it)
                WaterLevelR = (Alpha[idx] * np.power(Qnew[idx], Beta[idx])) / Bw[idx]
                Pr = Bw[idx] + (2.0 * WaterLevelR)
                Alpha[idx] = AlpTermR[idx] * np.power(Pr, AlpPow[idx])
                Qold[idx]= Qnew[idx]
    # remove last value from array and reshape to original format
    return acc_flow[:-1].reshape(shape)
    #return Qnew[:-1].reshape(shape)


@jit(nopython=True)
def kinematic_wave_ssf(ssf_in, ssf_old, zi_old, r, Ks_hor, Ks, slope, neff, f, D, dt, dx, w, ssf_max):
    
    epsilon = 1e-3
    MAX_ITERS = 3000
        
    if ((ssf_in+ssf_old) == 0. and (r <= 0)):
    #if (max(ssf_in+ssf_old+r,0.) == 0.):
        return 0., D, 0.
    else:
        #initial estimate
        ssf_n = (ssf_old + ssf_in)/2.0
        count = 0
        
        # Estimate zi on the basis of the relation between subsurfacel flow and zi           
        zi = math.log((f*ssf_n)/(w*Ks_hor*Ks*slope) + math.exp(-f*D))/-f
        # Reciprocal of derivative delta Q/ delta z_i, constrained w.r.t. neff on the basis of the continuity equation)
        Cn = (Ks_hor*Ks*np.exp(-f*zi)*slope)/neff
        # Term of the continuity equation for Newton-Raphson iteration for iteration 1
        # because celerity Cn is depending on zi, the increase or decrease of zi is moved to the recharge term of the continuity equation
        # then (1./Cn)*ssf_old can be replaced with (1./Cn)*ssf_n, and thus celerity and lateral flow rate ssf_n are then in line 
        c = (dt/dx)*ssf_in + (1./Cn)*ssf_n + dt*(r-(zi_old-zi)*neff*w)
   

        # Continuity equation of which solution should be zero
        fQ = (dt/dx)*ssf_n + (1./Cn)*ssf_n - c
        # Derivative of the continuity equation w.r.t. Q_out for iteration 1
        dfQ = (dt/dx) +  1/Cn
        # Update lateral outflow estimate ssf_n (Q_out) for iteration 1
        ssf_n =  ssf_n - (fQ/dfQ)

        # Check whether ssf_n is reasonable
        if math.isnan(ssf_n):
            ssf_n = 0.
        ssf_n   = max(ssf_n, 1e-30)
        
        # Start while loop of Newton-Raphson iteration m until continuity equation approaches zero        
        while abs(fQ) > epsilon and count < MAX_ITERS:
            # Estimate zi on the basis of the relation between lateral flow rate and groundwater level           
            zi = math.log((f*ssf_n)/(w*Ks_hor*Ks*slope) + math.exp(-f*D))/-f
            # Reciprocal of derivative delta Q/ delta z_i, constrained w.r.t. neff on the basis of the continuity equation
            Cn = (Ks_hor*Ks*np.exp(-f*zi)*slope)/neff

            # Term of the continuity equation for given Newton-Raphson iteration m
            # because celerity Cn is depending on zi, the increase or decrease of zi is moved to the recharge term of the continuity equation
            # then (1./Cn)*ssf_old can be replaced with (1./Cn)*ssf_n, and thus celerity and lateral flow rate ssf_n are then in line 
            c = (dt/dx)*ssf_in + (1./Cn)*ssf_n + dt*(r-(zi_old-zi)*neff*w)
    
            # Continuity equation of which solution should be zero
            fQ = (dt/dx)*ssf_n + (1./Cn)*ssf_n - c
            # Derivative of the continuity equation w.r.t. Q_out for iteration m+1
            dfQ = (dt/dx) +  1./Cn         
            # Update lateral outflow estimate ssf_n (Q_out) for iteration m+1
            ssf_n =  ssf_n - (fQ/dfQ)
            
            # Check whether ssf_n is reasonable
            if math.isnan(ssf_n):
                ssf_n = 0.
            ssf_n   = max(ssf_n, 1e-30)
            # Upload count
            count = count + 1

        # Constrain the lateral flow rate ssf_n
        ssf_n =  min(ssf_n,(ssf_max*w))
        # On the basis of the lateral flow rate, estimate the amount of groundwater level above surface (saturation excess conditions), then rest = negative
        rest = zi_old - (ssf_in + r*dx - ssf_n)/(w*dx)/neff
        # In case the groundwater level lies above surface (saturation excess conditions, rest = negative), calculate the exfiltration rate and set groundwater back to zero.     
        exfilt = min(rest,0.0) * -neff
        zi = max(0,zi)
        
        return ssf_n, zi, exfilt   


###############################################################################
# Rainfall interception modules                                               #
###############################################################################

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


@jit(nopython=True)
def propagate_downstream(rnodes, rnodes_up, material):
    shape = material.shape
    material = material.flatten()
    for i in range(len(rnodes)):
        for j in range(len(rnodes[i])):
            idx_ds = rnodes[i][j]
            idxs_us = rnodes_up[i][j] # NOTE: has nodata (-1) values
            v = 0
            for idx_us in idxs_us:
                if idx_us == -1: break
                v += material[idx_us]
            material[idx_ds] += v
    return material.reshape(shape)



###############################################################################
# Snow and glaciers modules                                                   #
###############################################################################

def SnowPackHBV(Snow, SnowWater, Precipitation, Temperature, TTI, TT, TTM, Cfmax, WHC):
    """
    HBV Type snowpack modelling using a Temperature degree factor. All correction
    factors (RFCF and SFCF) are set to 1. The refreezing efficiency factor is set to 0.05.

    :param Snow:
    :param SnowWater:
    :param Precipitation:
    :param Temperature:
    :param TTI:
    :param TT:
    :param TTM:
    :param Cfmax:
    :param WHC:
    :return: Snow,SnowMelt,Precipitation
    """

    RFCF = 1.0  # correction factor for rainfall
    CFR = 0.05000  # refreeing efficiency constant in refreezing of freewater in snow
    SFCF = 1.0  # correction factor for snowfall

    RainFrac = pcr.ifthenelse(
        1.0 * TTI == 0.0,
        pcr.ifthenelse(Temperature <= TT, pcr.scalar(0.0), pcr.scalar(1.0)),
        pcr.min((Temperature - (TT - TTI / 2)) / TTI, pcr.scalar(1.0)),
    )
    RainFrac = pcr.max(
        RainFrac, pcr.scalar(0.0)
    )  # fraction of precipitation which falls as rain
    SnowFrac = 1 - RainFrac  # fraction of precipitation which falls as snow
    Precipitation = (
        SFCF * SnowFrac * Precipitation + RFCF * RainFrac * Precipitation
    )  # different correction for rainfall and snowfall

    SnowFall = SnowFrac * Precipitation  # snowfall depth
    RainFall = RainFrac * Precipitation  # rainfall depth
    PotSnowMelt = pcr.ifthenelse(
        Temperature > TTM, Cfmax * (Temperature - TTM), pcr.scalar(0.0)
    )  # Potential snow melt, based on temperature
    PotRefreezing = pcr.ifthenelse(
        Temperature < TTM, Cfmax * CFR * (TTM - Temperature), 0.0
    )  # Potential refreezing, based on temperature
    Refreezing = pcr.ifthenelse(
        Temperature < TTM, pcr.min(PotRefreezing, SnowWater), 0.0
    )  # actual refreezing
    # No landuse correction here
    SnowMelt = pcr.min(PotSnowMelt, Snow)  # actual snow melt
    Snow = Snow + SnowFall + Refreezing - SnowMelt  # dry snow content
    SnowWater = SnowWater - Refreezing  # free water content in snow
    MaxSnowWater = Snow * WHC  # Max water in the snow
    SnowWater = (
        SnowWater + SnowMelt + RainFall
    )  # Add all water and potentially supersaturate the snowpack
    RainFall = pcr.max(SnowWater - MaxSnowWater, 0.0)  # rain + surpluss snowwater
    SnowWater = SnowWater - RainFall

    return Snow, SnowWater, SnowMelt, RainFall, SnowFall


def glacierHBV(GlacierFrac, 
                GlacierStore, 
                Snow, 
                Temperature, 
                TT, 
                Cfmax, 
                G_SIfrac,
                timestepsecs,
                basetimestep):
    """
    Run Glacier module and add the snowpack on-top of it.
    First, a fraction of the snowpack is converted into ice using the HBV-light
    model (fraction between 0.001-0.005 per day).
    Glacier melting is modelled using a Temperature degree factor and only
    occurs if the snow cover < 10 mm.


    :ivar GlacierFrac: Fraction of wflow cell covered by glaciers
    :ivar GlacierStore: Volume of the galcier in the cell in mm w.e.
    :ivar Snow: Snow pack on top of Glacier
    :ivar Temperature: Air temperature
    :ivar TT: Temperature threshold for ice melting
    :ivar Cfmax: Ice degree-day factor in mm/(°C/day)
    :ivar G_SIfrac: Fraction of the snow part turned into ice each timestep
    :ivar timestepsecs: Model timestep in seconds
    :ivar basetimestep: Model base timestep (86 400 seconds)

    :returns: Snow,Snow2Glacier,GlacierStore,GlacierMelt,
    """
    
    #Fraction of the snow transformed into ice (HBV-light model)
    Snow2Glacier = G_SIfrac * Snow

    Snow2Glacier = pcr.ifthenelse(
        GlacierFrac > 0.0, Snow2Glacier, pcr.scalar(0.0)
    )
    # Max conversion to 8mm/day
    Snow2Glacier = (
        pcr.min(Snow2Glacier, 8.0 * (timestepsecs / basetimestep))
    )

    Snow = Snow - (Snow2Glacier * GlacierFrac)
    GlacierStore = GlacierStore + Snow2Glacier

    PotMelt = pcr.ifthenelse(
        Temperature > TT, Cfmax * (Temperature - TT), pcr.scalar(0.0)
    )  # Potential snow melt, based on temperature

    GlacierMelt = pcr.ifthenelse(
        Snow < 10.0, pcr.min(PotMelt, GlacierStore), pcr.scalar(0.0)
    )  # actual Glacier melt
    GlacierMelt = pcr.cover(GlacierMelt, pcr.scalar(0.0))
    GlacierStore = GlacierStore - GlacierMelt  # dry snow content

    return Snow, Snow2Glacier, GlacierStore, GlacierMelt



###############################################################################
# Lake and reservoirmodules                                                   #
###############################################################################
    
def sCurve(X, a=0.0, b=1.0, c=1.0):
    """
    sCurve function:

    Input:
        - X input map
        - C determines the steepness or "stepwiseness" of the curve.
          The higher C the sharper the function. A negative C reverses the function.
        - b determines the amplitude of the curve
        - a determines the centre level (default = 0)

    Output:
        - result
    """
    try:
        s = 1.0 / (b + pcr.exp(-c * (X - a)))
    except:
        s = 1.0 / (b + pcr.exp(-c * (X - a)))
    return s

    
def lookupResRegMatr(ReserVoirLocs, values, hq, JDOY):

    np_res_ids = pcr.pcr2numpy(ReserVoirLocs, 0)
    npvalues = pcr.pcr2numpy(values, 0)
    out = np.copy(npvalues) * 0.0

    if len(hq) > 0:
        for key in hq:
            value = npvalues[np.where(np_res_ids == key)]

            val = np.interp(value, hq[key][:, 0], hq[key][:, JDOY])

            out[np.where(np_res_ids == key)] = val

    return pcr.numpy2pcr(pcr.Scalar, out, 0)


def lookupResFunc(ReserVoirLocs, values, sh, dirLookup):

    np_res_ids = pcr.pcr2numpy(ReserVoirLocs, 0)
    npvalues = pcr.pcr2numpy(values, 0)
    out = np.copy(npvalues) * 0.0

    if len(sh) > 0:
        for key in sh:
            value = npvalues[np.where(np_res_ids == key)]

            if dirLookup == "0-1":
                val = np.interp(value, sh[key][:, 0], sh[key][:, 1])
            if dirLookup == "1-0":
                val = np.interp(value, sh[key][:, 1], sh[key][:, 0])

            out[np.where(np_res_ids == key)] = val

    return pcr.numpy2pcr(pcr.Scalar, out, 0)


def naturalLake(
    waterlevel,
    LakeLocs,
    LinkedLakeLocs,
    LakeArea,
    LakeThreshold,
    LakeStorFunc,
    LakeOutflowFunc,
    sh,
    hq,
    lake_b,
    lake_e,
    inflow,
    precip,
    pet,
    LakeAreasMap,
    JDOY,
    timestepsecs=86400,
):
    
    """
    Run Natural Lake module to compute the new waterlevel and outflow.
    Solves lake water balance with linearisation and iteration procedure,
    for any rating and storage curve.
    For the case where storage curve is S = AH and Q=b(H-Ho)^2, uses the direct
    solution from the Modified Puls Approach (LISFLOOD).


    :ivar waterlevel: water level H in the lake
    :ivar LakeLocs: location of lake's outlet
    :ivar LinkedLakeLocs: ID of linked lakes
    :ivar LakeArea: total lake area
    :ivar LakeThreshold: water level threshold Ho under which outflow is zero
    :ivar LakeStorFunc: type of lake storage curve
                        1: S = AH
                        2: S = f(H) from lake data and interpolation
    :ivar LakeOutflowFunc: type of lake rating curve
                           1: Q = f(H) from lake data and interpolation
                           2: General Q = b(H - Ho)^e
                           3: Case of Puls Approach Q = b(H - Ho)^2
    :ivar sh: data for storage curve
    :ivar hq: data for rating curve
    :ivar lake_b: rating curve coefficient
    :ivar lake_e: rating curve exponent
    :ivar inflow: inflow to the lake (surface runoff + river discharge + seepage)
    :ivar precip: precipitation map
    :ivar pet: PET map
    :ivar LakeAreasMap: lake extent map (for filtering P and PET)
    :ivar JDOY: Julian Day of Year to read storage/rating curve from data
    :ivar timestepsecs: model timestep in seconds

    :returns: waterlevel, outflow, prec_av, pet_av, storage
    """

    mv = -999.0
    LakeZeros = LakeArea * 0.0
    
    waterlevel_start = waterlevel

    inflow = pcr.ifthen(pcr.boolean(LakeLocs), inflow)

    prec_av = pcr.ifthen(
        pcr.boolean(LakeLocs), pcr.areaaverage(precip, LakeAreasMap)
    )
    pet_av = pcr.ifthen(
        pcr.boolean(LakeLocs), pcr.areaaverage(pet, LakeAreasMap)
    )
    
    
    ### Modified Puls Approach (Burek et al., 2013, LISFLOOD) ###
    #ResOutflowFunc = 3 
    
    #Calculate lake factor and SI parameter
    LakeFactor = pcr.ifthenelse(
            LakeOutflowFunc == 3,
            LakeArea / (timestepsecs * (lake_b) ** 0.5),
            mv
            )
    
    storage_start = pcr.ifthenelse(
            LakeStorFunc == 1,
            LakeArea * waterlevel_start,
            lookupResFunc(LakeLocs, waterlevel_start, sh, "0-1"),
            )

    SIFactor = pcr.ifthenelse(
            LakeOutflowFunc == 3,
            ((storage_start + (prec_av-pet_av)*LakeArea/1000.0) / timestepsecs 
             + inflow),
            mv
            )
    #Adjust SIFactor for ResThreshold != 0
    SIFactorAdj = SIFactor - LakeArea * LakeThreshold / timestepsecs
    
    #Calculate the new lake outflow/waterlevel/storage
    outflow = pcr.ifthenelse(
            LakeOutflowFunc == 3,
            pcr.ifthenelse(
                    SIFactorAdj > 0.0,
                    (-LakeFactor + (LakeFactor**2 + 2*SIFactorAdj) ** 0.5) ** 2,
                    0.0),
            LakeZeros
            )
    storage = pcr.ifthenelse(
            LakeOutflowFunc == 3,
            (SIFactor - outflow) * timestepsecs,
            LakeZeros
            )
    waterlevel = pcr.ifthenelse(
            LakeOutflowFunc == 3,
            storage / LakeArea,
            LakeZeros
            )
    
    ### Linearisation and iteration for specific storage/rating curves ###
    np_lakeoutflowfunc = pcr.pcr2numpy(LakeOutflowFunc, 0.0)
    if ((bool(np.isin(1, np.unique(np_lakeoutflowfunc)))) or 
        (bool(np.isin(2, np.unique(np_lakeoutflowfunc))))):
        
        np_lakelocs = pcr.pcr2numpy(LakeLocs, 0.0)
        np_linkedlakelocs = pcr.pcr2numpy(LinkedLakeLocs, 0.0)
        waterlevel_loop = waterlevel_start
    
        _outflow = []
        nr_loop = np.max([int(timestepsecs / 21600), 1])
        for n in range(0, nr_loop):
            np_waterlevel = pcr.pcr2numpy(waterlevel_loop, np.nan)
            np_waterlevel_lower = np_waterlevel.copy()
    
            for val in np.unique(np_linkedlakelocs):
                if val > 0:
                    np_waterlevel_lower[np_linkedlakelocs == val] = np_waterlevel[
                        np.where(np_lakelocs == val)
                    ]
    
            diff_wl = np_waterlevel - np_waterlevel_lower
            diff_wl[np.isnan(diff_wl)] = mv
            np_waterlevel_lower[np.isnan(np_waterlevel_lower)] = mv
    
            pcr_diff_wl = pcr.numpy2pcr(pcr.Scalar, diff_wl, mv)
            pcr_wl_lower = pcr.numpy2pcr(pcr.Scalar, np_waterlevel_lower, mv)
    
            storage_start_loop = pcr.ifthenelse(
                LakeStorFunc == 1,
                LakeArea * waterlevel_loop,
                lookupResFunc(LakeLocs, waterlevel_loop, sh, "0-1"),
            )
    
            outflow_loop = pcr.ifthenelse(
                LakeOutflowFunc == 1,
                lookupResRegMatr(LakeLocs, waterlevel_loop, hq, JDOY),
                pcr.ifthenelse(
                    pcr_diff_wl >= 0,
                    pcr.max(lake_b * (waterlevel_loop - LakeThreshold) ** lake_e, 0),
                    pcr.min(-1 * lake_b * (pcr_wl_lower - LakeThreshold) ** lake_e, 0),
                ),
            )
    
            np_outflow = pcr.pcr2numpy(outflow_loop, np.nan)
            np_outflow_linked = np_lakelocs * 0.0
    
            with np.errstate(invalid="ignore"):
                if np_outflow[np_outflow < 0] is not None:
                    np_outflow_linked[
                        np.in1d(np_lakelocs, np_linkedlakelocs[np_outflow < 0]).reshape(
                            np_linkedlakelocs.shape
                        )
                    ] = np_outflow[np_outflow < 0]
    
            outflow_linked = pcr.numpy2pcr(pcr.Scalar, np_outflow_linked, 0.0)
    
            fl_nr_loop = float(nr_loop)
            storage_loop = (
                storage_start_loop
                + (inflow * timestepsecs / fl_nr_loop)
                + (prec_av / fl_nr_loop / 1000.0) * LakeArea
                - (pet_av / fl_nr_loop / 1000.0) * LakeArea
                - (pcr.cover(outflow_loop, 0.0) * timestepsecs / fl_nr_loop)
                + (pcr.cover(outflow_linked, 0.0) * timestepsecs / fl_nr_loop)
            )
    
            waterlevel_loop = pcr.ifthenelse(
                LakeStorFunc == 1,
                waterlevel_loop + (storage_loop - storage_start_loop) / LakeArea,
                lookupResFunc(LakeLocs, storage_loop, sh, "1-0"),
            )
    
            np_outflow_nz = np_outflow * 0.0
            with np.errstate(invalid="ignore"):
                np_outflow_nz[np_outflow > 0] = np_outflow[np_outflow > 0]
            _outflow.append(np_outflow_nz)
    
        outflow_av_temp = np.average(_outflow, 0)
        outflow_av_temp[np.isnan(outflow_av_temp)] = mv
        outflow_av = pcr.numpy2pcr(pcr.Scalar, outflow_av_temp, mv)
        
        #Add the discharge/waterlevel/storage from the loop to the one from puls approach
        outflow = pcr.ifthenelse(
                LakeOutflowFunc == 3,
                outflow,
                outflow_av
                )
        waterlevel = pcr.ifthenelse(
                LakeOutflowFunc == 3,
                waterlevel,
                waterlevel_loop
                )
        storage = pcr.ifthenelse(
                LakeOutflowFunc == 3,
                storage,
                storage_loop
                )

    return waterlevel, outflow, prec_av, pet_av, storage


def simplereservoir(
    storage,
    inflow,
    ResArea,
    maxstorage,
    target_perc_full,
    maximum_Q,
    demand,
    minimum_full_perc,
    ReserVoirLocs,
    precip,
    pet,
    ReservoirSimpleAreas,
    timestepsecs=86400,
):
    """

    :param storage: initial storage m^3
    :param inflow: inflow m^3/s
    :param maxstorage: maximum storage (above which water is spilled) m^3
    :param target_perc_full: target fraction full (of max storage) -
    :param maximum_Q: maximum Q to release m^3/s if below spillway
    :param demand: water demand (all combined) m^3/s
    :param minimum_full_perc: target minimum full fraction (of max storage) -
    :param ReserVoirLocs: map with reservoir locations
    :param timestepsecs: timestep of the model in seconds (default = 86400)
    :return: storage (m^3), outflow (m^3/s), PercentageFull (0-1), Release (m^3/sec)
    """

    inflow = pcr.ifthen(pcr.boolean(ReserVoirLocs), inflow)

    prec_av = pcr.cover(
        pcr.ifthen(
            pcr.boolean(ReserVoirLocs), pcr.areaaverage(precip, ReservoirSimpleAreas)
        ),
        pcr.scalar(0.0),
    )
    pet_av = pcr.cover(
        pcr.ifthen(
            pcr.boolean(ReserVoirLocs), pcr.areaaverage(pet, ReservoirSimpleAreas)
        ),
        pcr.scalar(0.0),
    )
    
    _outflow = 0
    _demandRelease = 0
    
    nr_loop = np.max([int(timestepsecs / 21600), 1])
    for n in range(0, nr_loop):
        
        fl_nr_loop = float(nr_loop)
        
        storage = (
            storage
            + (inflow * timestepsecs / fl_nr_loop)
            + (prec_av / fl_nr_loop / 1000.0) * ResArea
            - (pet_av / fl_nr_loop / 1000.0) * ResArea
        )

        percfull = storage / maxstorage
        # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
        fac = sCurve(percfull, a=minimum_full_perc, c=30.0)
        demandRelease = pcr.min(fac * demand * timestepsecs / fl_nr_loop, storage)
        storage = storage - demandRelease

        wantrel = pcr.max(0.0, storage - (maxstorage * target_perc_full))
        # Assume extra maximum Q if spilling
        overflowQ = pcr.max((storage - maxstorage), 0.0)
        torelease = pcr.min(wantrel, overflowQ + maximum_Q * timestepsecs / fl_nr_loop - demandRelease)
        storage = storage - torelease
        outflow = torelease + demandRelease
        percfull = storage / maxstorage
        
        _outflow = _outflow + outflow
        _demandRelease = _demandRelease + demandRelease

    return storage, _outflow / timestepsecs, percfull, prec_av, pet_av, _demandRelease / timestepsecs