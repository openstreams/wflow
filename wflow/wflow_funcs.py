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

from numba import jit
import math
import numpy as np
import pcraster as pcr


@jit(nopython=True)
def _up_nb(ldd_f, idx0, shape, _ldd_us, sr=1):
    """returns a numpy array with 1d indices of upstream neighbors on a ldd
    """
    nrow, ncol = shape
    r = idx0 // ncol
    c = idx0 % ncol
    wdw_idx = list()
    i = 0
    for dr in range(-sr, sr+1):
        row = r + dr
        if row >= 0 and row < nrow:
            for dc in range(-sr, sr+1):
                col = c + dc
                if col >= 0 and col < ncol:
                    idx = idx0 + dc + dr*ncol
                    if ldd_f[idx] == _ldd_us[i]:
                        wdw_idx.append(idx)
                i += 1
        else:
            i += sr*2+1
    return np.array(wdw_idx, dtype=np.int32)


@jit(nopython=True)
def set_dd(ldd, _ldd_us, river, pit_value=5):
    """set drainage direction network from downstream to upstream
    """
    shape = ldd.shape
    ldd = ldd.flatten()
    river = np.concatenate((river, np.array([0], dtype=river.dtype)))
    nodes = list()
    nodes_up = list()
    rnodes = list()
    rnodes_up = list()
    
    idx_ds_ = np.where(ldd==np.array(pit_value).astype(ldd.dtype))[0].astype(np.int32)
    
    for c, idx_ in enumerate(idx_ds_):
        idx_ds = np.array([idx_])
        
        # move upstream
        while True:
            nodes.append(idx_ds)
            idx_r_ds = idx_ds[np.where(river[idx_ds])]
            if idx_r_ds.size > 0:
                rnodes.append(idx_r_ds)
                r_nbs_up = np.ones((idx_r_ds.size, 8), dtype=np.int32)*-1
            idx_next = list()
            nbs_up = np.ones((idx_ds.size, 8), dtype=np.int32)*-1
            j = 0
            for i, idx in enumerate(idx_ds):
                idx_up = _up_nb(ldd, idx, shape, _ldd_us)
                if np.any(idx_up):
                    idx_next.extend(idx_up)
                    nbs_up[i, :idx_up.size] = idx_up
                    if river[idx]:
                        idx_r_up = idx_up[np.where(river[idx_up])]
                        r_nbs_up[j, :idx_r_up.size] = idx_r_up
                        j = j + 1
            nodes_up.append(nbs_up)
            if idx_r_ds.size > 0:
                rnodes_up.append(r_nbs_up)
            if len(idx_next) == 0:
                break
            idx_ds = np.array(idx_next, dtype=np.int32)
    return nodes[::-1], nodes_up[::-1], rnodes[::-1], rnodes_up[::-1]


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
def kinematic_wave_ssf(ssf_in, ssf_old, zi_old, r, Ks_hor, Ks ,slope ,neff, f, D, dt, dx, w, ssf_max):
    
    epsilon = 1e-6
    MAX_ITERS = 3000
        
    if (max(ssf_in+ssf_old+r,0.) == 0.):
        return 0., D, 0.
    else:
        #initial estimate
        ssf_n = (ssf_in + ssf_old)/2.
        count = 0        
        
        zi = np.log(f*ssf_n/(w*Ks_hor*Ks*slope) + np.exp(-f*D))/-f
        Cn = (Ks_hor*Ks*slope)/neff * np.exp(-f*zi)     
        c = (dt/dx)*ssf_in + 1./Cn*ssf_n + dt*(r-(zi_old-zi)*neff*w)
    
        fQ = (dt/dx)*ssf_n + 1./Cn*ssf_n - c      
        dfQ = (dt/dx) + 1./Cn        
        ssf_n =  ssf_n - (fQ/dfQ)
        
        if math.isnan(ssf_n):
            ssf_n = 0.
        ssf_n   = max(ssf_n, 1e-30)
                
        while abs(fQ) > epsilon and count < MAX_ITERS:
            
            zi = np.log(f*ssf_n/(w*Ks_hor*Ks*slope) + np.exp(-f*D))/-f
            Cn = (Ks_hor*Ks*slope)/neff * np.exp(-f*zi) 
            c = (dt/dx)*ssf_in + 1./Cn*ssf_n + dt*(r-(zi_old-zi)*neff*w)
            
            fQ = (dt/dx)*ssf_n + 1./Cn*ssf_n - c
            dfQ = (dt/dx) + 1./Cn          
            ssf_n =  ssf_n - (fQ/dfQ)
            
            if math.isnan(ssf_n):
                ssf_n = 0.
            ssf_n   = max(ssf_n, 1e-30)

            count = count + 1
        
        ssf_n =  min(ssf_n,(ssf_max*w))
        #exfilt = min(0,zi) * -neff
        exfilt = min(zi_old - (ssf_in + r*dx - ssf_n)/(w*dx)/neff,0.0) * -neff
        zi = max(0,zi)
        
        return ssf_n, zi, exfilt   


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