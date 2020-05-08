#!/usr/bin/python

# Wflow is Free software, see below:
#
# Copyright (c) J. Schellekens/Deltares 2005-2014
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
Run the wflow_sbm hydrological model..

usage

::

    wflow_sbm [-h][-v level][-L logfile][-C casename][-R runId]
          [-c configfile][-T last_step][-S first_step][-s seconds]
          [-P parameter multiplication][-X][-f][-I][-i tbl_dir][-x subcatchId]
          [-p inputparameter multiplication][-l loglevel]


    -X: save state at the end of the run over the initial conditions at the start

    -f: Force overwrite of existing results

    -T: Set end time of the run: yyyy-mm-dd hh:mm:ss

    -S: Set start time of the run: yyyy-mm-dd hh:mm:ss

    -s: Set the model timesteps in seconds

    -I: re-initialize the initial model conditions with default

    -i: Set input table directory (default is intbl)

    -x: Apply multipliers (-P/-p ) for subcatchment only (e.g. -x 1)

    -C: set the name  of the case (directory) to run

    -R: set the name runId within the current case

    -L: set the logfile

    -c: name of wflow configuration file (default: Casename/wflow_sbm.ini).

    -h: print usage information

    -P: set parameter change string (e.g: -P "self.FC = self.FC * 1.6") for non-dynamic variables

    -p: set parameter change string (e.g: -P "self.Precipitation = self.Precipitation * 1.11") for
        dynamic variables

    -l: loglevel (most be one of DEBUG, WARNING, ERROR)


"""

import os.path

import numpy as np
from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *
from wflow.wflow_funcs import *
import pcraster.framework
import pcraster as pcr

import pdb
import math
from numba import jit

from numba.errors import NumbaPendingDeprecationWarning
import warnings

warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

wflow = "wflow_sbm: "

updateCols = []

def usage(*args):
    sys.stdout = sys.stderr
    """Way"""
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)



@jit(nopython=True)
def _sCurve(X, a=0.0, b=1.0, c=1.0):
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
    s = 1.0 / (b + np.exp(-c * (X - a)))

    return s


@jit(nopython=True)
def actTransp_unsat_SBM(
    RootingDepth,
    UStoreLayerDepth,
    sumLayer,
    RestPotEvap,
    sumActEvapUStore,
    c,
    L,
    thetaS,
    thetaR,
    hb,
    ust=0):
    
    """
    Actual transpiration function for unsaturated zone:

      if ust is True, all ustore is available for transpiration

    Input:

        - RootingDepth, UStoreLayerDepth, sumLayer (depth (z) of upper boundary unsaturated layer),
          RestPotEvap (remaining evaporation), sumActEvapUStore (cumulative actual transpiration (more than one UStore layers))
          c (Brooks-Corey coefficient), L (thickness of unsaturated zone), thetaS, thetaR, hb (air entry pressure), ust

    Output:

        - UStoreLayerDepth,  sumActEvapUStore, ActEvapUStore
    """

    # AvailCap is fraction of unsat zone containing roots
    if ust >= 1:
        AvailCap = UStoreLayerDepth * 0.99
    else:
        if L > 0:
            AvailCap = min(1.0, max(0.0, (RootingDepth - sumLayer) / L))
        else:
            AvailCap = 0.0

    MaxExtr = AvailCap * UStoreLayerDepth

    # Next step is to make use of the Feddes curve in order to decrease ActEvapUstore when soil moisture values
    # occur above or below ideal plant growing conditions (see also Feddes et al., 1978). h1-h4 values are
    # actually negative, but all values are made positive for simplicity.
    h1 = hb  # cm (air entry pressure)
    h2 = 100  # cm (pF 2 for field capacity)
    h3 = 400  # cm (pF 3, critical pF value)
    h4 = 15849  # cm (pF 4.2, wilting point)

    # According to Brooks-Corey
    par_lambda = 2 / (c - 3)
    if L > 0.0:
        vwc = UStoreLayerDepth / L
    else:
        vwc = 0.0
    vwc = max(vwc, 0.0000001)
    head = hb / (
        ((vwc) / (thetaS - thetaR)) ** (1 / par_lambda)
    )  # Note that in the original formula, thetaR is extracted from vwc, but thetaR is not part of the numerical vwc calculation
    head = max(head,hb)

    # Transform h to a reduction coefficient value according to Feddes et al. (1978).
    # For now: no reduction for head < h2 until following improvement (todo):
    #       - reduction only applied to crops
    if(head <= h1):
        alpha = 1
    elif(head >= h4):
        alpha = 0
    elif((head < h2) & (head > h1)):
        alpha = 1
    elif((head > h3) & (head < h4)):
        alpha = 1 - (head - h3) / (h4 - h3)
    else:
        alpha = 1

    ActEvapUStore = (min(MaxExtr, RestPotEvap, UStoreLayerDepth)) * alpha

    UStoreLayerDepth = UStoreLayerDepth - ActEvapUStore

    RestPotEvap = RestPotEvap - ActEvapUStore
    sumActEvapUStore = ActEvapUStore + sumActEvapUStore

    return UStoreLayerDepth, sumActEvapUStore, RestPotEvap


@jit(nopython=True)    
def infiltration(AvailableForInfiltration, PathFrac, cf_soil, TSoil,InfiltCapSoil,InfiltCapPath, UStoreCapacity, modelSnow, soilInfReduction):
    # First determine if the soil infiltration capacity can deal with the amount of water
    # split between infiltration in undisturbed soil and compacted areas (paths)
    SoilInf = AvailableForInfiltration  * (1 - PathFrac)
    PathInf = AvailableForInfiltration * PathFrac
    if modelSnow & soilInfReduction:
        bb = 1.0 / (1.0 - cf_soil)
        soilInfRedu = _sCurve(TSoil, a=0.0, b=bb, c=8.0) + cf_soil
    else:
        soilInfRedu = 1.0
    MaxInfiltSoil = min(InfiltCapSoil * soilInfRedu, SoilInf)
    MaxInfiltPath = min(InfiltCapPath * soilInfRedu, PathInf)      
    InfiltSoilPath = min(MaxInfiltPath + MaxInfiltSoil, max(0.0, UStoreCapacity))
    
    if (MaxInfiltPath + MaxInfiltSoil) > 0.0:
        InfiltSoil = MaxInfiltSoil * min(1.0, max(0.0, UStoreCapacity)/(MaxInfiltPath + MaxInfiltSoil))
        InfiltPath = MaxInfiltPath * min(1.0, max(0.0, UStoreCapacity)/(MaxInfiltPath + MaxInfiltSoil))
    else:
        InfiltSoil = 0.0
        InfiltPath = 0.0
        
    InfiltExcess = (SoilInf - MaxInfiltSoil) + (PathInf - MaxInfiltPath)
    
    return InfiltSoilPath, InfiltSoil, InfiltPath, SoilInf, PathInf, InfiltExcess

@jit(nopython=True)
def unsatzone_flow_layer(USd,KsatVerFrac,KsatVer,f,z,L_sat,c):

    sum_ast = 0
    #first transfer soil water > maximum soil water capacity layer (iteration is not required because of steady theta (USd))
    st = KsatVerFrac * KsatVer * np.exp(-f * z) * min((USd/L_sat)**c,1.0)
    st_sat = max(0, USd-L_sat)
    USd = USd - min(st, st_sat)
    sum_ast = sum_ast + min(st, st_sat)
    ast = min(st - min(st, st_sat),USd)
    #number of iterations (to reduce "overshooting") based on fixed maximum change in soil water per iteration step (0.02 x maximum soil water capacity layer)
    its = max(int(math.ceil(ast/(L_sat*0.02))),1)
    k = KsatVerFrac * KsatVer/its * np.exp(-f * z)
    for i in range(0,its):
        st = k * min((USd/L_sat)**c,1.0)
        ast = min(st,USd)
        USd = USd - ast
        sum_ast = sum_ast + ast
    
    return USd, sum_ast


@jit(nopython=True)  
def unsatzone_flow(UStoreLayerDepth, InfiltSoilPath, L, z, KsatVerFrac, c, KsatVer, f, thetaS, thetaR, SoilWaterCapacity, SWDold, shape_layer, TransferMethod):
    
    m = 0
    UStoreLayerDepth[m] = UStoreLayerDepth[m] + InfiltSoilPath
    
    if L[m] > 0.0:
        #sbm option for vertical transfer (only for 1 layer)
        if (TransferMethod == 1 and shape_layer == 1):
            Sd = SoilWaterCapacity - SWDold
            if Sd <= 0.00001:
                st = 0.0
            else:
                st = KsatVerFrac[m] * KsatVer * (min(UStoreLayerDepth[m],L[m]*(thetaS-thetaR))/Sd)
                ast = min(st,UStoreLayerDepth[m])
                UStoreLayerDepth[m] = UStoreLayerDepth[m] - ast
        else:
            L_sat = L[m] * (thetaS-thetaR)
            UStoreLayerDepth[m], ast = unsatzone_flow_layer(UStoreLayerDepth[m],KsatVerFrac[m],KsatVer,f,z[m],L_sat,c[m])
    else: 
        ast = 0.0
            
                
    for m in range(1,len(L)):
        UStoreLayerDepth[m] = UStoreLayerDepth[m] + ast
        if L[m] > 0.0:
            L_sat = L[m] * (thetaS-thetaR)
            UStoreLayerDepth[m], ast = unsatzone_flow_layer(UStoreLayerDepth[m],KsatVerFrac[m],KsatVer,f,z[m],L_sat,c[m])
        else:
            ast = 0.0
            
    return ast, UStoreLayerDepth
    
    
@jit(nopython=True)    
def sbm_cell(nodes, nodes_up, ldd, layer, static, dyn, modelSnow, soilInfReduction, timestepsecs, basetimestep, deltaT, nrpaddyirri, shape, TransferMethod, it_kinL=1, ust=0):
        
    shape_layer = layer['UStoreLayerThickness'].shape
    
    # flat new state
    ssf_new = np.zeros(dyn['ssf'].size, dtype=dyn['ssf'].dtype)

    qo_new = np.zeros(dyn['LandRunoff'].size, dtype=dyn['LandRunoff'].dtype)
    qo_new = np.concatenate((qo_new, np.array([0], dtype=dyn['LandRunoff'].dtype)))

    
    # append zero to end to deal with nodata (-1) in indices
    ssf_new = np.concatenate((ssf_new, np.array([0], dtype=dyn['ssf'].dtype)))
    ldd_ = np.concatenate((ldd, np.array([0], dtype=ldd.dtype)))
    slope_ = np.concatenate((static['landSlope'], np.array([0], dtype=static['landSlope'].dtype)))
    
    SWDold = np.zeros(dyn['ssf'].size, dtype=dyn['ssf'].dtype)
    sumUSold = np.zeros(dyn['ssf'].size, dtype=dyn['ssf'].dtype)    
    
    for i in range(len(nodes)):
        for j in range(len(nodes[i])):
            idx = nodes[i][j]
            nbs = nodes_up[i][j]
                        
            sumlayer = np.unique(layer['UStoreLayerThickness'][:,idx].cumsum())
            sumlayer_0 = np.concatenate((np.array([0.0]), sumlayer))
            SWDold[idx] = dyn['SatWaterDepth'][idx] 
            sumUSold[idx] = layer['UStoreLayerDepth'][:,idx].sum()

            # Identify which layers contain the unsaturated zone
            n = np.where(dyn['zi'][idx] > sumlayer_0)[0]
            if len(n) > 1:     
                L = np.concatenate((layer['UStoreLayerThickness'][n[0:-1],idx], np.array([dyn['zi'][idx] - sumlayer_0[n[-1]]]))).astype(np.float64)
            else:
                L = np.array([dyn['zi'][idx]]).astype(np.float64)
            z = L.cumsum()
                
            dyn['ActEvapUStore'][idx] = 0.0
            
            if static['River'][idx]:
                ind = np.where(ldd_[nbs] != ldd_[idx])
                chanperc = np.zeros(ldd_[nbs].size)
                chanperc[ind] = slope_[nbs][ind]/(slope_[idx]+slope_[nbs][ind]) 
            
                ssf_in = np.sum((1-chanperc)*ssf_new[nbs])
                dyn['ssf_toriver'][idx] = np.sum((chanperc)*ssf_new[nbs])/(1000*1000*1000)/timestepsecs
                   
            else:
                ssf_in = np.sum(ssf_new[nbs])
            
            dyn['CellInFlow'][idx] = ssf_in
            
            # Calculate the initial capacity of the unsaturated store U 
            UStoreCapacity = static['SoilWaterCapacity'][idx] - dyn['SatWaterDepth'][idx] - layer['UStoreLayerDepth'][n,idx].sum()

            # Calculate the infiltration flux into the soil column
            InfiltSoilPath, InfiltSoil, InfiltPath, SoilInf, PathInf, InfiltExcess = infiltration(dyn['AvailableForInfiltration'][idx], static['PathFrac'][idx], static['cf_soil'][idx], 
                                          dyn['TSoil'][idx],static['InfiltCapSoil'][idx],static['InfiltCapPath'][idx],UStoreCapacity, modelSnow, soilInfReduction)
            
            dyn['InfiltSoilPath'][idx] = InfiltSoilPath
            dyn['InfiltExcess'][idx] = InfiltExcess
            dyn['SoilInf'][idx] = SoilInf
            dyn['PathInf'][idx] = PathInf
            
            # Using the surface infiltration rate, calculate the flow rate between the different soil layers that contain unsaturated storage assuming gravity based flow only,
            # estimate the gravity based flux rate to the saturated zone (ast) and the updated unsaturated storage for each soil layer.
            ast, layer['UStoreLayerDepth'][:,idx] = unsatzone_flow(layer['UStoreLayerDepth'][:,idx], InfiltSoilPath, L, z, layer['KsatVerFrac'][:,idx], layer['c'][:,idx], static['KsatVer'][idx], static['f'][idx],
                                         static['thetaS'][idx], static['thetaR'][idx], static['SoilWaterCapacity'][idx], SWDold[idx], shape_layer[0], TransferMethod)
            dyn['Transfer'][idx] = ast
            
            # then evapotranspiration from layers
            for k in range(len(L)):
                # For the most upper layer
                if k==0:
                    # Calculate saturation deficity (SC = SWC - S) 
                    SaturationDeficit = static['SoilWaterCapacity'][idx] - dyn['SatWaterDepth'][idx]

                    # First calculate the evaporation of unsaturated storage into the atmosphere from the layer closest to the surface
                    # In case of a s single soil layer
                    if shape_layer[0] == 1:
                        soilevapunsat = dyn['PotSoilEvap'][idx] * min(1.0, SaturationDeficit / static['SoilWaterCapacity'][idx])
                    else:
                        #In case only the most upper soil layer contains unsaturated storage
                        if len(L) == 1:
                            # Check if groundwater level lies below the surface
                            if dyn['zi'][idx] > 0:
                                soilevapunsat = dyn['PotSoilEvap'][idx] * min(1.0, layer['UStoreLayerDepth'][k,idx]/(dyn['zi'][idx]*(static['thetaS'][idx]-static['thetaR'][idx])))
                            else:
                               soilevapunsat = 0.0 
                        else:
                            # In case first layer contains no saturated storage
                            soilevapunsat = dyn['PotSoilEvap'][idx] * min(1.0, layer['UStoreLayerDepth'][k,idx]/(layer['UStoreLayerThickness'][k,idx]*(static['thetaS'][idx]-static['thetaR'][idx])))
                            
                    # Ensure that the unsaturated evaporation rate does not exceed the available unsaturated moisture
                    soilevapunsat = min(soilevapunsat, layer['UStoreLayerDepth'][k,idx])
                    # Update the additional atmospheric demand
                    dyn['PotSoilEvap'][idx] = dyn['PotSoilEvap'][idx] - soilevapunsat                        
                    layer['UStoreLayerDepth'][k,idx] = layer['UStoreLayerDepth'][k,idx] - soilevapunsat
                    
                    if shape_layer[0] == 1:
                        soilevapsat = 0.0
                    else:
                        if len(L) == 0 or len(L) == 1:
                            soilevapsat = dyn['PotSoilEvap'][idx] * min(1.0, (layer['UStoreLayerThickness'][k,idx] - dyn['zi'][idx])/ layer['UStoreLayerThickness'][k,idx])
                            soilevapsat = min(soilevapsat, (layer['UStoreLayerThickness'][k,idx] - dyn['zi'][idx]) * (static['thetaS'][idx] - static['thetaR'][idx]))
                        else:
                           soilevapsat = 0.0 
                    
                    
                    dyn['soilevap'][idx] = soilevapunsat + soilevapsat
                    dyn['soilevapunsat'][idx] = soilevapunsat
                    dyn['soilevapsat'][idx] = soilevapsat
                    dyn['SatWaterDepth'][idx] = dyn['SatWaterDepth'][idx] - soilevapsat
                    
                    # transpiration from saturated store
                    wetroots = _sCurve(dyn['zi'][idx], a=static['ActRootingDepth'][idx], c=static['rootdistpar'][idx])
                    dyn['ActEvapSat'][idx] = min(dyn['PotTrans'][idx] * wetroots, dyn['SatWaterDepth'][idx])
                    dyn['SatWaterDepth'][idx] = dyn['SatWaterDepth'][idx] - dyn['ActEvapSat'][idx]
                    RestPotTrans = dyn['PotTrans'][idx] - dyn['ActEvapSat'][idx]

                    # actual transpiration from UStore
                    layer['UStoreLayerDepth'][k,idx], dyn['ActEvapUStore'][idx], RestPotTrans = actTransp_unsat_SBM(static['ActRootingDepth'][idx], layer['UStoreLayerDepth'][k,idx], 
                                                                          sumlayer_0[k], RestPotTrans, dyn['ActEvapUStore'][idx], layer['c'][k,idx], L[k], static['thetaS'][idx], static['thetaR'][idx], static['AirEntryPressure'][idx], ust)
                # For the layers below  
                else:
                    # actual transpiration from UStore
                    layer['UStoreLayerDepth'][k,idx], dyn['ActEvapUStore'][idx], RestPotTrans = actTransp_unsat_SBM(static['ActRootingDepth'][idx], layer['UStoreLayerDepth'][k,idx], 
                                                                          sumlayer_0[k], RestPotTrans, dyn['ActEvapUStore'][idx], layer['c'][k,idx], L[k], static['thetaS'][idx], static['thetaR'][idx], static['AirEntryPressure'][idx], ust)

            #check soil moisture balance per layer
            du = 0.0
            for k in range(L.size-1,-1,-1):
                du = max(0,layer['UStoreLayerDepth'][k,idx] - L[k]*(static['thetaS'][idx]-static['thetaR'][idx]))
                layer['UStoreLayerDepth'][k,idx] = layer['UStoreLayerDepth'][k,idx] - du
                if k > 0:
                    layer['UStoreLayerDepth'][k-1,idx] = layer['UStoreLayerDepth'][k-1,idx] + du
            
            Ksat = layer['KsatVerFrac'][len(L)-1,idx] * static['KsatVer'][idx] * np.exp(-static['f'][idx] * dyn['zi'][idx])  
            
            UStoreCapacity = static['SoilWaterCapacity'][idx] - dyn['SatWaterDepth'][idx] - layer['UStoreLayerDepth'][n,idx].sum()
            
            MaxCapFlux = max(0.0, min(Ksat, dyn['ActEvapUStore'][idx], UStoreCapacity, dyn['SatWaterDepth'][idx]))
            
            if dyn['zi'][idx] > static['ActRootingDepth'][idx]:
                CapFluxScale = static['CapScale'][idx] / (static['CapScale'][idx] + dyn['zi'][idx] - static['ActRootingDepth'][idx]) * timestepsecs / basetimestep
            else:
                CapFluxScale = 0.0
                    
            CapFlux = MaxCapFlux * CapFluxScale
            
            netCapflux = CapFlux
            actCapFlux = 0.0
            for k in range(L.size-1,-1,-1):
                toadd = min(netCapflux, max(L[k]*(static['thetaS'][idx]-static['thetaR'][idx]) - layer['UStoreLayerDepth'][k,idx], 0.0))
                layer['UStoreLayerDepth'][k,idx] = layer['UStoreLayerDepth'][k,idx] + toadd
                netCapflux = netCapflux - toadd
                actCapFlux = actCapFlux + toadd
            
            dyn['CapFlux'][idx] = actCapFlux
                        
            DeepKsat = static['KsatVer'][idx] * np.exp(-static['f'][idx] * static['SoilThickness'][idx])
            DeepTransfer = min(dyn['SatWaterDepth'][idx], DeepKsat)
            dyn['ActLeakage'][idx] = max(0.0, min(static['MaxLeakage'][idx], DeepTransfer))
                                                  
            r = (ast - actCapFlux - dyn['ActLeakage'][idx] - dyn['ActEvapSat'][idx] - soilevapsat) * static['DW'][idx]*1000
            ssf_new[idx], dyn['zi'][idx], ExfiltSatWater = kinematic_wave_ssf(ssf_in, dyn['ssf'][idx], dyn['zi'][idx], r, static['KsatHorFrac'][idx], 
                                              static['KsatVer'][idx], static['landSlope'][idx], static['neff'][idx], static['f'][idx], 
                                              static['SoilThickness'][idx], deltaT, static['DL'][idx]*1000, static['DW'][idx]*1000, static['ssfmax'][idx])
            
            dyn['zi'][idx] = min(dyn['zi'][idx], static['SoilThickness'][idx])
            dyn['SatWaterDepth'][idx] =  (static['SoilThickness'][idx] - dyn['zi'][idx]) * (static['thetaS'][idx] - static['thetaR'][idx])         
                        
            n_new = np.where(dyn['zi'][idx] > sumlayer_0)[0]           
            if len(n_new) > 1:     
                L_new = np.concatenate((layer['UStoreLayerThickness'][n_new[0:-1],idx], np.array([dyn['zi'][idx] - sumlayer_0[n_new[-1]]]))).astype(np.float64)
            else:
                L_new = np.array([dyn['zi'][idx]]).astype(np.float64)                
                   
            ExfiltFromUstore = 0.0
            for k in range(L.size-1,-1,-1):
                if (np.where(n_new == k))[0].size > 0:
                    ExfiltFromUstore = max(0,layer['UStoreLayerDepth'][k,idx] - L_new[k]*(static['thetaS'][idx]-static['thetaR'][idx]))
                else:
                    ExfiltFromUstore = layer['UStoreLayerDepth'][k,idx]                
                layer['UStoreLayerDepth'][k,idx] = layer['UStoreLayerDepth'][k,idx] - ExfiltFromUstore
                if k > 0:
                    layer['UStoreLayerDepth'][k-1,idx] = layer['UStoreLayerDepth'][k-1,idx] + ExfiltFromUstore

            dyn['ExfiltWater'][idx] = ExfiltSatWater + ExfiltFromUstore
            dyn['ExfiltSatWater'][idx] = ExfiltSatWater
            dyn['ExfiltFromUstore'][idx] = ExfiltFromUstore
            
            dyn['ExcessWater'][idx] = dyn['AvailableForInfiltration'][idx] - InfiltSoilPath - InfiltExcess + du    
            dyn['ActInfilt'][idx] = InfiltSoilPath - du
            
            #Separation between compacted and non compacted areas (correction with the satflow du)
            if (InfiltSoil + InfiltPath) > 0.0:
                dyn['ActInfiltSoil'][idx] = InfiltSoil - du * InfiltSoil / (InfiltPath + InfiltSoil)
                dyn['ActInfiltPath'][idx] = InfiltPath - du * InfiltPath / (InfiltPath + InfiltSoil)
            else:
                dyn['ActInfiltSoil'][idx] = 0.0
                dyn['ActInfiltPath'][idx] = 0.0
                
            dyn['ExcessWaterSoil'][idx] = max(SoilInf - dyn['ActInfiltSoil'][idx], 0.0)
            dyn['ExcessWaterPath'][idx] = max(PathInf - dyn['ActInfiltPath'][idx], 0.0)
            
            ponding_add = 0
            if nrpaddyirri > 0:
                if static['h_p'][idx] > 0:
                    ponding_add = min(dyn['ExfiltWater'][idx] + dyn['ExcessWater'][idx] + dyn['InfiltExcess'][idx], static['h_p'][idx] - dyn['PondingDepth'][idx])
                    dyn['PondingDepth'][idx] = dyn['PondingDepth'][idx] + ponding_add
            
            dyn['InwaterO'][idx] = max(dyn['ExfiltWater'][idx] + dyn['ExcessWater'][idx] + dyn['InfiltExcess'][idx] + dyn['RunoffLandCells'][idx] - dyn['ActEvapOpenWaterLand'][idx] - ponding_add, 0.0) * (static['xl'][idx] * static['yl'][idx]) * 0.001 / timestepsecs
            
            dyn['sumUStoreLayerDepth'][idx] = layer['UStoreLayerDepth'][:,idx].sum()
            
            # volumetric water contents per soil layer and root zone            
            for k in range(layer['UStoreLayerThickness'][:,idx].size):
                if (np.where(n_new == k))[0].size > 0:
                    if layer['UStoreLayerThickness'][k,idx] > 0:
                        layer['vwc'][k,idx] =  (layer['UStoreLayerDepth'][k,idx] + (layer['UStoreLayerThickness'][k,idx] - L_new[k]) * (static['thetaS'][idx] - static['thetaR'][idx])) / layer['UStoreLayerThickness'][k,idx] + static['thetaR'][idx]
                else:
                    layer['vwc'][k,idx] = static['thetaS'][idx]
                
                layer['vwc_perc'][k,idx] = (layer['vwc'][k,idx]/static['thetaS'][idx]) * 100.0
                
                
            rootStore_unsat = 0
            for k in range(L_new.size):
                if L_new[k] > 0:
                    rootStore_unsat =  rootStore_unsat + (max(0.0, static['ActRootingDepth'][idx] - sumlayer_0[k])/L_new[k]) * layer['UStoreLayerDepth'][k,idx]

            dyn['RootStore_unsat'][idx] = rootStore_unsat 
            
            
    
    acc_flow = np.zeros(dyn['LandRunoff'].size, dtype=dyn['LandRunoff'].dtype)
    acc_flow = np.concatenate((acc_flow, np.array([0], dtype=dyn['LandRunoff'].dtype)))
    qo_toriver_acc = np.copy(acc_flow)

    q = dyn['InwaterO'] / static['DL']
    for v in range(0,it_kinL):
	
        qo_new = np.zeros(dyn['LandRunoff'].size, dtype=dyn['LandRunoff'].dtype)
        qo_new = np.concatenate((qo_new, np.array([0], dtype=dyn['LandRunoff'].dtype)))
        
        for i in range(len(nodes)):
            for j in range(len(nodes[i])):
                idx = nodes[i][j]
                nbs = nodes_up[i][j]
                
                if static['River'][idx]:
                    ind = np.where(ldd_[nbs] != ldd_[idx])
                    chanperc = np.zeros(ldd_[nbs].size)
                    chanperc[ind] = slope_[nbs][ind]/(slope_[idx]+slope_[nbs][ind]) 
 
                    if static['SW'][idx] > 0.0:
                        qo_in = np.sum((1-chanperc)*qo_new[nbs])
                        qo_toriver_vol = np.sum(chanperc*qo_new[nbs]) * (timestepsecs/it_kinL)
                    else:
                        qo_in = 0.0
                        qo_toriver_vol = np.sum(qo_new[nbs]) * (timestepsecs/it_kinL)
                else:
                    qo_in = np.sum(qo_new[nbs])
                    qo_toriver_vol = 0.0

                    
                qo_new[idx] = kinematic_wave(qo_in, dyn['LandRunoff'][idx], q[idx], dyn['AlphaL'][idx], static['Beta'][idx], timestepsecs/it_kinL, static['DL'][idx])
                                
                acc_flow[idx] = acc_flow[idx] + qo_new[idx] * (timestepsecs/it_kinL)
                dyn['Qo_in'][idx] = dyn['Qo_in'][idx] + qo_in * (timestepsecs/it_kinL)
                qo_toriver_acc[idx] = qo_toriver_acc[idx] + qo_toriver_vol
                if static['SW'][idx] > 0:
                    WaterLevelL = (dyn['AlphaL'][idx] * np.power(qo_new[idx], static['Beta'][idx])) / static['SW'][idx]
                Pl = static['SW'][idx] + (2.0 * WaterLevelL)
                dyn['AlphaL'][idx] = static['AlpTermL'][idx] * np.power(Pl, static['AlpPow'][idx])
                dyn['LandRunoff'][idx]= qo_new[idx]
    qo_new = acc_flow/timestepsecs
    dyn['qo_toriver'][:] = qo_toriver_acc[:-1]/timestepsecs
    dyn['Qo_in'][:] = dyn['Qo_in'][:] / timestepsecs
    
    dyn['SoilWatbal'][:] = (dyn['ActInfilt'][:] - ((dyn['SatWaterDepth'][:] + dyn['sumUStoreLayerDepth'][:]) - (sumUSold[:] + SWDold[:])) +
                  (dyn['CellInFlow'][:]-ssf_new[:-1])/(static['DW'][:]*static['DL'][:]*1000*1000) - dyn['ExfiltWater'][:] - dyn['soilevap'][:] - dyn['ActEvapUStore'][:] -
                  dyn['ActEvapSat'][:] - dyn['ActLeakage'][:])
    
    return ssf_new[:-1], qo_new[:-1], dyn, layer



class WflowModel(pcraster.framework.DynamicModel):
    """
    .. versionchanged:: 0.91
        - Calculation of GWScale moved to resume() to allow fitting.

    .. versionadded:: 0.91
        - added S-curve for freezing soil infiltration reduction calculations

    .. todo::
        - add slope based quick-runoff -> less percolation on hillslopes...
  """

    def __init__(self, cloneMap, Dir, RunDir, configfile):
        pcraster.framework.DynamicModel.__init__(self)

        self.UStoreLayerDepth = []
        self.caseName = os.path.abspath(Dir)
        self.clonemappath = os.path.join(os.path.abspath(Dir), "staticmaps", cloneMap)
        pcr.setclone(self.clonemappath)
        self.runId = RunDir
        self.Dir = os.path.abspath(Dir)
        self.configfile = configfile
        self.SaveDir = os.path.join(self.Dir, self.runId)

    def irrigationdemand(self, pottrans, acttrans, irareas):
        """
        Determine irrigation water demand from the difference bewteen potential
        transpiration and actual transpiration.

        :param pottrans: potential transpiration (epot minus interception and soil/open water evaporation)
        :param acttrans: actual transpiration
        :param ir_areas: maps of irrigation areas

        :return: demand
        """

        Et_diff = pcr.areaaverage(pottrans - acttrans, pcr.nominal(irareas))
        # Now determine demand in m^3/s for each area
        sqmarea = pcr.areatotal(self.reallength * self.reallength, pcr.nominal(irareas))
        m3sec = Et_diff * sqmarea / 1000.0 / self.timestepsecs

        return Et_diff, m3sec


    def updateRunOff(self):
        """
      Updates the kinematic wave reservoir. Should be run after updates to Q
      """
        self.WaterLevelR = (self.AlphaR * pow(self.RiverRunoff, self.Beta)) / self.Bw
        # wetted perimeter (m)
        Pr = self.Bw + (2 * self.WaterLevelR)
        # Alpha River
        self.AlphaR = self.AlpTermR * pow(Pr, self.AlpPow)
        self.OldKinWaveVolumeR = self.KinWaveVolumeR
        self.KinWaveVolumeR = self.WaterLevelR * self.Bw * self.DCL
        
        self.dyn['AlphaR'] = pcr.pcr2numpy(self.AlphaR,self.mv).ravel()
        
        self.WaterLevelL = pcr.ifthenelse( self.SW > 0, (self.AlphaL * pow(self.LandRunoff, self.Beta)) / self.SW, 0.0)
        
        Pl = self.SW + (2 * self.WaterLevelL)
        # Alpha Land
        self.AlphaL = self.AlpTermL * pow(Pl, self.AlpPow)
        self.OldKinWaveVolumeL = self.KinWaveVolumeL
        self.KinWaveVolumeL = self.WaterLevelL * self.SW * self.DL
        
        self.dyn['AlphaL'] = pcr.pcr2numpy(self.AlphaL,self.mv).ravel()       
        
        
    def stateVariables(self):
        """
        returns a list of state variables that are essential to the model.
        This list is essential for the resume and suspend functions to work.

        This function is specific for each model and **must** be present.

       :var self.RiverRunoff: Surface runoff in the kin-wave resrvoir [m^3/s]
       :var self.LandRunoff: Surface runoff in the kin-wave resrvoir [m^3/s]
       :var self.SurfaceRunoffDyn: Surface runoff in the dyn-wave resrvoir [m^3/s]
       :var self.WaterLevelR: Water level in the river kin-wave reservoir [m]
       :var self.WaterLevelL: Water level in the land kin-wave reservoir [m]
       :var self.WaterLevelDyn: Water level in the dyn-wave resrvoir [m^]
       :var self.Snow: Snow pack [mm]
       :var self.SnowWater: Snow pack water [mm]
       :var self.TSoil: Top soil temperature [oC]
       :var self.UStoreLayerDepth: Water in the Unsaturated Store, per layer [mm]
       :var self.SatWaterDepth: Water in the saturated store [mm]
       :var self.CanopyStorage: Amount of water on the Canopy [mm]
       :var self.ReservoirVolume: Volume of each reservoir [m^3]
       :var self.LakeWaterLevel: Water level in natural/uncontrolled lakes [m]
       :var self.GlacierStore: Thickness of the Glacier in a gridcell [mm]
       """
        states = [
            "RiverRunoff",
            "WaterLevelR",
            "LandRunoff",
            "WaterLevelL",
            "SatWaterDepth",
            "Snow",
            "TSoil",
            "UStoreLayerDepth",
            "SnowWater",
            "CanopyStorage",
            "SubsurfaceFlow",
        ]
        if hasattr(self, "GlacierFrac"):
            states.append("GlacierStore")

        if hasattr(self, "ReserVoirSimpleLocs"):
            states.append("ReservoirVolume")

        if hasattr(self, "LakeLocs"):
            states.append("LakeWaterLevel")

        if hasattr(self, "nrpaddyirri"):
            if self.nrpaddyirri > 0:
                states.append("PondingDepth")
        return states

    def supplyCurrentTime(self):
        """
      gets the current time in seconds after the start of the run
      """
        return self.currentTimeStep() * self.timestepsecs

    def suspend(self):

        self.logger.info("Saving initial conditions...")
        self.wf_suspend(os.path.join(self.SaveDir, "outstate"))

        if self.OverWriteInit:
            self.logger.info("Saving initial conditions over start conditions...")
            self.wf_suspend(self.Dir + "/instate/")

    def parameters(self):
        """
        Define all model parameters here that the framework should handle for the model
        See wf_updateparameters and the parameters section of the ini file
        If you use this make sure to all wf_updateparameters at the start of the dynamic section
        and at the start/end of the initial section
        """
        modelparameters = []

        # Static model parameters e.g.
        # modelparameters.append(self.ParamType(name="RunoffGeneratingGWPerc",stack="intbl/RunoffGeneratingGWPerc.tbl",type="static",default=0.1))
        # 3: Input time series ###################################################
        self.P_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Precipitation", "/inmaps/P"
        )  # timeseries for rainfall
        self.PET_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "EvapoTranspiration", "/inmaps/PET"
        )  # timeseries for rainfall"/inmaps/PET"          # potential evapotranspiration
        self.TEMP_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Temperature", "/inmaps/TEMP"
        )  # timeseries for rainfall "/inmaps/TEMP"          # global radiation
        self.Inflow_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Inflow", "/inmaps/IF"
        )  # timeseries for rainfall "/inmaps/IF" # in/outflow locations (abstractions)

        # Meteo and other forcing
        modelparameters.append(
            self.ParamType(
                name="Precipitation",
                stack=self.P_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="PotenEvap",
                stack=self.PET_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="Temperature",
                stack=self.TEMP_mapstack,
                type="timeseries",
                default=10.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="Inflow",
                stack=self.Inflow_mapstack,
                type="timeseries",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )

        modelparameters.append(
            self.ParamType(
                name="IrrigationAreas",
                stack="staticmaps/wflow_irrigationareas.map",
                type="staticmap",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="IrrigationSurfaceIntakes",
                stack="staticmaps/wflow_irrisurfaceintake.map",
                type="staticmap",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="IrrigationPaddyAreas",
                stack="staticmaps/wflow_irrigationpaddyareas.map",
                type="staticmap",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="IrrigationSurfaceReturn",
                stack="staticmaps/wflow_irrisurfacereturns.map",
                type="staticmap",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )

        modelparameters.append(
            self.ParamType(
                name="h_max",
                stack="staticmaps/wflow_hmax.map",
                type="staticmap",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="h_min",
                stack="staticmaps/wflow_hmin.map",
                type="staticmap",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="h_p",
                stack="staticmaps/wflow_hp.map",
                type="staticmap",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )

        return modelparameters

    def initial(self):
        """
    Initial part of the model, executed only once. Reads all static data from disk


    *Soil*

    :var M.tbl: M parameter in the SBM model. Governs the decay of Ksat with depth [-]
    :var thetaR.tbl: Residual water content [mm/mm]
    :var thetaS.tbl: Saturated water content (porosity) [mm/mm]
    :var KsatVer.tbl: Saturated conductivity [mm/d]
    :var PathFrac.tbl: Fraction of compacted area per grid cell [-]
    :var InfiltCapSoil.tbl: Soil infiltration capacity [m/d]
    :var InfiltCapPath.tbl: Infiltration capacity of the compacted areas [mm/d]
    :var SoilMinThickness.tbl: Minimum wdepth of the soil [mm]
    :var SoilThickness.tbl: Maximum depth of the soil [m]
    :var RootingDepth.tbl: Depth of the roots [mm]
    :var MaxLeakage.tbl: Maximum leakage out of the soil profile [mm/d]
    :var CapScale.tbl: Scaling factor in the Capilary rise calculations (100) [mm/d]
    :var RunoffGeneratingGWPerc: Fraction of the soil depth that contributes to subcell runoff (0.1) [-]
    :var rootdistpar.tbl: Determine how roots are linked to water table. The number
        should be negative. A more negative  number means that all roots are wet if the water
        table is above the lowest part of the roots.
        A less negative number smooths this. [mm] (default = -80000)



    *Canopy*

    :var CanopyGapFraction.tbl: Fraction of precipitation that does not hit the canopy directly [-]
    :var MaxCanopyStorage.tbl: Canopy interception storage capacity [mm]
    :var EoverR.tbl: Ratio of average wet canopy evaporation rate over rainfall rate [-]

    *Surface water*

    :var N.tbl: Manning's N parameter
    :var N_river.tbl: Manning's N parameter for cells marked as river


    *Snow and frozen soil modelling parameters*

    :var cf_soil.tbl: Soil infiltration reduction factor when soil is frozen [-] (< 1.0)
    :var TTI.tbl: critical temperature for snowmelt and refreezing  (1.000) [oC]
    :var TT.tbl: defines interval in which precipitation falls as rainfall and snowfall (-1.41934) [oC]
    :var Cfmax.tbl: meltconstant in temperature-index ( 3.75653) [-]
    :var WHC.tbl: fraction of Snowvolume that can store water (0.1) [-]
    :var w_soil.tbl: Soil temperature smooth factor. Given for daily timesteps. (0.1125) [-] Wigmosta, M. S., L. J. Lane, J. D. Tagestad, and A. M. Coleman (2009).

    """
        global statistics
        global multpars
        global updateCols

        self.thestep = pcr.scalar(0)
        self.basetimestep = 86400
        self.SSSF = False
        pcr.setglobaloption("unittrue")
        
        self.mv = -999
        self.count = 0

        self.logger.info("running for " + str(self.nrTimeSteps()) + " timesteps")

        # Set and get defaults from ConfigFile here ###################################

        self.Tslice = int(configget(self.config, "model", "Tslice", "1"))
        self.reinit = int(configget(self.config, "run", "reinit", "0"))
        self.OverWriteInit = int(configget(self.config, "model", "OverWriteInit", "0"))
        self.updating = int(configget(self.config, "model", "updating", "0"))
        self.updateFile = configget(self.config, "model", "updateFile", "no_set")
        self.TransferMethod = int(
            configget(self.config, "model", "transfermethod", "0")
        )
        self.maxitsupply = int(configget(self.config, "model", "maxitsupply", "5"))
        self.UST = int(configget(self.config, "model", "Whole_UST_Avail", "0"))
        self.NRiverMethod = int(configget(self.config, "model", "nrivermethod", "1"))
        self.kinwaveIters = int(configget(self.config, "model", "kinwaveIters", "0"))
        self.kinwaveRiverTstep = int(configget(self.config, "model", "kinwaveRiverTstep", "0"))
        self.kinwaveLandTstep = int(configget(self.config, "model", "kinwaveLandTstep", "0"))        
        if self.kinwaveIters == 1:
            self.logger.info(
                "Using sub timestep for kinematic wave (iterate)"
            )
            if self.kinwaveRiverTstep > 0:
                self.logger.info(
                    "Using a fixed timestep (seconds) for kinematic wave river flow: " + str(self.kinwaveRiverTstep)
                )
            if self.kinwaveLandTstep > 0:
                self.logger.info(
                    "Using a fixed timestep (seconds) for kinematic wave overland flow: " + str(self.kinwaveLandTstep)
                )                
        if self.TransferMethod == 1:
            self.logger.info(
                "Applying the original topog_sbm vertical transfer formulation"
            )
        self.sCatch = int(configget(self.config, "model", "sCatch", "0"))
        self.intbl = configget(self.config, "model", "intbl", "intbl")

        self.modelSnow = int(configget(self.config, "model", "ModelSnow", "1"))
        self.soilInfReduction = int(configget(self.config, "model", "soilInfRedu", "0"))
        sizeinmetres = int(configget(self.config, "layout", "sizeinmetres", "0"))
        alf = float(configget(self.config, "model", "Alpha", "60"))
        # TODO: make this into a list for all gauges or a map
        Qmax = float(configget(self.config, "model", "AnnualDischarge", "300"))
        self.UpdMaxDist = float(configget(self.config, "model", "UpdMaxDist", "100"))

        self.MaxUpdMult = float(configget(self.config, "model", "MaxUpdMult", "1.3"))
        self.MinUpdMult = float(configget(self.config, "model", "MinUpdMult", "0.7"))
        self.UpFrac = float(configget(self.config, "model", "UpFrac", "0.8"))

        # self.ExternalQbase=int(configget(self.config,'model','ExternalQbase','0'))
        self.waterdem = int(configget(self.config, "model", "waterdem", "0"))
        WIMaxScale = float(configget(self.config, "model", "WIMaxScale", "0.8"))
        self.MassWasting = int(configget(self.config, "model", "MassWasting", "0"))
        self.nrLayers = int(configget(self.config, "model", "nrLayers", "1"))

        # static maps to use (normally default)
        wflow_subcatch = configget(
            self.config, "model", "wflow_subcatch", "staticmaps/wflow_subcatch.map"
        )
        wflow_dem = configget(
            self.config, "model", "wflow_dem", "staticmaps/wflow_dem.map"
        )
        wflow_ldd = configget(
            self.config, "model", "wflow_ldd", "staticmaps/wflow_ldd.map"
        )
        wflow_river = configget(
            self.config, "model", "wflow_river", "staticmaps/wflow_river.map"
        )
        wflow_riverlength = configget(
            self.config,
            "model",
            "wflow_riverlength",
            "staticmaps/wflow_riverlength.map",
        )
        wflow_riverlength_fact = configget(
            self.config,
            "model",
            "wflow_riverlength_fact",
            "staticmaps/wflow_riverlength_fact.map",
        )
        wflow_landuse = configget(
            self.config, "model", "wflow_landuse", "staticmaps/wflow_landuse.map"
        )
        wflow_soil = configget(
            self.config, "model", "wflow_soil", "staticmaps/wflow_soil.map"
        )
        wflow_gauges = configget(
            self.config, "model", "wflow_gauges", "staticmaps/wflow_gauges.map"
        )
        wflow_inflow = configget(
            self.config, "model", "wflow_inflow", "staticmaps/wflow_inflow.map"
        )
        wflow_riverwidth = configget(
            self.config, "model", "wflow_riverwidth", "staticmaps/wflow_riverwidth.map"
        )
        wflow_streamorder = configget(
            self.config,
            "model",
            "wflow_streamorder",
            "staticmaps/wflow_streamorder.map",
        )

        # 2: Input base maps ########################################################
        subcatch = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_subcatch), 0.0, fail=True)
        )  # Determines the area of calculations (all cells > 0)
        subcatch = pcr.ifthen(subcatch > 0, subcatch)

        self.Altitude = self.wf_readmap(
            os.path.join(self.Dir, wflow_dem), 0.0, fail=True
        )  # * pcr.scalar(pcr.defined(subcatch)) # DEM
        self.TopoLdd = pcr.ldd(
            self.wf_readmap(os.path.join(self.Dir, wflow_ldd), 0.0, fail=True)
        )  # Local
        self.TopoId = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_subcatch), 0.0, fail=True)
        )  # area map
        self.River = pcr.cover(
            pcr.boolean(
                self.wf_readmap(os.path.join(self.Dir, wflow_river), 0.0, fail=True)
            ),
            0,
        )

        self.RiverLength = pcr.cover(
            self.wf_readmap(os.path.join(self.Dir, wflow_riverlength), 0.0), 0.0
        )
        # Factor to multiply riverlength with (defaults to 1.0)
        self.RiverLengthFac = self.wf_readmap(
            os.path.join(self.Dir, wflow_riverlength_fact), 1.0
        )

        # read landuse and soilmap and make sure there are no missing points related to the
        # subcatchment map. Currently sets the lu and soil type  type to 1
        self.LandUse = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_landuse), 0.0, fail=True)
        )
        self.LandUse = pcr.cover(self.LandUse, pcr.ordinal(subcatch > 0))
        self.Soil = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_soil), 0.0, fail=True)
        )
        self.Soil = pcr.cover(self.Soil, pcr.ordinal(subcatch > 0))
        self.OutputLoc = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_gauges), 0.0, fail=True)
        )  # location of output gauge(s)
        self.InflowLoc = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_inflow), 0.0)
        )  # location abstractions/inflows.
        self.RiverWidth = self.wf_readmap(os.path.join(self.Dir, wflow_riverwidth), 0.0)
        
        self.SubCatchFlowOnly = int(
            configget(self.config, "model", "SubCatchFlowOnly", "0")
        )
        self.OutputId = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_subcatch), 0.0, fail=True)
        )  # location of subcatchment
        # Temperature correction poer cell to add

        self.TempCor = self.wf_readmap(
            self.Dir
            + "\\"
            + configget(
                self.config,
                "model",
                "TemperatureCorrectionMap",
                "staticmaps/wflow_tempcor.map",
            ),
            0.0,
        )

        self.ZeroMap = 0.0 * pcr.scalar(subcatch)  # map with only zero's

        # Set static initial values here #########################################
        self.pi = 3.1416
        self.e = 2.7183
        self.SScale = 100.0
        self.Latitude = pcr.ycoordinate(pcr.boolean(self.Altitude))
        self.Longitude = pcr.xcoordinate(pcr.boolean(self.Altitude))

        # Read parameters NEW Method
        self.logger.info("Linking parameters to landuse, catchment and soil...")
        self.wf_updateparameters()
        
        if hasattr(self, "LAI"):
            # Sl must also be defined
            if not hasattr(self, "Sl"):
                logging.error(
                    "Sl (specific leaf storage) not defined! Needed becausee LAI is defined."
                )
                logging.error("Please add it to the modelparameters section. e.g.:")
                logging.error(
                    "Sl=inmaps/clim/LCtoSpecificLeafStorage.tbl,tbl,0.5,1,inmaps/clim/LC.map"
                )
            if not hasattr(self, "Kext"):
                logging.error(
                    "Kext (canopy extinction coefficient) not defined! Needed becausee LAI is defined."
                )
                logging.error("Please add it to the modelparameters section. e.g.:")
                logging.error(
                    "Kext=inmaps/clim/LCtoExtinctionCoefficient.tbl,tbl,0.5,1,inmaps/clim/LC.map"
                )
            if not hasattr(self, "Swood"):
                logging.error(
                    "Swood wood (branches, trunks) canopy storage not defined! Needed becausee LAI is defined."
                )
                logging.error("Please add it to the modelparameters section. e.g.:")
                logging.error(
                    "Swood=inmaps/clim/LCtoBranchTrunkStorage.tbl,tbl,0.5,1,inmaps/clim/LC.map"
                )

            self.Cmax = self.Sl * self.LAI + self.Swood
            self.CanopyGapFraction = pcr.exp(-self.Kext * self.LAI)
            self.np_CanopyGapFraction = pcr.pcr2numpy(self.CanopyGapFraction,self.mv)
            # TODO: Add MAXLAI and CWf lookup
        else:
            self.Cmax = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/MaxCanopyStorage.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                1.0,
            )
            self.CanopyGapFraction = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/CanopyGapFraction.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.1,
            )
            
            self.EoverR = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/EoverR.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.1,
            )

        if not hasattr(self, "DemandReturnFlowFraction"):
            self.DemandReturnFlowFraction = self.ZeroMap

        self.RootingDepth = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/RootingDepth.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            750.0,
        )  # rooting depth
        
        self.AirEntryPressure = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/AirEntryPressure.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            10.0,
        )  # air entry pressure (cm) of soil (Brooks-Corey)
        
        #rootdistpar determine how roots are linked to water table.
        self.rootdistpar = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/rootdistpar.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            -8000,
        )  # rootdistpar

        # Soil parameters
        # infiltration capacity if the soil [mm/day]
        self.InfiltCapSoil = (
            self.readtblDefault(
                self.Dir + "/" + self.intbl + "/InfiltCapSoil.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                100.0,
            )
            * self.timestepsecs
            / self.basetimestep
        )
        self.CapScale = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/CapScale.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            100.0,
        )  #

        # infiltration capacity of the compacted
        self.InfiltCapPath = (
            self.readtblDefault(
                self.Dir + "/" + self.intbl + "/InfiltCapPath.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                10.0,
            )
            * self.timestepsecs
            / self.basetimestep
        )
        self.MaxLeakage = (
            self.readtblDefault(
                self.Dir + "/" + self.intbl + "/MaxLeakage.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.0,
            )
            * self.timestepsecs
            / self.basetimestep
        )
        self.MaxPercolation = (
            self.readtblDefault(
                self.Dir + "/" + self.intbl + "/MaxPercolation.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.0,
            )
            * self.timestepsecs
            / self.basetimestep
        )

        # areas (paths) in [mm/day]
        # Fraction area with compacted soil (Paths etc.)
        self.PathFrac = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/PathFrac.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.01,
        )
        # thickness of the soil
        self.SoilThickness = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/SoilThickness.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            2000.0,
        )
                
        self.thetaR = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/thetaR.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.01,
        )
        
        self.thetaS = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/thetaS.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.6,
        )
        
        
        # minimum thickness of soild
        self.SoilMinThickness = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/SoilMinThickness.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            500.0,
        )

        # KsatVer = $2\inmaps\KsatVer.map
        self.KsatVer = (
            self.readtblDefault(
                self.Dir + "/" + self.intbl + "/KsatVer.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                3000.0,
            )
            * self.timestepsecs
            / self.basetimestep
        )
        self.MporeFrac = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/MporeFrac.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.0,
        )

        self.KsatHorFrac = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/KsatHorFrac.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            1.0,
        )

        # Check if we have irrigation areas
        tt = pcr.pcr2numpy(self.IrrigationAreas, 0.0)
        self.nrirri = tt.max()
        # Check of we have paddy irrigation areas
        tt = pcr.pcr2numpy(self.IrrigationPaddyAreas, 0.0)
        self.nrpaddyirri = tt.max()
            

        self.Beta = pcr.scalar(0.6)  # For sheetflow

        self.M = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/M.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            300.0,
        )  # Decay parameter in Topog_sbm
        self.N = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/N.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.072,
        )  # Manning overland flow
        if self.NRiverMethod == 1:
            self.NRiver = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/N_River.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.036,
            )  # Manning river
        if self.NRiverMethod == 2:
            self.NRiver = self.readtblFlexDefault(
                self.Dir + "/" + self.intbl + "/N_River.tbl", 0.036, wflow_streamorder
            )

        self.WaterFrac = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/WaterFrac.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.0,
        )  # Fraction Open water
        self.et_RefToPot = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/et_reftopot.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            1.0,
        )  # Fraction Open water

        if self.modelSnow:
            # HBV Snow parameters
            # critical temperature for snowmelt and refreezing:  TTI= 1.000
            self.TTI = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/TTI.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                1.0,
            )
            # TT = -1.41934 # defines interval in which precipitation falls as rainfall and snowfall
            self.TT = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/TT.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                -1.41934,
            )
            self.TTM = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/TTM.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                -1.41934,
            )
            # Cfmax = 3.75653 # meltconstant in temperature-index
            self.Cfmax = (self.readtblDefault(
                self.Dir + "/" + self.intbl + "/Cfmax.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                3.75653,
                )
                * self.timestepsecs
                / self.basetimestep
            )
            # WHC= 0.10000        # fraction of Snowvolume that can store water
            self.WHC = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/WHC.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.1,
            )
            # Wigmosta, M. S., L. J. Lane, J. D. Tagestad, and A. M. Coleman (2009).
            self.w_soil = (
                self.readtblDefault(
                    self.Dir + "/" + self.intbl + "/w_soil.tbl",
                    self.LandUse,
                    subcatch,
                    self.Soil,
                    0.9 * 3.0 / 24.0,
                )
                * self.timestepsecs
                / self.basetimestep
            )
            if self.soilInfReduction:
                self.cf_soil = pcr.min(
                    0.99,
                    self.readtblDefault(
                        self.Dir + "/" + self.intbl + "/cf_soil.tbl",
                        self.LandUse,
                        subcatch,
                        self.Soil,
                        0.038,
                    ),
                )  # Soil infiltration reduction factor for frozen soil
                
        # Determine real slope and cell length

        self.xl, self.yl, self.reallength = pcrut.detRealCellLength(
            self.ZeroMap, sizeinmetres
        )
        self.Slope = pcr.slope(self.Altitude)
        self.Slope = pcr.max(0.00001, self.Slope * pcr.celllength() / self.reallength)

        #self.N = pcr.ifthenelse(self.River, self.NRiver, self.N)

        if hasattr(self, "ReserVoirSimpleLocs") or hasattr(
            self, "LakeLocs"
        ):
            self.ReserVoirLocs = self.ZeroMap
            self.filter_P_PET = self.ZeroMap + 1.0

        if hasattr(self, "ReserVoirSimpleLocs"):
            # Check if we have simple and or complex reservoirs
            self.ReserVoirSimpleLocs = pcr.nominal(self.ReserVoirSimpleLocs)
            self.ReservoirSimpleAreas = pcr.nominal(self.ReservoirSimpleAreas)
            tt_simple = pcr.pcr2numpy(self.ReserVoirSimpleLocs, 0.0)
            self.nrresSimple = np.size(np.where(tt_simple > 0.0)[0])
            self.ReserVoirLocs = self.ReserVoirLocs + pcr.cover(
                pcr.scalar(self.ReserVoirSimpleLocs), 0.0
            )
            
            res_area = pcr.cover(pcr.scalar(self.ReservoirSimpleAreas), 0.0)
            self.filter_P_PET = pcr.ifthenelse(
                pcr.boolean(pcr.cover(res_area, pcr.scalar(0.0))),
                res_area * 0.0,
                self.filter_P_PET,
            )
        else:
            self.nrresSimple = 0

        if hasattr(self, "LakeLocs"):
            #add parameter for lake threshold estimation
            self.estimatelakethresh = int(configget(self.config, "model", "estimatelakethresh", "0"))
            self.LakeAreasMap = pcr.nominal(self.LakeAreasMap)
            self.LakeLocs = pcr.nominal(self.LakeLocs)
            tt_lake = pcr.pcr2numpy(self.LakeLocs, 0.0)
            #self.nrlake = tt_lake.max()
            self.nrlake = np.size(np.where(tt_lake > 0.0)[0])
            self.ReserVoirLocs = self.ReserVoirLocs + pcr.cover(
                pcr.scalar(self.LakeLocs), 0.0
            )
            lake_area = pcr.cover(pcr.scalar(self.LakeAreasMap), 0.0)
            self.filter_P_PET = pcr.ifthenelse(
                lake_area > 0, lake_area * 0.0, self.filter_P_PET
            )

            # read files
            self.sh = {}
            lake_ids = pcr.ifthen(self.LakeStorFunc == 2, self.LakeLocs)
            np_lake_ids = pcr.pcr2numpy(lake_ids, 0)
            np_lake_ids_u = np.unique(np_lake_ids[np.nonzero(np_lake_ids)])
            if np.size(np_lake_ids_u) > 0:
                for item in np.nditer(np_lake_ids_u):
                    self.sh[int(item)] = np.loadtxt(
                        self.Dir
                        + "/"
                        + self.intbl
                        + "/Lake_SH_"
                        + str(item)
                        + ".tbl"
                    )
            self.hq = {}
            lake_ids = pcr.ifthen(self.LakeOutflowFunc == 1, self.LakeLocs)
            np_lake_ids = pcr.pcr2numpy(lake_ids, 0)
            np_lake_ids_u = np.unique(np_lake_ids[np.nonzero(np_lake_ids)])
            if np.size(np_lake_ids_u) > 0:
                for item in np.nditer(np_lake_ids_u):
                    self.hq[int(item)] = np.loadtxt(
                        self.Dir
                        + "/"
                        + self.intbl
                        + "/Lake_HQ_"
                        + str(item)
                        + ".tbl",
                        skiprows=3,
                    )
                    
            #Ini for the Modified Puls Approach (Burek et al., 2013, LISFLOOD)
            #Check which lakes uses the puls approach (ResOutflowFunc = 3)
            #And if the corresponding res_e=2 and ResStorFunc=2
            
            #Update Lake_b in ini if ResThreshold different from zero
            
            np_lakeoutflowfunc_old = pcr.pcr2numpy(self.LakeOutflowFunc, 0)
            self.LakeOutflowFunc = pcr.ifthenelse(
                    pcr.pcrand(self.LakeOutflowFunc == 3, self.LakeStorFunc != 1),
                    2,
                    self.LakeOutflowFunc
                    )
            self.LakeOutflowFunc = pcr.ifthenelse(
                    pcr.pcrand(self.LakeOutflowFunc == 3, self.Lake_e != 2.0),
                    2,
                    self.LakeOutflowFunc
                    )
            np_lakeoutflowfunc = pcr.pcr2numpy(self.LakeOutflowFunc, 0)
            if np_lakeoutflowfunc_old.sum() != np_lakeoutflowfunc.sum():
                self.logger.warning("Lake outflow modelling using the modified puls approach selected "+ 
                                    "but found contradictory arguments for LakeStorFunc/Lake_e: "+
                                    "using the general iteration method instead")
            

        else:
            self.nrlake = 0

        if (self.nrresSimple + self.nrlake) > 0:
            self.ReserVoirLocs = pcr.ordinal(self.ReserVoirLocs)
            self.logger.info(
                "A total of "
                + str(self.nrresSimple)
                + " simple reservoirs and "
                + str(self.nrlake)
                + " lakes found."
            )
            self.ReserVoirDownstreamLocs = pcr.downstream(
                self.TopoLdd, self.ReserVoirLocs
            )
            self.TopoLddOrg = self.TopoLdd
            self.TopoLdd = pcr.lddrepair(
                pcr.cover(
                    pcr.ifthen(pcr.boolean(self.ReserVoirLocs), pcr.ldd(5)),
                    self.TopoLdd,
                )
            )

            tt_filter = pcr.pcr2numpy(self.filter_P_PET, 1.0)
            self.filterResArea = tt_filter.min()

        # Slope land surface either provided as map, or based om DEM (default)
        self.landSlope = pcr.max(0.00001, self.wf_readmap(os.path.join(self.Dir, "staticmaps/Slope.map"), self.Slope))            

        # soil thickness based on topographic wetness index (see Environmental modelling: finding simplicity in complexity)
        # 1: calculate wetness index
        # 2: Scale the soil thickness (now actually a max) based on the index, also apply a minimum soil thickness
        WI = pcr.ln(
            pcr.accuflux(self.TopoLdd, 1) / self.landSlope
        )  # Topographic wetnesss index. Scale WI by zone/subcatchment assuming these ara also geological units
        WIMax = pcr.areamaximum(WI, self.TopoId) * WIMaxScale
        self.SoilThickness = pcr.max(
            pcr.min(self.SoilThickness, (WI / WIMax) * self.SoilThickness),
            self.SoilMinThickness,
        )

        self.SoilWaterCapacity = self.SoilThickness * (self.thetaS - self.thetaR)

        # determine number of layers based on total soil thickness
        UStoreLayerThickness = configget(
            self.config, "model", "UStoreLayerThickness", "0"
        )
        if UStoreLayerThickness != "0":
            self.nrLayers = len(UStoreLayerThickness.split(","))
            self.maxLayers = self.nrLayers + 1
        else:
            UStoreLayerThickness = self.SoilThickness
            self.nrLayers = 1
            self.maxLayers = self.nrLayers

        self.KsatVerFrac = []
        self.c = []
        for n in range(self.maxLayers):
            self.KsatVerFrac.append(
                self.readtblLayersDefault(
                    self.Dir + "/" + self.intbl + "/KsatVerFrac.tbl",
                    self.LandUse,
                    subcatch,
                    self.Soil,
                    n,
                    1.0,
                )
            )
            self.c.append(
                self.readtblLayersDefault(
                    self.Dir + "/" + self.intbl + "/c.tbl",
                    self.LandUse,
                    subcatch,
                    self.Soil,
                    n,
                    10.0,
                )
            )

        # limit roots to top 99% of first zone
        self.RootingDepth = pcr.min(self.SoilThickness * 0.99, self.RootingDepth)

        # Which columns/gauges to use/ignore in updating
        self.UpdateMap = self.ZeroMap

        if self.updating:
            _tmp = pcr.pcr2numpy(self.OutputLoc, 0.0)
            gaugear = _tmp
            touse = np.zeros(gaugear.shape, dtype="int")

            for thecol in updateCols:
                idx = (gaugear == thecol).nonzero()
                touse[idx] = thecol

            self.UpdateMap = pcr.numpy2pcr(pcr.Nominal, touse, 0.0)
            # Calculate distance to updating points (upstream) annd use to scale the correction
            # ldddist returns zero for cell at the gauges so add 1.0 tp result
            self.DistToUpdPt = pcr.cover(
                pcr.min(
                    pcr.ldddist(
                        self.TopoLdd, pcr.boolean(pcr.cover(self.UpdateMap, 0)), 1
                    )
                    * self.reallength
                    / pcr.celllength(),
                    self.UpdMaxDist,
                ),
                self.UpdMaxDist,
            )

        # Initializing of variables
        self.logger.info("Initializing of model variables..")
        self.TopoLdd = pcr.lddmask(self.TopoLdd, pcr.boolean(self.TopoId))
        catchmentcells = pcr.maptotal(pcr.scalar(self.TopoId))

        # Limit lateral flow per subcatchment (make pits at all subcatch boundaries)
        # This is very handy for Ribasim etc...
        if self.SubCatchFlowOnly > 0:
            self.logger.info("Creating subcatchment-only drainage network (ldd)")
            ds = pcr.downstream(self.TopoLdd, self.TopoId)
            usid = pcr.ifthenelse(ds != self.TopoId, self.TopoId, 0)
            self.TopoLdd = pcr.lddrepair(
                pcr.ifthenelse(pcr.boolean(usid), pcr.ldd(5), self.TopoLdd)
            )

        # Used to seperate output per LandUse/management classes
        OutZones = self.LandUse

        self.QMMConv = self.timestepsecs / (
            self.reallength * self.reallength * 0.001
        )  # m3/s --> actial mm of water over the cell
        # self.QMMConvUp = 1000.0 * self.timestepsecs / ( pcr.catchmenttotal(pcr.cover(1.0), self.TopoLdd) * self.reallength * self.reallength)  #m3/s --> mm over upstreams
        temp = (
            pcr.catchmenttotal(pcr.cover(1.0), self.TopoLdd)
            * self.reallength
            * 0.001
            * 0.001
            * self.reallength
        )
        self.QMMConvUp = pcr.cover(self.timestepsecs * 0.001) / temp
        self.ToCubic = (
            self.reallength * self.reallength * 0.001
        ) / self.timestepsecs  # m3/s
        self.KinWaveVolumeR = self.ZeroMap
        self.OldKinWaveVolumeR = self.ZeroMap
        self.KinWaveVolumeL = self.ZeroMap
        self.OldKinWaveVolumeL = self.ZeroMap
        
        self.sumprecip = self.ZeroMap  # accumulated rainfall for water balance
        self.sumevap = self.ZeroMap  # accumulated evaporation for water balance
        self.sumrunoff = self.ZeroMap  # accumulated runoff for water balance
        self.sumint = self.ZeroMap  # accumulated interception for water balance
        self.sumleakage = self.ZeroMap
        self.sumoutflow = self.ZeroMap
        self.sumsnowmelt = self.ZeroMap
        self.CumRad = self.ZeroMap
        self.SnowMelt = self.ZeroMap
        self.CumPrec = self.ZeroMap
        self.CumInwaterMM = self.ZeroMap
        self.CumInfiltExcess = self.ZeroMap
        self.CumExfiltWater = self.ZeroMap
        self.CumSurfaceWater = self.ZeroMap
        self.watbal = self.ZeroMap
        self.CumEvap = self.ZeroMap
        self.CumPotenEvap = self.ZeroMap
        self.CumPotenTrans = self.ZeroMap
        self.CumInt = self.ZeroMap
        self.CumRad = self.ZeroMap
        self.CumLeakage = self.ZeroMap
        self.CumPrecPol = self.ZeroMap
        self.SatWaterFlux = self.ZeroMap
        self.SumCellWatBal = self.ZeroMap
        self.PathInfiltExceeded = self.ZeroMap
        self.SoilInfiltExceeded = self.ZeroMap
        self.CumOutFlow = self.ZeroMap
        self.CumCellInFlow = self.ZeroMap
        self.CumIF = self.ZeroMap
        self.CumActInfilt = self.ZeroMap
        self.IRSupplymm = self.ZeroMap
        self.Aspect = pcr.scalar(pcr.aspect(self.Altitude))  # aspect [deg]
        self.Aspect = pcr.ifthenelse(self.Aspect <= 0.0, pcr.scalar(0.001), self.Aspect)
        # On Flat areas the Aspect function fails, fill in with average...
        self.Aspect = pcr.ifthenelse(
            pcr.defined(self.Aspect),
            self.Aspect,
            pcr.areaaverage(self.Aspect, self.TopoId),
        )
        
        # Factor on river length (self.RiverLengthFac) only used in combination with
        # calculated (by wflow_sbm) slope    
        # Set DCL to riverlength if that is longer that the basic length calculated from grid
        drainlength = detdrainlength(self.TopoLdd, self.xl, self.yl)
        # Multiply with Factor (taken from upscaling operation, defaults to 1.0 if no map is supplied)
        self.DCL = drainlength * pcr.max(1.0, self.RiverLengthFac)
        # Correct slope for extra length of the river in a gridcel
        riverslopecor = drainlength / self.DCL
        self.riverSlope = self.Slope * riverslopecor
        
        # If river slope available as map, also provide river length 
        self.riverSlope = pcr.max(0.00001, self.wf_readmap(os.path.join(self.Dir, "staticmaps/RiverSlope.map"), self.riverSlope))
        if os.path.isfile(os.path.join(self.Dir, wflow_riverlength)):
            self.DCL = self.RiverLength # m


        # Determine river width from DEM, upstream area and yearly average discharge
        # Scale yearly average Q at outlet with upstream are to get Q over whole catchment
        # Alf ranges from 5 to > 60. 5 for hardrock. large values for sediments
        # "Noah J. Finnegan et al 2005 Controls on the channel width of rivers:
        # Implications for modeling fluvial incision of bedrock"
        if (self.nrresSimple + self.nrlake) > 0:
            upstr = pcr.catchmenttotal(1, self.TopoLddOrg)
        else:
            upstr = pcr.catchmenttotal(1, self.TopoLdd)
        Qscale = upstr / pcr.mapmaximum(upstr) * Qmax
        W = (
            (alf * (alf + 2.0) ** (0.6666666667)) ** (0.375)
            * Qscale ** (0.375)
            * (pcr.max(0.0001, pcr.windowaverage(self.riverSlope, pcr.celllength() * 4.0)))
            ** (-0.1875)
            * self.NRiver ** (0.375)
        )
        # Use supplied riverwidth if possible, else calulate
        self.RiverWidth = pcr.ifthenelse(self.RiverWidth <= 0.0, W, self.RiverWidth)
        #Use W instead of RiverWidth for reservoirs and lake cells
        if self.nrresSimple > 0:
            self.RiverWidth = pcr.ifthenelse(
                    pcr.cover(pcr.scalar(self.ReservoirSimpleAreas), 0.0) > 0.0,
                    W,
                    self.RiverWidth
                    )
        if self.nrlake > 0:
            self.RiverWidth = pcr.ifthenelse(
                    pcr.cover(pcr.scalar(self.LakeAreasMap), 0.0) > 0.0,
                    W,
                    self.RiverWidth
                    )
        
        # water depth (m)
        # set width for kinematic wave to cell width for all cells
        self.Bw = detdrainwidth(self.TopoLdd, self.xl, self.yl)
        # However, in the main river we have real flow so set the width to the
        # width of the river

        self.Bw = pcr.ifthenelse(self.River, self.RiverWidth, self.Bw)

        # Add rivers to the WaterFrac, but check with waterfrac map and correct
        self.RiverFrac = pcr.min(
            1.0,
            pcr.ifthenelse(
                self.River, (self.RiverWidth * self.DCL) / (self.xl * self.yl), 0
            ),
        )
        
        self.WaterFrac = pcr.max(self.WaterFrac - self.RiverFrac, 0)

        # term for Alpha (River)
        self.AlpTermR = pow((self.NRiver / (pcr.sqrt(self.riverSlope))), self.Beta)
        
        # power for Alpha
        self.AlpPow = (2.0 / 3.0) * self.Beta
        

        
        # term for Alpha (Land)
        self.AlpTermL = pow((self.N / (pcr.sqrt(self.landSlope))), self.Beta)
        
        #Estimate LakeThreshold depending on outlet characteristics
        if (self.nrlake > 0 and self.estimatelakethresh == 1):
            #initial waterLevel
            level_map_path = self.Dir + "/instate/WaterLevelR.map"
            if os.path.exists(level_map_path):
                level_map = pcr.readmap(level_map_path)
            else:
                level_map = self.ZeroMap
            alphaR = self.AlpTermR * (pcr.downstream(self.TopoLdd, (self.Bw+2*level_map))) ** self.AlpPow
            outletRivLevel = alphaR * self.LakeAvgOut ** self.Beta
            #Lake Threshold = Lake Level - 130% River Level
            #130% = Adjustment linked to possible uncertainties in River Level estimation
            self.LakeThreshold = pcr.ifthenelse(
                    self.LakeThreshold > 0.0,
                    self.LakeThreshold,
                    pcr.max(self.LakeAvgLevel - 1.3 * outletRivLevel, 0.0)
                    )
            
            #Reupdate rating curve coefficient
            self.Lake_b = pcr.ifthenelse(
                    self.LakeOutflowFunc == 3,
                    self.LakeAvgOut / (self.LakeAvgLevel - self.LakeThreshold) ** 2,
                    self.Lake_b
                    )
        
        # calculate catchmentsize
        self.upsize = pcr.catchmenttotal(self.xl * self.yl, self.TopoLdd)
        self.csize = pcr.areamaximum(self.upsize, self.TopoId)
        
        self.wf_multparameters()
        
        # determine flow network and upstream nodes
        self.np_ldd = pcr.pcr2numpy(self.TopoLdd, self.mv)

        # convert pcr objects to numpy for kinemativ wave surface water
        np_zeros = pcr.pcr2numpy(self.ZeroMap, self.mv).ravel()
        np_2d_zeros = pcr.pcr2numpy(self.ZeroMap, self.mv)
        
        self.neff = (self.thetaS - self.thetaR)
        self.f = pcr.abs((self.thetaS - self.thetaR) /self.M)
        
        self.DL = detdrainlength(self.TopoLdd, self.xl, self.yl)
        self.DW = (self.xl * self.yl)/self.DL
        
        # width for overland kinematic reservoir        
        self.SW = pcr.ifthenelse(self.River, pcr.max(self.DW - self.RiverWidth,0), self.DW)
        
        layer_dtype = np.dtype(
                [('c', np.float64),
                 ('UStoreLayerDepth', np.float64),
                 ('st', np.float64),
                 ('KsatVerFrac', np.float64),
                 ('ActEvapUstore', np.float64),
                 ('vwc', np.float64),
                 ('vwc_perc', np.float64),
                 ('UStoreLayerThickness', np.float64),
                 ('UStest', np.float64)
                ])        
        
        self.layer = np.zeros((self.maxLayers,np_zeros.size), dtype=layer_dtype)       
        
        for n in range(self.maxLayers):
            self.layer['c'][n] = pcr.pcr2numpy(self.c[n], self.mv).ravel()
            self.layer['KsatVerFrac'][n] = pcr.pcr2numpy(self.KsatVerFrac[n], self.mv).ravel()
        
        static_dtype = np.dtype(
                [('KsatVer', np.float64),
                 ('KsatHorFrac', np.float64),
                 ('f', np.float64),
                 ('thetaS', np.float64),
                 ('thetaR', np.float64),
                 ('SoilThickness', np.float64),
                 ('InfiltCapSoil', np.float64),
                 ('PathFrac', np.float64),
                 ('InfiltCapPath', np.float64),
                 ('cf_soil', np.float64),
                 ('neff', np.float64),
                 ('landSlope', np.float64),
                 ('riverSlope', np.float64),
                 ('SoilWaterCapacity', np.float64),
                 ('MaxLeakage', np.float64),
                 ('CanopyGapFraction', np.float64),
                 ('CapScale', np.float64),
                 ('rootdistpar', np.float64),
                 ('River', np.float64),
                 ('DL', np.float64),
                 ('DW', np.float64),
                 ('SW', np.float64),
                 ('Beta', np.float64),
                 ('DCL', np.float64),
                 ('reallength', np.float64),
                 ('RootingDepth', np.float64),
                 ('AirEntryPressure', np.float64),
                 ('ActRootingDepth', np.float64),
                 ('ssfmax', np.float64),
                 ('xl', np.float64),
                 ('yl', np.float64),
                 ('h_p', np.float64),
                 ('Bw', np.float64),
                 ('AlpPow', np.float64),
                 ('AlpTermR', np.float64),
                 ('AlpTermL', np.float64)
                 ])
        
            
        self.static = np.zeros(np_zeros.size, dtype=static_dtype)   
        self.static['KsatVer'] = pcr.pcr2numpy(self.KsatVer, self.mv).ravel()
        self.static['KsatHorFrac'] = pcr.pcr2numpy(self.KsatHorFrac, self.mv).ravel()
        self.static['f'] = pcr.pcr2numpy(self.f, self.mv).ravel()
        self.static['thetaS'] = pcr.pcr2numpy(self.thetaS, self.mv).ravel()
        self.static['thetaR'] = pcr.pcr2numpy(self.thetaR, self.mv).ravel()
        self.static['SoilThickness'] = pcr.pcr2numpy(self.SoilThickness, self.mv).ravel()
        self.static['InfiltCapSoil'] = pcr.pcr2numpy(self.InfiltCapSoil, self.mv).ravel()
        self.static['InfiltCapPath'] = pcr.pcr2numpy(self.InfiltCapPath, self.mv).ravel()
        self.static['PathFrac'] = pcr.pcr2numpy(self.PathFrac, self.mv).ravel()
        if self.modelSnow & self.soilInfReduction:
            self.static['cf_soil'] = pcr.pcr2numpy(self.cf_soil, self.mv).ravel()
        self.static['neff'] = pcr.pcr2numpy(self.neff, self.mv).ravel()
        self.static['landSlope'] = pcr.pcr2numpy(self.landSlope, self.mv).ravel()
        self.static['riverSlope'] = pcr.pcr2numpy(self.riverSlope, self.mv).ravel()
        self.static['MaxLeakage'] = pcr.pcr2numpy(self.MaxLeakage, self.mv).ravel()
        self.static['SoilWaterCapacity'] = pcr.pcr2numpy(self.SoilWaterCapacity, self.mv).ravel()
        self.static['CanopyGapFraction'] = pcr.pcr2numpy(self.CanopyGapFraction, self.mv).ravel()
        self.static['CapScale'] = pcr.pcr2numpy(self.CapScale, self.mv).ravel()
        self.static['rootdistpar'] = pcr.pcr2numpy(self.rootdistpar, self.mv).ravel()
        self.static['River'] = pcr.pcr2numpy(self.River, self.mv).ravel()
        self.static['DL'] = pcr.pcr2numpy(self.DL, self.mv).ravel()
        self.static['DW'] = pcr.pcr2numpy(self.DW, self.mv).ravel()
        self.static['SW'] = pcr.pcr2numpy(self.SW, self.mv).ravel()
        self.static['Beta'] = pcr.pcr2numpy(self.Beta, self.mv).ravel()
        self.static['DCL'] = pcr.pcr2numpy(self.DCL, self.mv).ravel()
        self.static['reallength'] = pcr.pcr2numpy(self.reallength, self.mv).ravel()
        self.static['xl'] = pcr.pcr2numpy(self.xl, self.mv).ravel()
        self.static['yl'] = pcr.pcr2numpy(self.yl, self.mv).ravel()
        self.static['RootingDepth'] = pcr.pcr2numpy(self.RootingDepth, self.mv).ravel()
        self.static['AirEntryPressure'] = pcr.pcr2numpy(self.AirEntryPressure, self.mv).ravel()
        self.static['h_p'] = pcr.pcr2numpy(self.h_p, self.mv).ravel()
        self.static['Bw'] = pcr.pcr2numpy(self.Bw, self.mv).ravel()
        self.static['AlpPow'] = pcr.pcr2numpy(self.AlpPow, self.mv).ravel()
        self.static['AlpTermR'] = pcr.pcr2numpy(self.AlpTermR, self.mv).ravel()
        self.static['AlpTermL'] = pcr.pcr2numpy(self.AlpTermL, self.mv).ravel()
        self.static['ssfmax'] = ((self.static['KsatHorFrac'] * self.static['KsatVer'] * self.static['landSlope']) / self.static['f'] * (np.exp(0)-np.exp(-self.static['f'] * self.static['SoilThickness'])))

        # implement layers as numpy arrays
        np_SumThickness = np_zeros
        nrLayersMap = np_zeros  
        for n in np.arange(0, self.maxLayers):
            np_sumL = np_SumThickness
            if self.nrLayers > 1 and n < self.nrLayers:
                UstoreThick_temp = (
                    float(UStoreLayerThickness.split(",")[n]) + np_zeros
                )
                UstoreThick = np.minimum(UstoreThick_temp,np.maximum(self.static['SoilThickness']-np_sumL,0.0))
            else:
                UstoreThick_temp = np.maximum(self.static['SoilThickness'] - np_sumL,0.0)
                UstoreThick = np.maximum(self.static['SoilThickness']-np_sumL, 0.0)                   
            
            np_SumThickness = UstoreThick_temp + np_SumThickness
            
            nrLayersMap = np.where(self.static['SoilThickness'] - np_sumL > 0.0, nrLayersMap + 1.0, nrLayersMap)
            
            self.layer['UStoreLayerThickness'][n] = UstoreThick
        
        dyn_dtype = np.dtype(
                [('AvailableForInfiltration', np.float64),
                 ('SoilInf', np.float64),
                 ('PathInf', np.float64),
                 ('TSoil', np.float64),
                 ('zi', np.float64),
                 ('SatWaterDepth', np.float64),
                 ('ssf', np.float64),
                 ('LandRunoff', np.float64),
                 ('PotTrans', np.float64),
                 ('PotSoilEvap', np.float64),
                 ('ActEvapOpenWaterLand', np.float64),
                 ('ActEvapOpenWaterRiver', np.float64),
                 ('RiverRunoff', np.float64),
                 ('AlphaL', np.float64),
                 ('AlphaR', np.float64),
                 ('InfiltExcess', np.float64),
                 ('ExcessWater', np.float64),
                 ('ExcessWaterSoil', np.float64),
                 ('ExcessWaterPath', np.float64),
                 ('ExfiltWater', np.float64),
                 ('ExfiltSatWater', np.float64),
                 ('ExfiltFromUstore', np.float64),
                 ('soilevap', np.float64),
                 ('soilevapunsat', np.float64),
                 ('soilevapsat', np.float64),
                 ('Transfer', np.float64),
                 ('CapFlux', np.float64),
                 ('qo_toriver', np.float64),
                 ('ssf_toriver', np.float64),
                 ('ActEvapUStore', np.float64),
                 ('ActEvapSat', np.float64),
                 ('ActInfilt', np.float64),
                 ('ActInfiltSoil', np.float64),
                 ('ActInfiltPath', np.float64),
                 ('RunoffLandCells', np.float64),
                 ('ActLeakage', np.float64),
                 ('PondingDepth', np.float64),
                 ('RootStore_unsat', np.float64),
                 ('CellInFlow', np.float64),
                 ('InwaterO', np.float64),
                 ('SoilWatbal', np.float64),
                 ('Qo_in', np.float64),
                 ('WaterLevelL', np.float64),
                 ('InfiltSoilPath',np.float64),
                 ('sumUStoreLayerDepth', np.float64),
                 ('SatWaterDepthOld', np.float64),
                 ('sumUStoreLayerDepthOld', np.float64),
                 ('InfiltWater', np.float64)
                 ])        
        
        self.dyn = np.zeros(np_zeros.size, dtype=dyn_dtype)
        
        self.shape = np_2d_zeros.shape
        
        # initialize us-ds network for all cells
        self.nodes, self.nodes_up = set_dd(self.np_ldd)
        # initialize us-ds network for river cells only
        np_ldd_riv = self.np_ldd.copy()
        np_ldd_riv.flat[self.static['River']==0] = self.mv
        self.rnodes, self.rnodes_up = set_dd(np_ldd_riv)
                
        # Save some summary maps
        self.logger.info("Saving summary maps...")

        self.logger.info("End of initial section")

    def default_summarymaps(self):
        """
          Returns a list of default summary-maps at the end of a run.
          This is model specific. You can also add them to the [summary]section of the ini file but stuff
          you think is crucial to the model should be listed here
          """
        lst = [
            "self.RiverWidth",
            "self.Cmax",
            "self.csize",
            "self.upsize",
            "self.EoverR",
            "self.RootingDepth",
            "self.CanopyGapFraction",
            "self.InfiltCapSoil",
            "self.InfiltCapPath",
            "self.PathFrac",
            "self.thetaR",
            "self.thetaS",
            "self.SoilMinThickness",
            "self.SoilThickness",
            "self.nrLayersMap",
            "self.KsatVer",
            "self.M",
            "self.SoilWaterCapacity",
            "self.et_RefToPot",
            "self.riverSlope",
            "self.landSlope",
            "self.N",
            "self.RiverFrac",
            "self.WaterFrac",
            "self.xl",
            "self.yl",
            "self.reallength",
            "self.DCL",
            "self.Bw",
            "self.SW",
            "self.PathInfiltExceeded",
            "self.SoilInfiltExceeded",
            "self.DW",
            "self.DL",
        ]

        return lst

    def resume(self):

        if self.reinit == 1:
            self.logger.info("Setting initial conditions to default")
            self.SatWaterDepth = self.SoilWaterCapacity * 0.85

            self.zi = pcr.max(
                    0.0, self.SoilThickness - self.SatWaterDepth / (self.thetaS - self.thetaR)
                    )
            
            self.SubsurfaceFlow = ((self.KsatHorFrac * self.KsatVer * self.landSlope)/ pcr.abs(self.f))* (pcr.exp(-pcr.abs(self.f)*(self.zi)) - pcr.exp(-pcr.abs(self.f)*(self.SoilThickness))) * self.DW * 1000

            for n in range(self.maxLayers):
                self.UStoreLayerDepth.append(self.ZeroMap)

            self.WaterLevelR = self.ZeroMap
            self.WaterLevelL = self.ZeroMap
            self.RiverRunoff = self.ZeroMap
            self.LandRunoff = self.ZeroMap
            
            self.Snow = self.ZeroMap
            self.SnowWater = self.ZeroMap
            self.TSoil = self.ZeroMap + 10.0
            self.CanopyStorage = self.ZeroMap
            if hasattr(self, "ReserVoirSimpleLocs"):
                self.ReservoirVolume = self.ResMaxVolume * self.ResTargetFullFrac
            if hasattr(self, "LakeLocs"):
                self.LakeWaterLevel = self.LakeAvgLevel
            if hasattr(self, "GlacierFrac"):
                self.GlacierStore = self.wf_readmap(
                    os.path.join(self.Dir, "staticmaps", "wflow_glacierstore.map"),
                    55.0 * 100
                        
                )
            if self.nrpaddyirri > 0:
                self.PondingDepth = self.ZeroMap
 
        else:
            self.logger.info("Setting initial conditions from state files")
            self.wf_resume(os.path.join(self.Dir, "instate"))

        Pr = self.Bw + (2.0 * self.WaterLevelR)
        self.AlphaR = self.AlpTermR * pow(Pr, self.AlpPow)
        self.dyn['AlphaR'] = pcr.pcr2numpy(self.AlphaR, self.mv).ravel()
        
        Pl = self.SW + (2.0 * self.WaterLevelL)
        self.AlphaL = self.AlpTermL * pow(Pl, self.AlpPow)
        self.dyn['AlphaL'] = pcr.pcr2numpy(self.AlphaL, self.mv).ravel()
        
        self.OldLandRunoff = self.LandRunoff
        self.LandRunoffMM = self.LandRunoff * self.QMMConv
        # Determine initial kinematic wave volume overland
        self.KinWaveVolumeL = self.WaterLevelL * self.SW * self.DL
        self.OldKinWaveVolumeL = self.KinWaveVolumeL
        
        self.OldRiverRunoff = self.RiverRunoff

        self.RiverRunoffMM = self.RiverRunoff * self.QMMConv
        # Determine initial kinematic wave volume
        self.KinWaveVolumeR = self.WaterLevelR * self.Bw * self.DCL
        self.OldKinWaveVolumeR = self.KinWaveVolumeR

        self.QCatchmentMM = self.RiverRunoff * self.QMMConvUp
       
        self.InitialStorage = (
            self.SatWaterDepth
            + sum_list_cover(self.UStoreLayerDepth, self.ZeroMap)
            + self.CanopyStorage
        )
        self.CellStorage = self.SatWaterDepth + sum_list_cover(
            self.UStoreLayerDepth, self.ZeroMap
        )

        # Determine actual water depth
        self.zi = pcr.max(
            0.0, self.SoilThickness - self.SatWaterDepth / (self.thetaS - self.thetaR)
        )
        # TOPOG_SBM type soil stuff
        self.f = (self.thetaS - self.thetaR) / self.M


    def dynamic(self):
        """
        Stuff that is done for each timestep of the model

        Below a list of variables that can be save to disk as maps or as
        timeseries (see ini file for syntax):

        Dynamic variables
        ~~~~~~~~~~~~~~~~~

        All these can be saved per timestep if needed (see the config file [outputmaps] section).

        :var self.Precipitation: Gross precipitation [mm]
        :var self.Temperature: Air temperature [oC]
        :var self.PotenEvap: Potential evapotranspiration [mm]
        :var self.PotTransSoil: Potential Transpiration/Openwater and soil evap (after subtracting Interception from PotenEvap) [mm]
        :var self.Transpiration: plant/tree transpiration [mm]
        :var self.ActEvapOpenWater: actual open water evaporation [mm]
        :var self.soilevap: base soil evaporation [mm]
        :var self.Interception: Actual rainfall interception [mm]
        :var self.ActEvap: Actual evaporation (transpiration + Soil evap + open water evap) [mm]
        :var self.RiverRunoff: Surface runoff in the kinematic wave [m^3/s]
        :var self.SurfaceRunoffDyn: Surface runoff in the dyn-wave resrvoir [m^3/s]
        :var self.SurfaceRunoffCatchmentMM: Surface runoff in the dyn-wave reservoir expressed in mm over the upstream (catchment) area
        :var self.WaterLevelDyn: Water level in the dyn-wave resrvoir [m^]
        :var self.ActEvap: Actual EvapoTranspiration [mm] (minus interception losses)
        :var self.ExcessWater: Water that cannot infiltrate due to saturated soil [mm]
        :var self.InfiltExcess: Infiltration excess water [mm]
        :var self.WaterLevel: Water level in the kinematic wave [m] (above the bottom)
        :var self.ActInfilt: Actual infiltration into the unsaturated zone [mm]
        :var self.CanopyStorage: actual canopystorage (only for subdaily timesteps) [mm]
        :var self.SatWaterDepth: Amount of water in the saturated store [mm]
        :var self.UStoreLayerDepth: Amount of water in the unsaturated store, per layer [mm]
        :var self.zi: depth of the water table in mm below the surface [mm]
        :var self.Snow: Snow depth [mm]
        :var self.SnowWater: water content of the snow [mm]
        :var self.TSoil: Top soil temperature [oC]
        :var self.SatWaterDepth: amount of available water in the saturated part of the soil [mm]
        :var self.UStoreDepth: amount of available water in the unsaturated zone [mm]
        :var self.Transfer: downward flux from unsaturated to saturated zone [mm]
        :var self.CapFlux: capilary flux from saturated to unsaturated zone [mm]
        :var self.CanopyStorage: Amount of water on the Canopy [mm]
        :var self.RunoffCoeff: Runoff coefficient (Q/P) for each cell taking into account the whole upstream area [-]
        :var self.SurfaceWaterSupply: the negative Inflow (water demand) that could be met from the surfacewater [m^3/s]
        :var self.Recharge: simple recharge to groundwater (transfer - capillary flux - saturated act evap) [mm]


        Static variables
        ~~~~~~~~~~~~~~~~

        :var self.Altitude: The altitude of each cell [m]
        :var self.Bw: Width of the river [m]
        :var self.River: booolean map indicating the presence of a river [-]
        :var self.DLC: length of the river within a cell [m]
        :var self.ToCubic: Mutiplier to convert mm to m^3/s for fluxes
        """

        # Read forcing data and dynamic parameters

        self.wf_updateparameters()
        self.Precipitation = pcr.max(0.0, self.Precipitation)

        # NB This may interfere with lintul link
        if hasattr(self, "LAI"):
            # Sl must also be defined
            ##TODO: add MAXLAI and CWf
            self.Cmax = self.Sl * self.LAI + self.Swood
            self.CanopyGapFraction = pcr.exp(-self.Kext * self.LAI)
            
            self.np_CanopyGapFraction = pcr.pcr2numpy(self.CanopyGapFraction,self.mv)
            
            self.Ewet = (1 - pcr.exp(-self.Kext * self.LAI)) * self.PotenEvap
            self.EoverR = pcr.ifthenelse(
                self.Precipitation > 0.0,
                pcr.min(
                    0.25,
                    pcr.cover(self.Ewet / pcr.max(0.0001, self.Precipitation), 0.0),
                ),
                0.0,
            )
            if hasattr(self, "MAXLAI") and hasattr(self, "CWf"):
                # Adjust rootinggdepth
                self.ActRootingDepth = self.CWf * (
                    self.RootingDepth * self.LAI / pcr.max(0.001, self.MAXLAI)
                ) + ((1 - self.CWf) * self.RootingDepth)
            else:
                self.ActRootingDepth = self.RootingDepth
        else:
            self.ActRootingDepth = self.RootingDepth

        # Apply forcing data corrections
        self.PotenEvap = self.PotenEvap * self.et_RefToPot
        if self.modelSnow:
            self.Temperature = self.Temperature + self.TempCor

        self.wf_multparameters()

        if (self.nrresSimple + self.nrlake) > 0 and self.filterResArea == 0:
            self.ReserVoirPotEvap = self.PotenEvap
            self.ReserVoirPrecip = self.Precipitation

            self.PotenEvap = self.filter_P_PET * self.PotenEvap
            self.Precipitation = self.filter_P_PET * self.Precipitation

        self.OrgStorage = (
            sum_list_cover(self.UStoreLayerDepth, self.ZeroMap) + self.SatWaterDepth
        )
        self.OldCanopyStorage = self.CanopyStorage
        if self.nrpaddyirri > 0:
            self.OldPondingDepth = self.PondingDepth
        self.PotEvap = self.PotenEvap  #

        ##########################################################################
        # Interception according to a modified Gash model
        ##########################################################################
        if self.timestepsecs >= (23 * 3600):
            self.ThroughFall, self.Interception, self.StemFlow, self.CanopyStorage = rainfall_interception_gash(
                self.Cmax,
                self.EoverR,
                self.CanopyGapFraction,
                self.Precipitation,
                self.CanopyStorage,
                maxevap=self.PotEvap,
            )

            self.PotTransSoil = pcr.cover(
                pcr.max(0.0, self.PotEvap - self.Interception), 0.0
            )  # now in mm

        else:
            NetInterception, self.ThroughFall, self.StemFlow, LeftOver, Interception, self.CanopyStorage = rainfall_interception_modrut(
                self.Precipitation,
                self.PotEvap,
                self.CanopyStorage,
                self.CanopyGapFraction,
                self.Cmax,
            )
            self.PotTransSoil = pcr.cover(pcr.max(0.0, LeftOver), 0.0)  # now in mm
            self.Interception = NetInterception

        self.EffectivePrecipitation = self.ThroughFall + self.StemFlow

        if self.modelSnow:
            self.OldSnow = self.Snow
            self.TSoil = self.TSoil + self.w_soil * (self.Temperature - self.TSoil)
            # return Snow,SnowWater,SnowMelt,RainFall
            self.Snow, self.SnowWater, self.SnowMelt, self.RainFallPlusMelt, self.SnowFall = SnowPackHBV(
                self.Snow,
                self.SnowWater,
                self.EffectivePrecipitation,
                self.Temperature,
                self.TTI,
                self.TT,
                self.TTM,
                self.Cfmax,
                self.WHC,
            )
            MaxSnowPack = 10000.0
            if self.MassWasting:
                # Masswasting of dry snow
                # 5.67 = tan 80 graden
                SnowFluxFrac = pcr.min(0.5, self.landSlope / 5.67) * pcr.min(
                    1.0, self.Snow / MaxSnowPack
                )
                MaxFlux = SnowFluxFrac * self.Snow
                self.Snow = pcr.accucapacitystate(self.TopoLdd, self.Snow, MaxFlux)
            else:
                SnowFluxFrac = self.ZeroMap
                MaxFlux = self.ZeroMap

            self.SnowCover = pcr.ifthenelse(self.Snow > 0, pcr.scalar(1), pcr.scalar(0))
            self.NrCell = pcr.areatotal(self.SnowCover, self.TopoId)

            if hasattr(self, "GlacierFrac"):
                """
                Run Glacier module and add the snowpack on-top of it.
                Estimate the fraction of snow turned into ice (HBV-light).
                Estimate glacier melt.
                glacierHBV function in wflow_lib.py
                """

                self.Snow, self.Snow2Glacier, self.GlacierStore, self.GlacierMelt = glacierHBV(
                    self.GlacierFrac,
                    self.GlacierStore,
                    self.Snow,
                    self.Temperature,
                    self.G_TT,
                    self.G_Cfmax,
                    self.G_SIfrac,
                    self.timestepsecs,
                    self.basetimestep
                )
                # Convert to mm per grid cell and add to snowmelt
                self.GlacierMelt = self.GlacierMelt * self.GlacierFrac
                self.RainFallPlusMelt = (
                    self.RainFallPlusMelt + self.GlacierMelt
                )
        else:
            self.RainFallPlusMelt = self.EffectivePrecipitation



        # Start with the soil calculations
        # --------------------------------
        # Code to be able to force zi from the outside
        #
        self.SatWaterDepth = (self.thetaS - self.thetaR) * (
            self.SoilThickness - self.zi
        )

        self.AvailableForInfiltration = (
            self.RainFallPlusMelt + self.IRSupplymm
        )
        self.oldIRSupplymm = self.IRSupplymm
        
        self.UStoreDepth = sum_list_cover(self.UStoreLayerDepth, self.ZeroMap)

        UStoreCapacity = (
            self.SoilWaterCapacity
            - self.SatWaterDepth
            - self.UStoreDepth
        )

        # Runoff from water bodies and river network
        self.RunoffRiverCells = (
            pcr.min(1.0, self.RiverFrac)
            * self.AvailableForInfiltration
        )
        self.RunoffLandCells = (
            pcr.min(1.0, self.WaterFrac)
            * self.AvailableForInfiltration
        )
        
        self.dyn['RunoffLandCells'] = pcr.pcr2numpy(self.RunoffLandCells,self.mv).ravel()
        
        # self.RunoffOpenWater = self.ZeroMap
        self.AvailableForInfiltration = pcr.max(
            self.AvailableForInfiltration - self.RunoffRiverCells - self.RunoffLandCells, 0.0
        )


        # Limit rootingdepth (if set externally)
        self.ActRootingDepth = pcr.min(self.SoilThickness * 0.99, self.ActRootingDepth)       
        self.static['ActRootingDepth'] = pcr.pcr2numpy(self.ActRootingDepth,self.mv).ravel()

        # Determine Open Water EVAP based on riverfrac and waterfrac. Later subtract this from water that
        # enters the Kinematic wave
        self.ActEvapOpenWaterRiver = pcr.min(
            self.WaterLevelR * 1000.0 * self.RiverFrac,
            self.RiverFrac * self.PotTransSoil,
        )
        
        self.ActEvapOpenWaterLand = pcr.min(
            self.WaterLevelL * 1000.0 * self.WaterFrac,
            self.WaterFrac * self.PotTransSoil,
        )
            
        self.dyn['ActEvapOpenWaterLand'] = pcr.pcr2numpy(self.ActEvapOpenWaterLand,self.mv).ravel()
        self.dyn['ActEvapOpenWaterRiver'] = pcr.pcr2numpy(self.ActEvapOpenWaterRiver,self.mv).ravel()

        self.RestEvap = self.PotTransSoil - self.ActEvapOpenWaterRiver - self.ActEvapOpenWaterLand

        self.ActEvapPond = self.ZeroMap
        if self.nrpaddyirri > 0:
            self.ActEvapPond = pcr.min(self.PondingDepth, self.RestEvap)
            self.PondingDepth = self.PondingDepth - self.ActEvapPond
            self.RestEvap = self.RestEvap - self.ActEvapPond

        # evap available for soil evaporation and transpiration
        self.PotSoilEvap = self.RestEvap * self.CanopyGapFraction
        self.PotTrans = self.RestEvap * (1-self.CanopyGapFraction)
        
        self.Inflow = pcr.cover(self.Inflow, self.ZeroMap)
        # only run the reservoir module if needed
        if self.nrresSimple > 0:
            self.ReservoirVolume, self.OutflowSR, self.ResPercFull, self.ResPrecipSR, self.ResEvapSR, self.DemandRelease = simplereservoir(
                self.ReservoirVolume,
                self.RiverRunoff + self.LandRunoff + self.SubsurfaceFlow/1000/1000/1000/self.timestepsecs,
                self.ResSimpleArea,
                self.ResMaxVolume,
                self.ResTargetFullFrac,
                self.ResMaxRelease,
                self.ResDemand,
                self.ResTargetMinFrac,
                self.ReserVoirSimpleLocs,
                self.ReserVoirPrecip,
                self.ReserVoirPotEvap,
                self.ReservoirSimpleAreas,
                timestepsecs=self.timestepsecs,
            )
            self.OutflowDwn = pcr.upstream(
                self.TopoLddOrg, pcr.cover(self.OutflowSR, pcr.scalar(0.0))
            )
            self.Inflow = self.OutflowDwn + self.Inflow
        if self.nrlake > 0:
            self.LakeWaterLevel, self.LakeOutflow, self.LakePrecip, self.LakeEvap, self.LakeVolume = naturalLake(
                self.LakeWaterLevel,
                self.LakeLocs,
                self.LinkedLakeLocs,
                self.LakeArea,
                self.LakeThreshold,
                self.LakeStorFunc,
                self.LakeOutflowFunc,
                self.sh,
                self.hq,
                self.Lake_b,
                self.Lake_e,
                self.RiverRunoff + self.LandRunoff + self.SubsurfaceFlow/1000/1000/1000/self.timestepsecs,
                self.ReserVoirPrecip,
                self.ReserVoirPotEvap,
                self.LakeAreasMap,
                self.wf_supplyJulianDOY(),
                timestepsecs=self.timestepsecs,
            )
            self.OutflowDwn = pcr.upstream(
                self.TopoLddOrg, pcr.cover(self.LakeOutflow, pcr.scalar(0.0))
            )
            self.Inflow = self.OutflowDwn + self.Inflow
        
        
        # convert to numpy for numba        
        self.dyn['AvailableForInfiltration'] = pcr.pcr2numpy(self.AvailableForInfiltration, self.mv).ravel()
        self.dyn['zi'] = pcr.pcr2numpy(self.zi,self.mv).ravel()
        self.dyn['SatWaterDepth'] = pcr.pcr2numpy(self.SatWaterDepth,self.mv).ravel()
        self.dyn['PotSoilEvap'] = pcr.pcr2numpy(self.PotSoilEvap,self.mv).ravel()
        self.dyn['PotTrans'] = pcr.pcr2numpy(self.PotTrans,self.mv).ravel()
        self.dyn['TSoil'] = pcr.pcr2numpy(self.TSoil, self.mv).ravel()
        self.dyn['ssf'] = pcr.pcr2numpy(self.SubsurfaceFlow,self.mv).ravel()
        self.dyn['WaterLevelL'] = pcr.pcr2numpy(self.WaterLevelL,self.mv).ravel()
        self.dyn['Qo_in'] = self.dyn['Qo_in'] * 0.0
        
        if self.nrpaddyirri > 0:
            self.dyn['PondingDepth'] = pcr.pcr2numpy(self.PondingDepth,self.mv).ravel()
    
        for i in range(self.maxLayers):
            self.layer['UStoreLayerDepth'][i] = pcr.pcr2numpy(self.UStoreLayerDepth[i],self.mv).ravel()
        
        self.dyn['LandRunoff'] = pcr.pcr2numpy(self.LandRunoff,self.mv).ravel()
        self.dyn['sumUStoreLayerDepth'] = pcr.pcr2numpy(self.UStoreDepth,self.mv).ravel()

        it_kinL = 1
        if self.kinwaveIters == 1:
            if self.kinwaveLandTstep == 0:
                it_kinL = estimate_iterations_kin_wave(self.LandRunoff, self.Beta, self.AlphaL, self.timestepsecs, self.DL, self.mv)
            else:
                it_kinL = int(np.ceil(self.timestepsecs/self.kinwaveLandTstep))
                
        ssf, qo, self.dyn, self.layer  = sbm_cell(self.nodes, 
                                             self.nodes_up,
                                             self.np_ldd.ravel(),
                                             self.layer,
                                             self.static,
                                             self.dyn,
                                             self.modelSnow,
                                             self.soilInfReduction,
                                             self.timestepsecs, 
                                             self.basetimestep,
                                             1.0,
                                             self.nrpaddyirri,
                                             self.shape,
                                             self.TransferMethod,
                                             it_kinL,
                                             self.UST
                                             )

        self.SubsurfaceFlow = pcr.numpy2pcr(pcr.Scalar,ssf.reshape(self.shape),self.mv)
        self.LandRunoff = pcr.numpy2pcr(pcr.Scalar,qo.reshape(self.shape),self.mv)
        self.zi = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['zi'].reshape(self.shape)),self.mv)
        self.SatWaterDepth = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['SatWaterDepth'].reshape(self.shape)),self.mv)
        self.UStoreDepth = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['sumUStoreLayerDepth'].reshape(self.shape)),self.mv)
        self.CapFlux = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['CapFlux'].reshape(self.shape)),self.mv)
        self.Transfer = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['Transfer'].reshape(self.shape)),self.mv)
        
        for i in range(self.maxLayers):
            self.UStoreLayerDepth[i] = pcr.numpy2pcr(pcr.Scalar,np.copy(self.layer['UStoreLayerDepth'][i].reshape(self.shape)),self.mv)
                    
        if self.nrpaddyirri > 0:
            self.PondingDepth = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['PondingDepth'].reshape(self.shape)),self.mv)
        
        # Determine transpiration
        self.ActEvapUStore = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ActEvapUStore'].reshape(self.shape)),self.mv)
        self.ActEvapSat = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ActEvapSat'].reshape(self.shape)),self.mv)
        self.Transpiration = self.ActEvapUStore + self.ActEvapSat
        
        self.soilevap = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['soilevap'].reshape(self.shape)),self.mv)
        
        self.ActEvap = (
                self.soilevap
                + self.Transpiration
                + self.ActEvapOpenWaterRiver
                + self.ActEvapOpenWaterLand
                + self.ActEvapPond
                )
        
        self.Recharge = self.Transfer - self.CapFlux - self.ActEvapSat

        # Run only if we have irrigation areas or an externally given demand, determine irrigation demand based on potrans and acttrans
        if self.nrirri > 0 or hasattr(self, "IrriDemandExternal"):
            if not hasattr(self, "IrriDemandExternal"):  # if not given
                self.IrriDemand, self.IrriDemandm3 = self.irrigationdemand(
                    self.PotTrans, self.Transpiration, self.IrrigationAreas
                )
                IRDemand = (
                    idtoid(
                        self.IrrigationAreas,
                        self.IrrigationSurfaceIntakes,
                        self.IrriDemandm3,
                    )
                    * -1.0
                )
            else:
                IRDemand = self.IrriDemandExternal
            # loop over irrigation areas and assign Q to linked river extraction points
            self.Inflow = pcr.cover(IRDemand, self.Inflow)
        
                
        if self.nrpaddyirri > 0:
            irr_depth = (
                pcr.ifthenelse(
                    self.PondingDepth < self.h_min, self.h_max - self.PondingDepth, 0.0
                )
                * self.CRPST
            )
            sqmarea = pcr.areatotal(
                self.reallength * self.reallength, self.IrrigationPaddyAreas
            )
            self.IrriDemandm3 = pcr.cover((irr_depth / 1000.0) * sqmarea, 0)
            IRDemand = idtoid(
                self.IrrigationPaddyAreas,
                self.IrrigationSurfaceIntakes,
                self.IrriDemandm3,
            ) * (-1.0 / self.timestepsecs)

            self.IRDemand = IRDemand
            self.Inflow = pcr.cover(IRDemand, self.Inflow)
            self.irr_depth = irr_depth


        SurfaceWater = self.WaterLevelR / 1000.0  # SurfaceWater (mm)
        self.CumSurfaceWater = self.CumSurfaceWater + SurfaceWater
        
        self.InwaterMM = self.RunoffRiverCells - self.ActEvapOpenWaterRiver
        
        self.qo_toriver = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['qo_toriver'].reshape(self.shape)),self.mv)
        self.ssf_toriver = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ssf_toriver'].reshape(self.shape)),self.mv)     
        self.Inwater = self.InwaterMM * self.ToCubic + self.qo_toriver + self.ssf_toriver # m3/s

        self.ExfiltWater = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ExfiltWater'].reshape(self.shape)),self.mv)
        self.InfiltExcess = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['InfiltExcess'].reshape(self.shape)),self.mv)
        self.ExcessWater = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ExcessWater'].reshape(self.shape)),self.mv)
        self.ActInfilt = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ActInfilt'].reshape(self.shape)),self.mv)
        self.ActLeakage = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ActLeakage'].reshape(self.shape)),self.mv)
        self.soilevap = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['soilevap'].reshape(self.shape)),self.mv)
        
        self.SoilInf = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['SoilInf'].reshape(self.shape)),self.mv)
        self.PathInf = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['PathInf'].reshape(self.shape)),self.mv)
        self.InfiltExcessSoil = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ExcessWaterSoil'].reshape(self.shape)),self.mv)
        self.InfiltExcessPath = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ExcessWaterPath'].reshape(self.shape)),self.mv)
        self.ActInfiltSoil = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ActInfiltSoil'].reshape(self.shape)),self.mv)
        self.ActInfiltPath = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ActInfiltPath'].reshape(self.shape)),self.mv)
        
        self.ExfiltFromUstore = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ExfiltFromUstore'].reshape(self.shape)),self.mv)
        self.ExfiltFromSat = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ExfiltSatWater'].reshape(self.shape)),self.mv)
        self.SumEvapSat = (pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['soilevapsat'].reshape(self.shape)),self.mv) +
                           pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ActEvapSat'].reshape(self.shape)),self.mv))
        self.SumEvapUstore = (pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['soilevapunsat'].reshape(self.shape)),self.mv) +
                           pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['ActEvapUStore'].reshape(self.shape)),self.mv))
          
        self.ExfiltWaterCubic = self.ExfiltWater * self.ToCubic
        self.InfiltExcessCubic = self.InfiltExcess * self.ToCubic


        # Check if we do not try to abstract more runoff then present
        self.InflowKinWaveCell = pcr.upstream(
            self.TopoLdd, self.RiverRunoff
        )  # NG The extraction should be equal to the discharge upstream cell. You should not make the abstraction depended on the downstream cell, because they are correlated. During a stationary sum they will get equal to each other.
        MaxExtract = self.InflowKinWaveCell + self.Inwater  # NG
        # MaxExtract = self.SurfaceRunoff + self.Inwater
        self.SurfaceWaterSupply = pcr.ifthenelse(
            self.Inflow < 0.0, pcr.min(MaxExtract, -1.0 * self.Inflow), self.ZeroMap
        )
        self.OldRiverRunoff = self.RiverRunoff  # NG Store for iteration
        OldRiverRunoff = pcr.pcr2numpy(self.OldRiverRunoff,self.mv).ravel()
        self.OldInwater = self.Inwater
        self.Inwater = self.Inwater + pcr.ifthenelse(
            self.SurfaceWaterSupply > 0, -1.0 * self.SurfaceWaterSupply, self.Inflow
        )

        ##########################################################################
        # Runoff calculation via Kinematic wave ##################################
        ##########################################################################
        
        qr = self.Inwater/self.DCL     
        qr_np =  pcr.pcr2numpy(qr,self.mv).ravel()
        
        RiverRunoff = pcr.pcr2numpy(self.RiverRunoff,self.mv).ravel()
        
        it_kinR = 1
        if self.kinwaveIters == 1:
            if self.kinwaveRiverTstep == 0:
                it_kinR = estimate_iterations_kin_wave(self.RiverRunoff, self.Beta, self.AlphaR, self.timestepsecs, self.DCL, self.mv)
            else:
                it_kinR = int(np.ceil(self.timestepsecs/self.kinwaveRiverTstep))
                
        acc_flow = kin_wave(
                self.rnodes,
                self.rnodes_up,
                RiverRunoff,
                qr_np,
                self.dyn['AlphaR'],
                self.static['Beta'],
                self.static['DCL'],
                self.static['River'],
                self.static['Bw'],
                self.static['AlpTermR'],
                self.static['AlpPow'],
                self.timestepsecs,
                it_kinR)
            
        Qriver = acc_flow/self.timestepsecs
        self.RiverRunoff = pcr.numpy2pcr(pcr.Scalar, np.copy(Qriver).reshape(self.shape),self.mv)


        # If inflow is negative we have abstractions. Check if demand can be met (by looking
        # at the flow in the upstream cell) and iterate if needed
        self.nrit = 0
        self.breakoff = 0.0001
        if float(pcr.mapminimum(pcr.spatial(self.Inflow))) < 0.0:
            while True:
                self.nrit += 1
                oldsup = self.SurfaceWaterSupply
                self.InflowKinWaveCell = pcr.upstream(self.TopoLdd, self.RiverRunoff)
                ##########################################################################
                # Iterate to make a better estimation for the supply #####################
                # (Runoff calculation via Kinematic wave) ################################
                ##########################################################################
                MaxExtract = self.InflowKinWaveCell + self.OldInwater
                self.SurfaceWaterSupply = pcr.ifthenelse(
                    self.Inflow < 0.0,
                    pcr.min(MaxExtract, -1.0 * self.Inflow),
                    self.ZeroMap,
                )
                # Fraction of demand that is not used but flows back into the river get fracttion and move to return locations
                self.DemandReturnFlow = pcr.cover(
                    idtoid(
                        self.IrrigationSurfaceIntakes,
                        self.IrrigationSurfaceReturn,
                        self.DemandReturnFlowFraction * self.SurfaceWaterSupply,
                    ),
                    0.0,
                )

                self.Inwater = (
                    self.OldInwater
                    + pcr.ifthenelse(
                        self.SurfaceWaterSupply > 0,
                        -1.0 * self.SurfaceWaterSupply,
                        self.Inflow,
                    )
                    + self.DemandReturnFlow
                )
                # per distance along stream
                q = self.Inwater / self.DCL
                q_np = pcr.pcr2numpy(q, self.mv).ravel()

                if self.kinwaveIters == 1:
                    if self.kinwaveRiverTstep == 0:
                        it_kinR = estimate_iterations_kin_wave(self.OldRiverRunoff, self.Beta, self.AlphaR, self.timestepsecs, self.DCL, self.mv)
                    else:
                        it_kinR = int(np.ceil(self.timestepsecs/self.kinwaveRiverTstep))
                                
                acc_flow = kin_wave(
                        self.rnodes,
                        self.rnodes_up,
                        OldRiverRunoff,
                        q_np,
                        self.dyn['AlphaR'],
                        self.static['Beta'],
                        self.static['DCL'],
                        self.static['River'],
                        self.static['Bw'],
                        self.static['AlpTermR'],
                        self.static['AlpPow'],
                        self.timestepsecs,
                        it_kinR)
                    
                Qriver = acc_flow/self.timestepsecs
                self.RiverRunoff = pcr.numpy2pcr(pcr.Scalar, np.copy(Qriver).reshape(self.shape),self.mv)
                                
                self.RiverRunoffMM = (
                    self.RiverRunoff * self.QMMConv
                )  # RiverRunoffMM (mm) from RiverRunoff (m3/s)

                self.InflowKinWaveCell = pcr.upstream(
                    self.TopoLdd, self.OldRiverRunoff
                )
                deltasup = float(pcr.mapmaximum(pcr.abs(oldsup - self.SurfaceWaterSupply)))

                if deltasup < self.breakoff or self.nrit >= self.maxitsupply:
                    break

            self.InflowKinWaveCell = pcr.upstream(self.TopoLdd, self.RiverRunoff)
            self.updateRunOff()
        else:
            self.RiverRunoffMM = (
                self.RiverRunoff * self.QMMConv
            )  # RiverRunoffMM (mm) from RiverRunoff (m3/s)
            self.updateRunOff()
            
        self.LandRunoffMM = (
            self.LandRunoff * self.QMMConv
        )  # LandRunoffMM (mm) from LandRunoff (m3/s)

        # Now add the supply that is linked to irrigation areas to extra precip
        if self.nrirri > 0:
            # loop over irrigation areas and spread-out the supply over the area
            IRSupplymm = idtoid(
                self.IrrigationSurfaceIntakes,
                self.IrrigationAreas,
                self.SurfaceWaterSupply * (1 - self.DemandReturnFlowFraction),
            )
            sqmarea = pcr.areatotal(
                self.reallength * self.reallength, pcr.nominal(self.IrrigationAreas)
            )

            self.IRSupplymm = pcr.cover(
                IRSupplymm / (sqmarea / 1000.0 / self.timestepsecs), 0.0
            )

        if self.nrpaddyirri > 0:
            # loop over irrigation areas and spread-out the supply over the area
            IRSupplymm = idtoid(
                self.IrrigationSurfaceIntakes,
                pcr.ifthen(self.IrriDemandm3 > 0, self.IrrigationPaddyAreas),
                self.SurfaceWaterSupply,
            )
            sqmarea = pcr.areatotal(
                self.reallength * self.reallength,
                pcr.nominal(
                    pcr.ifthen(self.IrriDemandm3 > 0, self.IrrigationPaddyAreas)
                ),
            )

            self.IRSupplymm = pcr.cover(
                ((IRSupplymm * self.timestepsecs * 1000) / sqmarea), 0.0
            )
            
        #Inflow in kin wave is updated with the new surface runoff of the current timestep
        self.InflowKinWaveCell = pcr.upstream(
            self.TopoLdd, self.RiverRunoff
        )


        self.MassBalKinWaveR = (
            (-self.KinWaveVolumeR + self.OldKinWaveVolumeR) / self.timestepsecs
            + self.InflowKinWaveCell
            + self.Inwater
            - self.RiverRunoff
        )
        

        self.InwaterL = pcr.numpy2pcr(pcr.Scalar, np.copy(self.dyn['InwaterO'].reshape(self.shape)), self.mv)
        self.InflowKinWaveCellLand = pcr.numpy2pcr(pcr.Scalar, np.copy(self.dyn['Qo_in'].reshape(self.shape)), self.mv)
                                
        self.MassBalKinWaveL = (
            (-self.KinWaveVolumeL + self.OldKinWaveVolumeL) / self.timestepsecs
            + self.InflowKinWaveCellLand
            + self.InwaterL
            - self.LandRunoff
        )

        Runoff = self.RiverRunoff

        # Updating
        # --------
        # Assume a tss file with as many columns as outputlocs. Start updating for each non-missing value and start with the
        # first column (nr 1). Assumes that outputloc and columns match!

        if self.updating:
            self.QM = (
                pcr.timeinputscalar(self.updateFile, self.UpdateMap) * self.QMMConv
            )

            # Now update the state. Just add to the Ustore
            # self.UStoreDepth =  result
            # No determine multiplication ratio for each gauge influence area.
            # For missing gauges 1.0 is assumed (no change).
            # UpDiff = pcr.areamaximum(QM,  self.UpdateMap) - pcr.areamaximum(self.SurfaceRunoffMM, self.UpdateMap)
            UpRatio = pcr.areamaximum(self.QM, self.UpdateMap) / pcr.areamaximum(
                self.RiverRunoffMM, self.UpdateMap
            )

            UpRatio = pcr.cover(pcr.areaaverage(UpRatio, self.TopoId), 1.0)
            # Now split between Soil and Kyn  wave
            self.UpRatioKyn = pcr.min(
                self.MaxUpdMult,
                pcr.max(self.MinUpdMult, (UpRatio - 1.0) * self.UpFrac + 1.0),
            )
            UpRatioSoil = pcr.min(
                self.MaxUpdMult,
                pcr.max(self.MinUpdMult, (UpRatio - 1.0) * (1.0 - self.UpFrac) + 1.0),
            )

            # update/nudge self.UStoreDepth for the whole upstream area,
            # not sure how much this helps or worsens things
            # TODO: FIx this for multiple layers
            UpdSoil = True
            if UpdSoil:
                toadd = (self.UStoreDepth * UpRatioSoil) - self.UStoreDepth
                self.UStoreDepth = self.UStoreDepth + toadd

            # Update the kinematic wave reservoir up to a maximum upstream distance
            MM = (1.0 - self.UpRatioKyn) / self.UpdMaxDist
            self.UpRatioKyn = MM * self.DistToUpdPt + self.UpRatioKyn
            self.RiverRunoff = self.RiverRunoff * self.UpRatioKyn
            self.RiverRunoffMM = (
                self.RiverRunoff * self.QMMConv
            )  # SurfaceRunoffMM (mm) from SurfaceRunoff (m3/s)
            self.updateRunOff()
            Runoff = self.RiverRunoff

        # Determine Soil moisture profile
        # self.vwc, self.vwcRoot: volumetric water content [m3/m3] per soil layer and root zone (including thetaR and saturated store)
        # self.vwc_perc, self.vwc_percRoot: volumetric water content [%] per soil layer and root zone (including thetaR and saturated store)
        # self.RootStore_sat: root water storage [mm] in saturated store (excluding thetaR)
        # self.RootStore_unsat: root water storage [mm] in unsaturated store (excluding thetaR)
        # self.RootStore: total root water storage [mm] (excluding thetaR)
        
        self.vwc = []
        self.vwc_perc = []
        self.RootStore_sat = pcr.max(0.0, self.ActRootingDepth - self.zi) * (
            self.thetaS - self.thetaR
        )
        self.RootStore_unsat = pcr.numpy2pcr(pcr.Scalar, np.copy(self.dyn['RootStore_unsat'].reshape(self.shape)), self.mv)
         
        for i in range(self.maxLayers):
            self.vwc.append(pcr.numpy2pcr(pcr.Scalar,np.copy(self.layer['vwc'][i].reshape(self.shape)),self.mv))
            self.vwc_perc.append(pcr.numpy2pcr(pcr.Scalar,np.copy(self.layer['vwc_perc'][i].reshape(self.shape)),self.mv))
            
        self.RootStore = self.RootStore_sat + self.RootStore_unsat
        self.vwcRoot = self.RootStore / self.ActRootingDepth + self.thetaR
        self.vwc_percRoot = (self.vwcRoot / self.thetaS) * 100.0
        

        # 2:
        ##########################################################################
        # water balance ###########################################
        ##########################################################################

        self.QCatchmentMM = self.RiverRunoff * self.QMMConvUp
        self.RunoffCoeff = (
            self.QCatchmentMM
            / pcr.catchmenttotal(self.RainFallPlusMelt, self.TopoLdd)
            / pcr.catchmenttotal(pcr.cover(1.0), self.TopoLdd)
        )
        # Single cell based water budget. snow not included yet.

        self.CellStorage = (
            sum_list_cover(self.UStoreLayerDepth, self.ZeroMap) + self.SatWaterDepth
        )

        self.sumUstore = sum_list_cover(self.UStoreLayerDepth, self.ZeroMap)

        self.DeltaStorage = self.CellStorage - self.OrgStorage
        self.SatWaterFlux = self.SubsurfaceFlow /(self.xl*1000*self.yl*1000)
        OutFlow = self.SatWaterFlux
            
        CellInFlow = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['CellInFlow'].reshape(self.shape)),self.mv)
        self.CellInFlow = CellInFlow

        self.CumOutFlow = self.CumOutFlow + OutFlow
        self.CumActInfilt = self.CumActInfilt + self.ActInfilt
        self.CumCellInFlow = self.CumCellInFlow + CellInFlow
        self.CumPrec = self.CumPrec + self.Precipitation
        self.CumEvap = self.CumEvap + self.ActEvap
        self.CumPotenEvap = self.CumPotenEvap + self.PotenEvap

        self.CumInt = self.CumInt + self.Interception

        self.SnowCover = pcr.ifthenelse(
            self.Snow > 0.0, self.ZeroMap + 1.0, self.ZeroMap
        )
        self.CumLeakage = self.CumLeakage + self.ActLeakage
        self.CumExfiltWater = self.CumExfiltWater + self.ExfiltWater
        
        
        self.SoilWatbal = pcr.numpy2pcr(pcr.Scalar,np.copy(self.dyn['SoilWatbal'].reshape(self.shape)),self.mv)
        
        self.InterceptionWatBal = (
            self.Precipitation
            - self.Interception
            - self.StemFlow
            - self.ThroughFall
            - (self.OldCanopyStorage - self.CanopyStorage)
        )
        
        
        self.SurfaceWatbal = (
            self.RainFallPlusMelt
            + self.oldIRSupplymm
            - self.InfiltExcess
            - self.ExcessWater
            - self.RunoffRiverCells 
            - self.RunoffLandCells
            - self.ActInfilt
            - (self.OldCanopyStorage - self.CanopyStorage)
        )

        self.watbal = self.SoilWatbal + self.SurfaceWatbal
        
        self.count = self.count + 1
        

def main(argv=None):
    """
    Perform command line execution of the model.
    """
    caseName = "default_sbm"
    global multpars
    runId = "run_default"
    configfile = "wflow_sbm.ini"
    _lastTimeStep = 0
    _firstTimeStep = 0
    LogFileName = "wflow.log"

    runinfoFile = "runinfo.xml"
    timestepsecs = 86400
    wflow_cloneMap = "wflow_subcatch.map"
    _NoOverWrite = 1
    global updateCols
    loglevel = logging.DEBUG

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return
    ########################################################################
    ## Process command-line options                                        #
    ########################################################################
    try:
        opts, args = getopt.getopt(argv, "XL:hC:Ii:v:S:T:WR:u:s:EP:p:Xx:U:fOc:l:")
    except getopt.error as msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == "-C":
            caseName = a
        if o == "-R":
            runId = a
        if o == "-c":
            configfile = a
        if o == "-L":
            LogFileName = a
        if o == "-s":
            timestepsecs = int(a)
        if o == "-h":
            usage()
        if o == "-f":
            _NoOverWrite = 0
        if o == "-l":
            exec("loglevel = logging." + a)

    starttime = dt.datetime(1990, 1, 1)

    if _lastTimeStep < _firstTimeStep:
        print(
            "The starttimestep ("
            + str(_firstTimeStep)
            + ") is smaller than the last timestep ("
            + str(_lastTimeStep)
            + ")"
        )
        usage()

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(
        myModel, _lastTimeStep, firstTimestep=_firstTimeStep, datetimestart=starttime
    )
    dynModelFw.createRunId(
        NoOverWrite=_NoOverWrite,
        level=loglevel,
        logfname=LogFileName,
        model="wflow_sbm",
        doSetupFramework=False,
    )

    for o, a in opts:
        if o == "-P":
            left = a.split("=")[0]
            right = a.split("=")[1]
            configset(
                myModel.config, "variable_change_once", left, right, overwrite=True
            )
        if o == "-p":
            left = a.split("=")[0]
            right = a.split("=")[1]
            configset(
                myModel.config, "variable_change_timestep", left, right, overwrite=True
            )
        if o == "-X":
            configset(myModel.config, "model", "OverWriteInit", "1", overwrite=True)
        if o == "-I":
            configset(myModel.config, "run", "reinit", "1", overwrite=True)
        if o == "-i":
            configset(myModel.config, "model", "intbl", a, overwrite=True)
        if o == "-s":
            configset(myModel.config, "run", "timestepsecs", a, overwrite=True)
        if o == "-x":
            configset(myModel.config, "model", "sCatch", a, overwrite=True)
        if o == "-c":
            configset(myModel.config, "model", "configfile", a, overwrite=True)
        if o == "-M":
            configset(myModel.config, "model", "MassWasting", "0", overwrite=True)
        if o == "-R":
            runId = a
        if o == "-T":
            configset(myModel.config, "run", "endtime", a, overwrite=True)
        if o == "-S":
            configset(myModel.config, "run", "starttime", a, overwrite=True)

    dynModelFw.setupFramework()
    dynModelFw._runInitial()
    dynModelFw._runResume()
    # dynModelFw._runDynamic(0, 0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
