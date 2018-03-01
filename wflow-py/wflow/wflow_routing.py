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
Run the wflow_routing model..

usage

::

    wflow_routing [-h][-v level][-F runinfofile][-L logfile][-C casename][-R runId]
          [-c configfile][-T last_step][-S first_step][-s seconds][-l loglevel]


    -X: save state at the end of the run over the initial conditions at the start

    -T: Set end time of the run: yyyy-mm-dd hh:mm:ss

    -S: Set start time of the run: yyyy-mm-dd hh:mm:ss

    -s: Set the model timesteps in seconds

    -I: re-initialize the initial model conditions with default

    -C: set the name  of the case (directory) to run

    -R: set the name runId within the current case

    -L: set the logfile

    -c: name of wflow the configuration file (default: Casename/wflow_routing.ini).

    -h: print usage information

    -P: set parameter change string (e.g: -P 'self.FC = self.FC * 1.6') for non-dynamic variables

    -p: set parameter change string (e.g: -P 'self.Precipitation = self.Precipitation * 1.11') for
        dynamic variables

    -l: loglevel (must be one of DEBUG, WARNING, ERROR)


"""

import numpy
import os
import os.path
import shutil, glob
import getopt


from wflow.wf_DynamicFramework import *
from wflow.wflow_funcs import *
from wflow.wflow_adapt import *
import ConfigParser


wflow = "wflow_routing: "

updateCols = []

def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)


class WflowModel(DynamicModel):
    """
    .. versionchanged:: 0.1
        - initial version.

  """


    def __init__(self, cloneMap, Dir, RunDir, configfile):
        DynamicModel.__init__(self)

        self.caseName = os.path.abspath(Dir)
        self.clonemappath = os.path.join(os.path.abspath(Dir),"staticmaps",cloneMap)
        setclone(self.clonemappath)
        self.runId = RunDir
        self.Dir = os.path.abspath(Dir)
        self.configfile = configfile
        self.SaveDir = os.path.join(self.Dir,self.runId)


    def wetPerimiterFP(self,Waterlevel, floodplainwidth,threshold=0.0,sharpness=0.5):
        """

        :param Waterlevel:
        :param bottomwidth:
        :param bankfull:
        :param floodplainwidth:
        :return P: wetted perimiter
        """
        a = threshold
        b = 1.0
        c = sharpness  # not very sharp

        floodplainfact = max(0.001,sCurve(Waterlevel, a=a, c=c,b=b) -0.5)
        floodplainperimiter = min(1.0,2.0 * floodplainfact) * floodplainwidth

        return floodplainperimiter


    def wetPerimiterCH(self,Waterlevel,channelWidth):
        """

        :param Waterlevel:
        :param bottomwidth:
        :param bankfull:
        :param floodplainwidth:
        :return P: wetted perimiter
        """

        channelperimiter = 2.0 * Waterlevel + channelWidth

        return  channelperimiter



    def updateRunOff(self):
        """
        Updates the kinematic wave reservoir water level. Should be run after updates to Q

        WL = A * pow(Q,B)/Bw
        WL/A * Bw = pow(Q,B)
        Q = pow(WL/A * Bw,1/B)
        """

        self.Qbankfull = pow(self.bankFull/self.AlphaCh * self.Bw,1.0/self.Beta)
        self.Qchannel = min(self.SurfaceRunoff,self.Qbankfull)
        self.floodcells  = boolean(ifthenelse(self.WaterLevelCH > self.bankFull, boolean(1), boolean(0)))
        self.Qfloodplain = max(0.0,self.SurfaceRunoff - self.Qbankfull)

        self.WaterLevelCH = self.AlphaCh * pow(self.Qchannel, self.Beta) / (self.Bw)
        self.WaterLevelFP = ifthenelse(self.River,self.AlphaFP * pow(self.Qfloodplain, self.Beta) / (self.Bw + self.Pfp),0.0)
        self.WaterLevel = self.WaterLevelCH + self.WaterLevelFP

        # Determine Qtot as a check
        self.Qtot = pow(self.WaterLevelCH/self.AlphaCh * self.Bw,1.0/self.Beta) + pow(self.WaterLevelFP/self.AlphaFP * (self.Pfp + self.Bw),1.0/self.Beta)
        # wetted perimeter (m)
        self.Pch = self.wetPerimiterCH(self.WaterLevelCH,self.Bw)
        self.Pfp = ifthenelse(self.River,self.wetPerimiterFP(self.WaterLevelFP,self.floodPlainWidth,sharpness=self.floodPlainDist),0.0)

        # Alpha
        self.WetPComb = self.Pch + self.Pfp
        self.Ncombined = (self.Pch/self.WetPComb*self.N**1.5 + self.Pfp/self.WetPComb*self.NFloodPlain**1.5)**(2./3.)
        self.AlpTermFP = pow((self.NFloodPlain / (sqrt(self.SlopeDCL))), self.Beta)
        self.AlpTermComb = pow((self.Ncombined / (sqrt(self.SlopeDCL))), self.Beta)
        self.AlphaFP = self.AlpTermFP * pow(self.Pfp, self.AlpPow)
        self.AlphaCh = self.AlpTerm * pow(self.Pch, self.AlpPow)
        self.Alpha = ifthenelse(self.River,self.AlpTermComb * pow(self.Pch + self.Pfp, self.AlpPow),self.AlphaCh)
        self.OldKinWaveVolume = self.KinWaveVolume
        self.KinWaveVolume = (self.WaterLevelCH * self.Bw * self.DCL) + (self.WaterLevelFP * (self.Pfp + self.Bw) * self.DCL)


    def stateVariables(self):
        """
        returns a list of state variables that are essential to the model.
        This list is essential for the resume and suspend functions to work.

        This function is specific for each model and **must** be present.

        :var self.SurfaceRunoff: Surface runoff in the kin-wave resrvoir [m^3/s]
        :var self.WaterLevel: Water level in the kin-wave resrvoir [m]
        """
        states = ['SurfaceRunoff', 'WaterLevelCH','WaterLevelFP','ReservoirVolume']

        return states

    def supplyCurrentTime(self):
        """
        gets the current time in seconds after the start of the run
        """
        return self.currentTimeStep() * self.timestepsecs


    def suspend(self):

        self.logger.info("Saving initial conditions...")
        self.wf_suspend(os.path.join(self.SaveDir,"outstate"))

        if self.OverWriteInit:
            self.logger.info("Saving initial conditions over start conditions...")
            self.wf_suspend(self.SaveDir + "/instate/")





    def initial(self):
        """
    Initial part of the model, executed only once. Reads all static data from disk


    *Surface water*

    :var N.tbl: Manning's N parameter
    :var N_river.tbl: Manning's N parameter fro cells marked as river

    """
        global statistics
        global multpars
        global updateCols

        self.thestep = scalar(0)
        self.basetimestep = 86400
        self.SSSF = False
        setglobaloption("unittrue")

        self.inflowTss = "/intss/Inflow.tss"
        self.logger.info("running for " + str(self.nrTimeSteps()) + " timesteps")

        # Set and get defaults from ConfigFile here ###################################
        self.maxitsupply = int(configget(self.config, "model", "maxitsupply", "5"))
        # max number of iteration in abstraction calculations
        self.reinit = int(configget(self.config, "run", "reinit", "0"))
        self.OverWriteInit = int(configget(self.config, "model", "OverWriteInit", "0"))
        self.updating = int(configget(self.config, "model", "updating", "0"))
        self.updateFile = configget(self.config, "model", "updateFile", "no_set")
        self.Tslice = int(configget(self.config, "model", "Tslice", "1"))
        self.sCatch = int(configget(self.config, "model", "sCatch", "0"))
        self.intbl = configget(self.config, "model", "intbl", "intbl")
        self.timestepsecs = int(configget(self.config, "model", "timestepsecs", "86400"))
        sizeinmetres = int(configget(self.config, "layout", "sizeinmetres", "0"))
        alf = float(configget(self.config, "model", "Alpha", "60"))
        Qmax = float(configget(self.config, "model", "AnnualDischarge", "300"))
        self.UpdMaxDist = float(configget(self.config, "model", "UpdMaxDist", "100"))

        self.MaxUpdMult = float(configget(self.config, "model", "MaxUpdMult", "1.3"))
        self.MinUpdMult = float(configget(self.config, "model", "MinUpdMult", "0.7"))
        self.UpFrac = float(configget(self.config, "model", "UpFrac", "0.8"))
        self.SubCatchFlowOnly = int(configget(self.config, 'model', 'SubCatchFlowOnly', '0'))

        WIMaxScale = float(configget(self.config, 'model', 'WIMaxScale', '0.8'))

        # static maps to use (normally default)
        wflow_subcatch = configget(self.config, "model", "wflow_subcatch", "staticmaps/wflow_subcatch.map")
        wflow_dem = configget(self.config, "model", "wflow_dem", "staticmaps/wflow_dem.map")
        wflow_ldd = configget(self.config, "model", "wflow_ldd", "staticmaps/wflow_ldd.map")
        wflow_river = configget(self.config, "model", "wflow_river", "staticmaps/wflow_river.map")
        wflow_riverlength = configget(self.config, "model", "wflow_riverlength", "staticmaps/wflow_riverlength.map")
        wflow_riverlength_fact = configget(self.config, "model", "wflow_riverlength_fact",
                                           "staticmaps/wflow_riverlength_fact.map")
        wflow_gauges = configget(self.config, "model", "wflow_gauges", "staticmaps/wflow_gauges.map")
        wflow_inflow = configget(self.config, "model", "wflow_inflow", "staticmaps/wflow_inflow.map")
        wflow_riverwidth = configget(self.config, "model", "wflow_riverwidth", "staticmaps/wflow_riverwidth.map")
        wflow_floodplainwidth = configget(self.config, "model", "wflow_floodplainwidth", "staticmaps/wflow_floodplainwidth.map")
        wflow_bankfulldepth = configget(self.config, "model", "wflow_bankfulldepth", "staticmaps/wflow_bankfulldepth.map")
        wflow_floodplaindist = configget(self.config, "model", "wflow_floodplaindist", "staticmaps/wflow_floodplaindist.map")

        wflow_landuse = configget(self.config, "model", "wflow_landuse", "staticmaps/wflow_landuse.map")
        wflow_soil = configget(self.config, "model", "wflow_soil", "staticmaps/wflow_soil.map")

        # 2: Input base maps ########################################################
        self.instate = configget(self.config,"model","instate","instate")

        subcatch = ordinal(self.wf_readmap(os.path.join(self.Dir,wflow_subcatch),0.0,fail=True))  # Determines the area of calculations (all cells > 0)
        subcatch = ifthen(subcatch > 0, subcatch)

        self.Altitude = self.wf_readmap(os.path.join(self.Dir,wflow_dem),0.0,fail=True)  # * scalar(defined(subcatch)) # DEM
        self.TopoLdd = self.wf_readmap(os.path.join(self.Dir,wflow_ldd),0.0,fail=True)  # Local
        self.TopoId = self.wf_readmap(os.path.join(self.Dir,wflow_subcatch),0.0,fail=True)  # area map
        self.River = cover(boolean(self.wf_readmap(os.path.join(self.Dir,wflow_river),0.0,fail=True)), 0)

        self.RiverLength = cover(self.wf_readmap(os.path.join(self.Dir,wflow_riverlength), 0.0), 0.0)
        # Factor to multiply riverlength with (defaults to 1.0)
        self.RiverLengthFac = self.wf_readmap(os.path.join(self.Dir,wflow_riverlength_fact), 1.0)

        # read landuse and soilmap and make sure there are no missing points related to the
        # subcatchment map. Currently sets the lu and soil type  type to 1
        self.LandUse = ordinal(self.wf_readmap(os.path.join(self.Dir,wflow_landuse),0.0,fail=True))
        self.LandUse = cover(self.LandUse, ordinal(subcatch > 0))
        self.Soil = ordinal(self.wf_readmap(os.path.join(self.Dir,wflow_soil),0.0,fail=True))
        self.Soil = cover(self.Soil, ordinal(subcatch > 0))

        self.OutputLoc = ordinal(self.wf_readmap(os.path.join(self.Dir,wflow_gauges),0.0,fail=True))  # location of output gauge(s)
        self.InflowLoc = ordinal(self.wf_readmap(os.path.join(self.Dir,wflow_inflow), 0.0))  # location abstractions/inflows.
        self.RiverWidth = self.wf_readmap(os.path.join(self.Dir,wflow_riverwidth), 0.0)
        self.bankFull = self.wf_readmap(os.path.join(self.Dir,wflow_bankfulldepth), 999999.0)
        self.floodPlainWidth = self.wf_readmap(os.path.join(self.Dir,wflow_floodplainwidth), 8000.0)
        self.floodPlainDist = self.wf_readmap(os.path.join(self.Dir,wflow_floodplaindist), 0.5)

        self.OutputId = ordinal(self.wf_readmap(os.path.join(self.Dir,wflow_subcatch),0.0,fail=True))  # location of subcatchment
        self.ZeroMap = 0.0 * scalar(subcatch)  #map with only zero's


        self.Latitude = ycoordinate(boolean(self.Altitude))
        self.Longitude = xcoordinate(boolean(self.Altitude))

        self.logger.info("Linking parameters to landuse, catchment and soil...")
        self.wf_updateparameters()

        # Check if we have reservoirs
        tt = pcr2numpy(self.ReserVoirLocs,0.0)
        self.nrres = tt.max()
        if self.nrres > 0:
            self.logger.info("A total of " +str(self.nrres) + " reservoirs found.")
            self.ReserVoirDownstreamLocs = downstream(self.TopoLdd,self.ReserVoirLocs)
            self.TopoLddOrg = self.TopoLdd
            self.TopoLdd = lddrepair(cover(ifthen(boolean(self.ReserVoirLocs),ldd(5)), self.TopoLdd))

        # Check if we have irrigation areas
        tt = pcr2numpy(self.IrrigationAreas, 0.0)
        self.nrirri = tt.max()

        self.Beta = scalar(0.6)  # For sheetflow

        self.N = self.readtblDefault(self.Dir + "/" + self.intbl + "/N.tbl", self.LandUse, subcatch, self.Soil,
                                     0.072)  # Manning overland flow
        self.NRiver = self.readtblDefault(self.Dir + "/" + self.intbl + "/N_River.tbl", self.LandUse, subcatch,
                                          self.Soil, 0.036)  # Manning river
        self.NFloodPlain = self.readtblDefault(self.Dir + "/" + self.intbl + "/N_FloodPlain.tbl", self.LandUse, subcatch,
                                          self.Soil, self.NRiver * 2.0)  # Manning river
        self.xl, self.yl, self.reallength = pcrut.detRealCellLength(self.ZeroMap, sizeinmetres)
        self.Slope = slope(self.Altitude)
        #self.Slope=ifthen(boolean(self.TopoId),max(0.001,self.Slope*celllength()/self.reallength))
        self.Slope = max(0.00001, self.Slope * celllength() / self.reallength)
        Terrain_angle = scalar(atan(self.Slope))


        self.wf_multparameters()
        self.N = ifthenelse(self.River, self.NRiver, self.N)

        # Determine river width from DEM, upstream area and yearly average discharge
        # Scale yearly average Q at outlet with upstream are to get Q over whole catchment
        # Alf ranges from 5 to > 60. 5 for hardrock. large values for sediments
        # "Noah J. Finnegan et al 2005 Controls on the channel width of rivers:
        # Implications for modeling fluvial incision of bedrock"

        upstr = catchmenttotal(1, self.TopoLdd)
        Qscale = upstr / mapmaximum(upstr) * Qmax
        W = (alf * (alf + 2.0) ** (0.6666666667)) ** (0.375) * Qscale ** (0.375) * (
            max(0.0001, windowaverage(self.Slope, celllength() * 4.0))) ** (-0.1875) * self.N ** (0.375)
        # Use supplied riverwidth if possible, else calulate
        self.RiverWidth = ifthenelse(self.RiverWidth <= 0.0, W, self.RiverWidth)

        self.UpdateMap = self.ZeroMap

        if self.updating:
            _tmp = pcr2numpy(self.OutputLoc, 0.0)
            gaugear = _tmp
            touse = numpy.zeros(gaugear.shape, dtype='int')

            for thecol in updateCols:
                idx = (gaugear == thecol).nonzero()
                touse[idx] = thecol

            self.UpdateMap = numpy2pcr(Nominal, touse, 0.0)
            # Calculate distance to updating points (upstream) annd use to scale the correction
            # ldddist returns zero for cell at the gauges so add 1.0 tp result
            self.DistToUpdPt = cover(
                min(ldddist(self.TopoLdd, boolean(cover(self.UpdateMap, 0)), 1) * self.reallength / celllength(),
                    self.UpdMaxDist), self.UpdMaxDist)
            #self.DistToUpdPt = ldddist(self.TopoLdd,boolean(cover(self.OutputId,0.0)),1)
            #* self.reallength/celllength()

        # Initializing of variables
        self.logger.info("Initializing of model variables..")
        self.TopoLdd = lddmask(self.TopoLdd, boolean(self.TopoId))
        catchmentcells = maptotal(scalar(self.TopoId))

        # Limit lateral flow per subcatchment (make pits at all subcatch boundaries)
        # This is very handy for Ribasim etc...
        if self.SubCatchFlowOnly > 0:
            self.logger.info("Creating subcatchment-only drainage network (ldd)")
            ds = downstream(self.TopoLdd,self.TopoId)
            usid = ifthenelse(ds != self.TopoId,self.TopoId,0)
            self.TopoLdd = lddrepair(ifthenelse(boolean(usid),ldd(5),self.TopoLdd))

        self.QMMConv = self.timestepsecs / (self.reallength * self.reallength * 0.001)  #m3/s --> mm
        self.ToCubic = (self.reallength * self.reallength * 0.001) / self.timestepsecs  # m3/s
        self.KinWaveVolume = self.ZeroMap
        self.OldKinWaveVolume = self.ZeroMap

        self.Aspect = scalar(aspect(self.Altitude))  # aspect [deg]
        self.Aspect = ifthenelse(self.Aspect <= 0.0, scalar(0.001), self.Aspect)

        # On Flat areas the Aspect function fails, fill in with average...
        self.Aspect = ifthenelse(defined(self.Aspect), self.Aspect, areaaverage(self.Aspect, self.TopoId))

        # Set DCL to riverlength if that is longer that the basic length calculated from grid
        drainlength = detdrainlength(self.TopoLdd, self.xl, self.yl)
        # Multiply with Factor (taken from upscaling operation, defaults to 1.0 if no map is supplied
        self.DCL = drainlength * max(1.0, self.RiverLengthFac)
        self.DCL = max(self.DCL, self.RiverLength)  # m

        self.SlopeDCL = self.Slope * self.reallength/self.DCL

        # water depth (m)
        # set width for kinematic wave to cell width for all cells
        self.Bw = detdrainwidth(self.TopoLdd, self.xl, self.yl)
        # However, in the main river we have real flow so set the width to the
        # width of the river
        self.Bw = ifthenelse(self.River, self.RiverWidth, self.Bw)

        #riverslopecor = drainlength / self.DCL
        #report(riverslopecor,"cor.map")
        #report(self.Slope * riverslopecor,"slope.map")
        self.AlpTerm = pow((self.N / (sqrt(self.SlopeDCL))), self.Beta)
        # power for Alpha
        self.AlpPow = (2.0 / 3.0) * self.Beta
        # initial approximation for Alpha

        self.logger.info("Saving summary maps...")

        if self.updating:
            report(self.DistToUpdPt, self.Dir + "/" + self.runId + "/outsum/DistToUpdPt.map")

        #self.IF = self.ZeroMap
        self.logger.info("End of initial section")

    def default_summarymaps(self):
          """
          Returns a list of default summary-maps at the end of a run.
          This is model specific. You can also add them to the [summary]section of the ini file but stuff
          you think is crucial to the model should be listed here


          """
          lst = ['self.N',
                'self.NRiver',
                'self.NFloodPlain',
                'self.xl',
                'self.yl',
                'self.RiverWidth',
                'self.Bw',
                'self.RiverLength',
                'self.RiverLengthFac',
                'self.DCL',
                'self.Slope',
                'self.SlopeDCL',
                'self.bankFull',
                'self.floodPlainWidth',
                'self.floodPlainDist']

          return lst


    def parameters(self):
        """
        Define all model parameters here that the framework should handle for the model
        See wf_updateparameters and the parameters section of the ini file
        If you use this make sure to all wf_updateparameters at the start of the dynamic section
        and at the start/end of the initial section
        """
        modelparameters = []

        #Static model parameters e.g.
        #modelparameters.append(self.ParamType(name="RunoffGeneratingGWPerc",stack="intbl/RunoffGeneratingGWPerc.tbl",type="static",default=0.1))
        # 3: Input time series ###################################################
        self.IW_mapstack = self.Dir + configget(self.config, "inputmapstacks", "Inwater",
                                               "/inmaps/IW")  # timeseries for specific runoff

        self.Inflow_mapstack = self.Dir + configget(self.config, "inputmapstacks", "Inflow",
                                                    "/inmaps/IF")  # timeseries for rainfall "/inmaps/IF" # in/outflow locations (abstractions)

        # Meteo and other forcing
        modelparameters.append(self.ParamType(name="InwaterForcing",stack=self.IW_mapstack ,type="timeseries",default=0.0,verbose=True,lookupmaps=[]))
        modelparameters.append(self.ParamType(name="Inflow",stack=self.Inflow_mapstack,type="timeseries",default=0.0,verbose=False,lookupmaps=[]))
        modelparameters.append(self.ParamType(name="ReserVoirLocs",stack='staticmaps/wflow_reservoirlocs.map',type="staticmap",default=0.0,verbose=False,lookupmaps=[]))
        modelparameters.append(self.ParamType(name="ResTargetFullFrac",stack='intbl/ResTargetFullFrac.tbl',type="tblsparse",default=0.8,verbose=False,lookupmaps=['staticmaps/wflow_reservoirlocs.map']))
        modelparameters.append(self.ParamType(name="ResTargetMinFrac",stack='intbl/ResTargetMinFrac.tbl',type="tblsparse",default=0.4,verbose=False,lookupmaps=['staticmaps/wflow_reservoirlocs.map']))
        modelparameters.append(self.ParamType(name="ResMaxVolume",stack='intbl/ResMaxVolume.tbl',type="tblsparse",default=0.0,verbose=False,lookupmaps=['staticmaps/wflow_reservoirlocs.map']))
        modelparameters.append(self.ParamType(name="ResMaxRelease",stack='intbl/ResMaxRelease.tbl',type="tblsparse",default=1.0,verbose=False,lookupmaps=['staticmaps/wflow_reservoirlocs.map']))
        modelparameters.append(self.ParamType(name="ResDemand",stack='intbl/ResDemand.tbl',type="tblsparse",default=1.0,verbose=False,lookupmaps=['staticmaps/wflow_reservoirlocs.map']))
        modelparameters.append(self.ParamType(name="IrrigationAreas", stack='staticmaps/wflow_irrigationareas.map',
                                              type="staticmap", default=0.0, verbose=False, lookupmaps=[]))
        modelparameters.append(self.ParamType(name="IrrigationSurfaceIntakes", stack='staticmaps/wflow_irrisurfaceintakes.map',
                       type="staticmap", default=0.0, verbose=False, lookupmaps=[]))
        modelparameters.append(self.ParamType(name="IrrigationSurfaceReturn", stack='staticmaps/wflow_irrisurfacereturns.map',
                       type="staticmap", default=0.0, verbose=False, lookupmaps=[]))

        return modelparameters

    def resume(self):

        if self.reinit == 1:
            self.logger.info("Setting initial conditions to default")
            self.WaterLevelCH = self.ZeroMap
            self.WaterLevelFP   = self.ZeroMap
            self.SurfaceRunoff = self.ZeroMap
            self.WaterLevel = self.WaterLevelCH + self.WaterLevelFP
            self.ReservoirVolume = self.ResMaxVolume * self.ResTargetFullFrac
        else:
            self.logger.info("Setting initial conditions from state files")
            self.wf_resume(os.path.join(self.Dir, self.instate))

        self.Pch = self.wetPerimiterCH(self.WaterLevelCH,self.Bw)
        self.Pfp =  ifthenelse(self.River,self.wetPerimiterFP(self.WaterLevelFP,self.floodPlainWidth,sharpness=self.floodPlainDist),0.0)
        self.WetPComb = self.Pch + self.Pfp
        self.Ncombined = (self.Pch/self.WetPComb*self.N**1.5 + self.Pfp/self.WetPComb*self.NFloodPlain**1.5)**(2./3.)


        self.AlpTermFP = pow((self.NFloodPlain / (sqrt(self.SlopeDCL))), self.Beta)
        self.AlpTermComb = pow((self.Ncombined / (sqrt(self.SlopeDCL))), self.Beta)

        self.AlphaFP = self.AlpTermFP * pow(self.Pfp, self.AlpPow)
        self.AlphaCh = self.AlpTerm * pow(self.Pch, self.AlpPow)
        self.Alpha = ifthenelse(self.River,self.AlpTermComb * pow(self.Pch + self.Pfp, self.AlpPow),self.AlphaCh)
        self.OldSurfaceRunoff = self.SurfaceRunoff

        self.SurfaceRunoffMM = self.SurfaceRunoff * self.QMMConv
        # Determine initial kinematic wave volume
        self.KinWaveVolume = (self.WaterLevelCH * self.Bw * self.DCL) + (self.WaterLevelFP * (self.Pfp + self.Bw) * self.DCL)
        self.OldKinWaveVolume = self.KinWaveVolume
        self.SurfaceRunoffMM = self.SurfaceRunoff * self.QMMConv

    def dynamic(self):
        """
        Stuff that is done for each timestep

        Below a list of variables that can be save to disk as maps or as
        timeseries (see ini file for syntax):

        *Dynamic variables*

        :var self.SurfaceRunoff: Total Surface runoff in the kinematic wave [m^3/s]
        :var self.Qbankfull: Discharge at bankfull level [m^3/s]
        :var self.Qchannel: Discharge in the channel (maxed at bankfull) [m^3/s]
        :var self.floodcells: All cells with an active floodplain [-]
        :var self.Qfloodplain: Discharge over the floodplain [m^3/s]
        :var self.WaterLevelCH: Water level in the channel [m] Cannot go above bankfull above the bottom
        :var self.WaterLevelFP: Water level on the floodplain [m] above the floodplain level (bottom + bankfull)
        :var self.WaterLevel: Total water level in the kinematic wave [m] (above the bottom)
        :var self.Pfp: Actual wetted perimiter of the floodplain [m]
        :var self.Pch: Actual wetted perimiter of the channel [m]

        *Static variables*

        :var self.Altitude: The altitude of each cell [m]
        :var self.Bw: Width of the river [m]
        :var self.River: boolean map indicating the presence of a river [-]
        :var self.DLC: length of the river within a cell [m]
        :var self.ToCubic: Mutiplier to convert mm to m^3/s for fluxes
        """

        self.thestep = self.thestep + 1
        self.wf_updateparameters()

        self.wf_multparameters()
        # The MAx here may lead to watbal error. However, if inwaterMMM becomes < 0, the kinematic wave becomes very slow......
        self.InwaterMM = max(0.0,self.InwaterForcing)
        self.Inwater = self.InwaterMM * self.ToCubic  # m3/s

        #only run the reservoir module if needed
        if self.nrres > 0:
            self.ReservoirVolume, self.Outflow, self.ResPercFull,\
            self.DemandRelease = simplereservoir(self.ReservoirVolume, self.SurfaceRunoff,\
                                                 self.ResMaxVolume, self.ResTargetFullFrac,
                                                 self.ResMaxRelease, self.ResDemand,
                                                 self.ResTargetMinFrac, self.ReserVoirLocs,
                                                 timestepsecs=self.timestepsecs)
            self.OutflowDwn = upstream(self.TopoLddOrg,cover(self.Outflow,scalar(0.0)))
            self.Inflow = self.OutflowDwn + cover(self.Inflow,self.ZeroMap)
        else:
            self.Inflow= cover(self.Inflow,self.ZeroMap)


        # Run only if we have irrigation areas or an externally given demand, determine irrigation demand based on potrans and acttrans
        if self.nrirri > 0 or hasattr(self,"IrriDemandExternal"):
            if not hasattr(self,"IrriDemandExternal"): # if not given
                self.IrriDemand, self.IrriDemandm3 = self.irrigationdemand(self.PotTrans,self.Transpiration,self.IrrigationAreas)
                IRDemand = idtoid(self.IrrigationAreas, self.IrrigationSurfaceIntakes, self.IrriDemandm3)  * -1.0
            else:
                IRDemand = self.IrriDemandExternal
            # loop over irrigation areas and assign Q to linked river extraction points
            self.Inflow = cover(IRDemand,self.Inflow)


        # Check if we do not try to abstract more runoff then present
        self.InflowKinWaveCell = upstream(self.TopoLdd, self.SurfaceRunoff)
        # The extraction should be equal to the discharge upstream cell.
        # You should not make the abstraction depended on the downstream cell, because they are correlated.
        # During a stationary sum they will get equal to each other.
        MaxExtract = self.InflowKinWaveCell + self.Inwater
        self.SurfaceWaterSupply = ifthenelse (self.Inflow < 0.0 , min(MaxExtract,-1.0 * self.Inflow), self.ZeroMap)
        self.OldSurfaceRunoff=self.SurfaceRunoff
        self.OldInwater=self.Inwater
        self.Inwater = self.Inwater + ifthenelse(self.SurfaceWaterSupply> 0, -1.0 * self.SurfaceWaterSupply,self.Inflow)

        ##########################################################################
        # Runoff calculation via Kinematic wave ##################################
        ##########################################################################
        # per distance along stream
        q = self.Inwater / self.DCL
        # discharge (m3/s)
        self.SurfaceRunoff = kinematic(self.TopoLdd, self.SurfaceRunoff, q, self.Alpha, self.Beta, self.Tslice,
                                       self.timestepsecs, self.DCL)  # m3/s

        self.InflowKinWaveCell = upstream(self.TopoLdd, self.SurfaceRunoff)

        # If inflow is negative we have abstractions. Check if demand can be met (by looking
        # at the flow in the upstream cell) and iterate if needed
        self.nrit = 0
        self.breakoff = 0.0001
        if float(mapminimum(spatial(self.Inflow))) < 0.0:
            while True:
                self.nrit += 1
                oldsup = self.SurfaceWaterSupply
                self.InflowKinWaveCell = upstream(self.TopoLdd, self.SurfaceRunoff)
                ##########################################################################
                # Iterate to make a better estimation for the supply #####################
                # (Runoff calculation via Kinematic wave) ################################
                ##########################################################################
                MaxExtract = self.InflowKinWaveCell + self.OldInwater
                self.SurfaceWaterSupply = ifthenelse(self.Inflow < 0.0, min(MaxExtract, -1.0 * self.Inflow),\
                                                      self.ZeroMap)
                # Fraction of demand that is not used but flows back into the river get fracttion and move to return locations
                self.DemandReturnFlow = cover(idtoid(self.IrrigationSurfaceIntakes,self.IrrigationSurfaceReturn,
                                               self.DemandReturnFlowFraction * self.SurfaceWaterSupply),0.0)

                self.Inwater = self.OldInwater + ifthenelse(self.SurfaceWaterSupply> 0, -1.0 * self.SurfaceWaterSupply,\
                                                            self.Inflow) + self.DemandReturnFlow
                # per distance along stream
                q = self.Inwater / self.DCL
                # discharge (m3/s)
                self.SurfaceRunoff = kinematic(self.TopoLdd, self.OldSurfaceRunoff, q, self.Alpha, self.Beta, self.Tslice,
                                               self.timestepsecs, self.DCL)  # m3/s
                self.SurfaceRunoffMM = self.SurfaceRunoff * self.QMMConv  # SurfaceRunoffMM (mm) from SurfaceRunoff (m3/s)

                self.InflowKinWaveCell = upstream(self.TopoLdd, self.OldSurfaceRunoff)
                deltasup = float(mapmaximum(abs(oldsup - self.SurfaceWaterSupply)))



                if deltasup < self.breakoff or self.nrit >= self.maxitsupply:
                    break

            self.InflowKinWaveCell = upstream(self.TopoLdd, self.SurfaceRunoff)
            self.updateRunOff()
        else:
            self.SurfaceRunoffMM = self.SurfaceRunoff * self.QMMConv  # SurfaceRunoffMM (mm) from SurfaceRunoff (m3/s)
            self.updateRunOff()

        # Now add the supply that is linked to irrigation areas to extra precip
        if self.nrirri > 0:
            # loop over irrigation areas and spread-out the supply over the area
            IRSupplymm = idtoid(self.IrrigationSurfaceIntakes, self.IrrigationAreas,
                                self.SurfaceWaterSupply * (1 - self.DemandReturnFlowFraction))
            sqmarea = areatotal(self.reallength * self.reallength, nominal(self.IrrigationAreas))

            self.IRSupplymm = cover(IRSupplymm/ (sqmarea / 1000.0 / self.timestepsecs),0.0)

        self.MassBalKinWave = (-self.KinWaveVolume + self.OldKinWaveVolume) / self.timestepsecs +\
                                self.InflowKinWaveCell + self.Inwater - self.SurfaceRunoff



        Runoff = self.SurfaceRunoff

        # Updating
        # --------
        # Assume a tss file with as many columns as outputlocs. Start updating for each non-missing value and start with the
        # first column (nr 1). Assumes that outputloc and columns match!

        if self.updating:
            self.QM = timeinputscalar(self.updateFile, self.UpdateMap) * self.QMMConv

            # Now update the state. Just add to the Ustore
            # self.UStoreDepth =  result
            # No determine multiplication ratio for each gauge influence area.
            # For missing gauges 1.0 is assumed (no change).
            # UpDiff = areamaximum(QM,  self.UpdateMap) - areamaximum(self.SurfaceRunoffMM, self.UpdateMap)
            UpRatio = areamaximum(self.QM, self.UpdateMap) / areamaximum(self.SurfaceRunoffMM, self.UpdateMap)

            UpRatio = cover(areaaverage(UpRatio, self.TopoId), 1.0)
            # Now split between Soil and Kyn  wave
            self.UpRatioKyn = min(self.MaxUpdMult, max(self.MinUpdMult, (UpRatio - 1.0) * self.UpFrac + 1.0))
            UpRatioSoil = min(self.MaxUpdMult, max(self.MinUpdMult, (UpRatio - 1.0) * (1.0 - self.UpFrac) + 1.0))

            # Update the kinematic wave reservoir up to a maximum upstream distance
            MM = (1.0 - self.UpRatioKyn) / self.UpdMaxDist
            self.UpRatioKyn = MM * self.DistToUpdPt + self.UpRatioKyn
            self.SurfaceRunoff = self.SurfaceRunoff * self.UpRatioKyn
            self.SurfaceRunoffMM = self.SurfaceRunoff * self.QMMConv  # SurfaceRunoffMM (mm) from SurfaceRunoff (m3/s)
            self.updateRunOff()
            Runoff = self.SurfaceRunoff

        ##########################################################################
        # water balance ###########################################
        ##########################################################################

        # Single cell based water budget. snow not included yet.


def main(argv=None):
    """
    Perform command line execution of the model.
    """
    caseName = "default_routing"
    global multpars
    runId = "run_default"
    configfile = "wflow_routing.ini"
    _lastTimeStep = 0
    _firstTimeStep = 0
    LogFileName = "wflow_routing.log"
    runinfoFile = "runinfo.xml"
    timestepsecs = 86400
    wflow_cloneMap = 'wflow_subcatch.map'
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
        opts, args = getopt.getopt(argv, 'F:L:hC:Ii:v:S:T:WR:u:s:EP:p:Xx:U:fOc:l:g:')
    except getopt.error, msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == '-C': caseName = a
        if o == '-R': runId = a
        if o == '-c': configfile = a
        if o == '-L': LogFileName = a
        if o == '-s': timestepsecs = int(a)
        if o == '-h': usage()
        if o == '-f': _NoOverWrite = 0
        if o == '-l': exec "loglevel = logging." + a



    starttime = dt.datetime(1990,01,01)

    if _lastTimeStep < _firstTimeStep:
        print "The starttimestep (" + str(_firstTimeStep) + ") is smaller than the last timestep (" + str(
            _lastTimeStep) + ")"
        usage()

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(myModel, _lastTimeStep, firstTimestep=_firstTimeStep,datetimestart=starttime)
    dynModelFw.createRunId(NoOverWrite=_NoOverWrite, level=loglevel, logfname=LogFileName,doSetupFramework=False)

    for o, a in opts:
        if o == '-X': configset(myModel.config, 'model', 'OverWriteInit', '1', overwrite=True)
        if o == '-I': configset(myModel.config, 'run', 'reinit', '1', overwrite=True)
        if o == '-i': configset(myModel.config, 'model', 'intbl', a, overwrite=True)
        if o == '-s': configset(myModel.config, 'model', 'timestepsecs', a, overwrite=True)
        if o == '-x': configset(myModel.config, 'model', 'sCatch', a, overwrite=True)
        if o == '-c': configset(myModel.config, 'model', 'configfile', a, overwrite=True)
        if o == '-g': configset(myModel.config,'model','instate',a,overwrite=True)

        if o == '-U':
            configset(myModel.config, 'model', 'updateFile', a, overwrite=True)
            configset(myModel.config, 'model', 'updating', "1", overwrite=True)
        if o == '-u':
            exec "zz =" + a
            updateCols = zz
        if o == '-P':
            left = a.split('=')[0]
            right = a.split('=')[1]
            configset(myModel.config,'variable_change_once',left,right,overwrite=True)
        if o == '-p':
            left = a.split('=')[0]
            right = a.split('=')[1]
            configset(myModel.config,'variable_change_timestep',left,right,overwrite=True)
        if o == '-T':
            configset(myModel.config, 'run', 'endtime', a, overwrite=True)
        if o == '-S':
            configset(myModel.config, 'run', 'starttime', a, overwrite=True)

    dynModelFw.setupFramework()
    dynModelFw._runInitial()
    dynModelFw._runResume()
    #dynModelFw._runDynamic(0, 0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
