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
          
    -F: if set wflow is expected to be run by FEWS. It will determine
        the timesteps from the runinfo.xml file and save the output initial
        conditions to an alternate location. Also set fewsrun=1 in the .ini file!
        
    -X: save state at the end of the run over the initial conditions at the start        

    -T: Set last timestep
    
    -S: Set the start timestep (default = 1)
    
    -s: Set the model timesteps in seconds
    
    -I: re-initialize the initial model conditions with default

    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -L: set the logfile
    
    -c: name of wflow the configuration file (default: Casename/wflow_routing.ini).
    
    -h: print usage information

    -l: loglevel (most be one of DEBUG, WARNING, ERROR)


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
wflowVersion = "$Revision: 900 $  $Date: 2014-01-09 18:41:06 +0100 (Thu, 09 Jan 2014) $"

updateCols = []

multpars = {}
multdynapars = {}


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

    def _initAPIVars(self):
        """
        Sets vars in the API that are forcing variables to the model
        """
        apivars = self.wf_supplyVariableNamesAndRoles()

        for var in apivars:
            exec "self."+ var[0] + " = self.ZeroMap"

    def updateRunOff(self):
        """
      Updates the kinematic wave reservoir. Should be run after updates to Q
      """
        self.WaterLevel = (self.Alpha * pow(self.SurfaceRunoff, self.Beta)) / self.Bw
        # wetted perimeter (m)
        P = self.Bw + (2 * self.WaterLevel)
        # Alpha
        self.Alpha = self.AlpTerm * pow(P, self.AlpPow)
        self.OldKinWaveVolume = self.KinWaveVolume
        self.KinWaveVolume = self.WaterLevel * self.Bw * self.DCL

    def stateVariables(self):
        """
        returns a list of state variables that are essential to the model.
        This list is essential for the resume and suspend functions to work.

        This function is specific for each model and **must** be present.

        :var self.SurfaceRunoff: Surface runoff in the kin-wave resrvoir [m^3/s]
        :var self.WaterLevel: Water level in the kin-wave resrvoir [m]
        """
        states = ['SurfaceRunoff', 'WaterLevel']

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

        if self.fewsrun:
            self.logger.info("Saving initial conditions for FEWS...")
            self.wf_suspend(self.Dir + "/outstate/")

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
        self.reinit = int(configget(self.config, "model", "reinit", "0"))
        self.fewsrun = int(configget(self.config, "model", "fewsrun", "0"))
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
        wflow_landuse = configget(self.config, "model", "wflow_landuse", "staticmaps/wflow_landuse.map")
        wflow_soil = configget(self.config, "model", "wflow_soil", "staticmaps/wflow_soil.map")

        # 2: Input base maps ########################################################
        subcatch = ordinal(readmap(os.path.join(self.Dir,wflow_subcatch)))  # Determines the area of calculations (all cells > 0)
        subcatch = ifthen(subcatch > 0, subcatch)

        self.Altitude = readmap(os.path.join(self.Dir,wflow_dem))  # * scalar(defined(subcatch)) # DEM
        self.TopoLdd = readmap(os.path.join(self.Dir,wflow_ldd))  # Local
        self.TopoId = readmap(os.path.join(self.Dir,wflow_subcatch))  # area map
        self.River = cover(boolean(readmap(os.path.join(self.Dir,wflow_river))), 0)

        self.RiverLength = cover(pcrut.readmapSave(os.path.join(self.Dir,wflow_riverlength), 0.0), 0.0)
        # Factor to multiply riverlength with (defaults to 1.0)
        self.RiverLengthFac = pcrut.readmapSave(os.path.join(self.Dir,wflow_riverlength_fact), 1.0)

        # read landuse and soilmap and make sure there are no missing points related to the
        # subcatchment map. Currently sets the lu and soil type  type to 1
        self.LandUse = readmap(os.path.join(self.Dir,wflow_landuse))
        self.LandUse = cover(self.LandUse, nominal(ordinal(subcatch) > 0))
        self.Soil = readmap(os.path.join(self.Dir,wflow_soil))
        self.Soil = cover(self.Soil, nominal(ordinal(subcatch) > 0))

        self.OutputLoc = readmap(os.path.join(self.Dir,wflow_gauges))  # location of output gauge(s)
        self.InflowLoc = pcrut.readmapSave(os.path.join(self.Dir,wflow_inflow), 0.0)  # location abstractions/inflows.
        self.RiverWidth = pcrut.readmapSave(os.path.join(self.Dir,wflow_riverwidth), 0.0)
        self.OutputId = readmap(os.path.join(self.Dir,wflow_subcatch))  # location of subcatchment
        self.ZeroMap = 0.0 * scalar(subcatch)  #map with only zero's

        self.IW_mapstack = self.Dir + configget(self.config, "inputmapstacks", "Inwater",
                                               "/inmaps/IW")  # timeseries for specific runoff

        self.Inflow_mapstack = self.Dir + configget(self.config, "inputmapstacks", "Inflow",
                                                    "/inmaps/IF")  # timeseries for rainfall "/inmaps/IF" # in/outflow locations (abstractions)

        self.Latitude = ycoordinate(boolean(self.Altitude))
        self.Longitude = xcoordinate(boolean(self.Altitude))

        self.logger.info("Linking parameters to landuse, catchment and soil...")
        self.Beta = scalar(0.6)  # For sheetflow

        self.N = self.readtblDefault(self.Dir + "/" + self.intbl + "/N.tbl", self.LandUse, subcatch, self.Soil,
                                     0.072)  # Manning overland flow
        self.NRiver = self.readtblDefault(self.Dir + "/" + self.intbl + "/N_River.tbl", self.LandUse, subcatch,
                                          self.Soil, 0.036)  # Manning river

        self.xl, self.yl, self.reallength = pcrut.detRealCellLength(self.ZeroMap, sizeinmetres)
        self.Slope = slope(self.Altitude)
        #self.Slope=ifthen(boolean(self.TopoId),max(0.001,self.Slope*celllength()/self.reallength))
        self.Slope = max(0.00001, self.Slope * celllength() / self.reallength)
        Terrain_angle = scalar(atan(self.Slope))

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
        self._initAPIVars()
        self.logger.info("End of initial section")

    def default_summarymaps(self):
          """
          Returns a list of default summary-maps at the end of a run.
          This is model specific. You can also add them to the [summary]section of the ini file but stuff
          you think is crucial to the model should be listed here


          """
          lst = ['self.RiverWidth',
                'self.N',
                'self.xl',
                'self.yl',
                'self.DCL',
                'self.Bw',
                'self.Slope',
                'self.SlopeDCL']

          return lst

    def resume(self):

        if self.reinit == 1:
            self.logger.info("Setting initial conditions to default")
            self.WaterLevel = self.ZeroMap
            self.SurfaceRunoff = self.ZeroMap
        else:
            self.logger.info("Setting initial conditions from state files")
            self.wf_resume(os.path.join(self.Dir,"instate"))

        P = self.Bw + (2.0 * self.WaterLevel)
        self.Alpha = self.AlpTerm * pow(P, self.AlpPow)
        self.OldSurfaceRunoff = self.SurfaceRunoff

        self.SurfaceRunoffMM = self.SurfaceRunoff * self.QMMConv
        # Determine initial kinematic wave volume
        self.KinWaveVolume = self.WaterLevel * self.Bw * self.DCL
        self.OldKinWaveVolume = self.KinWaveVolume
        self.SurfaceRunoffMM = self.SurfaceRunoff * self.QMMConv

    def dynamic(self):
        """
        Stuf that is done for each timestep

        Below a list of variables that can be save to disk as maps or as
        timeseries (see ini file for syntax):

        *Dynamic variables*

        :var self.SurfaceRunoff: Surface runoff in the kinematic wave [m^3/s]
        :var self.WaterLevel: Water level in the kinematic wave [m] (above the bottom)


        *Static variables*

        :var self.Altitude: The altitude of each cell [m]
        :var self.Bw: Width of the river [m]
        :var self.River: booolean map indicating the presence of a river [-]
        :var self.DLC: length of the river within a cell [m]
        :var self.ToCubic: Mutiplier to convert mm to m^3/s for fluxes
        """

        self.logger.debug("Step: " + str(int(self.currentStep)) + "/" + str(int(self._d_nrTimeSteps)))
        self.thestep = self.thestep + 1
        self.InwaterForcing = cover(self.wf_readmap(self.IW_mapstack, 1.0), scalar(0.0))

        if self.thestep > 28:
            self.InwaterForcing = cover(0.0)
        #self.Inflow=cover(self.wf_readmap(self.Inflow),0)
        if (os.path.exists(self.caseName + self.inflowTss)):
            self.Inflow = cover(timeinputscalar(self.caseName + self.inflowTss, nominal(self.InflowLoc)), 0)
        else:
            self.Inflow = cover(self.wf_readmap(self.Inflow_mapstack, 0.0,verbose=False),0)

        # The MAx here may lead to watbal error. Howevere, if inwaterMMM becomes < 0, the kinematic wave becomes very slow......
        self.InwaterMM = max(0.0,self.InwaterForcing)
        self.Inwater = self.InwaterMM * self.ToCubic  # m3/s

        self.Inwater = self.Inwater + self.Inflow  # Add abstractions/inflows in m^3/sec

        ##########################################################################
        # Runoff calculation via Kinematic wave ##################################
        ##########################################################################
        # per distance along stream
        q = self.Inwater / self.DCL
        # discharge (m3/s)
        self.SurfaceRunoff = kinematic(self.TopoLdd, self.SurfaceRunoff, q, self.Alpha, self.Beta, self.Tslice,
                                       self.timestepsecs, self.DCL)  # m3/s
        self.SurfaceRunoffMM = self.SurfaceRunoff * self.QMMConv  # SurfaceRunoffMM (mm) from SurfaceRunoff (m3/s)
        self.updateRunOff()
        self.InflowKinWaveCell = upstream(self.TopoLdd, self.SurfaceRunoff)
        self.MassBalKinWave = (self.KinWaveVolume - self.OldKinWaveVolume) / self.timestepsecs + self.InflowKinWaveCell + self.Inwater - self.SurfaceRunoff

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
    LogFileName = "wflow.log"
    fewsrun = False
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
        opts, args = getopt.getopt(argv, 'F:L:hC:Ii:v:S:T:WR:u:s:EP:p:Xx:U:fOc:l:')
    except getopt.error, msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == '-F':
            runinfoFile = a
            fewsrun = True
        if o == '-C': caseName = a
        if o == '-R': runId = a
        if o == '-c': configfile = a
        if o == '-L': LogFileName = a
        if o == '-s': timestepsecs = int(a)
        if o == '-T': _lastTimeStep = int(a)
        if o == '-S': _firstTimeStep = int(a)
        if o == '-h': usage()
        if o == '-f': _NoOverWrite = 0
        if o == '-l': exec "loglevel = logging." + a

    if fewsrun:
        ts = getTimeStepsfromRuninfo(runinfoFile, timestepsecs)
        starttime = getStartTimefromRuninfo(runinfoFile)
        if (ts):
            _lastTimeStep = ts
            _firstTimeStep = 1
        else:
            print "Failed to get timesteps from runinfo file: " + runinfoFile
            exit(2)
    else:
        starttime = dt.datetime(1990,01,01)
        
    if _lastTimeStep < _firstTimeStep:
        print "The starttimestep (" + str(_firstTimeStep) + ") is smaller than the last timestep (" + str(
            _lastTimeStep) + ")"
        usage()

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(myModel, _lastTimeStep, firstTimestep=_firstTimeStep,datetimestart=starttime)
    dynModelFw.createRunId(NoOverWrite=_NoOverWrite, level=loglevel, logfname=LogFileName)

    for o, a in opts:
        if o == '-X': configset(myModel.config, 'model', 'OverWriteInit', '1', overwrite=True)
        if o == '-I': configset(myModel.config, 'model', 'reinit', '1', overwrite=True)
        if o == '-i': configset(myModel.config, 'model', 'intbl', a, overwrite=True)
        if o == '-s': configset(myModel.config, 'model', 'timestepsecs', a, overwrite=True)
        if o == '-x': configset(myModel.config, 'model', 'sCatch', a, overwrite=True)
        if o == '-c': configset(myModel.config, 'model', 'configfile', a, overwrite=True)
        if o == '-U':
            configset(myModel.config, 'model', 'updateFile', a, overwrite=True)
            configset(myModel.config, 'model', 'updating', "1", overwrite=True)
        if o == '-u':
            exec "zz =" + a
            updateCols = zz

    dynModelFw._runInitial()
    dynModelFw._runResume()
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
