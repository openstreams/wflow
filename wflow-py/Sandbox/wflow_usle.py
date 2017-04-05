#!/usr/bin/python
# coding=utf-8

"""
Definition of the wflow_usle model.
-----------------------------------


Usage:
wflow_usle  -C case -R Runid -c inifile

    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -c name of the config file (in the case directory)
    


This script assumes monthly timesteps.

usle_r, Rainfall erosivity is determined according to:
Beskow, S., Mello, C.R., Norton, L.D., Curi, N., Viola, M.R., Avanzi, J.C., 2009.
Soil erosion prediction in the Grande River Basin, Brazil using distributed modeling.
CATENA 79, 49–59. doi:10.1016/j.catena.2009.05.010

Determine m using:
------------------

Slope gradient Value of m
1%              0.2
1–3%            0.3
3–4.5%          0.4
4.5% or more    0.5

Determien L using the cell-length

self.usle_l = (self.Lambda/22.13)**self.m # The L factor
# Now determime S using: S=(0.43+0.30s+0.043s^2)/6.613 # s = slope
# Wischmeier and Smith (1978),
self.usle_s = (0.43+ 0.30*(self.Slope*100) + 0.043 * (self.Slope * 100)**2)/6.613


"""

import numpy
import os
import os.path
import shutil, glob
import getopt

from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *




def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)

class WflowModel(DynamicModel):
  """
  The user defined model class. This is your work!
  """

  def __init__(self, cloneMap,Dir,RunDir,configfile):
      """
      *Required*

      The init function **must** contain what is shown below. Other functionality
      may be added by you if needed.

      """
      DynamicModel.__init__(self)
      setclone(Dir + "/staticmaps/" + cloneMap)
      self.runId=RunDir
      self.caseName=Dir
      self.Dir = Dir
      self.configfile = configfile


  def parameters(self):
      """
      List all the parameters (both static and forcing here). Use the wf_updateparameters()
      function to update them in the initial section (static) and the dynamic section for
      dynamic parameters and forcing date.

      Possible parameter types are:

      + staticmap: Read at startup from map
      + statictbl: Read at startup from tbl, fallback to map (need Landuse, Soil and TopoId (subcatch) maps!
      + timeseries: read map for each timestep
      + monthlyclim: read a map corresponding to the current month (12 maps in total)
      + dailyclim: read a map corresponding to the current day of the year
      + hourlyclim: read a map corresponding to the current hour of the day (24 in total)


      :return: List of modelparameters
      """
      modelparameters = []

      #Static map model parameters
      modelparameters.append(self.ParamType(name="Altitude",stack="staticmaps/wflow_dem.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[]))
      modelparameters.append(self.ParamType(name="LandUse",stack="staticmaps/wflow_landuse.map",type="staticmap",default=1,verbose=False,lookupmaps=[]))
      modelparameters.append(self.ParamType(name="Soil",stack="staticmaps/wflow_soil.map",type="staticmap",default=1,verbose=False,lookupmaps=[]))
      modelparameters.append(self.ParamType(name="TopoId",stack="staticmaps/wflow_subcatch.map",type="staticmap",default=1,verbose=False,lookupmaps=[]))
      modelparameters.append(self.ParamType(name="TopoLdd",stack="staticmaps/wflow_ldd.map",type="staticmap",default=1,verbose=False,lookupmaps=[]))

      # These should be linked to soil type
      modelparameters.append(self.ParamType(name="percent_clay",stack="intbl/percent_clay.tbl",type="statictbl",default=0.1, verbose=False,lookupmaps=[]))
      modelparameters.append(self.ParamType(name="percent_silt",stack="intbl/percent_silt.tbl",type="statictbl",default=0.1, verbose=False,lookupmaps=[]))
      modelparameters.append(self.ParamType(name="percent_oc",stack="intbl/percent_oc.tbl",type="statictbl",default=0.01, verbose=False,lookupmaps=[]))

      # Sediment delivery ratio
#      modelparameters.append(self.ParamType(name="dratio",stack="intbl/dratio.tbl",type="statictbl",default=1.0, verbose=False,lookupmaps=[]))
      #modelparameters.append(self.ParamType(name="usle_k", stack="intbl/usle_k.tbl", type="statictbl", default=1.0, verbose=False,lookupmaps=[]))
      modelparameters.append(self.ParamType(name="usle_c", stack="intbl/usle_c.tbl", type="statictbl", default=1.0, verbose=False,lookupmaps=[]))
      modelparameters.append(self.ParamType(name="usle_p", stack="intbl/usle_p.tbl", type="statictbl", default=1.0, verbose=False,lookupmaps=[]))

      #Climatology
      modelparameters.append(self.ParamType(name="LAI",stack="inmaps/climatology/LAI",type="monthlyclim",default=0.9, verbose=False,lookupmaps=[]))

      # Meteo and other forcing
      modelparameters.append(self.ParamType(name="Precipitation",stack="inmaps/P",type="timeseries",default=0.0,verbose=True,lookupmaps=[]))


      return modelparameters

  def stateVariables(self):
      """
      *Required*

      Returns a list of state variables that are essential to the model.
      This list is essential for the resume and suspend functions to work.

      This function is specific for each model and **must** be present. This is
      where you specify the state variables of you model. If your model is stateless
      this function must return and empty array (states = [])

      In the simple example here the TSoil variable is a state
      for the model.

      :var TSoil: Temperature of the soil [oC]
      """
      states = []

      return states


  def supplyCurrentTime(self):
      """
      *Optional*

      Supplies the current time in seconds after the start of the run
      This function is optional. If it is not set the framework assumes
      the model runs with daily timesteps.

      Output:

          - time in seconds since the start of the model run

      """

      return self.currentTimeStep() * int(configget(self.config,'model','timestepsecs','86400'))

  def suspend(self):
    """
      *Required*
      
      Suspends the model to disk. All variables needed to restart the model
      are saved to disk as pcraster maps. Use resume() to re-read them
      
      This function is required. 
      
    """

    self.logger.info("Saving initial conditions...")
    #: It is advised to use the wf_suspend() function
    #: here which will suspend the variables that are given by stateVariables
    #: function.
    self.wf_suspend(self.Dir + "/outstate/")


  def initial(self):

    """
    *Required*
    
    Initial part of the model, executed only once. It reads all static model
    information (parameters) and sets-up the variables used in modelling.
    
    This function is required. The contents is free. However, in order to
    easily connect to other models it is advised to adhere to the directory
    structure used in the other models.

    :var percent_clay: Map or tbl
    :var percent_silt: Map or tbl
    :var wflow_dem.map: elevation model
    :var wflow_ldd.map: Drainage network

    
    """
    #: pcraster option to calculate with units or cells. Not really an issue
    #: in this model but always good to keep in mind.
    setglobaloption("unittrue")


    self.timestepsecs = int(configget(self.config,'run','timestepsecs','86400'))
    sizeinmetres = int(configget(self.config, "layout", "sizeinmetres", "0"))
    self.basetimestep=86400
    # Reads all parameter from disk
    self.wf_updateparameters()
    self.ZeroMap = 0.0 * self.Altitude
    self.SedStore = self.ZeroMap

    self.xl, self.yl, self.reallength = pcrut.detRealCellLength(self.ZeroMap, sizeinmetres)


    # Calulate slope taking into account that x,y may be in lat,lon
    self.Slope = slope(self.Altitude)
    self.Slope = max(0.00001, self.Slope * celllength() / self.reallength)
    # limit slope
    self.DUSlope = min(1.0,max(0.005,self.Slope))

    """
    First determine m exponent based on Slope (https://www.researchgate.net/publication/226655635_Estimation_of_Soil_Erosion_for_a_Himalayan_Watershed_Using_GIS_Technique)

    """
    self.m = ifthenelse(self.Slope <= scalar(0.01), scalar(0.2),
                        ifthenelse(self.Slope <= 0.03, scalar(0.03),
                                   ifthenelse(self.Slope <= 0.045, scalar(0.5), scalar(0.5))))
    # Sel Lambda equal to grid size
    self.Lambda = self.reallength
    self.usle_l = (self.Lambda/22.13)**self.m # The L factor
    # Now determime S using: S=(0.43+0.30s+0.043s^2)/6.613 # s = slope
    self.usle_s = (0.43+ 0.30*(self.Slope*100) + 0.043 * (self.Slope * 100)**2)/6.613

    self.percent_sand = 100-self.percent_clay-self.percent_silt
    """
    Calculation of USLE K factor  based on:
    Maeda et al. (2010) - 'Potential impacts of agricultural expansion and climate change on soil erosion in the Eastern Arc Mountains of Kenya', doi:10.1016/j.geomorph.2010.07.019
    """

    self.SN = 1 - self.percent_sand / 100
    self.usle_k = (0.2 + 0.3*exp(-0.0256*self.percent_sand * (1 - self.percent_silt/100))) \
        * (self.percent_silt/(max(0.01,self.percent_clay+self.percent_silt)))**0.3 \
        * (1 - (0.25*self.percent_oc)/(self.percent_oc + exp(3.72 - 2.95*self.percent_oc))) \
        * (1 - (0.7*self.SN) / (self.SN + exp(-5.51 + 22.9*self.SN)))

    """
    Calculation of sediment delivery ratio based on channel:
    https://peerj.com/preprints/2227.pdf
    http://data.naturalcapitalproject.org/nightly-build/invest-users-guide/html/sdr.html
    """
    self.unitareaupstr = catchmenttotal(1, self.TopoLdd)
    self.Dup = catchmenttotal(max(0.001,self.usle_c),self.TopoLdd)/self.unitareaupstr * catchmenttotal(self.DUSlope,self.TopoLdd)/self.unitareaupstr *\
        sqrt(catchmenttotal(self.reallength,self.TopoLdd))

    self.drainlength = detdrainlength(self.TopoLdd, self.xl, self.yl)
    self.Ddn = self.drainlength/(max(0.001,max(0.001,self.usle_c)) * self.DUSlope)
    self.IC = log10(self.Dup/self.Ddn)
    self.SDRMax = 0.8 #  (Vigiak et al., 2012)
    self.IC0 = 0.5
    self.k = 2.0
    expfact = exp((self.IC0 - self.IC)/self.k)
    self.SDR = self.SDRMax/(1 + expfact)

    self.SDR_area = 0.472 * catchmenttotal(self.reallength/1000.0,self.TopoLdd)**(-0.125)



    self.logger.info("Starting Dynamic run...")
    self.thestep = 0

  def resume(self):
    """
    *Required*

    This function is required. Read initial state maps (they are output of a 
    previous call to suspend()). The implementation shown here is the most basic
    setup needed.
    
    """
    self.logger.info("Reading initial conditions...")
    #: It is advised to use the wf_resume() function
    #: here which pick up the variable save by a call to wf_suspend()
    try:
        self.wf_resume(self.Dir + "/instate/")
    except:
        self.logger.warn("Cannot load initial states, setting to default")
        for s in self.stateVariables():
            exec "self." + s + " = cover(1.0)"


  def default_summarymaps(self):
      """
      *Optional*

      Return a default list of variables to report as summary maps in the outsum dir.
      The ini file has more options, including average and sum
      """
      return ['self.Altitude']

  def dynamic(self):
    """
    *Required*

    This is where all the time dependent functions are executed. Time dependent
    output should also be saved here.
    """


    self.thestep = self.thestep + 1

    self.wf_updateparameters() # Read forcing and dynamic variables

    rintnsty = self.Precipitation/self.timestepsecs

    # Weird assumption for now, should be a lookuptabel of LAI and landuse type...
    self.canopy_cover = min(1.0,self.LAI)
    # Determine erosivity from monthly rainfall and average yearly sum
    # taken from Beskow, S., Mello, C.R., Norton, L.D., Curi, N., Viola, M.R., Avanzi, J.C., 2009.
    # Soil erosion prediction in the Grande River Basin, Brazil using distributed modeling.
    # CATENA 79, 49–59. doi:10.1016/j.catena.2009.05.010

    self.usle_r = (125.92 * (self.Precipitation/self.Pmean)**0.603 + \
        111.173 * (self.Precipitation/self.Pmean) ** 0.691 + \
        68.73 * (self.Precipitation/self.Pmean)** 0.841) / 3.0


    self.SoilLoss = self.usle_l * self.usle_s * self.usle_k * self.usle_r *self.usle_c * self.usle_p
    self.SoilLossAvgUpstr = catchmenttotal(self.SoilLoss,self.TopoLdd)/catchmenttotal(1,self.TopoLdd)
    self.SoilLossTotalUpstr = self.SoilLossAvgUpstr * catchmenttotal((self.reallength * self.reallength)/10000.0,self.TopoLdd)
    self.SedimentYieldUpstrAvg = self.SDR * self.SoilLossAvgUpstr # Average upstream of each pixel
    self.SedimentYieldUpstrTot = self.SDR * self.SoilLossTotalUpstr # Total upstream of each pixel
    self.SedimentYield = self.SDR * self.SoilLoss









# The main function is used to run the program from the command line

def main(argv=None):
    """
    *Optional but needed it you want to run the model from the command line*
    
    Perform command line execution of the model. This example uses the getopt
    module to parse the command line options.
    
    The user can set the caseName, the runDir, the timestep and the configfile.
    """
    global multpars
    caseName = "default"
    runId = "run_default"
    configfile="wflow_usle.ini"
    _lastTimeStep = 10
    _firstTimeStep = 1
    timestepsecs=86400
    wflow_cloneMap = 'wflow_dem.map'

    # This allows us to use the model both on the command line and to call 
    # the model usinge main function from another python script.

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    opts, args = getopt.getopt(argv, 'C:S:T:c:s:R:')

    for o, a in opts:
        if o == '-C': caseName = a
        if o == '-R': runId = a
        if o == '-c': configfile = a
        if o == '-s': timestepsecs = int(a)
        if o == '-T': _lastTimeStep=int(a)
        if o == '-S': _firstTimeStep=int(a)

    if (len(opts) <=1):
        usage()

    myModel = WflowModel(wflow_cloneMap, caseName,runId,configfile)
    dynModelFw = wf_DynamicFramework(myModel, _lastTimeStep,firstTimestep=_firstTimeStep)
    dynModelFw.createRunId(NoOverWrite=False,level=logging.DEBUG)
    dynModelFw._runInitial()
    dynModelFw._runResume()
    dynModelFw._runDynamic(_firstTimeStep,_lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
