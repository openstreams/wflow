#!/usr/bin/python

from math import cos, sin, asin, sqrt
import math
import numpy as np
import os
import os.path
import shutil, glob
import getopt

from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *
from time import strftime
#import scipy

#from datetime import date
#from pcse_util import * 
#from math import exp, expm1

"""
Definition of the wflow_lintul model.

***************************************************************************
*   LINTUL3 simulates potential rice production using weather data as     *
*   forcing variables. It can also simulate water-limited production; in  *
*   that case, soil data are required. LINTUL3 is an extended version of  *
*   LINTUL1 (the version of LINTUL for optimal growth conditions) and     *
*   LINTUL2 includes a simple water balance for studying effects of       *
*   drought. LINTUL2-N includes N-limitation on crop growth. The latter   *
*   program is called LINTUL3.                                            *
***************************************************************************

Usage:
wflow_sceleton  -C case -R Runid -c inifile

    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -c name of the config file (in the case directory)
    
$Author: SanderCdeVries $
$Id: lintul1.py 2015-11-20 16:08:06Z 
$Rev: 001 $
"""

runinfoFile           = "runinfo.xml"
runinfoFile_Lintul    = os.path.join(os.path.abspath("wflow_lintul"), "inmaps", runinfoFile)
starttime             = getStartTimefromRuninfo(runinfoFile_Lintul)
DOY_wflow             = starttime.strftime('%j')

def NOTNUL(matrix):
    """    
    NOTNUL is a FST Fortran Simulation Translator intrinsic function.
	If the input value is positive: NOTNUL will just return the value as is. 
	If the input value equals zero: NOTNUL will return a value of 1 instead.
    Sander de Vries, 01-12-2015
    """
    b = np.shape(matrix)
    nullen = np.zeros(b)
    NOTNUL_check = np.equal(nullen, matrix)
    matrix += NOTNUL_check[:] 
    return matrix

def astro2(DAY, LAT):
    """    
* ---------------------------------------------------------------------*
*   SUBROUTINE ASTRO                                                    *
*   Purpose: This subroutine calculates the astronomic daylength,       *
*            based on Goudriaan and van Laar 1994, around p.30.         *
* ---------------------------------------------------------------------*

    Originally: Daniel van Kraalingen, April 1991
    This version: Sander de Vries, 30-10-2015
    """

   # Constants
    PI    = 3.1415926
    
   # SINE AND COSINE OF LATITUDE
    SINLAT = sin(PI * LAT / 180.)
    COSLAT = cos(PI * LAT / 180.)

   # MAXIMAL SINE OF DECLINATION
    SINDCM = sin(PI * 23.45 / 180.)

   # SINE AND COSINE OF DECLINATION  (EQUATIONS 3.4, 3.5)

   # SINDEC = -SINDCM * cos(2.* PI * (DAY+10.)/365.)
   # The '9' below (instead of 10) keeps things in sync with FST...
   # Todo: try to understand this at a certain point... for now it works perfectly. 
    SINDEC = -SINDCM * cos(2. * PI * (DAY + 11.) / 365.)
    COSDEC = sqrt(1.- SINDEC * SINDEC)

   # THE TERMS A AND B ACCORDING TO EQUATION 3.3

    A     = SINLAT * SINDEC
    B     = COSLAT * COSDEC

   # DAYLENGTH ACCORDING TO EQUATION 3.6. Make sure not to use asin from PCRaster...
    DAYL  = 12. * (1. + (2. / PI) * math.asin(A / B))

    return DAYL
    END


class Afgen2(object):
    """Emulates the AFGEN function in TTUTIL with Numpy.interp. Author: Allard de Wit
	   Adaptation for wflow/PCRaster: Sander de Vries, 11-2015
    """
    
    def __init__(self, tbl_xy, unit=None):
    
        x = tbl_xy[0::2]
        y = tbl_xy[1::2]
        
        # Determine if there are empty pairs in tbl_xy by searching
        # for the point where x stops monotonically increasing
        xpos = np.arange(1, len(x))
        ibreak = False
        for i in xpos:
            if x[i] <= x[i-1]:
                ibreak = True
                break
        
        if ibreak is True:
            x = x[0:i]
            y = y[0:i]
        
        self.x = x
        self.y = y
        self.unit = unit
    
    def __call__(self, p):
        
        v = np.interp(p, self.x, self.y)
        
        # if a unum unit is defined, multiply with a unit
        if self.unit is not None:
            v *= self.unit
        return v
    
    def __str__(self):
        msg = "AFGEN interpolation over (X,Y) pairs:\n"
        for (x,y) in zip(self.x, self.y):
            msg += ("(%f,%f)\n " % (x,y))
        msg += "\n"
        if self.unit is not None:
             msg += "Return value as unit: %s" % self.unit
        return msg


#-------------------------------------------------------------------------------






###############################################################################
# Some general information for setting up the runs:
StartDay    = 8.          # TODO later via kaart/forcing inlezen?
DOYEM       = 25.         # Day of Emergence (transplanting in this case, actually). 
LAT         = 3.16        # the latitude of Bah Lias, North Sumatra, for now.
DELT        = 1	          # Time step = 1 day.
TSUMI       = 362.        # TSUM at transplanting, kind of crop specific, actually.
TINY        = 1e-6
iDOY_wflow             = int(starttime.strftime('%j'))

# This needs to be set according to the geographic extent (map dimensions) of your study area/catchment:
np_Zero     = numpy.zeros((219,286))
np_One      = numpy.ones ((219,286))

	 
# Crop specific coefficients for rice:
K           = 0.6         # light extinction coefficient
LUE         = 3.          # Light use efficiency.
SLAC        = 0.02	      # Specific leaf area constant.
TSUMAN      = 1420.
TSUMMT      = 580.
TBASE       = 8.          # Base temperature for spring wheat crop.
RGRL        = 0.009       # Relative growth rate of LAI at the exponential growth phase
WLVGI       = 0.86
WSTI        = 0.71
WRTLI       = 1.58
WSOI        = 0.00	 
RDRNS       = 0.03
DVSDR       = 0.8         # Development stage above which death of leaves and roots start.
RDRRT       = 0.03        # Relative death rate of roots.
TTSUM       = 2000. 
RDRSHM      = 0.03        # and the maximum relative death rate of leaves due to shading.
LAICR       = 4.0         # (oC d)-1, critical LAI above which mutual shading of leaves occurs.
ROOTDM      = 1.0         # Maximum rooting depth (m) for a rice crop.
RRDMAX      = 0.010       # maximum rate of increase in rooting depth (m d-1) for a rice crop.
ROOTDI      = 0.05        # initial rooting depth (m)
TRANCO      = 3.          # Transpiration constant (mm/day) indicating the level of drought tolerance of the rice crop.


	
# Fixing some values here for now: to be simulated or read from maps later on.
#TRANRF      = np_One[:] # If enabled: switches off the included water balance and removes water limitation effects!
                        # TRANRF equation further below needs to be switched off if this is enabled.
NNI         = 1.
WCI         = 0.60
WCAD        = 0.01
WCWP        = 0.21
WCFC        = 0.47
WCWET       = 0.50
#WA = WCST   = 0.55
WCST   = 0.55 # aangepast voor BMI sdv 01-06-2015
WCWP        = 0.20
DRATE       = 30.   
IRRIGF      = 0.
IRRIG       = np_Zero[:]

     
# Nitrogen-related coefficients, specific for rice.  
NPART       = 1.0         # Coefficient for the effect of N stress on leaf biomass reduction
NSLA        = 1.0         # Coefficient for the effect of N stress on SLA reduction
NLAI        = 1.0         # Coefficient for the effect of N stress on LAI reduction(during juvenile phase)
     
# ----------------------- Interpolation functions---------------------#
	
# Relative death rate of leaves as a function of Developmental stage 
# Note: from FST Lintul 3-rice, 29-10-15 (adopted from ORYZA2000), SdV
RDRTB = [0.,0.00,
         0.6,0.00,
	     1.0,.015,
	     1.6,0.025,
	     2.1,0.05]

# Note: from FST Lintul 3-rice, 30-10-15, SdV		
PHOTTB = [0.0,0.0,
          8.,0.0, 
	      10.,1.0,
	      12.,1.0, 
     	  13.,0.8, 
	      14.,0.6,
	      18.,0.0]		
		
		
# Leaf area correction function as a function of development stage, DVS. 
# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
SLACF = [0.0 ,1.72,
         0.21,1.72,
         0.24,1.72,
         0.33,1.32,
         0.7 ,1.20,
         1.01,1.00,
         2.0 ,0.75, 
         2.1 ,0.75]

# Maximum N concentration in the leaves, from which the N-conc.values of the
# stem and roots are derived, as a function of development stage.
# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
NMXLV = [0.0,0.05,
         0.4,0.05,
		 0.7,0.04,
	     1.0,0.03,
	     2.0,0.02,
	     2.1,0.02]

# ********** Partitioning coefficients ***********************************
# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
FRTTB = [0.0 ,0.30,
         0.48,0.30,
	     0.62,0.12,
	     0.69,0.11,
	     0.84,0.11,
	     0.92,0.10,
	     1.00,0.08,
	     1.38,0.00,
	     2.10,0.0]
      
# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
FLVTB = [0.0 ,0.315,
         0.48,0.35,
	     0.62,0.44, 
	     0.69,0.463,
	     0.84,0.463,
	     0.92,0.45, 
	     1.00,0.00,
	     1.38,0.00,
	     2.10, 0.0]

# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
FSTTB = [0.0 ,0.385,
         0.48,0.35,
	     0.62,0.44,
	     0.69,0.427,
	     0.84,0.427,
	     0.92,0.27,
	     1.00,0.00,
	     1.38,0.00,
	     2.1, 0.0]

# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
FSOTB = [0.0 ,0.00,
         0.48,0.00,
         0.62,0.00,
	     0.69,0.00,
	     0.84,0.00, 
	     0.92,0.18,
	     1.00,0.92, 
	     1.38,1.00,
	     2.1, 1.00]

# Calculate initial LAI
np_TSUMI       = TSUMI * np_One[:]
DVSI           = np_TSUMI[:]/TSUMAN
Interpol_SLACF = Afgen2(SLACF)       
SLACFI         = Interpol_SLACF(DVSI)
ISLA           = SLAC  * SLACFI
LAII           = WLVGI * ISLA

RAINSUM = 0.

# Initial soil moisture		 
#np_WAI                = 1000. * ROOTDI * WCI * np_One[:]		 
		 
def supplyCurrentTime(self):
    """
    *Optional*
      
      Supplies the current time in seconds after the start of the run
      This function is optional. If it is not set the framework assumes
      the model runs with daily timesteps.
      
      Output:
      
          - time in seconds since the start of the model run
          
      """
      
    return self.currentTimeStep(self) * int(configget(self.config,'model','timestepsecs','86400'))

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
      --------copied from wflow_sbm: Define all model parameters here that the framework should handle for the model
      See wf_updateparameters and the parameters section of the ini file
      If you use this make sure to all wf_updateparameters at the start of the dynamic section
      and at the start/end of the initial section 
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

      #Static model parameters
      #e.g.: modelparameters.append(self.ParamType(name="RunoffGeneratingGWPerc",stack="intbl/RunoffGeneratingGWPerc.tbl",type="static",default=0.1)) 
      #modelparameters.append(self.ParamType(name="Altitude",stack="staticmaps/wflow_dem.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),
      #modelparameters.append(self.ParamType(name="WCAD",stack="staticmaps/wcad.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),
      #modelparameters.append(self.ParamType(name="WCWP",stack="staticmaps/wcwp.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),	  
      #modelparameters.append(self.ParamType(name="WCFC",stack="staticmaps/wcfc.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),	  
      #modelparameters.append(self.ParamType(name="WCWET",stack="staticmaps/wcwet.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),	  
      #modelparameters.append(self.ParamType(name="WCST",stack="staticmaps/wcst.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),
      #modelparameters.append(self.ParamType(name="DRATE",stack="staticmaps/drate.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[]))
  
      # Meteo and other forcing
      #modelparameters.append(self.ParamType(name="Temperature",stack="inmaps/TEMP",type="timeseries",default=10.0,verbose=False,lookupmaps=[])),
      modelparameters.append(self.ParamType(name="IRRAD",stack="inmaps/IRRAD",type="timeseries",default=11.0,verbose=False,lookupmaps=[])),
      modelparameters.append(self.ParamType(name="TMIN",stack="inmaps/TMIN",type="timeseries",default=10.0,verbose=False,lookupmaps=[])),
      modelparameters.append(self.ParamType(name="TMAX",stack="inmaps/TMAX",type="timeseries",default=10.0,verbose=False,lookupmaps=[])),
      modelparameters.append(self.ParamType(name="T",stack="inmaps/T",type="timeseries",default=10.0,verbose=False,lookupmaps=[])),
      modelparameters.append(self.ParamType(name="VAP",stack="inmaps/VAP",type="timeseries",default=10.0,verbose=False,lookupmaps=[])),
      modelparameters.append(self.ParamType(name="WIND",stack="inmaps/WIND",type="timeseries",default=2.0,verbose=False,lookupmaps=[])),
      modelparameters.append(self.ParamType(name="RAIN",stack="inmaps/RAIN",type="timeseries",default=2.0,verbose=False,lookupmaps=[])),  	  
      return modelparameters

  def stateVariables(self):
      """ 
      *Required*
      
      Returns a list of state variables that are essential to the model. 
      This list is essential for the resume and suspend functions to work.
      
      This function is specific for each model and **must** be present. This is
      where you specify the state variables of you model. If your model is stateless
      this function must return and empty array (states = [])
      
      """
      #states = ['Test', 'LAI', 'WLVG', 'WLVD', 'WST', 'WSO', 'WRT', 'ROOTD', 'WDRT', 'TSUM', 'DAY', 'DVS', 'WA']
      states = ['Test', 'LAI', 'WLVG', 'WLVD', 'WST', 'WSO', 'WRT', 'ROOTD', 'ROOTD_mm', 'WDRT', 'TSUM', 'DAY', 'DVS']
	  
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
        
    self.logger.info("Saving all variables...")
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
    
    """
    #: pcraster option to calculate with units or cells. Not really an issue
    #: in this model but always good to keep in mind.
    setglobaloption("unittrue")

    
    self.timestepsecs = int(configget(self.config,'model','timestepsecs','86400'))
    self.basetimestep=86400
	
    # Reads all parameter from disk
    self.wf_updateparameters()
    self.logger.info("Starting LINTUL Dynamic Crop Growth Simulation...")
		 
    

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
	  
    ########################################################################################

  def dynamic(self):
     """
     *Required*
      
     This is where all the time dependent functions are executed. Time dependent
     output should also be saved here.
     """
     self.wf_updateparameters() # read the temperature map for each step (see parameters())
     
	# Create numpy array from states (make pointers): 
     np_Day                = pcr_as_numpy(self.DAY)
     np_DVS                = pcr_as_numpy(self.DVS)
     np_LAI                = pcr_as_numpy(self.LAI)
     np_WA                 = pcr_as_numpy(self.WA)
     np_TSUM               = pcr_as_numpy(self.TSUM)
     np_WRT                = pcr_as_numpy(self.WRT)
     np_WSO                = pcr_as_numpy(self.WSO)
     np_WST                = pcr_as_numpy(self.WST)
     np_WDRT               = pcr_as_numpy(self.WDRT)
     np_Test               = pcr_as_numpy(self.Test)
     np_WLVG               = pcr_as_numpy(self.WLVG)
     np_WLVD               = pcr_as_numpy(self.WLVD)
     np_ROOTD              = pcr_as_numpy(self.ROOTD)

    # Make numpy arrays from soil characteristics. Todo: later read them from PCRaster maps, then pcr_as_numpy. 
     np_WCAD               = WCAD  * np_One[:]
     np_WCWP               = WCWP  * np_One[:]
     np_WCFC               = WCFC  * np_One[:]
     np_WCWET              = WCWET * np_One[:]
     np_WCST               = WCST  * np_One[:]

#     runinfoFile           = "runinfo.xml"
#     runinfoFile_Lintul    = os.path.join(os.path.abspath("wflow_lintul"), "inmaps", runinfoFile)
#     starttime             = getStartTimefromRuninfo(runinfoFile_Lintul)
#     DOY_wflow             = starttime.strftime('%j')
     
	# Implement forcing data, coefficients and some handy numbers spatially, as numpy arrays:
     np_StartDay           = StartDay * np_One[:]
     np_StartNow           = np.equal(np_Day[:], np_StartDay[:] - np_One[:])    
     np_TMIN               = pcr_as_numpy(self.TMIN)
     np_TMAX               = pcr_as_numpy(self.TMAX)
     np_T                  = pcr_as_numpy(self.T)
     np_DAVTMP             = np_T #0.5     * (np_TMIN + np_TMAX)
     np_VAP                = pcr_as_numpy(self.VAP)
     np_WIND               = pcr_as_numpy(self.WIND)
     np_RAIN               = pcr_as_numpy(self.RAIN)
     np_TSUMAN             = TSUMAN  * np_One[:]
     np_TTSUM              = TTSUM   * np_One[:]
     np_LAICR              = LAICR   * np_One[:]
     np_NNI                = NNI     * np_One[:]
     np_WCWP               = WCWP    * np_One[:]     
     np_DVSDR              = DVSDR   * np_One[:]
     np_DELT               = DELT    * np_One[:]
     np_RGRL               = RGRL    * np_One[:]
     	 
    # Define tests for conditions that influence the behaviour of the crop: 
     DAVTMP                = numpy2pcr(Scalar, np_DAVTMP, -99) 
     Warm_Enough           = DAVTMP >= TBASE
     np_Warm_Enough        = np.greater_equal(np_DAVTMP[:], TBASE * np_One[:]) #!
     DegreeDay             = DAVTMP -  TBASE
     DTEFF                 = ifthenelse(Warm_Enough, DegreeDay, 0.)
     np_DTEFF              = pcr_as_numpy(DTEFF) 
	
     #np_Enough_water       = np.greater(np_WA[:], (WCWP * np_ROOTD[:] * 1000.))
     np_Enough_water       = np_One[:]
     Enough_water          = numpy2pcr(Boolean, np_Enough_water, -99)
     Leaves_Present        = self.LAI > 0.
     #CanGrowDownward       = self.ROOTD <= ROOTDM 
     #np_CanGrowDownward    = np.less_equal(np_ROOTD[:], ROOTDM * np_One[:])
     
	# Check when certain important decision moments are reached, in chronological order: 
     SimStart              = self.DAY == 0. #todo: do we need ROOTDI before emergence?
     np_Simstart           = pcr_as_numpy(SimStart)
     Just_Before_Start     = self.DAY == StartDay - 3.
     np_Just_Before_Start  = pcr_as_numpy(Just_Before_Start)
     DevStarted            = self.DAY >= StartDay - 1. # a bit mysterious, but the '-1.' keeps things in perfect sync with the FST version.
     np_DevStarted         = pcr_as_numpy(DevStarted)  
     BeforeAnthesis        = self.TSUM < TSUMAN
     np_BeforeAnthesis     = pcr_as_numpy(BeforeAnthesis)
     UntilAnthesis         = self.TSUM <= TSUMAN
     AtAndAfterAnthesis    = self.TSUM >= TSUMAN
     np_AtAndAfterAnthesis = pcr_as_numpy(AtAndAfterAnthesis)
     AfterAnthesis         = self.TSUM > TSUMAN
     Roots_Dying           = self.DVS >= DVSDR
     np_Roots_Dying        = pcr_as_numpy(Roots_Dying)
     
     Vegetative            = DevStarted & UntilAnthesis
     Generative            = DevStarted & AfterAnthesis
     EarlyStages           = self.DVS   <   0.2
     LaterStages           = self.DVS   >=  0.2
     SmallLeaves           = self.LAI   <   0.75
     BiggerLeaves          = self.LAI   >=  0.75
     Juvenile              = EarlyStages & SmallLeaves
     Adult                 = LaterStages | BiggerLeaves
     np_Juvenile           = pcr_as_numpy(Juvenile)
     np_Adult              = pcr_as_numpy(Adult)
	 
     TSUM_not_Finished     = self.TSUM <= TTSUM
     DVS_not_Finished      = self.DVS  <= 2.01
     Not_Finished          = TSUM_not_Finished & DVS_not_Finished
     np_Not_Finished       = pcr_as_numpy(Not_Finished)
  
     TSUMINIT              = numpy2pcr(Scalar, np_TSUMI[:], -99)
     One                   = numpy2pcr(Scalar, np_One[:],   -99)
     Zero                  = numpy2pcr(Scalar, np_Zero,     -99)
	 
    # Specific Leaf area(m2/g)
     Interpol_SLACF        = Afgen2(SLACF)       
     SLA                   = SLAC * np_One[:] * Interpol_SLACF(np_DVS[:]) #* exp(-NSLA * (1.-NNI))     should become numpy.exp
     
	# Calculate daylength, based on latitude and day of year:
	# Daylength assumed similar throughout the catchment area -> scalar, no array - SdV	 
     DAYL                  = astro2(np_Day.item((1,1)), LAT) 	 
     
     Interpol_FRTTB        = Afgen2(FRTTB)
     Interpol_FLVTB        = Afgen2(FLVTB)
     Interpol_FSTTB        = Afgen2(FSTTB)
     Interpol_FSOTB        = Afgen2(FSOTB)
     Interpol_RDRTB        = Afgen2(RDRTB)	
     Interpol_PHOTTB       = Afgen2(PHOTTB)
     
     FRTWET                = Interpol_FRTTB(np_DVS[:])
     FLVT                  = Interpol_FLVTB(np_DVS[:])
     FSTT                  = Interpol_FSTTB(np_DVS[:])
     FSOT                  = Interpol_FSOTB(np_DVS[:])
     RDRTMP                = Interpol_RDRTB(np_DVS[:]) 
     PHOTT                 = Interpol_PHOTTB(DAYL)     
	
     EMERG                 = DevStarted & Enough_water & Leaves_Present & Not_Finished
     np_EMERG              = pcr_as_numpy(EMERG)
     PHOTPF                = ifthenelse(BeforeAnthesis, PHOTT, One)
     RTSUMP                = DTEFF * PHOTPF
     self.TSUM             = self.TSUM + ifthenelse(Just_Before_Start, TSUMINIT, 0.) + ifthenelse(EMERG, RTSUMP, 0.) 
     
     TSUM_veg              = self.TSUM/TSUMAN
     TSUM_gen              = 1. + (self.TSUM - TSUMAN)/TSUMMT
     self.DVS              = ifthenelse(Vegetative, TSUM_veg, 0.) + ifthenelse(Generative, TSUM_gen, 0.)
	 # In TSUM_gen, '1' in fact stands for TSUM/TSUMAN, at the moment that TSUM = TSUMAN
	 
     np_IRRAD              = pcr_as_numpy(self.IRRAD) 

# Computation of the old PENMAN EQUATION, for testing/comparison purposes (to be replaced with P-M soon). SdV 30-11-2015
     #DTRJM2                = np_IRRAD[:] * 1000.
     #BOLTZM                = 5.668E-8
     #LHVAP                 = 2.4E6
     #PSYCH                 = 0.067 * np_One[:]

     #BBRAD                 = BOLTZM * (np_DAVTMP + 273. * np_One[:])**4 * 86400.
     #SVP                   = 0.611 * np.exp(17.4 * np_DAVTMP[:] / (np_DAVTMP[:] + 239. * np_One[:]))
     #SLOPE                 = 4158.6 * SVP[:] / (np_DAVTMP[:] + 239.)**2
     #RLWN                  = BBRAD * np.maximum(np_Zero[:],0.55*(1.- np_VAP[:]/SVP[:]))
     #NRADS                 = DTRJM2 * (1.-0.15) - RLWN[:]
     #NRADC                 = DTRJM2 * (1.-0.25) - RLWN[:]
     #PENMRS                = NRADS[:] * SLOPE[:]/(SLOPE[:]+PSYCH[:])
     #PENMRC                = NRADC[:] * SLOPE[:]/(SLOPE[:]+PSYCH[:])

     #WDF                   = 2.63 * (1.0 + 0.54 * np_WIND[:])
     #PENMD                 = LHVAP * WDF[:] * (SVP[:] - np_VAP[:]) * PSYCH[:]/(SLOPE[:]+PSYCH[:])
     #PEVAP                 = np.exp(-0.5 * np_LAI[:])  * (PENMRS[:] + PENMD[:]) / LHVAP
     #PTRAN                 = (np_One[:] - np.exp(-0.5 * np_LAI[:])) * (PENMRC[:] + PENMD[:]) / LHVAP
     #RNINTC                = np.minimum(np_RAIN[:], 0.25 * np_LAI[:])
     #PTRAN                 = np.maximum(np_Zero[:], PTRAN[:] - 0.5 * RNINTC[:])

	########################################################################################
    # Root depth growth:

     #Roots_can_grow        = Enough_water & BeforeAnthesis & EMERG & CanGrowDownward
     #RROOTD                = ifthenelse(Roots_can_grow, RRDMAX, Zero)
     np_CanGrowDownward    = np.less_equal(np_ROOTD[:], ROOTDM * np_One[:])
	 
     np_RROOTD             = np_Enough_water[:] * np_BeforeAnthesis[:] * np_EMERG[:] * np_CanGrowDownward[:] * RRDMAX 
     #np_ROOTD[:]           = np_Simstart[:] * ROOTDI + np_ROOTD[:] + RRDMAX * np_CanGrowDownward[:]
     #self.ROOTD           += ifthenelse(SimStart, ROOTDI, Zero) + RROOTD
     self.ROOTD            = self.ROOTD + ifthenelse(SimStart, ROOTDI, Zero) + numpy2pcr(Scalar, np_RROOTD[:], -99)
     self.ROOTD_mm         = self.ROOTD * 1000.
     np_EXPLOR             = 1000. * np_RROOTD[:] * np_WCFC[:] 

	# Subroutine EVAPTR
    # np_WC                 = 0.001 * np_WA[:]/ NOTNUL(np_ROOTD[:])
    # np_WAAD               = 1000. * np_WCAD[:] * np_ROOTD[:]
    # np_WAFC               = 1000. * np_WCFC[:] * np_ROOTD[:]
    # 
    # np_RelWatAvail4Evap   = (np_WC[:] - np_WCAD[:])/(np_WCFC[:] - np_WCAD[:])
    # np_EVAP               = PEVAP[:] * np.clip(np_RelWatAvail4Evap[:], np_Zero[:], np_One[:])
	# 
    # np_RelWatAvail4Tran   = (PTRAN[:]/(PTRAN[:] + TRANCO * np_One[:])) * (np_WCFC[:] - np_WCWP[:])
    # np_WCCR               = np_WCWP[:] + np.maximum(0.01, np_RelWatAvail4Tran)
    # 
    # WC_Above_WCCR         = np.greater(np_WC[:], np_WCCR[:])
    # WC_At_or_Below_WCCR   = np.less_equal(np_WC[:], np_WCCR[:])

    # np_RelMoist1          = (np_WCST[:] - np_WC[:])   / (np_WCST[:] - np_WCWET[:])
    # np_RelMoist2	       = (np_WC[:]   - np_WCWP[:]) / (np_WCCR[:] - np_WCWP[:])
    # FR                    = WC_Above_WCCR[:]       * np.clip(np_RelMoist1, np_Zero[:], np_One[:]) + \
    #                         WC_At_or_Below_WCCR[:] * np.clip(np_RelMoist2, np_Zero[:], np_One[:])
    # np_TRAN               = PTRAN[:] * FR[:]
    # np_EVAPOTRAN          = np_EVAP[:] + np_TRAN[:]
	# 
    # np_RelMoist3          = ((np_WA[:] - np_WAAD[:])/DELT) / (NOTNUL(np_EVAPOTRAN[:]))
    # np_AVAILF             = np.minimum(np_One[:], np_RelMoist3[:])
    # np_EVAP[:]           *= np_AVAILF[:]
    # np_TRAN[:]           *= np_AVAILF[:]
     
     #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! this enables the original FST water balance:
    # TRANRF                = np_TRAN[:] / NOTNUL(PTRAN[:])
    # TRANRF      = np_One[:]
	
	
	
	
	
	
	
     #BMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMI
	 #############################################################################################################
	 # Voor BMI nodig om eerst self.Transpiration, self.PotTrans om te zetten naar numpy.
	 # Niet zeker of "self." nodig is bij variabelen die worden uitgewisseld via de BMI.
	 
     np_Transpiration = pcr_as_numpy(self.Transpiration) 	 
     np_PotTrans      = pcr_as_numpy(self.PotTrans)
     TRANRF           = np_Transpiration[:] / NOTNUL(np_PotTrans[:])
     #TRANRF           = np_One[:]
	 
	 #############################################################################################################
	##BMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMIBMI
    
	
	
	
	
	
	 
	 
    # Water Limitation: effects on partitioning 
     FRTMOD                = np_One[:] #np.maximum(np_One[:], np_One[:]/(TRANRF[:]+ 0.5 * np_One[:]))
     FRTMOD                = np.maximum(np_One[:], np_One[:]/(TRANRF[:]+ 0.5 * np_One[:]))
     #FRT                   = FRTMOD[:] #FRTWET[:] * FRTMOD[:]
     FRT                   = FRTWET[:] * FRTMOD[:]
     FSHMOD                = (1.-FRT[:]) / (1.-(FRT[:]/FRTMOD[:]))
     FLV                   = FLVT[:] * FSHMOD
     FST                   = FSTT[:] * FSHMOD
     FSO                   = FSOT[:] * FSHMOD	 
     PartitionCheck        = FLV[:] + FST[:] + FSO[:] + FRT[:]


     
	
	# Subroutine DRUNIR:
     #WC                    = 0.001 * WA/(self.ROOTD + TINY)
     #np_WAFC               = 1000. * np_WCFC[:] * np_ROOTD[:]
     #np_WAST               = 1000. * np_WCST[:] * np_ROOTD[:]
     #np_DRAIN_tmp          = (np_WA[:] - np_WAFC[:])/DELT + (np_RAIN[:] - RNINTC[:] - np_EVAP[:] - np_TRAN[:])
     #np_DRAIN              = np.clip(np_DRAIN_tmp, 0., DRATE)  
	 
     #np_RUNOFF_tmp         = (np_WA[:] - np_WAST[:])/DELT + (np_RAIN[:] - RNINTC[:] - np_EVAP[:] - np_TRAN[:] - np_DRAIN[:])
     #np_RUNOFF             = np.maximum(np_Zero[:], np_RUNOFF_tmp)
	 # 
     # np_RWC                = np_RAIN[:] + np_EXPLOR[:] + IRRIG - RNINTC[:] - np_RUNOFF[:] - np_TRAN[:] - np_EVAP[:] - np_DRAIN[:]
     
     #np_WA[:]             += np_Simstart[:] * np_WAI[:] + np_RWC[:]
     #np_WA[:]             += np_Simstart[:] * np_WAI[:] 
     
	 
    # todo: incorporate effects of N stress - sdv 30-11-15
     PARINT                = 0.5 * np_IRRAD[:] * 0.001 * (1.- np.exp(-K * np_LAI[:])) * np_Not_Finished[:]
     GTOTAL                = LUE * PARINT[:] * TRANRF[:]
    #FRT, FLV, FST, FSO    = dryMatterPartitioningFractions(self, NPART, TRANRF, NNI, FRTWET[:], FLVT[:], FSTT[:], FSOT[:])

    
    
	 
    #----------------------------------------------------------------------
	# Rel. Death rate due to ageing:
     RDRDV                 = np_AtAndAfterAnthesis[:] * RDRTMP[:]
     RDRSH                 = np.maximum(np_Zero[:], RDRSHM * (np_LAI[:] - np_LAICR[:])/ np_LAICR[:])
     RDR                   = np.maximum(RDRDV[:], RDRSH[:]) * np_Not_Finished[:]
	 
	# (Abs.) impact of leaf dying on leaf weight - todo
     N_Limitation          = np.less(np_NNI[:], np_One[:])
     DLVNS                 = np_DevStarted[:] * N_Limitation[:] * np_WLVG[:] * RDRNS * (1. - np_NNI[:])
     DLVS                  = np_WLVG[:]    * RDR[:]
     DLV                   = (DLVS[:]        + DLVNS[:]) * np_Not_Finished[:]
	 
     RWLVG                 = GTOTAL[:] * FLV[:] - DLV[:]
     np_WLVG[:]           += np_StartNow[:] * WLVGI + RWLVG[:] 
     np_WLVD[:]           += DLV[:]
	 
    #----------------------------------------------------------------------
	# Leaf totalGrowthRate and LAI.
     GLV                   = FLV[:] * GTOTAL[:]
     GLV_pcr               = numpy2pcr(Scalar, GLV, -99)
	 
     GLAI = np_Adult[:] * SLA[:] * GLV[:] + \
            np_Juvenile[:]     * (np_LAI[:] * (np.exp(np_RGRL[:] * np_DTEFF[:] * np_DELT[:]) - np_One[:])/ np_DELT[:] )* TRANRF[:] * np.exp(-NLAI* (1.0 - NNI)) + \
	        np.equal(np_LAI[:], np_Zero[:]) * np_Enough_water[:] * np_DevStarted[:] * LAII / DELT  
     
	# (Abs.) impact of leaf dying on LAI
	# Death of leaves due to ageing and shading:
     DLAIS                 = np_LAI[:]     * RDR[:]
     DLAINS                = np_DevStarted[:] * N_Limitation[:] * DLVNS[:] * SLA[:]
     DLAI                  = (DLAIS         + DLAINS) * np_Not_Finished[:]
     np_LAI[:]            += GLAI[:] - DLAI[:] 
     

	#----------------------------------------------------------------------	 
	#DRRT = np.greater_equal(np_DVS[:], np_DVSDR[:]) * np_WRT[:] * RDRRT
     DRRT                  = np_Roots_Dying[:] * np_WRT[:] * RDRRT
     RWRT                  = GTOTAL[:] * FRT[:] - DRRT[:] * np_Not_Finished[:]
     np_WRT[:]            += np_StartNow[:] * WRTLI + RWRT[:] 	
     np_WDRT[:]           += DRRT[:]
	 
    #----------------------------------------------------------------------
     RWSO                  = GTOTAL[:] * FSO[:] 
     np_WSO[:]            += np_StartNow[:] * WSOI + RWSO[:]
     np_WSOTHA             = np_WSO[:]/100.
	 
    #----------------------------------------------------------------------
     RWST                  = GTOTAL[:] * FST[:] 
     np_WST[:]            += np_StartNow[:] * WSTI + RWST[:]
 
    #----------------------------------------------------------------------
     np_WLV                = np_WLVG[:] + np_WLVD[:]
     TAGBM                 = np_WLV[:] + np_WST[:] + np_WSO[:]
    
	#----------------------------------------------------------------------
     
     self.Test            = numpy2pcr(Scalar, np_PotTrans[:], -99)
     bla                  = self.timestepsecs/self.basetimestep

     self.DAY              += 1. 
     DOY_wflow = iDOY_wflow + self.currentTimeStep() - 1

     np_ROOTD_mm = pcr_as_numpy(self.ROOTD_mm)
     
	 
     print "********************************************"
     print self.currentTimeStep(), DOY_wflow, np_Day[1,1]
     print starttime
     #print "LAI      via Lintul",  "\n"
     #print np_LAI[:] , "\n"
     #print "ROOTD    via Lintul",  "\n"
     #print cellvalue(self.ROOTD, 1), "\n"
     #print np_ROOTD_mm[:], "\n"
     #print "PotTrans   via Lintul",  "\n"
     #print np_PotTrans[0,:5],  "\n"
     #print "np_EMERG via Lintul"
     #print cellvalue(EMERG, 1)#[1,1], np_Day[1,1]
     #print "WA       via Lintul",  "\n"
     #print np_WA, "\n"
	 
    #----------------------------------------------------------------------
  
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
    configfile="wflow_sceleton.ini"
    _lastTimeStep = 150
    _firstTimeStep = 1
    timestepsecs=86400
    wflow_cloneMap = 'wflow_subcatch.map'
    
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
