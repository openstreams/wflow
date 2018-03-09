#!/usr/bin/python
# 

from math import pi
import math
import numpy as np
import os
import os.path
import shutil, glob
import getopt

from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *
from time import strftime
import time

"""

wflow_lintul simulates potential or water-limited rice production using weather data as forcing variables. For water-limited production,  
soil data are required. wflow_lintul is a modified version of LINTUL1 (for simulating potential growth) and
LINTUL2 (for water-limited production). Rice-specific features were derived from LINTUL3 (N-limited production; Shibu et al., 2010).

Usage:
wflow_lintul  -C case -R Runid -c inifile

    -C: set the name  of the case (directory) to run

    -R: set the name runId within the current case

    -c name of the config file (in the case directory)

$Author: SanderCdeVries $
$Id: lintul1.py 2018-02-21 16:08:06Z
$Rev: 001 $
"""

# This needs to be set according to the geographic extent (map dimensions) of your study area/catchment:
np_Zero = numpy.zeros((219, 286))
np_One = numpy.ones((219, 286))

# Some last remaining hardcoded (& partly non-functional) parameters:
DELT  = 1.   # Time step (delta T) = 1 day; changing it is not recommended.
TINY  = 1e-6 # A tiny number (from the original LINTUL code:)
WCWP  = 0.21 # Volumetric soil water content at wilting point... how soil specific is this for puddled soils...? todo
WCFC  = 0.47 # Volumetric soil water content at field capacity... how soil specific is this for puddled soils...? todo
WCST  = 0.55 # Volumetric soil water content at saturation (normal condition for irrigated rice soil) ... how soil specific is this for puddled soils...? todo
NNI   = 1.   # Nitrogen Nutrition Index (non-functional, for future development)
NPART = 1.0  # Coefficient for the effect of N stress on leaf biomass reduction (presently non-functional, for future development)
NSLA  = 1.0  # Coefficient for the effect of N stress on SLA reduction (presently non-functional, for future development)
NLAI  = 1.0  # Coefficient for the effect of N stress on LAI reduction(during juvenile phase; presently non-functional, for future development)

Point_Output = open('Point_Output.csv', 'w')

def NOTNUL(matrix):
    """
    NOTNUL was originally a FST Fortran Simulation Translator intrinsic function.
    Here it is applied to arrays. If a value in the array is positive, NOTNUL will 
    just return the value as is. If it equals zero: NOTNUL will return a value of 1 instead.
    Sander de Vries, 01-12-2015
    """
    b = np.shape(matrix)
    nullen = np.zeros(b)
    NOTNUL_check = np.equal(nullen, matrix)
    matrix += NOTNUL_check[:]
    return matrix

def NOTNUL_pcr(pcr_map):
    """
    NOTNUL was originally a FST Fortran Simulation Translator intrinsic function.
    Here it is applied to arrays. If a value in the array is positive, NOTNUL will 
    just return the value as is. If it equals zero: NOTNUL will return a value of 1 instead.
    Sander de Vries, 01-12-2015
    """
    checkzeros = pcr_map == False
    map = ifthenelse(checkzeros, 1., scalar(pcr_map))
    return pcr_map

    
def astro2(DAY, LAT):
    """
* ---------------------------------------------------------------------*
*   SUBROUTINE ASTRO                                                    *
*   Purpose: This subroutine calculates the astronomic daylength,       *
*            based on Goudriaan and van Laar 1994, around p.30.         *
* ---------------------------------------------------------------------*

    Originally: Daniel van Kraalingen, April 1991
    Python version: Sander de Vries, 30-10-2015
    """

    # SINE AND COSINE OF LATITUDE
    sinLAT = math.sin(pi * LAT / 180.)
    cosLAT = math.cos(pi * LAT / 180.)

    # MAXIMAL SINE OF DECLINATION
    sinDCM = math.sin(math.pi * 23.45 / 180.)

    # SINE AND COSINE OF DECLINATION  (EQUATIONS 3.4, 3.5)

    # SINDEC = -SINDCM * cos(2.* PI * (DAY+10.)/365.)
    # The '9' below (instead of 10) keeps things in sync with FST...
    # Todo: try to understand this at a certain point... for now it works perfectly.
    sinDEC = -sinDCM * math.cos(2. * pi * (DAY + 11.) / 365.)
    cosDEC = math.sqrt(1. - sinDEC * sinDEC)

    # THE TERMS A AND B ACCORDING TO EQUATION 3.3

    A = sinLAT * sinDEC
    B = cosLAT * cosDEC

    # DAYLENGTH ACCORDING TO EQUATION 3.6. Make sure not to use asin from PCRaster...
    DAYL = 12. * (1. + (2. / pi) * math.asin(A / B))

    return DAYL        
    
    
class Interpol_Obj(object):
# not used in present version... possibly for future development.
    """
    Class to facilitate use of the 'lookuplinear' PCraster function. 
    Upon initialization of an interpolation object, a temporary file 
    containing x, y value pairs is created and saved in the case directory. 
    This file is accessed by (PCraster) lookuplinear if the lookup_linear method is called.
    """
    def __init__(self, name):
        self.data = name
        self.name = name[-1]
        self.filename = self.name + ".tmp"        
        
        temptablefile = open(self.filename, 'w')
        index = range(0, len(self.data)-1)
        for i in index:
            if i < (len(self.data)-1):
                if i > i + 1:
                    print "x values of lookuplinear table not sorted in strictly ascending order..."
            if i//2. - i/2. <> 0.:
                string = str(self.data[i]) + " "
            else: 
                string = '\n' + str(self.data[i]) + " " 
            temptablefile.write(string)
        temptablefile.close()
        #print "method done"
    def lookup_linear(self, x):
        y = lookuplinear(self.filename, x)
        return y
    
    
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
            if x[i] <= x[i - 1]:
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
        for (x, y) in zip(self.x, self.y):
            msg += ("(%f,%f)\n " % (x, y))
        msg += "\n"
        if self.unit is not None:
            msg += "Return value as unit: %s" % self.unit
        return msg


def supplyCurrentTime(self):
    """
    *Optional*

      Supplies the current time in seconds after the start of the run
      This function is optional. If it is not set the framework assumes
      the model runs with daily timesteps.

      Output:

          - time in seconds since the start of the model run

      """

    return self.currentTimeStep(self) * int(configget(self.config, 'model', 'timestepsecs', '86400'))


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)


class WflowModel(DynamicModel):
    """
    The user defined model class. This is your work!
    """

    def __init__(self, cloneMap, Dir, RunDir, configfile):
        """
        *Required*

        The init function **must** contain what is shown below. Other functionality
        may be added by you if needed.

        """
        DynamicModel.__init__(self)
        setclone(Dir + "/staticmaps/" + cloneMap)
        self.runId = RunDir
        self.caseName = Dir
        self.Dir = Dir
        self.configfile = configfile
        self.SaveDir = os.path.join(self.Dir,self.runId)

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
        
        self.RainSumStart_Month   = int(configget(self.config, "model", "RainSumStart_Month", "11"))        
        self.RainSumStart_Day     = float(configget(self.config, "model", "RainSumStart_Day", "1"))        
        self.Sim3rdSeason         = eval(configget(self.config, "model", "Sim3rdSeason", "False"))
        self.RainSumReq           = float(configget(self.config, "model", "RainSumReq", "200."))        
        self.Pause                = int(configget(self.config, "model", "Pause", "13"))        
        self.AutoStartStop        = eval(configget(self.config, "model", "AutoStartStop", "False")) #default changed to 'False', for running from 'crop profile' maps (CRPST.xxx) under DEWS. sdv 21-2-2018
        self.WATERLIMITED         = (configget(self.config, "model", "WATERLIMITED", "True"))
        self.CropStartDOY         = int(configget(self.config, "model", "CropStartDOY", "0")) - 1 # to keep things in sync with the original LINTUL version in FST
        self.HarvestDAP           = int(configget(self.config, "model", "HarvestDAP", "150"))
        self.LAT                  = float(configget(self.config, "model", "LAT", "3.16"))        
        self.TSUMI                = float(configget(self.config, "model", "TSUMI", "362."))        
        self.K                    = float(configget(self.config, "model", "K", "0.6"))
        self.LUE                  = float(configget(self.config, "model", "LUE", "2.47")) # The default value from Shibu et al. (2010) is 3.0; 2.47 was obtained by calibration for central Java. (sdv)
        self.SLAC                 = float(configget(self.config, "model", "SLAC", "0.02"))
        self.TSUMAN               = float(configget(self.config, "model", "TSUMAN", "1420."))
        self.TSUMMT               = float(configget(self.config, "model", "TSUMMT", "580."))
        self.TBASE                = float(configget(self.config, "model", "TBASE", "8."))
        self.RGRL                 = float(configget(self.config, "model", "RGRL", "0.009"))
        self.WLVGI                = float(configget(self.config, "model", "WLVGI", "0.86"))
        self.WSTI                 = float(configget(self.config, "model", "WSTI", "0.71"))
        self.WRTLI                = float(configget(self.config, "model", "WRTLI", "1.58"))
        self.WSOI                 = float(configget(self.config, "model", "WSOI", "0."))
        self.RDRNS                = float(configget(self.config, "model", "RDRNS", "0.03"))
        self.DVSDR                = float(configget(self.config, "model", "DVSDR", "0.8"))
        self.RDRRT                = float(configget(self.config, "model", "RDRRT", "0.03"))
        #self.TTSUM                = float(configget(self.config, "model", "TTSUM", "2000."))
        self.RDRSHM               = float(configget(self.config, "model", "RDRSHM", "0.03"))
        self.LAICR                = float(configget(self.config, "model", "LAICR", "4."))
        self.ROOTDM_mm            = float(configget(self.config, "model", "ROOTDM_mm", "1000."))
        self.RRDMAX_mm            = float(configget(self.config, "model", "RRDMAX_mm", "10."))
        self.ROOTDI_mm            = float(configget(self.config, "model", "ROOTDI_mm", "50."))
        self.NLAI                 = float(configget(self.config, "model", "NLAI", "1."))
        self.RDRTB                = eval(configget(self.config, "model", "RDRTB", "[0.0, 0.00 , 0.6 , 0.00, 1.0 , .015, 1.6 , 0.025, 2.1 , 0.05]"))
        self.PHOTTB               = eval(configget(self.config, "model", "PHOTTB", "[0.0, 0.0  , 8.  , 0.0 , 10. , 1.0 , 12. , 1.0  , 13. , 0.8  , 14.,  0.6 , 18. , 0.0]"))
        self.SLACF                = eval(configget(self.config, "model", "SLACF", "[0.0, 1.72 , 0.21, 1.72, 0.24, 1.72, 0.33, 1.32 , 0.7 , 1.20 , 1.01, 1.00, 2.0 , 0.75, 2.1 , 0.75]"))
        #self.SLACF2              = eval(configget(self.config, "model", "SLACF2", "[0.0, 1.72 , 0.21, 1.72, 0.24, 1.72, 0.33, 1.32 , 0.7 , 1.20 , 1.01, 1.00, 2.0 , 0.75, 2.1 , 0.75]")) #for testing/development
        self.NMXLV                = eval(configget(self.config, "model", "NMXLV", "[0.0, 0.05 , 0.4 , 0.05, 0.7 , 0.04, 1.0 , 0.03 , 2.0 , 0.02 , 2.1 , 0.02]"))
        self.FRTTB                = eval(configget(self.config, "model", "FRTTB", "[0.0, 0.300, 0.48, 0.30, 0.62, 0.12, 0.69, 0.11 , 0.84, 0.11 , 0.92, 0.10, 1.00, 0.08, 1.38, 0.00, 2.10, 0.0]"))
        #self.FRTTB2              = eval(configget(self.config, "model", "FRTTB2", "[0.0, 0.300, 0.48, 0.30, 0.62, 0.12, 0.69, 0.11 , 0.84, 0.11 , 0.92, 0.10, 1.00, 0.08, 1.38, 0.00, 2.10, 0.0]")) #for testing/development
        self.FLVTB                = eval(configget(self.config, "model", "FLVTB", "[0.0, 0.315, 0.48, 0.35, 0.62, 0.44, 0.69, 0.463, 0.84, 0.463, 0.92, 0.45, 1.00, 0.00, 1.38, 0.00, 2.10, 0.0]"))
        #self.FLVTB2              = eval(configget(self.config, "model", "FLVTB2", "[0.0, 0.315, 0.48, 0.35, 0.62, 0.44, 0.69, 0.463, 0.84, 0.463, 0.92, 0.45, 1.00, 0.00, 1.38, 0.00, 2.10, 0.0]")) #for testing/development
        self.FSTTB                = eval(configget(self.config, "model", "FSTTB", "[0.0, 0.385, 0.48, 0.35, 0.62, 0.44, 0.69, 0.427, 0.84, 0.427, 0.92, 0.27, 1.00, 0.00, 1.38, 0.00, 2.10, 0.0]"))
        self.FSOTB                = eval(configget(self.config, "model", "FSOTB", "[0.0, 0.00 , 0.48, 0.00, 0.62, 0.00, 0.69, 0.00 , 0.84, 0.00 , 0.92, 0.18, 1.00, 0.92, 1.38, 1.00, 2.10, 1.00]"))

        # Static model parameters
        # modelparameters.append(self.ParamType(name="RunoffGeneratingGWPerc",stack="intbl/RunoffGeneratingGWPerc.tbl",type="static",default=0.1))
        # modelparameters.append(self.ParamType(name="Altitude",stack="staticmaps/wflow_dem.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),
        # modelparameters.append(self.ParamType(name="WCAD",stack="staticmaps/wcad.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),
        # modelparameters.append(self.ParamType(name="WCWP",stack="staticmaps/wcwp.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),
        # modelparameters.append(self.ParamType(name="WCFC",stack="staticmaps/wcfc.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),
        # modelparameters.append(self.ParamType(name="WCWET",stack="staticmaps/wcwet.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),
        # modelparameters.append(self.ParamType(name="WCST",stack="staticmaps/wcst.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[])),
        # modelparameters.append(self.ParamType(name="DRATE",stack="staticmaps/drate.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[]))

        # Meteo and other forcing
        modelparameters.append(self.ParamType(name="IRRAD", stack="inmaps/IRRAD", type="timeseries", default=11.0, verbose=False,lookupmaps=[])),
        modelparameters.append(self.ParamType(name="T", stack="inmaps/T", type="timeseries", default=10.0, verbose=False, lookupmaps=[])),
        modelparameters.append(self.ParamType(name="TMIN",stack="inmaps/TMIN",type="timeseries",default=10.0,verbose=False,lookupmaps=[])),
        modelparameters.append(self.ParamType(name="TMAX",stack="inmaps/TMAX",type="timeseries",default=10.0,verbose=False,lookupmaps=[])),
        modelparameters.append(self.ParamType(name="RAIN", stack="inmaps/P", type="timeseries", default=0., verbose=False, lookupmaps=[])),
        modelparameters.append(self.ParamType(name="CRPST", stack="inmaps/CRPST", type="timeseries", default=11.0, verbose=False, lookupmaps=[])),
        return modelparameters

    def stateVariables(self):
        """
        *Required*

        Returns a list of state variables that are essential to the model.
        This list is essential for the resume and suspend functions to work.

        This function is specific for each model and **must** be present. This is
        where you specify the state variables of you model. If your model is stateless
        this function must return and empty array (states = [])
        
        'Test' is a handy variable to keep if you quickly want to get map output for a certain variable 
        (e.g. include a statement self.Test = OutputVariable; sdv)

        """

        states = ['Season', 'PSUM', 'Test', 'LAI', 'WLVG', 'WLVD', 'WST', 'WSO', 'WRT', 'ROOTD_mm', 'WDRT', 'TSUM',
                  'STARTED', 'DVS']

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

        return self.currentTimeStep() * int(configget(self.config, 'model', 'timestepsecs', '86400'))

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
        self.wf_suspend(self.SaveDir + "/outstate/")

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

        self.timestepsecs  = int(configget(self.config, 'model', 'timestepsecs', '86400'))
        self.basetimestep  = 86400

      # Reads all parameter from disk
        self.wf_updateparameters()
        self.logger.info("Starting LINTUL Dynamic Crop Growth Simulation...")
        
      # Read a static map of the rice area. To be replaced with real-time radar images of rice area in the future? Todo
      # Simulation is mostly restricted to the rice area (to be checked), which saves calculation time. Todo
        wflow_ricemask     = configget(self.config, "model", "wflow_ricemask", "staticmaps/wflow_ricemask.map")
        self.ricemask      = self.wf_readmap(os.path.join(self.Dir,wflow_ricemask),0.0,fail=True)
        self.ricemask_BOOL = boolean(self.ricemask)
        self.Pausedays     = self.Pause + 1
        
      # Calculate initial development stage (at the time of transplanting)
        self.DVSI               = self.TSUMI / self.TSUMAN
      # Calculate the initial leaf area correction function as a function of development stage, DVS. 
        Interpol_SLACF     = Afgen2(self.SLACF)
        SLACFI             = Interpol_SLACF(self.DVSI)
      # Multiply with specific leaf area constant => initial specific leaf area
        ISLA               = self.SLAC * SLACFI
      # Multiply with weight of green leaves to obtain initial LAI
        self.LAII          = self.WLVGI * ISLA
      # Calculate total temperature sum from transplanting to crop maturity:
        self.TTSUM = self.TSUMAN + self.TSUMMT      
      #For future development - Todo - sdv
        #self.SLACF2        = Interpol_Obj(self.SLACF2)
        #self.FLVTB2        = Interpol_Obj(self.FLVTB2)
        #self.FRTTB2        = Interpol_Obj(self.FRTTB2)
        #SLACFI2            = self.SLACF2.lookup_linear(DVSI)
        #report(SLACFI2, "SLACFI2")
        #self.SLAF2.naam()
        
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
        if self.reinit == 1:
            self.logger.info("Setting initial conditions to default")
            for s in self.stateVariables():
                exec "self." + s + " = cover(0.0)"
        else:
            self.wf_resume(self.Dir + "/instate/")


            # try:
            #    self.wf_resume(self.Dir + "/instate/")
            # except:
            #    self.logger.warn("Cannot load initial states, setting to default")
            #    for s in self.stateVariables():
            #        exec "self." + s + " = cover(1.0)"

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
        self.wf_updateparameters()
        
      # Get the date as a Python datetime object and as the day of the year (DOY): used for cropping calendar and daylength.
        self.date         = datetime.utcfromtimestamp(self.wf_supplyStartTime()) + dt.timedelta(self.currentTimeStep() - 1)      
        self.enddate     = datetime.utcfromtimestamp(self.wf_supplyEndTime()) #wf_supplyEndTime() in wflow_dyamicframework? todo
        DOY               = self.wf_supplyJulianDOY() 

      # Some Boolean PCRaster variables:
        TSUM_not_Finished = self.TSUM <= self.TTSUM
        DVS_not_Finished  = self.DVS <= 2.01
        Not_Finished      = TSUM_not_Finished #& DVS_not_Finished
        np_Not_Finished   = pcr_as_numpy(Not_Finished)
        np_STARTED        = pcr_as_numpy(self.STARTED)
        np_ricemask       = pcr_as_numpy(self.ricemask)

      # Start calculating the accumulated preciptation from November one on, 
      # to judge when there's enough water for crop establishment. 
        if (self.date.month == self.RainSumStart_Month and self.date.day == self.RainSumStart_Day):  
            Calc_RainSum     = True
            self.PSUM       += self.RAIN + TINY
        else:
            Calc_RainSum     = False
       
      # Check whether the precipitation sum is positive
        WeveGotRain       = self.PSUM > 0.
      # Check whether the precipitation sum is still below the threshhold for crop establishment
        NotEnoughRainYet  = self.PSUM <= self.RainSumReq
        EnoughRain        = self.PSUM >= self.RainSumReq
        KeepAddingRain    = WeveGotRain & NotEnoughRainYet
      # The first season is defined here as starting on November 1
        #self.Season      += ifthenelse(Calc_RainSum, self.ricemask, 0.)
        self.Season      += ifthenelse(EnoughRain, self.ricemask, 0.) # beware, this variable is also modified in another equation
        FirstSeason       = self.Season == 1
        SecondSeason      = self.Season == 2
        ThirdSeason       = self.Season == 3
      # Add rain when the precipitation sum is positive but still below the threshold for crop establishment, reset to 0. when this is no longer the case.
        self.PSUM         = (self.PSUM + ifthenelse(KeepAddingRain, self.RAIN, 0.)) * ifthenelse(KeepAddingRain, scalar(1.), 0.)
        #pcr_PrepareField_temp = scalar(self.Pausedays)
        #pcr_PrepareField      = ifthenelse(FirstSeason, pcr_PrepareField_temp, 0.)

      # Create numpy array from states:
        np_CRPST          = pcr_as_numpy(self.CRPST)
        np_IRRAD          = pcr_as_numpy(self.IRRAD)
        np_DVS            = pcr_as_numpy(self.DVS)
        np_LAI            = pcr_as_numpy(self.LAI)
        np_WA             = pcr_as_numpy(self.WA)
        np_TSUM           = pcr_as_numpy(self.TSUM)
        np_WRT            = pcr_as_numpy(self.WRT)
        np_WSO            = pcr_as_numpy(self.WSO)
        np_WST            = pcr_as_numpy(self.WST)
        np_WDRT           = pcr_as_numpy(self.WDRT)
        np_Test           = pcr_as_numpy(self.Test)
        np_PSUM           = pcr_as_numpy(self.PSUM)
        np_WLVG           = pcr_as_numpy(self.WLVG)
        np_WLVD           = pcr_as_numpy(self.WLVD)
        np_ROOTD_mm       = pcr_as_numpy(self.ROOTD_mm)
        np_Season         = pcr_as_numpy(self.Season)


        
      # Initializing crop harvest:
      # If a fixed planting and a fixed harvest date are forced for the whole catchment:
        if self.CropStartDOY   > -1:        
            if self.HarvestDAP > 0:
               #np_CropHarvNow     = np.greater_equal(DOY, self.CropStartDOY + self.HarvestDAP) * np_One[:]
                np_CropHarvNow     = np.greater_equal(DOY, self.CropStartDOY + self.HarvestDAP) * np_ricemask[:]
                CropHarvNow        = numpy2pcr(Boolean, np_CropHarvNow, -99)    
                print "Warning: harvest date read from ini file, not from Crop Profile map..."
            elif  self.HarvestDAP == 0:
                np_CropHarvNow     = (1 - np_Not_Finished)* np_ricemask[:]
                HarvNow            = Not_Finished == False
                CropHarvNow        = HarvNow & self.ricemask_BOOL
                print "Harvest date not specified; crop harvest at crop maturity"
            else: 
                print "Crop harvest not initialized, found strange values in ini file... CTRL + C to exit..."
                
            # Initializing crop growth, optionally from a single start day (CropStartDOY in the ini file),
            # but normally from a crop profile forcing variable.
            np_CropStartNow    = np.equal(DOY, self.CropStartDOY) * np_ricemask[:]
            CropStartNow       = numpy2pcr(Boolean, np_CropStartNow, -99)
            CropStartNow_scalar= scalar(CropStartNow)
            #np_CropStarted    = np.greater_equal(DOY, self.CropStartDOY) * np_ricemask[:]
            #CropStarted       = numpy2pcr(Boolean, np_CropStarted, -99)
            Started            = self.STARTED > 0
            CropStarted        = Started & self.ricemask_BOOL
            np_CropStarted     = pcr_as_numpy(CropStarted)
            self.STARTED       = (self.STARTED + CropStartNow_scalar + scalar(CropStarted)) * ifthenelse(CropHarvNow, scalar(0.), 1.) 
            print "Warning: using start date from ini file, not read from Crop Profile..."
                
      # If planting is initiated gridcell-by-gridcell when a pre-determined rainfall requirement is met OR based on a (remotely sensed) crop profile map,
      # and a fixed harvest date is used across all gridcells:
        elif self.CropStartDOY == -1:
            if self.HarvestDAP > 0:
                HarvNow            = self.STARTED == self.HarvestDAP
                CropHarvNow        = HarvNow & self.ricemask_BOOL
                np_CropHarvNow     = pcr_as_numpy(CropHarvNow)
                
      # If planting is initiated gridcell-by-gridcell based on accumulated rainfall OR a (remotely sensed) crop profile map and crop harvest occurs at crop maturity (may vary depending on ambient temperatures during the growing season):            
            elif self.CropStartDOY == -1 and self.HarvestDAP == 0:
                if self.AutoStartStop == False:
                    started_gt_zero    = self.STARTED > 0.
                    crpprfl_eq_zero    = self.CRPST == 0.
                    #CropHarvNow        = pcrand(started_gt_zero, crpprfl_eq_zero)
                    CropHarvNow        = started_gt_zero & crpprfl_eq_zero & self.ricemask_BOOL
                    np_CropHarvNow     = pcr_as_numpy(CropHarvNow)
            
                    print "Start date read from Crop Profile..."
                    # Two auxilliary variables:
                    np_CRPST_gt_0       = np.greater(np_CRPST[:], 0)
                    np_CRPST_eq_STARTED = np.equal(np_CRPST[:], np_STARTED[:])  # of course started has to become positive then.
                    np_CropStartNow     = np.logical_and(np_CRPST_gt_0[:], np_CRPST_eq_STARTED[:]) * np_ricemask[:]  ##!
                    CropStartNow        = numpy2pcr(Boolean, np_CropStartNow, -99)
                    Started             = self.STARTED > 0
                    CropStarted         = Started & self.ricemask_BOOL
                    np_CropStarted      = pcr_as_numpy(CropStarted) * np_ricemask[:]  ##!
                    self.STARTED        = (self.STARTED + self.CRPST) * ifthenelse(CropHarvNow, scalar(0.), 1.)  # - ifthenelse(CropHarvNow, self.STARTED, 0.)
                elif self.AutoStartStop == True:
                    np_CropHarvNow     = (1 - np_Not_Finished)* np_ricemask[:]
                    HarvNow            = Not_Finished == False
                    CropHarvNow        = HarvNow & self.ricemask_BOOL
                    #print "Transpl. date based on cumulative rain after November 1..."
                    # Two auxilliary variables:
                    Time2Plant1stCrop    = self.PSUM >= self.RainSumReq
                    StdMin1              = self.STARTED == -1
                    CropStartNow_Season1 = Time2Plant1stCrop & self.ricemask_BOOL
                    CropStartNow_Season2 = StdMin1 & self.ricemask_BOOL
                    CropStartNow         = CropStartNow_Season1 | CropStartNow_Season2
                    np_CropStartNow      = pcr_as_numpy(CropStartNow) 
                    CropStartNow_scalar  = scalar(CropStartNow)
                    if self.Sim3rdSeason == False:   
                        HarvSeason1_temp     = FirstSeason & CropHarvNow
                    
                        HarvSeasonOne        = HarvSeason1_temp & self.ricemask_BOOL
                        HarvSeason2_temp     = SecondSeason & CropHarvNow
                        HarvSeasonTwo        = HarvSeason2_temp & self.ricemask_BOOL
                        self.Season          = self.Season + ifthenelse(HarvSeasonOne, self.ricemask, 0.) - ifthenelse(HarvSeasonTwo, self.ricemask * 2., 0.) # beware, this variable is also modified in another equation
                        Started              = self.STARTED > 0
                        CropStarted          = Started & self.ricemask_BOOL
                        SeasonOneHarvd       = self.STARTED < 0
                        SeasonOneHarvd_Scalar= scalar(SeasonOneHarvd)
                        np_CropStarted       = pcr_as_numpy(CropStarted)
                        pcr_PrepareField_temp = scalar(self.Pausedays)
                        pcr_PrepareField      = ifthenelse(FirstSeason, pcr_PrepareField_temp, 0.)
                        self.STARTED         = (self.STARTED + CropStartNow_scalar + scalar(CropStarted)) * ifthenelse(CropHarvNow, scalar(0.), 1.) - ifthenelse(HarvSeasonOne, pcr_PrepareField, 0.) + SeasonOneHarvd_Scalar
                    elif self.Sim3rdSeason == True:   
                        HarvSeason12_temp    = FirstSeason | SecondSeason
                        HarvSeasonOneTwo     = HarvSeason12_temp & CropHarvNow
                        #HarvSeason3_temp     = SecondSeason & CropHarvNow
                        HarvSeasonThree      = ThirdSeason & CropHarvNow
                        self.Season          = self.Season + ifthenelse(HarvSeasonOneTwo, scalar(1.), 0.) - ifthenelse(HarvSeasonThree, scalar(3.), 0.) # beware, this variable is also modified in another equation
                        Started              = self.STARTED > 0
                        CropStarted          = Started & self.ricemask_BOOL
                        Season12Harvd       = self.STARTED < 0
                        Season12Harvd_Scalar= scalar(Season12Harvd)
                        np_CropStarted       = pcr_as_numpy(CropStarted)
                        pcr_PrepareField_temp = scalar(self.Pausedays)
                        FirstorSecondSeason     = FirstSeason | SecondSeason
                        pcr_PrepareField      = ifthenelse(FirstorSecondSeason, pcr_PrepareField_temp, 0.)
                        self.STARTED         = (self.STARTED + CropStartNow_scalar + scalar(CropStarted)) * ifthenelse(CropHarvNow, scalar(0.), 1.) - ifthenelse(HarvSeasonOneTwo, pcr_PrepareField, 0.) + Season12Harvd_Scalar
                    else: 
                        print self.Sim3rdSeason
                        time.sleep(10)
                    # self.started= 0 and season = 2
                    # change season directly after harvest season 1
                else:
                    np_CropStartNow      = np_Zero[:]
                    CropStartNow         = numpy2pcr(Boolean, np_CropStartNow, -99)
                    np_CropStarted       = np_Zero[:]
                    CropStarted          = numpy2pcr(Boolean, np_CropStarted, -99)
                    print "Crop growth and/or harvest not initializing, pls. check wflow_lintul.ini..."            
                    time.sleep(100)
                           
        
        if self.WATERLIMITED == "True":
            #pcr_TRANRF       = self.Transpiration/self.PotTrans # Via numpy, because the NOTNUL is essential here. 
            np_PotTrans       = pcr_as_numpy(self.PotTrans)
            np_Transpiration  = pcr_as_numpy(self.Transpiration)
            np_TRANRF         = np_Transpiration[:] / NOTNUL(np_PotTrans[:]) # Via numpy, because the NOTNUL is essential here (todo: implement NOTNUL for pcr). 
            pcr_TRANRF        = numpy2pcr(Scalar, np_TRANRF, -99)
            WAWP              = WCWP * self.ROOTD_mm 
            Enough_water      = ifthenelse(CropStartNow, True, self.WA > WAWP) # timestep delay...! todo 
            np_Enough_water   = pcr_as_numpy(Enough_water)
        else:
            print "Warning, run without water effects on crop growth..."
            np_TRANRF         = np_One[:]
            pcr_TRANRF        = numpy2pcr(Scalar, np_TRANRF, -99)
            np_Enough_water   = np_One[:] #Todo: check!
            Enough_water      = numpy2pcr(Boolean, np_Enough_water, -99)
            
        #self.T = (self.TMIN + self.TMAX)/2. # for testing with Wageningen weather files only - sdv
        np_T                     = pcr_as_numpy(self.T)
        np_DAVTMP                = np_T  
        np_NNI                   = NNI * np_One[:]

      # Define tests for conditions that influence the behaviour of the crop:
        DAVTMP         = numpy2pcr(Scalar, np_DAVTMP, -99)
        Warm_Enough    = DAVTMP >= self.TBASE
        np_Warm_Enough = np.greater_equal(np_DAVTMP[:], self.TBASE * np_One[:])  # !
        DegreeDay      = DAVTMP - self.TBASE
        DTEFF          = ifthenelse(Warm_Enough, DegreeDay, 0.)
        np_DTEFF       = pcr_as_numpy(DTEFF)

        #np_Enough_water = np.greater(np_WA[:], (WCWP * np_ROOTD_mm[:]))
        #np_Enough_water   = np_One[:] #Todo: check!
        #Enough_water      = numpy2pcr(Boolean, np_Enough_water, -99)
        Leaves_Present    = self.LAI > 0.
        np_Leaves_Present = pcr_as_numpy(Leaves_Present)

      # Check when certain important decision moments are reached, in chronological order:
        BeforeAnthesis        = self.TSUM < self.TSUMAN
        np_BeforeAnthesis     = pcr_as_numpy(BeforeAnthesis)
        UntilAnthesis         = self.TSUM <= self.TSUMAN
        AtAndAfterAnthesis    = self.TSUM >= self.TSUMAN
        np_AtAndAfterAnthesis = pcr_as_numpy(AtAndAfterAnthesis)
        AfterAnthesis         = self.TSUM > self.TSUMAN
        Roots_Dying           = self.DVS >= self.DVSDR
        np_Roots_Dying        = pcr_as_numpy(Roots_Dying)

        Vegetative   = CropStarted & UntilAnthesis
        Generative   = CropStarted & AfterAnthesis
        EarlyStages  = self.DVS < 0.2
        LaterStages  = self.DVS >= 0.2
        SmallLeaves  = self.LAI < 0.75
        BiggerLeaves = self.LAI >= 0.75
        Juvenile     = EarlyStages & SmallLeaves
        Adult        = LaterStages | BiggerLeaves
        np_Juvenile  = pcr_as_numpy(Juvenile)
        np_Adult     = pcr_as_numpy(Adult)

      # Calculate daylength, based on latitude and day of year (assumed similar throughout the catchment area -> scalar, no array) - SdV
        DAYL = astro2(DOY, self.LAT)
        
      # Interpolation functions
        Interpol_SLACF = Afgen2(self.SLACF)
        Interpol_FRTTB  = Afgen2(self.FRTTB)
        Interpol_FLVTB  = Afgen2(self.FLVTB)
        Interpol_FSTTB  = Afgen2(self.FSTTB)
        Interpol_FSOTB  = Afgen2(self.FSOTB)
        Interpol_RDRTB  = Afgen2(self.RDRTB)
        Interpol_PHOTTB = Afgen2(self.PHOTTB)
        
        SLA       = self.SLAC * Interpol_SLACF(np_DVS[:])  # * np.exp(-NSLA * (1.-NNI))     
        #SLA2     = self.SLAC * self.SLACF2.lookup_linear(self.DVS) 
        FRTWET    = Interpol_FRTTB(np_DVS[:])
        #FRTWET2  = self.FRTTB2.lookup_linear(self.DVS)
        FLVT      = Interpol_FLVTB(np_DVS[:])
        #FLVT2    = self.FLVTB2.lookup_linear(self.DVS)
        FSTT      = Interpol_FSTTB(np_DVS[:])
        FSOT      = Interpol_FSOTB(np_DVS[:])
        RDRTMP    = Interpol_RDRTB(np_DVS[:])
        PHOTT     = Interpol_PHOTTB(DAYL)

        EMERG     = CropStarted & Enough_water & Leaves_Present & Not_Finished
        np_EMERG  = pcr_as_numpy(EMERG)
        PHOTPF    = ifthenelse(BeforeAnthesis, PHOTT, scalar(1.))
        RTSUMP    = DTEFF * PHOTPF
        self.TSUM = (self.TSUM + ifthenelse(CropStartNow, scalar(self.TSUMI), 0.) + ifthenelse(EMERG, RTSUMP, 0.)) * ifthenelse(CropHarvNow, scalar(0.), 1.)  

        TSUM_veg  = self.TSUM / self.TSUMAN * ifthenelse(CropHarvNow, scalar(0.), 1.)
        TSUM_gen  = (1. + (self.TSUM - self.TSUMAN) / self.TSUMMT) * ifthenelse(CropHarvNow, scalar(0.), 1.)
        self.DVS  = ifthenelse(Vegetative, TSUM_veg, 0.) + ifthenelse(Generative, TSUM_gen, 0.)

      # Root depth growth:
        CanGrowDownward      = self.ROOTD_mm <= self.ROOTDM_mm
        np_CanGrowDownward   = pcr_as_numpy(CanGrowDownward)
        #np_CanGrowDownward  = np.less_equal(np_ROOTD_mm[:], self.ROOTDM_mm * np_One[:])
        np_RROOTD_mm         = np_Enough_water[:] * np_BeforeAnthesis[:] * np_EMERG[:] * np_CanGrowDownward[:] * self.RRDMAX_mm
        self.ROOTD_mm        = (self.ROOTD_mm + ifthenelse(CropStartNow, self.ROOTDI_mm, scalar(0.)) + numpy2pcr(Scalar, np_RROOTD_mm[:], -99)) * ifthenelse(CropHarvNow, scalar(0.), 1.)
      # np_EXPLOR            = np_RROOTD_mm[:] * WCFC # for ROOTD in mm (!)
        np_EXPLOR            = np_RROOTD_mm[:] * WCST # for ROOTD in mm (!)

        
        #############################################################################################################
      # Water Limitation: effects on partitioning
        #FRTMOD              = np_One[:]  
        #np_FRTMOD               = np.maximum(np_One[:], np_One[:] / (np_TRANRF[:] + 0.5 * np_One[:]))
        FRTMOD               = max(1., 1./(pcr_TRANRF + 0.5)) # was FRTMOD2
        np_FRTMOD            = pcr_as_numpy(FRTMOD)
        FRT                  = FRTWET[:] * np_FRTMOD[:]
        #FRT2                = FRTWET2 * FRTMOD2
        FSHMOD               = (1. - FRT[:]) / (1. - (FRT[:] / np_FRTMOD[:]))
        #FSHMOD2             = (1. -FRT2)/(1 - FRT2/FRTMOD2)
        FLV                  = FLVT[:] * FSHMOD
        #FLV2                = FLVT2 * FSHMOD2
        FST                  = FSTT[:] * FSHMOD
        FSO                  = FSOT[:] * FSHMOD
        PartitionCheck       = FLV[:] + FST[:] + FSO[:] + FRT[:]

      # todo: incorporate effects of N stress - sdv 30-11-15
        #PARINT               = 0.5 * np_IRRAD[:] * 0.001 * (1. - np.exp(-self.K * np_LAI[:])) * np_Not_Finished[:]
        PARINT2              = ifthenelse(Not_Finished, 0.5 * self.IRRAD * 0.001 * (1. - exp(-self.K * self.LAI)), 0.)
        np_PARINT            = pcr_as_numpy(PARINT2)
        GTOTAL               = self.LUE * np_PARINT[:] * np_TRANRF[:]
        GTOTAL2              = self.LUE * PARINT2 * pcr_TRANRF

      # Rel. Death rate due to ageing:
        np_RDRDV             = np_AtAndAfterAnthesis[:] * RDRTMP[:]
        RDRDV                = numpy2pcr(Scalar, np_RDRDV, -99)
        RDRSH                = max(0., self.RDRSHM * (self.LAI - self.LAICR)/self.LAICR)
        np_RDRSH             = pcr_as_numpy(RDRSH)
        #RDRSH               = np.maximum(np_Zero[:], self.RDRSHM * (np_LAI[:] - self.LAICR[:]) / self.LAICR[:])
        RDR                  = max(RDRDV, RDRSH)
        np_RDR               = pcr_as_numpy(RDR)
        #RDR                 = np.maximum(RDRDV[:], np_RDRSH[:]) * np_Not_Finished[:]

      # Impact of leaf dying on leaf weight - N limitation stuff not (yet) implemented
        N_Limitation         = np.less(np_NNI[:], 1.)
        DLVNS                = np_CropStarted[:] * N_Limitation[:] * np_WLVG[:] * self.RDRNS * (1. - np_NNI[:])
        DLVS                 = np_WLVG[:] * np_RDR[:]
        DLV                  = (DLVS[:] + DLVNS[:]) * np_Not_Finished[:]

        RWLVG                = np_EMERG[:] * (GTOTAL[:] * FLV[:] - DLV[:])
        np_WLVG[:]           = (np_WLVG[:] + np_CropStartNow[:] * self.WLVGI + RWLVG[:]) * (1. - np_CropHarvNow[:])
        np_WLVD[:]          += DLV[:]

      # Leaf totalGrowthRate and LAI.
        GLV                  = FLV[:] * GTOTAL[:]
        GLV_pcr              = numpy2pcr(Scalar, GLV, -99)
        #GLV2                = FLV2 * GTOTAL2
        Adt_or_Harv = pcror(Adult, CropHarvNow)
        Juv_or_Harv = pcror(Juvenile, CropHarvNow)
        NoLeavesYet = self.LAI == 0.
        LetsGo      = pcrand(Enough_water, CropStartNow)
        LetsGro     = pcrand(NoLeavesYet, LetsGo)
        GLAI        = (np_Adult[:] * SLA[:] * GLV[:] * (1. - np_CropHarvNow[:]) + 
                       np_Juvenile[:] * (np_LAI[:] * (np.exp(self.RGRL * np_DTEFF[:] * DELT) - 1.) / DELT) * np_TRANRF[:] * np.exp(-NLAI * (1. - NNI)) * (1. - np_CropHarvNow[:]) + 
                       np.equal(np_LAI[:], np_Zero[:]) * np_Enough_water[:] * np_CropStartNow[:] * self.LAII / DELT)
        #GLAI2      = ifthenelse(Adt_or_Harv, SLA2 * GLV2, scalar(0.))  + \
        #             ifthenelse(Juv_or_Harv, self.LAI * (exp(self.RGRL * DTEFF * DELT )- 1.)/DELT * pcr_TRANRF * exp(-self.LAI * (1.0-NNI)), 0.)  + \
        #             ifthenelse(LetsGro, self.LAII/DELT, scalar(0.))
        
        
      # For future development (sdv): 
        #pcr_GLAI    = numpy2pcr(Scalar, GLAI, -99)
        #np_GLA2    = pcr_as_numpy(GLAI2)
        #GLAI2name  = "GLAI2000."+str(self.currentTimeStep()) 
        #report(GLAI2, GLAI2name)
        #np_GLAI_check  = np.equal(GLAI, np_GLA2)
        #pcr_GLAI_check = pcr_GLAI == GLAI2
        #GLAI2checkname = "GLAICHK00."+str(self.currentTimeStep()) 
        #report(pcr_GLAI_check, GLAI2checkname)
               
      # (Abs.) impact of leaf dying on LAI
      # Death of leaves due to ageing and shading:
        DLAIS      = np_LAI[:] * np_RDR[:]
        DLAINS     = np_CropStarted[:] * N_Limitation[:] * DLVNS[:] * SLA[:]
        DLAI       = (DLAIS + DLAINS) * np_Not_Finished[:]
      # The initial LAI (LAII, transplanted rice) is added to GLAI at crop establishment, not in below state equation as done by Shibu et al. (2010).         
        np_LAI[:]  = (np_LAI[:] + GLAI[:] - DLAI[:]) * (1. - np_CropHarvNow[:]) 
        
        # DRRT = np.greater_equal(np_DVS[:], np_DVSDR[:]) * np_WRT[:] * RDRRT
        DRRT       = np_Roots_Dying[:] * np_WRT[:] * self.RDRRT
        #RWRT       = GTOTAL[:] * FRT[:] - DRRT[:] * np_Not_Finished[:]
        RWRT       = np_EMERG[:] * (GTOTAL[:] * FRT[:] - DRRT[:])
        np_WRT[:]  = (np_WRT[:] + np_CropStartNow[:] * self.WRTLI + RWRT[:]) * (1. - np_CropHarvNow[:])
        np_WDRT[:]+= DRRT[:]

        RWSO       = np_EMERG[:] * (GTOTAL[:] * FSO[:])
        np_WSO[:]  = (np_WSO[:] + np_CropStartNow[:] * self.WSOI + RWSO[:]) * (1. - np_CropHarvNow[:])
        np_WSOTHA  = np_WSO[:] / 100.

        RWST       = np_EMERG[:] * (GTOTAL[:] * FST[:])
        np_WST[:]  = (np_WST[:] + np_CropStartNow[:] * self.WSTI + RWST[:]) * (1. - np_CropHarvNow[:])

        np_WLV     = np_WLVG[:] + np_WLVD[:]
        TAGBM      = np_WLV[:] + np_WST[:] + np_WSO[:]
        
        
        #FRTTB2file = 'frttb2'
        self.Test = EMERG
        
        print '\n', cellvalue(ThirdSeason,100,100)[0], '\n'
        #print self.DVSI, "self.dvsi"
        #time.sleep(0.25)
       #print '\n', cellvalue(Enough_water,100,100)[0], cellvalue(self.WA,100,100)[0], cellvalue(self.TSUM,100,100)[0],cellvalue(CropStartNow,100,100)[0], self.currentTimeStep(), '\n'
        
       #For quickly getting point output (sdv). Works only with a wf_supplyEndTime() implemented in wf_dynamicframework... todo?
        Point_Output_Line = (str(cellvalue (self.LAI, 100,100)[0]) + "," + str(cellvalue (self.TSUM, 100,100)[0]) + "," + str(cellvalue (self.Test, 100,100)[0]) + "," +
                            str(cellvalue (CropStarted, 100,100)[0]) + "," + str(cellvalue (Enough_water, 100,100)[0])  + "," + str(cellvalue (Leaves_Present, 100,100)[0]) + ","
                             + str(cellvalue (HarvSeasonOne, 100,100)[0]) +"," + str(cellvalue (self.STARTED, 100,100)[0]) + "," +   str(cellvalue (self.Season, 100,100)[0]) + "," +
                            str(np_CropStartNow[100,100])+ "," + str(np_CropHarvNow[100,100]) + '\n')
        if self.date < self.enddate:
            Point_Output.write(Point_Output_Line)
        elif self.date == self.enddate:
            Point_Output.close()

#EMERG     = CropStarted & Enough_water & Leaves_Present & Not_Finished
# The main function is used to run the program from the command line

def main(argv=None):
    """
    *Optional but needed it you want to run the model from the command line*

    Perform command line execution of the model. This example uses the getopt
    module to parse the command line options.
x
    The user can set the caseName, the runDir, the timestep and the configfile.
    """
    global multpars
    caseName = "default_lintul"
    runId = "run_default"
    configfile = "wflow_lintul.ini"
    _lastTimeStep = 0
    _firstTimeStep = 0
    runinfoFile = "runinfo.xml"
    timestepsecs = 86400
    wflow_cloneMap = 'wflow_subcatch.map'
    _NoOverWrite = 1
    loglevel = logging.DEBUG

    # This allows us to use the model both on the command line and to call
    # the model usinge main function from another python script.

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return
    ########################################################################
    ## Process command-line options                                        #
    ########################################################################
    try:
        opts, args = getopt.getopt(argv, 'F:C:S:T:c:s:R:l')
    except getopt.error, msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == '-C': caseName = a
        if o == '-R': runId = a
        if o == '-c': configfile = a
        if o == '-s': timestepsecs = int(a)
        if o == '-T': _lastTimeStep = int(a)
        if o == '-S': _firstTimeStep = int(a)
        if o == '-f': _NoOverWrite = 0
        if o == '-l': exec "loglevel = logging." + a

    if (len(opts) <= 1):
        usage()

    # starttime = dt.datetime(1990,01,01)
    starttime = dt.datetime(1981, 9, 27)

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(myModel, _lastTimeStep, firstTimestep=_firstTimeStep, datetimestart=starttime)
    dynModelFw.createRunId(NoOverWrite=False, level=loglevel)
    #dynModelFw.createRunId(NoOverWrite=_NoOverWrite, level=loglevel, logfname=LogFileName,model="wflow_lintul",doSetupFramework=False)

    dynModelFw.setupFramework()
    dynModelFw._runInitial()
    dynModelFw._runResume()
    # dynModelFw._runDynamic(0,0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()

