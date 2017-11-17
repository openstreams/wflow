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
import time



# from datetime import date
# from pcse_util import *
# from math import exp, expm1

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


# runinfoFile           = "runinfo.xml"
# runinfoFile_Lintul    = os.path.join(os.path.abspath("wflow_lintul"), "inmaps", runinfoFile)
# starttime             = getStartTimefromRuninfo(runinfoFile_Lintul)
# DOY_wflow             = starttime.strftime('%j')

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
    PI = 3.1415926

    # SINE AND COSINE OF LATITUDE
    SINLAT = math.sin(PI * LAT / 180.)
    COSLAT = math.cos(PI * LAT / 180.)

    # MAXIMAL SINE OF DECLINATION
    SINDCM = math.sin(PI * 23.45 / 180.)

    # SINE AND COSINE OF DECLINATION  (EQUATIONS 3.4, 3.5)

    # SINDEC = -SINDCM * cos(2.* PI * (DAY+10.)/365.)
    # The '9' below (instead of 10) keeps things in sync with FST...
    # Todo: try to understand this at a certain point... for now it works perfectly.
    SINDEC = -SINDCM * math.cos(2. * PI * (DAY + 11.) / 365.)
    COSDEC = math.sqrt(1. - SINDEC * SINDEC)

    # THE TERMS A AND B ACCORDING TO EQUATION 3.3

    A = SINLAT * SINDEC
    B = COSLAT * COSDEC

    # DAYLENGTH ACCORDING TO EQUATION 3.6. Make sure not to use asin from PCRaster...
    DAYL = 12. * (1. + (2. / PI) * math.asin(A / B))

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


# -------------------------------------------------------------------------------



AutoStartStop = True
Pause = 13

###############################################################################
# Some general information for setting up the runs:
LAT = 3.16  # the latitude of Bah Lias, North Sumatra, for now.
DELT = 1  # Time step = 1 day.
TSUMI = 362.  # TSUM at transplanting, kind of crop specific, actually.
TINY = 1e-6

# This needs to be set according to the geographic extent (map dimensions) of your study area/catchment:
np_Zero = numpy.zeros((219, 286))
np_One = numpy.ones((219, 286))
Pausedays = np_One[:] * (Pause + 1)

# Crop specific coefficients for rice:
K = 0.6  # light extinction coefficient
#LUE = 3.  # Light use efficiency.
SLAC = 0.02  # Specific leaf area constant.
TSUMAN = 1420.
TSUMMT = 580.
TBASE = 8.  # Base temperature for spring wheat crop.
RGRL = 0.009  # Relative growth rate of LAI at the exponential growth phase
WLVGI = 0.86
WSTI = 0.71
WRTLI = 1.58
WSOI = 0.00
RDRNS = 0.03
DVSDR = 0.8  # Development stage above which death of leaves and roots start.
RDRRT = 0.03  # Relative death rate of roots.
TTSUM = 2000.
RDRSHM = 0.03  # and the maximum relative death rate of leaves due to shading.
LAICR = 4.0  # (oC d)-1, critical LAI above which mutual shading of leaves occurs.
ROOTDM_mm = 1000.  # Maximum rooting depth (mm!) for a rice crop.
RRDMAX_mm = 10.  # maximum rate of increase in rooting depth (mm! d-1) for a rice crop.
ROOTDI_mm = 50.  # initial rooting depth (mm!)
TRANCO = 3.  # Transpiration constant (mm/day) indicating the level of drought tolerance of the rice crop.

# Fixing some values here for now: to be simulated or read from maps later on.
# TRANRF      = np_One[:] # If enabled: switches off the included water balance and removes water limitation effects!
# TRANRF equation further below needs to be switched off if this is enabled.
NNI = 1.
WCI = 0.60
WCAD = 0.01
WCWP = 0.21
WCFC = 0.47
WCWET = 0.50
# WA = WCST   = 0.55
WCST = 0.55  # aangepast voor BMI sdv 01-06-2015
WCWP = 0.20
DRATE = 30.
IRRIGF = 0.
IRRIG = np_Zero[:]

# Nitrogen-related coefficients, specific for rice.
NPART = 1.0  # Coefficient for the effect of N stress on leaf biomass reduction
NSLA = 1.0  # Coefficient for the effect of N stress on SLA reduction
NLAI = 1.0  # Coefficient for the effect of N stress on LAI reduction(during juvenile phase)

# ----------------------- Interpolation functions---------------------#

# Relative death rate of leaves as a function of Developmental stage
# Note: from FST Lintul 3-rice, 29-10-15 (adopted from ORYZA2000), SdV
RDRTB = [0., 0.00,
         0.6, 0.00,
         1.0, .015,
         1.6, 0.025,
         2.1, 0.05]

# Note: from FST Lintul 3-rice, 30-10-15, SdV
PHOTTB = [0.0, 0.0,
          8., 0.0,
          10., 1.0,
          12., 1.0,
          13., 0.8,
          14., 0.6,
          18., 0.0]

# Leaf area correction function as a function of development stage, DVS.
# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
SLACF = [0.0, 1.72,
         0.21, 1.72,
         0.24, 1.72,
         0.33, 1.32,
         0.7, 1.20,
         1.01, 1.00,
         2.0, 0.75,
         2.1, 0.75]

# Maximum N concentration in the leaves, from which the N-conc.values of the
# stem and roots are derived, as a function of development stage.
# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
NMXLV = [0.0, 0.05,
         0.4, 0.05,
         0.7, 0.04,
         1.0, 0.03,
         2.0, 0.02,
         2.1, 0.02]

# ********** Partitioning coefficients ***********************************
# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
FRTTB = [0.0, 0.30,
         0.48, 0.30,
         0.62, 0.12,
         0.69, 0.11,
         0.84, 0.11,
         0.92, 0.10,
         1.00, 0.08,
         1.38, 0.00,
         2.10, 0.0]

# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
FLVTB = [0.0, 0.315,
         0.48, 0.35,
         0.62, 0.44,
         0.69, 0.463,
         0.84, 0.463,
         0.92, 0.45,
         1.00, 0.00,
         1.38, 0.00,
         2.10, 0.0]

# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
FSTTB = [0.0, 0.385,
         0.48, 0.35,
         0.62, 0.44,
         0.69, 0.427,
         0.84, 0.427,
         0.92, 0.27,
         1.00, 0.00,
         1.38, 0.00,
         2.1, 0.0]

# Note: taken from FST Lintul 3-rice, 29-10-15, SdV
FSOTB = [0.0, 0.00,
         0.48, 0.00,
         0.62, 0.00,
         0.69, 0.00,
         0.84, 0.00,
         0.92, 0.18,
         1.00, 0.92,
         1.38, 1.00,
         2.1, 1.00]

# Calculate initial LAI
np_TSUMI = TSUMI * np_One[:]
DVSI = np_TSUMI[:] / TSUMAN
Interpol_SLACF = Afgen2(SLACF)
SLACFI = Interpol_SLACF(DVSI)
ISLA = SLAC * SLACFI
LAII = WLVGI * ISLA


# Initial soil moisture
# np_WAI                = 1000. * ROOTDI * WCI * np_One[:]

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

        self.LUE = float(configget(self.config, "model", "LUE", "2.47"))
        self.AutoStartStop = eval(configget(self.config, "model", "AutoStartStop", "True"))
        self.BMI_RUN = configget(self.config, "model", "BMI_RUN", "True")
        self.WATERLIMITED = (configget(self.config, "model", "WATERLIMITED", "True"))
        self.CropStartDOY = int(configget(self.config, "model", "CropStartDOY", "0"))
        self.HarvestDAP = int(configget(self.config, "model", "HarvestDAP", "150"))
        #self.stdt = (configget(self.config, "run", "starttime", "1979-01-02 00:00:00")).rsplit('-')
        #self.startyr_lintul, self.startmo_lintul, self.startd_lintul = int(self.stdt[0]), int(self.stdt[1]), int(
        #    self.stdt[2].rsplit()[0])
        #self.startdate_lintul = datetime(self.startyr_lintul, self.startmo_lintul, self.startd_lintul)

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
        # modelparameters.append(self.ParamType(name="Temperature",stack="inmaps/TEMP",type="timeseries",default=10.0,verbose=False,lookupmaps=[])),
        modelparameters.append(self.ParamType(name="IRRAD", stack="inmaps/IRRAD", type="timeseries", default=11.0, verbose=False,lookupmaps=[])),
        # modelparameters.append(self.ParamType(name="TMIN",stack="inmaps/TMIN",type="timeseries",default=10.0,verbose=False,lookupmaps=[])),
        # modelparameters.append(self.ParamType(name="TMAX",stack="inmaps/TMAX",type="timeseries",default=10.0,verbose=False,lookupmaps=[])),
        modelparameters.append(self.ParamType(name="T", stack="inmaps/T", type="timeseries", default=10.0, verbose=False, lookupmaps=[])),
        # modelparameters.append(self.ParamType(name="VAP",stack="inmaps/VAP",type="timeseries",default=10.0,verbose=False,lookupmaps=[])),
        # modelparameters.append(self.ParamType(name="WIND",stack="inmaps/WIND",type="timeseries",default=2.0,verbose=False,lookupmaps=[])),
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

        """

        # states = ['Season', 'Test', 'LAI', 'WLVG', 'WLVD', 'WST', 'WSO', 'WRT', 'ROOTD', 'ROOTD_mm', 'WDRT', 'TSUM', 'STARTED', 'DVS']
        states = ['Season', 'PSUM', 'Test', 'LAI', 'WLVG', 'WLVD', 'WST', 'WSO', 'WRT', 'ROOTD_mm', 'WDRT', 'TSUM',
                  'STARTED', 'DVS']
        if self.BMI_RUN == "False":
            states.append('WA')
            # states.pop(-6)

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

        self.timestepsecs = int(configget(self.config, 'model', 'timestepsecs', '86400'))
        self.basetimestep = 86400

        # Reads all parameter from disk
        self.wf_updateparameters()
        self.logger.info("Starting LINTUL Dynamic Crop Growth Simulation...")
        
        wflow_ricemask = configget(self.config, "model", "wflow_ricemask", "staticmaps/wflow_ricemask.map")
        self.ricemask = self.wf_readmap(os.path.join(self.Dir,wflow_ricemask),0.0,fail=True)
        self.ricemask_BOOL = boolean(self.ricemask)

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
        

        #self.date = self.startdate_lintul + dt.timedelta(self.currentTimeStep() - 1)
        self.date = datetime.utcfromtimestamp(self.wf_supplyStartTime()) + dt.timedelta(self.currentTimeStep() - 1)      
        #DOY = int(self.date.strftime('%j'))
        DOY = self.wf_supplyJulianDOY() 

        One = numpy2pcr(Scalar, np_One[:], -99)
        Zero = numpy2pcr(Scalar, np_Zero, -99)
        TSUM_not_Finished = self.TSUM <= TTSUM
        DVS_not_Finished = self.DVS <= 2.01
        Not_Finished = TSUM_not_Finished & DVS_not_Finished
        np_Not_Finished = pcr_as_numpy(Not_Finished)
        np_STARTED = pcr_as_numpy(self.STARTED)
        #self.ricemask = readmap("D:\\Reference_run\\wflow\\wflow_lintul\\inmaps\\CRPST000.031")
        #self.ricemask_BOOL = boolean(self.ricemask)
        np_ricemask = pcr_as_numpy(self.ricemask)

        if (self.date.month == 11 and self.date.day == 1):  # and np.sum(np_RAIN) == 0.:
            november1 = True
            self.PSUM += self.RAIN + 0.001
        else:
            november1 = False
            # self.Season = self.Season + ifthenelse(self.ricemask_BOOL, scalar(1.), 0.)
            # self.Season = self.Season + self.ricemask
        self.Season += ifthenelse(november1, self.ricemask, 0.)
        self.Test += ifthenelse(november1, self.ricemask, 0.)
        # toto: check above 1000 new aguila batch needed
        # check with e.g. psum?

        # self.Season = self.Season + ifthenelse(HarvSeasonOne, self.ricemask, 0.)  - ifthenelse

        FirstSeason = self.Season == 1
        SecondSeason = self.Season == 2
        TestPos = self.PSUM > 0.
        Test200 = self.PSUM <= 200.
        TestRange = TestPos & Test200
        # self.PSUM = ifthenelse(TestRange, self.PSUM + self.RAIN, 0.)
        self.PSUM = (self.PSUM + ifthenelse(TestRange, self.RAIN, 0.)) * ifthenelse(TestRange, scalar(1.), 0.)

        pcr_PrepareField_temp = numpy2pcr(Scalar, Pausedays, -99)
        pcr_PrepareField = ifthenelse(FirstSeason, pcr_PrepareField_temp, 0.)

        # Create numpy array from states:
        np_CRPST = pcr_as_numpy(self.CRPST)
        np_DVS = pcr_as_numpy(self.DVS)
        np_LAI = pcr_as_numpy(self.LAI)
        np_WA = pcr_as_numpy(self.WA)
        np_TSUM = pcr_as_numpy(self.TSUM)
        np_WRT = pcr_as_numpy(self.WRT)
        np_WSO = pcr_as_numpy(self.WSO)
        np_WST = pcr_as_numpy(self.WST)
        np_WDRT = pcr_as_numpy(self.WDRT)
        np_Test = pcr_as_numpy(self.Test)
        np_PSUM = pcr_as_numpy(self.PSUM)
        np_WLVG = pcr_as_numpy(self.WLVG)
        np_WLVD = pcr_as_numpy(self.WLVD)
        np_ROOTD_mm = pcr_as_numpy(self.ROOTD_mm)
        np_Season = pcr_as_numpy(self.Season)

        # Make numpy arrays from soil characteristics. Todo: later read them from PCRaster maps, then pcr_as_numpy.
        np_WCAD = WCAD * np_One[:]
        np_WCWP = WCWP * np_One[:]
        np_WCFC = WCFC * np_One[:]
        np_WCWET = WCWET * np_One[:]
        np_WCST = WCST * np_One[:]
                

        # Initializing crop harvest:
        if self.CropStartDOY > 0 and self.HarvestDAP > 0:
            np_CropHarvNow = np.greater_equal(DOY, self.CropStartDOY + self.HarvestDAP) * np_One[:]
            CropHarvNow = numpy2pcr(Boolean, np_CropHarvNow, -99)
            print "Warning: harvest date read from ini file, not from Crop Profile map..."
        elif self.CropStartDOY > 0 and self.HarvestDAP == 0:
            np_CropHarvNow = 1 - np_Not_Finished
            CropHarvNow = Not_Finished == False
            print "Harvest date not specified; crop harvest at crop maturity"
        elif self.CropStartDOY == 0 and self.HarvestDAP > 0:
            # CropHarvNow            = self.STARTED -2 == self.HarvestDAP
            CropHarvNow = self.STARTED == self.HarvestDAP
            np_CropHarvNow = pcr_as_numpy(CropHarvNow)
        elif self.CropStartDOY == 0 and self.HarvestDAP == 0:
            started_gt_zero = self.STARTED > 0.
            crpprfl_eq_zero = self.CRPST == 0.
            CropHarvNow = pcrand(started_gt_zero, crpprfl_eq_zero)
            np_CropHarvNow = pcr_as_numpy(CropHarvNow)
        else:
            np_CropHarvNow = np_Zero[:]
            CropHarvNow = numpy2pcr(Boolean, np_CropHarvNow, -99)
            print "Crop harvest not initialized, found a strange values in ini file..."
            # Initializing crop growth, optionally from a single start day (CropStartDOY in the ini file),
        # but normally from a crop profile forcing variable.
        if self.CropStartDOY > 0.:  # this part works fine, sdv, 13 sept 2017

            # np_CropStartNow        = np.equal(DOY, self.CropStartDOY) * np_One[:]
            np_CropStartNow = np.equal(DOY, self.CropStartDOY) * np_ricemask[:]
            CropStartNow = numpy2pcr(Boolean, np_CropStartNow, -99)
            # np_CropStarted         = np.greater_equal(DOY, self.CropStartDOY) * np_One[:]
            np_CropStarted = np.greater_equal(DOY, self.CropStartDOY) * np_ricemask[:]
            CropStarted = numpy2pcr(Boolean, np_CropStarted, -99)
            print "Warning: using start date from ini file, not read from Crop Profile..."

        elif self.CropStartDOY == 0 and AutoStartStop == False:
            print "Start date read from Crop Profile..."
            # Two auxilliary variables:
            np_CRPST_gt_0 = np.greater(np_CRPST[:], 0)
            np_CRPST_eq_STARTED = np.equal(np_CRPST[:], np_STARTED[:])  # of course started has to become positive then.
            np_CropStartNow = np.logical_and(np_CRPST_gt_0[:], np_CRPST_eq_STARTED[:]) * np_ricemask[:]  ##!
            CropStartNow = numpy2pcr(Boolean, np_CropStartNow, -99)
            # CropStarted            = self.CRPST > 0
            CropStarted = self.STARTED > 0
            # np_CropStarted         = pcr_as_numpy(self.CRPST)
            np_CropStarted = pcr_as_numpy(CropStarted) * np_ricemask[:]  ##!
            self.STARTED = (self.STARTED + self.CRPST) * ifthenelse(CropHarvNow, Zero,
                                                                    1.)  # - ifthenelse(CropHarvNow, self.STARTED, 0.)

        elif self.CropStartDOY == 0 and AutoStartStop == True:
            #print "Transpl. date based on cumulative rain after November 1..."
            # Two auxilliary variables:
            # CropStartNow          = self.PSUM >= 200.
            PSUM200 = self.PSUM >= 200.
            StdMin1 = self.STARTED == -1
            CropStartNow_Season1 = pcrand(PSUM200, self.ricemask_BOOL)
            CropStartNow_Season2 = pcrand(StdMin1, self.ricemask_BOOL)
            CropStartNow = pcror(CropStartNow_Season1, CropStartNow_Season2)
            np_CropStartNow = pcr_as_numpy(CropStartNow) * np_ricemask[:]
            CropStartNow_scalar = scalar(CropStartNow)

            HarvSeason1_temp = pcrand(FirstSeason, CropHarvNow)
            HarvSeasonOne = pcrand(HarvSeason1_temp, self.ricemask_BOOL)
            HarvSeason2_temp = pcrand(SecondSeason, CropHarvNow)
            HarvSeasonTwo = pcrand(HarvSeason2_temp, self.ricemask_BOOL)
            # self.Season           += ifthenelse(HarvSeasonOne, scalar(1.), 0.) #HarvSeason1 is al genoeg?
            self.Season = self.Season + ifthenelse(HarvSeasonOne, self.ricemask, 0.) - ifthenelse(HarvSeasonTwo,
                                                                                                  self.ricemask * 2.,
                                                                                                  0.)

            # CropStarted           = self.CRPST > 0
            CropStarted = self.STARTED > 0
            SeasonOneHarvd = self.STARTED < 0
            SeasonOneHarvd_Scalar = scalar(SeasonOneHarvd)
            np_CropStarted = pcr_as_numpy(CropStarted)
            # ReallyStarted         = self.STARTED > 0.
            self.STARTED = (self.STARTED + CropStartNow_scalar + scalar(CropStarted)) * ifthenelse(CropHarvNow, Zero,
                                                                                                   1.) - \
                           ifthenelse(HarvSeasonOne, pcr_PrepareField, 0.) + SeasonOneHarvd_Scalar
            # note to self:          apply elsewhere too: self.STARTED + CropStartNow_scalar + scalar(CropStarted)
            # self.started= 0 and season = 2
            # change season directly after harvest season 1
        else:
            np_CropStartNow = np_Zero[:]
            CropStartNow = numpy2pcr(Boolean, np_CropStartNow, -99)
            np_CropStarted = np_Zero[:]
            CropStarted = numpy2pcr(Boolean, np_CropStarted, -99)
            print "Crop growth not initializing, pls. check wflow_lintul.ini..."

            # self.STARTED         += self.CRPST * ifthenelse(CropHarvNow, Zero, 1.) - ifthenelse(CropHarvNow, self.STARTED, 0.)
            # self.STARTED         = (self.STARTED + self.CRPST) * ifthenelse(CropHarvNow, Zero, 1.) #- ifthenelse(CropHarvNow, self.STARTED, 0.)
            # self.STARTED         += self.CRPST


            # Implement forcing data, coefficients and some handy numbers spatially, as numpy arrays:

            # np_TMIN               = pcr_as_numpy(self.TMIN)
            # np_TMAX               = pcr_as_numpy(self.TMAX)
        np_T = pcr_as_numpy(self.T)
        np_DAVTMP = np_T  # 0.5     * (np_TMIN + np_TMAX)
        # np_VAP                = pcr_as_numpy(self.VAP)
        # np_WIND               = pcr_as_numpy(self.WIND)
        np_RAIN = pcr_as_numpy(self.RAIN)
        np_TSUMAN = TSUMAN * np_One[:]
        np_TTSUM = TTSUM * np_One[:]
        np_LAICR = LAICR * np_One[:]
        np_NNI = NNI * np_One[:]
        np_WCWP = WCWP * np_One[:]
        np_DVSDR = DVSDR * np_One[:]
        np_DELT = DELT * np_One[:]
        np_RGRL = RGRL * np_One[:]

        # Define tests for conditions that influence the behaviour of the crop:
        DAVTMP = numpy2pcr(Scalar, np_DAVTMP, -99)
        Warm_Enough = DAVTMP >= TBASE
        np_Warm_Enough = np.greater_equal(np_DAVTMP[:], TBASE * np_One[:])  # !
        DegreeDay = DAVTMP - TBASE
        DTEFF = ifthenelse(Warm_Enough, DegreeDay, 0.)
        np_DTEFF = pcr_as_numpy(DTEFF)

        # np_Enough_water       = np.greater(np_WA[:], (WCWP * np_ROOTD[:] * 1000.))
        np_Enough_water = np_One[:]
        Enough_water = numpy2pcr(Boolean, np_Enough_water, -99)
        Leaves_Present = self.LAI > 0.
        np_Leaves_Present = pcr_as_numpy(Leaves_Present)

        # Check when certain important decision moments are reached, in chronological order:
        BeforeAnthesis = self.TSUM < TSUMAN
        np_BeforeAnthesis = pcr_as_numpy(BeforeAnthesis)
        UntilAnthesis = self.TSUM <= TSUMAN
        AtAndAfterAnthesis = self.TSUM >= TSUMAN
        np_AtAndAfterAnthesis = pcr_as_numpy(AtAndAfterAnthesis)
        AfterAnthesis = self.TSUM > TSUMAN
        Roots_Dying = self.DVS >= DVSDR
        np_Roots_Dying = pcr_as_numpy(Roots_Dying)

        Vegetative = CropStarted & UntilAnthesis
        Generative = CropStarted & AfterAnthesis
        EarlyStages = self.DVS < 0.2
        LaterStages = self.DVS >= 0.2
        SmallLeaves = self.LAI < 0.75
        BiggerLeaves = self.LAI >= 0.75
        Juvenile = EarlyStages & SmallLeaves
        Adult = LaterStages | BiggerLeaves
        np_Juvenile = pcr_as_numpy(Juvenile)
        np_Adult = pcr_as_numpy(Adult)

        TSUMINIT = numpy2pcr(Scalar, np_TSUMI[:], -99)

        # Specific Leaf area(m2/g)
        Interpol_SLACF = Afgen2(SLACF)
        SLA = SLAC * np_One[:] * Interpol_SLACF(np_DVS[:])  # * exp(-NSLA * (1.-NNI))     should become numpy.exp

        # Calculate daylength, based on latitude and day of year:
        # Daylength assumed similar throughout the catchment area -> scalar, no array - SdV
        DAYL = astro2(DOY, LAT)

        Interpol_FRTTB = Afgen2(FRTTB)
        Interpol_FLVTB = Afgen2(FLVTB)
        Interpol_FSTTB = Afgen2(FSTTB)
        Interpol_FSOTB = Afgen2(FSOTB)
        Interpol_RDRTB = Afgen2(RDRTB)
        Interpol_PHOTTB = Afgen2(PHOTTB)

        FRTWET = Interpol_FRTTB(np_DVS[:])
        FLVT = Interpol_FLVTB(np_DVS[:])
        FSTT = Interpol_FSTTB(np_DVS[:])
        FSOT = Interpol_FSOTB(np_DVS[:])
        RDRTMP = Interpol_RDRTB(np_DVS[:])
        PHOTT = Interpol_PHOTTB(DAYL)

        EMERG = CropStarted & Enough_water & Leaves_Present & Not_Finished
        np_EMERG = pcr_as_numpy(EMERG)
        PHOTPF = ifthenelse(BeforeAnthesis, PHOTT, One)
        RTSUMP = DTEFF * PHOTPF
        # self.TSUM             += ifthenelse(CropStartNow, TSUMINIT, 0.) + ifthenelse(EMERG, RTSUMP, 0.) * ifthenelse(CropHarvNow, Zero, 1.) - ifthenelse(CropHarvNow, self.TSUM ,0.)
        self.TSUM = (self.TSUM + ifthenelse(CropStartNow, TSUMINIT, 0.) + ifthenelse(EMERG, RTSUMP, 0.)) * ifthenelse(
            CropHarvNow, Zero, 1.)  # - ifthenelse(CropHarvNow, self.TSUM ,0.)

        TSUM_veg = self.TSUM / TSUMAN * ifthenelse(CropHarvNow, Zero, 1.)
        TSUM_gen = (1. + (self.TSUM - TSUMAN) / TSUMMT) * ifthenelse(CropHarvNow, Zero, 1.)
        self.DVS = ifthenelse(Vegetative, TSUM_veg, 0.) + ifthenelse(Generative, TSUM_gen, 0.)
        # In TSUM_gen, '1' in fact stands for TSUM/TSUMAN, at the moment that TSUM = TSUMAN

        np_IRRAD = pcr_as_numpy(self.IRRAD)

        ########################################################################################
        # Root depth growth:

        np_CanGrowDownward = np.less_equal(np_ROOTD_mm[:], ROOTDM_mm * np_One[:])
        np_RROOTD_mm = np_Enough_water[:] * np_BeforeAnthesis[:] * np_EMERG[:] * np_CanGrowDownward[:] * RRDMAX_mm
        # self.ROOTD_mm        += ifthenelse(CropStartNow, ROOTDI_mm, Zero) + numpy2pcr(Scalar, np_RROOTD_mm[:], -99) * ifthenelse(CropHarvNow, Zero, 1.)- ifthenelse(CropHarvNow, self.ROOTD_mm ,0.)
        self.ROOTD_mm = (self.ROOTD_mm + ifthenelse(CropStartNow, ROOTDI_mm, Zero) + numpy2pcr(Scalar, np_RROOTD_mm[:],
                                                                                               -99)) * ifthenelse(
            CropHarvNow, Zero, 1.)
        # self.ROOTD_mm         = self.ROOTD * 1000.
        np_EXPLOR = np_RROOTD_mm[:] * np_WCFC[:]  # now for ROOTD in mm => factor 1000 omitted (!)

        #############################################################################################################
        if self.BMI_RUN == "True" and self.WATERLIMITED == "True":
            np_Transpiration = pcr_as_numpy(self.Transpiration)
            np_PotTrans = pcr_as_numpy(self.PotTrans)
            TRANRF = np_Transpiration[:] / NOTNUL(np_PotTrans[:])
            print "BMI run with water limitation"
        elif self.BMI_RUN == "True" and self.WATERLIMITED == "False":
            TRANRF = np_One[:]
            print "BMI run w/o water limitation"
        else:
            print "BMI not engaged => no water reduction"
            TRANRF = np_One[:]

            #############################################################################################################
            # Water Limitation: effects on partitioning
        FRTMOD = np_One[:]  # np.maximum(np_One[:], np_One[:]/(TRANRF[:]+ 0.5 * np_One[:]))
        FRTMOD = np.maximum(np_One[:], np_One[:] / (TRANRF[:] + 0.5 * np_One[:]))
        # FRT                   = FRTMOD[:] #FRTWET[:] * FRTMOD[:]
        FRT = FRTWET[:] * FRTMOD[:]
        FSHMOD = (1. - FRT[:]) / (1. - (FRT[:] / FRTMOD[:]))
        FLV = FLVT[:] * FSHMOD
        FST = FSTT[:] * FSHMOD
        FSO = FSOT[:] * FSHMOD
        PartitionCheck = FLV[:] + FST[:] + FSO[:] + FRT[:]

        # todo: incorporate effects of N stress - sdv 30-11-15
        PARINT = 0.5 * np_IRRAD[:] * 0.001 * (1. - np.exp(-K * np_LAI[:])) * np_Not_Finished[:]
        # pcr_PARINT            = 0.5 * self.IRRAD * 0.001 * (1.- exp(-K * self.LAI)) * scalar(Not_Finished)
        # np_pcr_PARINT         = pcr_as_numpy(pcr_PARINT)
        # check = np.equal(PARINT, np_pcr_PARINT)
        GTOTAL = self.LUE * PARINT[:] * TRANRF[:]
        # FRT, FLV, FST, FSO    = dryMatterPartitioningFractions(self, NPART, TRANRF, NNI, FRTWET[:], FLVT[:], FSTT[:], FSOT[:])




        # ----------------------------------------------------------------------
        # Rel. Death rate due to ageing:
        RDRDV = np_AtAndAfterAnthesis[:] * RDRTMP[:]
        RDRSH = np.maximum(np_Zero[:], RDRSHM * (np_LAI[:] - np_LAICR[:]) / np_LAICR[:])
        RDR = np.maximum(RDRDV[:], RDRSH[:]) * np_Not_Finished[:]

        # (Abs.) impact of leaf dying on leaf weight - todo
        N_Limitation = np.less(np_NNI[:], np_One[:])
        DLVNS = np_CropStarted[:] * N_Limitation[:] * np_WLVG[:] * RDRNS * (1. - np_NNI[:])
        DLVS = np_WLVG[:] * RDR[:]
        DLV = (DLVS[:] + DLVNS[:]) * np_Not_Finished[:]

        RWLVG = GTOTAL[:] * FLV[:] - DLV[:]
        # np_WLVG[:]           += np_CropStartNow[:] * WLVGI + RWLVG[:] * (1. - np_CropHarvNow[:]) - np_CropHarvNow[:] * np_WLVG[:]
        np_WLVG[:] = (np_WLVG[:] + np_CropStartNow[:] * WLVGI + RWLVG[:]) * (1. - np_CropHarvNow[:])
        np_WLVD[:] += DLV[:]

        # ----------------------------------------------------------------------
        # Leaf totalGrowthRate and LAI.
        GLV = FLV[:] * GTOTAL[:]
        GLV_pcr = numpy2pcr(Scalar, GLV, -99)

        GLAI = np_Adult[:] * SLA[:] * GLV[:] * (1. - np_CropHarvNow[:]) + \
               np_Juvenile[:] * (
               np_LAI[:] * (np.exp(np_RGRL[:] * np_DTEFF[:] * np_DELT[:]) - np_One[:]) / np_DELT[:]) * TRANRF[
                                                                                                       :] * np.exp(
                   -NLAI * (1.0 - NNI)) * (1. - np_CropHarvNow[:]) + \
               np.equal(np_LAI[:], np_Zero[:]) * np_Enough_water[:] * np_CropStartNow[:] * LAII / DELT

        # (Abs.) impact of leaf dying on LAI
        # Death of leaves due to ageing and shading:
        DLAIS = np_LAI[:] * RDR[:]
        DLAINS = np_CropStarted[:] * N_Limitation[:] * DLVNS[:] * SLA[:]
        DLAI = (DLAIS + DLAINS) * np_Not_Finished[:]
        # np_LAI[:]            += GLAI[:] - DLAI[:] - np_CropHarvNow[:] * np_LAI[:]
        np_LAI[:] = (np_LAI[:] + GLAI[:] - DLAI[:]) * (1. - np_CropHarvNow[:])

        # ----------------------------------------------------------------------
        # DRRT = np.greater_equal(np_DVS[:], np_DVSDR[:]) * np_WRT[:] * RDRRT
        DRRT = np_Roots_Dying[:] * np_WRT[:] * RDRRT
        RWRT = GTOTAL[:] * FRT[:] - DRRT[:] * np_Not_Finished[:]
        # np_WRT[:]            += np_CropStartNow[:] * WRTLI + RWRT[:] * (1. - np_CropHarvNow[:]) - np_CropHarvNow[:] * np_WRT[:]
        np_WRT[:] = (np_WRT[:] + np_CropStartNow[:] * WRTLI + RWRT[:]) * (1. - np_CropHarvNow[:])
        np_WDRT[:] += DRRT[:]

        # ----------------------------------------------------------------------
        RWSO = GTOTAL[:] * FSO[:]
        # np_WSO[:]            += np_CropStartNow[:] * WSOI + RWSO[:] * (1. - np_CropHarvNow[:]) - np_CropHarvNow[:] * np_WSO[:]
        np_WSO[:] = (np_WSO[:] + np_CropStartNow[:] * WSOI + RWSO[:]) * (1. - np_CropHarvNow[:])
        np_WSOTHA = np_WSO[:] / 100.

        # ----------------------------------------------------------------------
        RWST = GTOTAL[:] * FST[:]
        # np_WST[:]            += np_CropStartNow[:] * WSTI + RWST[:] * (1. - np_CropHarvNow[:]) - np_CropHarvNow[:] * np_WST[:]
        np_WST[:] = (np_WST[:] + np_CropStartNow[:] * WSTI + RWST[:]) * (1. - np_CropHarvNow[:])

        # ----------------------------------------------------------------------
        np_WLV = np_WLVG[:] + np_WLVD[:]
        TAGBM = np_WLV[:] + np_WST[:] + np_WSO[:]

        # ----------------------------------------------------------------------

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
        if o == '-l': exec "loglevel = logging." + a

    if (len(opts) <= 1):
        usage()

    # starttime = dt.datetime(1990,01,01)
    starttime = dt.datetime(1981, 9, 27)

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(myModel, _lastTimeStep, firstTimestep=_firstTimeStep, datetimestart=starttime)
    dynModelFw.createRunId(NoOverWrite=False, level=loglevel)

    dynModelFw.setupFramework()
    dynModelFw._runInitial()
    dynModelFw._runResume()
    # dynModelFw._runDynamic(0,0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
