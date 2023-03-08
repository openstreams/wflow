#!/usr/bin/python
#

import os.path
import math

import pcraster.framework
import pcraster as pcr
from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *

"""

wflow_lintul simulates potential or water-limited rice production using weather data as forcing variables. For water-limited production,  
soil data are required. wflow_lintul is a modified version of LINTUL1 (for simulating potential growth) and
LINTUL2 (for water-limited production). Rice-specific features were derived from LINTUL3 (N-limited production; Shibu et al., 2010).

Usage:
wflow_lintul  -C case -R Runid -c inifile

    -C: set the name  of the case (directory) to run

    -R: set the name runId within the current case

    -c name of the config file (in the case directory)

$Author: SanderCdeVries, for Wageningen Plant Research, WUR, The Netherlands
$Id: lintul1.py 2018-02-21 16:08:06Z
$Rev: 001 $
"""

# This needs to be set according to the geographic extent (map dimensions) of your study area/catchment:
# np_Zero = numpy.zeros((219, 286))
# np_One = numpy.ones((219, 286))

# Some last remaining hardcoded (& partly non-functional) parameters:
DELT = 1.0  # Time step (delta T) = 1 day; changing it is not recommended.
TINY = 1e-6  # A tiny number (from the original LINTUL code:)
WCWP = 0.21  # Volumetric soil water content at wilting point... how soil specific is this for puddled soils...? todo
WCFC = 0.47  # Volumetric soil water content at field capacity... how soil specific is this for puddled soils...? todo
WCST = 0.55  # Volumetric soil water content at saturation (normal condition for irrigated rice soil) ... how soil specific is this for puddled soils...? todo
NNI = 1.0  # Nitrogen Nutrition Index (non-functional, for future development)
NPART = 1.0  # Coefficient for the effect of N stress on leaf biomass reduction (presently non-functional, for future development)
NSLA = 1.0  # Coefficient for the effect of N stress on SLA reduction (presently non-functional, for future development)
NLAI = 1.0  # Coefficient for the effect of N stress on LAI reduction(during juvenile phase; presently non-functional, for future development)


def NOTNUL_pcr(pcr_map):
    """
    NOTNUL was originally a FST Fortran Simulation Translator intrinsic function.
    Here it is applied to arrays. If a value in the array is positive, NOTNUL will
    just return the value as is. If it equals zero: NOTNUL will return a value of 1 instead.
    Sander de Vries, March 2018
    """
    checkzeros = pcr_map == 0.0
    checkzeros_scalar = pcr.scalar(checkzeros)
    pcr_map += checkzeros_scalar
    return pcr_map


def astro_py(DAY, LAT):
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
    sinLAT = math.sin(math.pi * LAT / 180.0)
    cosLAT = math.cos(math.pi * LAT / 180.0)

    # MAXIMAL SINE OF DECLINATION
    sinDCM = math.sin(math.pi * 23.45 / 180.0)

    # SINE AND COSINE OF DECLINATION  (EQUATIONS 3.4, 3.5)

    # SINDEC = -SINDCM * cos(2.* PI * (DAY+10.)/365.)
    # The '9' below (instead of 10) keeps things in sync with FST...
    # Todo: try to understand this at a certain point... for now it works perfectly.
    sinDEC = -sinDCM * math.cos(2.0 * math.pi * (DAY + 11.0) / 365.0)
    cosDEC = math.sqrt(1.0 - sinDEC * sinDEC)

    # THE TERMS A AND B ACCORDING TO EQUATION 3.3

    A = sinLAT * sinDEC
    B = cosLAT * cosDEC

    # DAYLENGTH ACCORDING TO EQUATION 3.6.
    DAYL = 12.0 * (1.0 + (2.0 / math.pi) * math.asin(A / B))

    return DAYL


class Interpol_Obj(object):
    """
    Class to facilitate use of the 'lookuplinear' PCraster function.
    Upon initialization of an interpolation object, a temporary file
    containing x, y value pairs is created and saved in the case directory.
    This file is accessed by (PCraster) lookuplinear if the lookup_linear method is called.

    Sander de Vries, March 2018
    """

    def __init__(self, name):
        self.data = name
        self.name = name[-1]
        self.filename = self.name + ".tmp"

        temptablefile = open(self.filename, "w")
        index = list(range(0, len(self.data) - 1))
        for i in index:
            if i < (len(self.data) - 1):
                if i >= i + 1:
                    print(
                        "x values of lookuplinear table not sorted in strictly ascending order..."
                    )
            if i // 2.0 - i / 2.0 != 0.0:
                string = str(self.data[i]) + " "
            else:
                string = "\n" + str(self.data[i]) + " "
            temptablefile.write(string)
        temptablefile.close()

    def lookup_linear(self, x):
        y = pcr.lookuplinear(self.filename, x)
        return y


def supplyCurrentTime(self):
    """
    *Optional*

      Supplies the current time in seconds after the start of the run
      This function is optional. If it is not set the framework assumes
      the model runs with daily timesteps.

      Output:

          - time in seconds since the start of the model run

    """

    return self.currentTimeStep(self) * int(
        configget(self.config, "model", "timestepsecs", "86400")
    )


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


class WflowModel(pcraster.framework.DynamicModel):
    """
    wflow_lintul
    """

    def __init__(self, cloneMap, Dir, RunDir, configfile):
        """
        *Required*

        The init function **must** contain what is shown below. Other functionality
        may be added by you if needed.

        """
        pcraster.framework.DynamicModel.__init__(self)
        pcr.setclone(Dir + "/staticmaps/" + cloneMap)
        self.runId = RunDir
        self.caseName = Dir
        self.Dir = Dir
        self.configfile = configfile
        self.SaveDir = os.path.join(self.Dir, self.runId)

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

        self.RainSumStart_Month = int(
            configget(self.config, "model", "RainSumStart_Month", "11")
        )
        self.RainSumStart_Day = float(
            configget(self.config, "model", "RainSumStart_Day", "1")
        )
        self.Sim3rdSeason = eval(
            configget(self.config, "model", "Sim3rdSeason", "False")
        )
        self.RainSumReq = float(configget(self.config, "model", "RainSumReq", "200."))
        self.Pause = int(configget(self.config, "model", "Pause", "13"))
        self.AutoStartStop = eval(
            configget(self.config, "model", "AutoStartStop", "False")
        )  # default changed to 'False', for running from 'crop profile' maps (CRPST.xxx) under DEWS. sdv 21-2-2018
        self.WATERLIMITED = configget(self.config, "model", "WATERLIMITED", "True")
        self.CropStartDOY = (
            int(configget(self.config, "model", "CropStartDOY", "0")) - 1
        )  # to keep things in sync with the original LINTUL version in FST
        self.HarvestDAP = int(configget(self.config, "model", "HarvestDAP", "150"))
        self.LAT = float(configget(self.config, "model", "LAT", "3.16"))
        self.TSUMI = float(configget(self.config, "model", "TSUMI", "362."))
        self.K = float(configget(self.config, "model", "K", "0.6"))
        self.LUE = float(
            configget(self.config, "model", "LUE", "2.47")
        )  # The default value from Shibu et al. (2010) is 3.0; 2.47 was obtained by calibration for central Java. (sdv)
        self.SLAC = float(configget(self.config, "model", "SLAC", "0.02"))
        self.TSUMAN = float(configget(self.config, "model", "TSUMAN", "1420."))
        self.TSUMMT = float(configget(self.config, "model", "TSUMMT", "580."))
        self.TBASE = float(configget(self.config, "model", "TBASE", "8."))
        self.RGRL = float(configget(self.config, "model", "RGRL", "0.009"))
        self.WLVGI = float(configget(self.config, "model", "WLVGI", "0.86"))
        self.WSTI = float(configget(self.config, "model", "WSTI", "0.71"))
        self.WRTLI = float(configget(self.config, "model", "WRTLI", "1.58"))
        self.WSOI = float(configget(self.config, "model", "WSOI", "0."))
        self.RDRNS = float(configget(self.config, "model", "RDRNS", "0.03"))
        self.DVSDR = float(configget(self.config, "model", "DVSDR", "0.8"))
        self.RDRRT = float(configget(self.config, "model", "RDRRT", "0.03"))
        self.RDRSHM = float(configget(self.config, "model", "RDRSHM", "0.03"))
        self.LAICR = float(configget(self.config, "model", "LAICR", "4."))
        self.ROOTDM_mm = float(configget(self.config, "model", "ROOTDM_mm", "1000."))
        self.RRDMAX_mm = float(configget(self.config, "model", "RRDMAX_mm", "10."))
        self.ROOTDI_mm = float(configget(self.config, "model", "ROOTDI_mm", "50."))
        self.NLAI = float(configget(self.config, "model", "NLAI", "1."))
        self.RDRTB = eval(
            configget(
                self.config,
                "model",
                "RDRTB",
                "[0.0, 0.00 , 0.6 , 0.00, 1.0 , .015, 1.6 , 0.025, 2.1 , 0.05, 'RDRTB']",
            )
        )
        self.PHOTTB = eval(
            configget(
                self.config,
                "model",
                "PHOTTB",
                "[0.0, 0.0  , 8.  , 0.0 , 10. , 1.0 , 12. , 1.0  , 13. , 0.8  , 14.,  0.6 , 18. , 0.0, 'PHOTTB']",
            )
        )
        self.SLACF = eval(
            configget(
                self.config,
                "model",
                "SLACF",
                "[0.0, 1.72 , 0.21, 1.72, 0.24, 1.72, 0.33, 1.32 , 0.7 , 1.20 , 1.01, 1.00, 2.0 , 0.75, 2.1 , 0.75, 'SLACF']",
            )
        )  # for testing/development
        self.NMXLV = eval(
            configget(
                self.config,
                "model",
                "NMXLV",
                "[0.0, 0.05 , 0.4 , 0.05, 0.7 , 0.04, 1.0 , 0.03 , 2.0 , 0.02 , 2.1 , 0.02]",
            )
        )
        self.FRTTB = eval(
            configget(
                self.config,
                "model",
                "FRTTB",
                "[0.0, 0.300, 0.48, 0.30, 0.62, 0.12, 0.69, 0.11 , 0.84, 0.11 , 0.92, 0.10, 1.00, 0.08, 1.38, 0.00, 2.10, 0.0, 'FRTTB']",
            )
        )  # for testing/development
        self.FLVTB = eval(
            configget(
                self.config,
                "model",
                "FLVTB",
                "[0.0, 0.315, 0.48, 0.35, 0.62, 0.44, 0.69, 0.463, 0.84, 0.463, 0.92, 0.45, 1.00, 0.00, 1.38, 0.00, 2.10, 0.0, 'FLVTB']",
            )
        )  # for testing/development
        self.FSTTB = eval(
            configget(
                self.config,
                "model",
                "FSTTB",
                "[0.0, 0.385, 0.48, 0.35, 0.62, 0.44, 0.69, 0.427, 0.84, 0.427, 0.92, 0.27, 1.00, 0.00, 1.38, 0.00, 2.10, 0.0, 'FSTTB']",
            )
        )
        self.FSOTB = eval(
            configget(
                self.config,
                "model",
                "FSOTB",
                "[0.0, 0.00 , 0.48, 0.00, 0.62, 0.00, 0.69, 0.00 , 0.84, 0.00 , 0.92, 0.18, 1.00, 0.92, 1.38, 1.00, 2.10, 1.00, 'FSOTB']",
            )
        )

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
        modelparameters.append(
            self.ParamType(
                name="IRRAD",
                stack="inmaps/IRRAD",
                type="timeseries",
                default=11.0,
                verbose=False,
                lookupmaps=[],
            )
        ),
        modelparameters.append(
            self.ParamType(
                name="T",
                stack="inmaps/T",
                type="timeseries",
                default=10.0,
                verbose=False,
                lookupmaps=[],
            )
        ),
        modelparameters.append(
            self.ParamType(
                name="TMIN",
                stack="inmaps/TMIN",
                type="timeseries",
                default=10.0,
                verbose=False,
                lookupmaps=[],
            )
        ),
        modelparameters.append(
            self.ParamType(
                name="TMAX",
                stack="inmaps/TMAX",
                type="timeseries",
                default=10.0,
                verbose=False,
                lookupmaps=[],
            )
        ),
        modelparameters.append(
            self.ParamType(
                name="RAIN",
                stack="inmaps/P",
                type="timeseries",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        ),
        modelparameters.append(
            self.ParamType(
                name="CRPST",
                stack="inmaps/CRPST",
                type="timeseries",
                default=11.0,
                verbose=False,
                lookupmaps=[],
            )
        ),
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

        states = [
            "Season",
            "PSUM",
            "Test",
            "LAI",
            "WLVG",
            "WLVD",
            "WST",
            "WSO",
            "WRT",
            "ROOTD_mm",
            "WDRT",
            "TSUM",
            "STARTED",
            "DVS",
        ]

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

        return self.currentTimeStep() * int(
            configget(self.config, "model", "timestepsecs", "86400")
        )

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
        pcr.setglobaloption("unittrue")

        self.timestepsecs = int(
            configget(self.config, "model", "timestepsecs", "86400")
        )
        self.basetimestep = 86400

        # Reads all parameter from disk
        self.wf_updateparameters()
        self.logger.info("Starting LINTUL Dynamic Crop Growth Simulation...")

        # Read a static map of the rice area. To be replaced with real-time radar images of rice area in the future? Todo
        # Simulation is mostly restricted to the rice area (to be checked), which saves calculation time. Todo
        wflow_ricemask = configget(
            self.config, "model", "wflow_ricemask", "staticmaps/wflow_ricemask.map"
        )
        self.ricemask = self.wf_readmap(
            os.path.join(self.Dir, wflow_ricemask), 0.0, fail=True
        )
        # Create a PCRaster boolean map too:
        self.ricemask_BOOL = pcr.boolean(self.ricemask)
        self.Pausedays = self.Pause + 1

        # Calculate initial development stage (at the time of transplanting)
        self.DVSI = self.TSUMI / self.TSUMAN

        # Turn all interpolation tables (model parameters) into instances of the Interpol_Obj class
        self.RDRTB = Interpol_Obj(self.RDRTB)
        self.PHOTTB = Interpol_Obj(self.PHOTTB)
        self.SLACF = Interpol_Obj(self.SLACF)
        self.FRTTB = Interpol_Obj(self.FRTTB)
        self.FLVTB = Interpol_Obj(self.FLVTB)
        self.FSTTB = Interpol_Obj(self.FSTTB)
        self.FSOTB = Interpol_Obj(self.FSOTB)

        # Calculate the initial leaf area correction function as a function of development stage, DVS.
        SLACFI = self.SLACF.lookup_linear(self.DVSI)
        # Multiply with specific leaf area constant => initial specific leaf area
        ISLA = self.SLAC * SLACFI
        # Multiply with weight of green leaves to obtain initial LAI
        self.LAII = self.WLVGI * ISLA
        # Calculate total temperature sum from transplanting to crop maturity:
        self.TTSUM = self.TSUMAN + self.TSUMMT

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
                exec("self." + s + " = pcr.cover(0.0)")
        else:
            self.wf_resume(self.Dir + "/instate/")

            # try:
            #    self.wf_resume(self.Dir + "/instate/")
            # except:
            #    self.logger.warning("Cannot load initial states, setting to default")
            #    for s in self.stateVariables():
            #        exec "self." + s + " = pcr.cover(1.0)"

    def default_summarymaps(self):
        """
        *Optional*

        Return a default list of variables to report as summary maps in the outsum dir.
        The ini file has more options, including average and sum
        """
        return ["self.Altitude"]

        ########################################################################################

    def dynamic(self):
        """
        *Required*

        This is where all the time dependent functions are executed. Time dependent
        output should also be saved here.
        """
        self.wf_updateparameters()

        # Get the date as a Python datetime object and as the day of the year (DOY): used for cropping calendar and daylength.
        self.date = datetime.utcfromtimestamp(self.wf_supplyStartTime()) + dt.timedelta(
            self.currentTimeStep() - 1
        )
        self.enddate = datetime.utcfromtimestamp(
            self.wf_supplyEndTime()
        )  # wf_supplyEndTime() in wflow_dyamicframework? todo
        DOY = self.wf_supplyJulianDOY()

        # Some Boolean PCRaster variables to check if crop is still developing, in terms of thermal time/phenology:
        TSUM_not_Finished = self.TSUM <= self.TTSUM
        DVS_not_Finished = self.DVS <= 2.01
        Not_Finished = TSUM_not_Finished

        # Start calculating the accumulated preciptation from a certain date on (defined by RainSumStart_Month, RainSumStart_Day),
        # to judge when there's enough water for rice crop establishment.
        if (
            self.date.month == self.RainSumStart_Month
            and self.date.day == self.RainSumStart_Day
        ):
            Calc_RainSum = True
            self.PSUM += self.RAIN + TINY
        else:
            Calc_RainSum = False

        # Check whether the precipitation sum is positive
        WeveGotRain = self.PSUM > 0.0
        # Check whether the precipitation sum is still below the threshhold for crop establishment (and hence calculation of the sum should still proceed):
        NotEnoughRainYet = self.PSUM <= self.RainSumReq
        EnoughRain = self.PSUM >= self.RainSumReq
        KeepAddingRain = WeveGotRain & NotEnoughRainYet

        # The first season is defined here as starting on November 1. The 2md and 3rd season are following the 1st with break periods of <self.Pausedays> days,
        # to account for the time that the farmer needs for rice harvesting and crop establishment.
        # self.Season      += pcr.ifthenelse(Calc_RainSum, self.ricemask, 0.)
        FirstSeason = self.Season == 1
        SecondSeason = self.Season == 2
        ThirdSeason = self.Season == 3
        self.Season += pcr.ifthenelse(
            EnoughRain, self.ricemask, 0.0
        )  # beware, this variable is also modified in another equation
        # Add rain when the precipitation sum is positive but still below the threshold for crop establishment, reset to 0. when this is no longer the case.
        self.PSUM = (
            self.PSUM + pcr.ifthenelse(KeepAddingRain, self.RAIN, 0.0)
        ) * pcr.ifthenelse(KeepAddingRain, pcr.scalar(1.0), 0.0)

        # Initializing crop harvest:
        # If a fixed planting and a fixed harvest date are forced for the whole catchment:
        if self.CropStartDOY > -1:

            if self.HarvestDAP > 0:
                HarvNow = DOY >= (self.CropStartDOY + self.HarvestDAP)
                print(
                    "Warning: harvest date read from ini file, not from Crop Profile map..."
                )
            elif self.HarvestDAP == 0:
                HarvNow = Not_Finished == False
                print("Harvest date not specified; crop harvest at crop maturity")
            else:
                print(
                    "Crop harvest not initialized, found strange values in ini file... CTRL + C to exit..."
                )
                time.sleep(100)
            CropHarvNow = HarvNow & self.ricemask_BOOL

            # Initializing crop growth, optionally from a single start day (CropStartDOY in the ini file),
            # but normally from a crop profile forcing variable.
            StartNow = DOY == self.CropStartDOY
            CropStartNow = StartNow & self.ricemask_BOOL
            CropStartNow_scalar = pcr.scalar(CropStartNow)
            Started = self.STARTED > 0
            CropStarted = Started & self.ricemask_BOOL
            self.STARTED = (
                self.STARTED + CropStartNow_scalar + pcr.scalar(CropStarted)
            ) * pcr.ifthenelse(CropHarvNow, pcr.scalar(0.0), 1.0)
            print(
                "Warning: using start date from ini file, not read from Crop Profile..."
            )

        elif self.CropStartDOY == -1:

            if self.AutoStartStop == False:
                Started = self.STARTED > 0.0
                if self.HarvestDAP == 0:
                    crpprfl_eq_zero = self.CRPST == 0.0
                    CropHarvNow = Started & crpprfl_eq_zero & self.ricemask_BOOL
                elif self.HarvestDAP > 0:
                    HarvNow = self.STARTED == self.HarvestDAP
                    CropHarvNow = HarvNow & self.ricemask_BOOL
                print("Start date read from Crop Profile...")
                # Two auxilliary variables:
                CRPST_gt_0 = self.CRPST > 0.0
                CRPST_eq_STARTED = self.CRPST == self.STARTED
                CropStartNow = CRPST_gt_0 & CRPST_eq_STARTED & self.ricemask_BOOL
                CropStarted = Started & self.ricemask_BOOL
                self.STARTED = (self.STARTED + self.CRPST) * pcr.ifthenelse(
                    CropHarvNow, pcr.scalar(0.0), 1.0
                )  # - pcr.ifthenelse(CropHarvNow, self.STARTED, 0.)

            elif self.AutoStartStop == True:
                if self.HarvestDAP == 0:
                    HarvNow = (Not_Finished == False) | Calc_RainSum
                    CropHarvNow = HarvNow & self.ricemask_BOOL
                elif self.HarvestDAP > 0:
                    HarvNow = self.STARTED == self.HarvestDAP
                    CropHarvNow = (HarvNow & self.ricemask_BOOL) | Calc_RainSum
                # Two auxilliary variables:
                Time2Plant1stCrop = self.PSUM >= self.RainSumReq
                StdMin1 = self.STARTED == -1
                CropStartNow_Season1 = Time2Plant1stCrop & self.ricemask_BOOL
                CropStartNow_Season2 = StdMin1 & self.ricemask_BOOL
                CropStartNow = CropStartNow_Season1 | CropStartNow_Season2
                CropStartNow_scalar = pcr.scalar(CropStartNow)
                if self.Sim3rdSeason == False:
                    HarvSeason1_temp = FirstSeason & CropHarvNow
                    HarvSeasonOne = HarvSeason1_temp & self.ricemask_BOOL
                    HarvSeason2_temp = SecondSeason & CropHarvNow
                    HarvSeasonTwo = HarvSeason2_temp & self.ricemask_BOOL
                    self.Season = (
                        self.Season
                        + pcr.ifthenelse(HarvSeasonOne, self.ricemask, 0.0)
                        - pcr.ifthenelse(HarvSeasonTwo, self.ricemask * 2.0, 0.0)
                    )  # beware, this variable is also modified in another equation
                    Started = self.STARTED > 0
                    CropStarted = Started & self.ricemask_BOOL
                    SeasonOneHarvd = self.STARTED < 0
                    SeasonOneHarvd_Scalar = pcr.scalar(SeasonOneHarvd)
                    PrepareField_temp = pcr.scalar(self.Pausedays)
                    PrepareField = pcr.ifthenelse(FirstSeason, PrepareField_temp, 0.0)
                    self.STARTED = (
                        (self.STARTED + CropStartNow_scalar + pcr.scalar(CropStarted))
                        * pcr.ifthenelse(CropHarvNow, pcr.scalar(0.0), 1.0)
                        - pcr.ifthenelse(HarvSeasonOne, PrepareField, 0.0)
                        + SeasonOneHarvd_Scalar
                    )
                elif self.Sim3rdSeason == True:
                    HarvSeason12_temp = FirstSeason | SecondSeason
                    HarvSeasonOneTwo = HarvSeason12_temp & CropHarvNow
                    HarvSeasonThree = (ThirdSeason & CropHarvNow) | (
                        ThirdSeason & Calc_RainSum
                    )
                    self.Season = (
                        self.Season
                        + pcr.ifthenelse(HarvSeasonOneTwo, pcr.scalar(1.0), 0.0)
                        - pcr.ifthenelse(HarvSeasonThree, pcr.scalar(3.0), 0.0)
                    )  # beware, this variable is also modified in another equation
                    Started = self.STARTED > 0
                    CropStarted = Started & self.ricemask_BOOL
                    Season12Harvd = self.STARTED < 0
                    Season12Harvd_Scalar = pcr.scalar(Season12Harvd)
                    PrepareField_temp = pcr.scalar(self.Pausedays)
                    FirstorSecondSeason = FirstSeason | SecondSeason
                    PrepareField = pcr.ifthenelse(
                        FirstorSecondSeason, PrepareField_temp, 0.0
                    )
                    self.STARTED = (
                        (self.STARTED + CropStartNow_scalar + pcr.scalar(CropStarted))
                        * pcr.ifthenelse(CropHarvNow, pcr.scalar(0.0), 1.0)
                        - pcr.ifthenelse(HarvSeasonOneTwo, PrepareField, 0.0)
                        + Season12Harvd_Scalar
                    )
                else:
                    print(self.Sim3rdSeason)
                    time.sleep(10)

            else:
                print(
                    "Strange value of variable AutoStartStop found... ctrl + c to exit..."
                )
                time.sleep(100)
        else:
            print(
                "Strange (negative?) value of variable CropStartDOY found... ctrl + c to exit..."
            )
            time.sleep(100)

        if self.WATERLIMITED == "True":
            TRANRF = self.Transpiration / NOTNUL_pcr(self.PotTrans)
            WAWP = WCWP * self.ROOTD_mm
            Enough_water = pcr.ifthenelse(
                CropStartNow, True, self.WA > WAWP
            )  # timestep delay...! todo
        else:
            print("Warning, run without water effects on crop growth...")
            TRANRF = pcr.scalar(1.0)
            Enough_water = True

        # self.T = (self.TMIN + self.TMAX)/2. # for testing with Wageningen weather files only - sdv
        # Calculate thermal time (for TSUM and DVS):
        Warm_Enough = self.T >= self.TBASE
        DegreeDay = self.T - self.TBASE
        DTEFF = pcr.ifthenelse(Warm_Enough, DegreeDay, 0.0)
        # Check if leaves are present:
        Leaves_Present = self.LAI > 0.0

        # Check whether certain critical moments, external circumstances or crop growth stages occur that influence crop growth and development:
        BeforeAnthesis = self.TSUM < self.TSUMAN
        UntilAnthesis = self.TSUM <= self.TSUMAN
        AtAndAfterAnthesis = self.TSUM >= self.TSUMAN
        AfterAnthesis = self.TSUM > self.TSUMAN
        Roots_Dying = self.DVS >= self.DVSDR

        Vegetative = CropStarted & UntilAnthesis
        Generative = CropStarted & AfterAnthesis
        EarlyStages = self.DVS < 0.2
        LaterStages = self.DVS >= 0.2
        SmallLeaves = self.LAI < 0.75
        BiggerLeaves = self.LAI >= 0.75
        Juvenile = EarlyStages & SmallLeaves
        Adult = LaterStages | BiggerLeaves

        # Calculate daylength (assumed similar throughout the catchment area -> scalar, no array), based on latitude and Day Of Year (DOY)
        DAYL = astro_py(DOY, self.LAT)

        # Calculate the specific leaf area (m2 (leaf) g−1 (leaf)) by interpolation of development stage in SLAF, multiplication with self.SLACF
        SLA = self.SLAC * self.SLACF.lookup_linear(self.DVS)
        # Obtain the fractions (-) of daily dry matter production allocated (in absence of water shortage) to, respectively, root growth (FRTWET), leaf growth (FLVT),
        # growth of stems (FSTT) and growth of storage organs (FSO, i.e. rice grains), as a function of development stage (DVS), by interpolation.
        FRTWET = self.FRTTB.lookup_linear(self.DVS)
        FLVT = self.FLVTB.lookup_linear(self.DVS)
        FSTT = self.FSTTB.lookup_linear(self.DVS)
        FSOT = self.FSOTB.lookup_linear(self.DVS)
        RDRTMP = self.RDRTB.lookup_linear(self.DVS)

        # Many growth processes can only occur when EMERG = TRUE; this is the case when crop phenological development has started, soil water content is above
        # permanent wilting point, the crop has leaves and is not yet harvested or growth has otherwise been terminated:
        EMERG = CropStarted & Enough_water & Leaves_Present & Not_Finished

        # Determine the influence of astronomical daylength on crop development (via thermal time - TSUM) by interpolation in the PHOTTB table
        PHOTT = self.PHOTTB.lookup_linear(DAYL)
        # Daylength only potentially has a  modifying influence on crop development (via thermal time, TSUM) before anthesis:
        PHOTPF = pcr.ifthenelse(BeforeAnthesis, PHOTT, pcr.scalar(1.0))
        # Influence (if any) of daylength results in a modified daily change in thermal time (TSUM).
        RTSUMP = DTEFF * PHOTPF
        # TSUM (state): at crop establishment TSUMI is added (the TSUM that was accumulated in the nursery in the case of transplanted rice);
        # during crop growth, the daily rate of change RTSUMP is added if EMERG = TRUE. Upon crop harvest, TSUM is reset (i.e. multiplied with 0.).
        self.TSUM = (
            self.TSUM
            + pcr.ifthenelse(CropStartNow, pcr.scalar(self.TSUMI), 0.0)
            + pcr.ifthenelse(EMERG, RTSUMP, 0.0)
        ) * pcr.ifthenelse(CropHarvNow, pcr.scalar(0.0), 1.0)

        # Calculation of DVS (state).
        # In LINTUL1 and LINTUL2, TSUM directly steered all processes influenced by crop phenological development.
        # However, Shibu et al. (2010) derived some code from ORYZA_2000 (Bouman et al., 2001), including the use DVS instead of TSUM.
        # Hence in LINTUL3, some processes are still controlled directly by TSUM and some are controlled by its derived variable DVS
        # – a somewhat confusing situation that offers scope for future improvement.
        # After anthesis DVS proceeds at a different rate (DVS_gen) than before (DVS_veg). Throughout crop development DVS is calculated as the DVS_veg + DVS_gen.
        DVS_veg = (
            self.TSUM / self.TSUMAN * pcr.ifthenelse(CropHarvNow, pcr.scalar(0.0), 1.0)
        )
        DVS_gen = (1.0 + (self.TSUM - self.TSUMAN) / self.TSUMMT) * pcr.ifthenelse(
            CropHarvNow, pcr.scalar(0.0), 1.0
        )
        self.DVS = pcr.ifthenelse(Vegetative, DVS_veg, 0.0) + pcr.ifthenelse(
            Generative, DVS_gen, 0.0
        )

        # Root depth growth can occur as long as the maximum rooting depth has not yet been achieved:
        CanGrowDownward = self.ROOTD_mm <= self.ROOTDM_mm
        # Root growth occurs before anthesis if there is crop growth (EMERG = TRUE), enough water (already in EMERG - todo) and the maximum rooting depth has not yet been reached.
        RootGrowth = Enough_water & BeforeAnthesis & EMERG & CanGrowDownward
        # If root growth occurs, it occurs at a fixed pace (mm/day):
        RROOTD_mm = pcr.ifthenelse(RootGrowth, self.RRDMAX_mm, pcr.scalar(0.0))
        # Rooting depth (state): at crop establishment ROOTDI_mm is added (the rooting depth at transplanting); during crop growth, the daily rate of change
        # self.ROOTDI_mm is added. Upon crop harvest, rooting depth is reset (i.e. multiplied with 0.).
        self.ROOTD_mm = (
            self.ROOTD_mm
            + pcr.ifthenelse(CropStartNow, self.ROOTDI_mm, pcr.scalar(0.0))
            + RROOTD_mm
        ) * pcr.ifthenelse(CropHarvNow, pcr.scalar(0.0), 1.0)
        # By depth growth, roots explore deeper layers of soil that contain previously untapped water supplies, the assumption is.
        # In the case of irrigated rice, it seems reasonable to assume that those layers are saturated with water (WCST = volumetric soil water content at saturation).
        # The volume of additional water that becomes available to the crop is then equal to EXPLOR:
        EXPLOR = RROOTD_mm * WCST

        #############################################################################################################
        # Water Limitation: effects on partitioning
        # If TRANRF falls below 0.5, root growth is accelerated:
        FRTMOD = pcr.max(1.0, 1.0 / (TRANRF + 0.5))
        FRT = FRTWET * FRTMOD
        # ... and shoot growth (i.e. growth of all aboveground parts) diminshed:
        FSHMOD = (1.0 - FRT) / (1 - FRT / FRTMOD)
        FLV = FLVT * FSHMOD
        FST = FSTT * FSHMOD
        FSO = FSOT * FSHMOD

        # Daily intercepted Photosynthetically Active Radiation (PAR), according to (Lambert-)Beer's law.
        # The factor 0.5 accounts for the fact that about 50% (in terms of energy) of the frequency spectrum of incident solar radiation
        # can be utilized for photosynthesis by green plants.
        PARINT = pcr.ifthenelse(
            Not_Finished,
            0.5 * self.IRRAD * 0.001 * (1.0 - pcr.exp(-self.K * self.LAI)),
            0.0,
        )
        # The total growth rate is proportional to the intercepted PAR with a fixed Light Use Efficiency (LUE) - the core of the LINTUL apporach.
        GTOTAL = self.LUE * PARINT * TRANRF

        # Leaf dying due to ageing occurs from anthesis on (actually that is already arranged in the interpolation table - double!), with a relative death rate RDRTMP:
        RDRDV = pcr.ifthenelse(AtAndAfterAnthesis, RDRTMP, pcr.scalar(0.0))
        # Leaf dying due to mutual shading occurs when LAI > LAICR:
        RDRSH = pcr.max(0.0, self.RDRSHM * (self.LAI - self.LAICR) / self.LAICR)
        # The largest of the two effects determines the relative death rate of foliage:
        RDR = pcr.max(RDRDV, RDRSH)

        # Impact of leaf dying on leaf weight - N limitation stuff not (yet) implemented
        N_Limitation = NNI < 1.0
        DLVNS = pcr.ifthenelse(CropStarted, pcr.scalar(1.0), 0.0) * pcr.ifthenelse(
            N_Limitation, self.WLVG * self.RDRNS * (1.0 - NNI), 0.0
        )
        DLVS = self.WLVG * RDR
        DLV = (DLVS + DLVNS) * pcr.scalar(Not_Finished)
        RWLVG = pcr.ifthenelse(EMERG, GTOTAL * FLV - DLV, pcr.scalar(0.0))
        self.WLVG = (
            self.WLVG
            + pcr.ifthenelse(CropStartNow, self.WLVGI, pcr.scalar(0.0))
            + RWLVG
        ) * (1.0 - pcr.scalar(CropHarvNow))
        self.WLVD = (self.WLVD + DLV) * (1.0 - pcr.scalar(CropHarvNow))

        # Growth of leaves in terms of mass (GLV) and in terms of LAI (GLAI).
        GLV = FLV * GTOTAL
        Adt_or_Harv = pcr.pcror(Adult, CropHarvNow)
        Juv_or_Harv = pcr.pcror(Juvenile, CropHarvNow)
        NoLeavesYet = self.LAI == 0.0
        LetsGo = pcr.pcrand(Enough_water, CropStartNow)
        LetsGro = pcr.pcrand(NoLeavesYet, LetsGo)

        GLAI = (
            pcr.ifthenelse(Adt_or_Harv, SLA * GLV, pcr.scalar(0.0))
            + pcr.ifthenelse(
                Juv_or_Harv,
                self.LAI
                * (pcr.exp(self.RGRL * DTEFF * DELT) - 1.0)
                / DELT
                * TRANRF
                * pcr.exp(-self.LAI * (1.0 - NNI)),
                0.0,
            )
            + pcr.ifthenelse(LetsGro, self.LAII / DELT, pcr.scalar(0.0))
        )

        # (Abs.) impact of leaf dying on LAI
        # Daily decrease in LAI due to dying of leaves (if any), due to aging and/or mutual shading:
        DLAIS = self.LAI * RDR
        # Daily decrease in LAI due to nitrogen shortage (presently non-functional):
        DLAINS = pcr.ifthenelse(CropStarted, pcr.scalar(1.0), 0.0) * pcr.ifthenelse(
            N_Limitation, DLVNS * SLA, 0.0
        )
        # Total daily decrease in LAI due to leaf death (aging, mutual shading, N shortage):
        DLAI = (DLAIS + DLAINS) * pcr.scalar(Not_Finished)
        # The initial LAI (LAII, transplanted rice) is added to GLAI at crop establishment, not in below state equation as done by Shibu et al. (2010).
        self.LAI = (self.LAI + GLAI - DLAI) * pcr.ifthenelse(
            CropHarvNow, pcr.scalar(0.0), 1.0
        )

        # Daily death rate of roots: if self.DVS >= self.DVSDR, a fraction self.DRRT of the roots is dying every day:
        DRRT = pcr.ifthenelse(Roots_Dying, self.WRT * self.RDRRT, pcr.scalar(0.0))
        # Net daily change in root weight: if there is crop growth, this is equal to daily weight increase (GTOTAL * FRT) minus the daily decrease due to dying roots (if any).
        RWRT = pcr.ifthenelse(EMERG, GTOTAL * FRT - DRRT, pcr.scalar(0.0))
        # Calculation of the root weight (state): when the crop is planted, the initial root weight WRTLI is added. After that, the net (daily) rate of change RWRT is added.
        # Upon crop harvest, RWRT is reset (i.e. multiplied with 0.)
        self.WRT = (
            self.WRT + pcr.ifthenelse(CropStartNow, self.WRTLI, pcr.scalar(0.0)) + RWRT
        ) * (1.0 - pcr.scalar(CropHarvNow))
        # WDRT (state) is the total quantity of leaves that has died (mostly relevant for mass balance checking purposes). Simply calculated by adding the daily dying roots (if any);
        # (the variable is reset upon crop harvest).
        self.WDRT = (self.WDRT + DRRT) * (1.0 - pcr.scalar(CropHarvNow))

        # Daily change in the weight of the storage organs (rice grains) is fraction FSO of the total growth rate (GTOTAL). FSO is determined by phenology and moisture stress
        # (which can modify the root/shoort ratio)
        RWSO = pcr.ifthenelse(EMERG, GTOTAL * FSO, pcr.scalar(0.0))
        # Weight of storage organs (state) is simply calculated by accumulating the daily growth RWSO. At crop harvest, it is reset to 0.
        self.WSO = (
            self.WSO + pcr.ifthenelse(CropStartNow, self.WSOI, pcr.scalar(0.0)) + RWSO
        ) * (1.0 - pcr.scalar(CropHarvNow))
        # WSO in tons/ha:
        WSOTHA = self.WSO / 100.0

        # Daily change in the weight of the stems is a fraction FST of the total growth rate (GTOTAL). FST is determined by phenology and moisture stress
        RWST = pcr.ifthenelse(EMERG, GTOTAL * FST, pcr.scalar(0.0))
        # Weight of storage organs (state) is simply calculated by accumulating the daily growth RWSO. At crop harvest, it is reset to 0.
        self.WST = (
            self.WST + pcr.ifthenelse(CropHarvNow, self.WSTI, pcr.scalar(0.0)) + RWST
        ) * (1.0 - pcr.scalar(CropHarvNow))

        # An additional state, handy for testing purposes:
        self.Test += 1.0

    # For quickly getting point output (sdv). Works only with a wf_supplyEndTime() implemented in wf_dynamicframework... todo?
    # Point_Output = open('Point_Output.csv', 'w')
    # Point_Output_Line = (str(cellvalue (self.LAI, 100,100)[0]) ) + '\n'
    #                    #str(np_CropStartNow[100,100])+ "," + str(np_CropHarvNow[100,100]) + '\n')
    # if self.date < self.enddate:
    #    Point_Output.write(Point_Output_Line)
    # elif self.date == self.enddate:
    #    Point_Output.close()


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
    wflow_cloneMap = "wflow_subcatch.map"
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
        opts, args = getopt.getopt(argv, "F:C:S:T:c:s:R:l")
    except getopt.error as msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == "-C":
            caseName = a
        if o == "-R":
            runId = a
        if o == "-c":
            configfile = a
        if o == "-s":
            timestepsecs = int(a)
        if o == "-T":
            _lastTimeStep = int(a)
        if o == "-S":
            _firstTimeStep = int(a)
        if o == "-f":
            _NoOverWrite = 0
        if o == "-l":
            exec("loglevel = logging." + a)

    if len(opts) <= 1:
        usage()

    # starttime = dt.datetime(1990, 1, 1)
    starttime = dt.datetime(1981, 9, 27)

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(
        myModel, _lastTimeStep, firstTimestep=_firstTimeStep, datetimestart=starttime
    )
    dynModelFw.createRunId(NoOverWrite=False, level=loglevel)
    # dynModelFw.createRunId(NoOverWrite=_NoOverWrite, level=loglevel, logfname=LogFileName,model="wflow_lintul",doSetupFramework=False)

    dynModelFw.setupFramework()
    dynModelFw._runInitial()
    dynModelFw._runResume()
    # dynModelFw._runDynamic(0,0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
