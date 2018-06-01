#!/usr/bin/python

"""
Definition of the wflow_W3RA model.
---------------------------------------
The model is modified from the Australian Water Resources Assessment Landscape (AWRA-L) model version 0.5
W3RA is documented in van Dijk et al. (2013), Water Resour. Res., 49, 2729-2746, doi:10.1002/wrcr.20251
URL: http://onlinelibrary.wiley.com/doi/10.1002/wrcr.20251/abstract
More comprehensive documentation of AWRA-L version 0.5 can be found in:
Van Dijk, A.I.J.M. (2010) The Australian water resources assessment system
(version 0.5), 3.0.5.Technical description of the landscape hydrology model
(AWRA-L). WIRADA Technical Report, CSIRO Water for a Healthy Country
Flagship, Canberra.
URL: http://www.clw.csiro.au/publications/waterforahealthycountry/2010/wfhc-aus-water-resources-assessment-system.pdf
The section references below refer to the sections in the AWRA-L report.
Changes compared to that code are indicated, e.g. by commenting out
redundant code.
Further question please contact albert.vandijk@anu.edu.au
Port to Python/PCRaster: Deltares
Usage:
wflow_W3RA  -C case -R Runid -c inifile
    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -c name of the config file (in the case directory)
    
$Author: schelle $
$Id: wflow_sceleton.py 898 2014-01-09 14:47:06Z schelle $
$Rev: 898 $
"""

import numpy
import os
import os.path
import shutil, glob
import getopt

from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *

# TODO: Make the script HRU independent (loop over the nr of HRU's)
# TODO:


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print(msg)
    print(__doc__)
    sys.exit(0)


def pcr_tanh(x):
    """
    define tanh for pcraster objects
    
    """
    return (exp(x) - exp(-x)) / (exp(x) + exp(-x))


class WflowModel(DynamicModel):
    """
  The user defined model class. T
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
        self.SaveDir = self.Dir + "/" + self.runId + "/"

    def stateVariables(self):
        """ 
      *Required*
      
      Returns a list of state variables that are essential to the model. 
      This list is essential for the resume and suspend functions to work.
      
      This function is specific for each model and **must** be present. This is
      where you specify the state variables of you model. If your model is stateless
      this function must return and empty array (states = [])
      """

        states = [
            "S0",
            "Ss",
            "Sd",
            "Mleaf",
            "FreeWater",
            "DrySnow",
            "Sg",
            "Sr",
            "OpenWaterFrac",
        ]

        return states

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
        self.wf_suspend(self.SaveDir + "/outstate/")

        if self.fewsrun:
            self.logger.info("Saving initial conditions for FEWS...")
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
        setglobaloption("radians")  # Needed as W3RA was originally written in matlab

        self.timestepsecs = int(
            configget(self.config, "model", "timestepsecs", "86400")
        )
        self.UseETPdata = int(
            configget(self.config, "model", "UseETPdata", "1")
        )  #  1: Use ETP data, 0: Compute ETP from meteorological variables
        self.logger.debug("use DATA: " + str(self.UseETPdata))
        self.basetimestep = 86400
        self.SaveMapDir = self.Dir + "/" + self.runId + "/outmaps"

        # Define here the W3RA mapstacks (best to read these via netcdf)

        self.TMAX_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "TMAX", "/inmaps/TMAX"
        )
        self.TMIN_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "TMIN", "/inmaps/TMIN"
        )
        self.TDAY_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "TDAY", "/inmaps/TDAY"
        )
        self.EPOT_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "EPOT", "/inmaps/EPOT"
        )
        self.PRECIP_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "PRECIP", "/inmaps/PRECIP"
        )
        self.RAD_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "RAD", "/inmaps/RAD"
        )
        # self.WINDSPEED_mapstack=self.Dir + configget(self.config,"inputmapstacks","WINDSPEED","/inmaps/ClimatologyMapFiles/WINDS/WNDSPEED")
        # self.AIRPRESS_mapstack=self.Dir + configget(self.config,"inputmapstacks","AIRPRESS","/inmaps/ClimatologyMapFiles/AIRPRESS/AIRPRESS")
        self.ALBEDO_mapstack = self.Dir + configget(
            self.config,
            "inputmapstacks",
            "ALBEDO",
            "/inmaps/ClimatologyMapFiles/ALBEDO/ALBEDO",
        )
        self.WINDSPEED_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "WINDSPEED", "/inmaps/WIND"
        )
        self.AIRPRESS_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "AIRPRESS", "/inmaps/PRES"
        )

        self.Altitude = readmap(self.Dir + "/staticmaps/wflow_dem")

        self.fewsrun = int(configget(self.config, "model", "fewsrun", "0"))

        self.latitude = ycoordinate(boolean(self.Altitude))

        # Add reading of parameters here
        self.ER_coef = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/ER_coef.map"), 0.0, fail=True
        )
        self.flmp = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/flmp.map"), 0.0, fail=True
        )
        self.fPotDeep = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fPotDeep.map"), 0.0, fail=True
        )
        self.FsoilEmax = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/FsoilEmax.map"), 0.0, fail=True
        )
        self.Gs_scalar = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/Gs_scalar.map"), 0.0, fail=True
        )
        self.hand_percentile0 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile0.map"), 0.0, fail=True
        )
        self.hand_percentile1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile1.map"), 0.0, fail=True
        )
        self.hand_percentile2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile2.map"), 0.0, fail=True
        )
        self.hand_percentile3 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile3.map"), 0.0, fail=True
        )
        self.hand_percentile4 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile4.map"), 0.0, fail=True
        )
        self.hand_percentile5 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile5.map"), 0.0, fail=True
        )
        self.hand_percentile6 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile6.map"), 0.0, fail=True
        )
        self.hand_percentile7 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile7.map"), 0.0, fail=True
        )
        self.hand_percentile8 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile8.map"), 0.0, fail=True
        )
        self.hand_percentile9 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile9.map"), 0.0, fail=True
        )
        self.hand_percentile10 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile10.map"), 0.0, fail=True
        )
        self.hand_percentile11 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile11.map"), 0.0, fail=True
        )
        self.hand_percentile12 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile12.map"), 0.0, fail=True
        )
        self.hand_percentile13 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile13.map"), 0.0, fail=True
        )
        self.hand_percentile14 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile14.map"), 0.0, fail=True
        )
        self.hand_percentile15 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile15.map"), 0.0, fail=True
        )
        self.hand_percentile16 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile16.map"), 0.0, fail=True
        )
        self.hand_percentile17 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile17.map"), 0.0, fail=True
        )
        self.hand_percentile18 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile18.map"), 0.0, fail=True
        )
        self.hand_percentile19 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hand_percentile19.map"), 0.0, fail=True
        )
        self.hveg = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hveg.map"), 0.0, fail=True
        )
        self.K_gw = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/K_gw.map"), 0.0, fail=True
        )
        self.k_s = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/k_s.map"), 0.0, fail=True
        )
        self.K0_scalar = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/K0_scalar.map"), 0.0, fail=True
        )
        self.Ksat_exp = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/Ksat_exp.map"), 0.0, fail=True
        )
        self.lambd = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/lambda.map"), 0.0, fail=True
        )
        self.OpenWaterFrac = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/OpenWaterFrac.map"), 0.0, fail=True
        )
        self.porosity = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/porosity.map"), 0.0, fail=True
        )
        self.Pref = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/pref.map"), 0.0, fail=True
        )
        self.psi_s = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/psi_s.map"), 0.0, fail=True
        )
        self.S_sls = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/S_sls.map"), 0.0, fail=True
        )
        self.slope = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/slope.map"), 0.0, fail=True
        )
        self.snow_Cfmax = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/snow_Cfmax.map"), 0.0, fail=True
        )
        self.snow_Cfr = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/snow_Cfr.map"), 0.0, fail=True
        )
        self.snow_TT = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/snow_TT.map"), 0.0, fail=True
        )
        self.snow_WHC = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/snow_WHC.map"), 0.0, fail=True
        )
        self.T_offset = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/T_offset.map"), 0.0, fail=True
        )
        self.theta_s = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/theta_s.map"), 0.0, fail=True
        )
        # self.K_rout = self.wf_readmap(os.path.join(self.Dir,  "staticmaps/k_rout.map"),0.0,fail=True)
        # self.Sgref = self.wf_readmap(os.path.join(self.Dir, "staticmaps/sgref.map"),0.0,fail=True)
        # self.alb_dry = self.wf_readmap(os.path.join(self.Dir, "staticmaps/alb_dry.map"),0.0,fail=True)
        # self.alb_wet = self.wf_readmap(os.path.join(self.Dir, "staticmaps/alb_wet.map"),0.0,fail=True)
        # self.beta = self.wf_readmap(os.path.join(self.Dir, "staticmaps/beta.map"),0.0,fail=True)
        # self.cGsmax = self.wf_readmap(os.path.join(self.Dir, "staticmaps/cgsmax.map"),0.0,fail=True)
        # self.ER_frac_ref = self.wf_readmap(os.path.join(self.Dir, "staticmaps/er_frac_ref.map"),0.0,fail=True)
        # self.Fhru = self.wf_readmap(os.path.join(self.Dir, "staticmaps/fHRU.map"),0.0,fail=True)
        # self.FdrainFC = self.wf_readmap(os.path.join(self.Dir, "staticmaps/fdrainfc.map"),0.0,fail=True)
        # self.Fgw_conn = self.wf_readmap(os.path.join(self.Dir, "staticmaps/fgw_conn.map"),0.0,fail=True)
        # self.SLA  = self.wf_readmap(os.path.join(self.Dir, "staticmaps/sla.map"),0.0,fail=True)
        # self.LAIref = self.wf_readmap(os.path.join(self.Dir, "staticmaps/lairef.map"),0.0,fail=True)
        # self.fvegref_G = self.wf_readmap(os.path.join(self.Dir, "staticmaps/fvegref_g.map"),0.0,fail=True)
        # self.FwaterE = self.wf_readmap(os.path.join(self.Dir, "staticmaps/fwatere.map"),0.0,fail=True)
        # self.Gfrac_max = self.wf_readmap(os.path.join(self.Dir, "staticmaps/gfrac_max.map"),0.0,fail=True)
        # self.InitLoss = self.wf_readmap(os.path.join(self.Dir, "staticmaps/initloss.map"),0.0,fail=True)
        # self.LAImax = self.wf_readmap(os.path.join(self.Dir, "staticmaps/laimax.map"),0.0,fail=True)
        # self.S0FC = self.wf_readmap(os.path.join(self.Dir, "staticmaps/s0fc.map"),0.0,fail=True)
        # self.SsFC = self.wf_readmap(os.path.join(self.Dir, "staticmaps/ssfc.map"),0.0,fail=True)
        # self.SdFC = self.wf_readmap(os.path.join(self.Dir, "staticmaps/sdfc.map"),0.0,fail=True)
        # self.Vc  = self.wf_readmap(os.path.join(self.Dir, "staticmaps/vc.map"),0.0,fail=True)
        # self.w0ref_alb = self.wf_readmap(os.path.join(self.Dir, "staticmaps/w0ref_alb.map"),0.0,fail=True)
        # self.Us0  = self.wf_readmap(os.path.join(self.Dir, "staticmaps/us0.map"),0.0,fail=True)
        # self.Ud0  = self.wf_readmap(os.path.join(self.Dir, "staticmaps/ud0.map"),0.0,fail=True)
        # self.wslimU = self.wf_readmap(os.path.join(self.Dir, "staticmaps/wslimu.map"),0.0,fail=True)
        # self.wdlimU = self.wf_readmap(os.path.join(self.Dir, "staticmaps/wdlimu.map"),0.0,fail=True)
        # self.w0limE = self.wf_readmap(os.path.join(self.Dir, "staticmaps/w0lime.map"),0.0,fail=True)
        # self.Tgrow = self.wf_readmap(os.path.join(self.Dir, "staticmaps/tgrow.map"),0.0,fail=True)
        # self.Tsenc = self.wf_readmap(os.path.join(self.Dir, "staticmaps/tsenc.map"),0.0,fail=True)

        self.wf_multparameters()
        # Static, for the computation of Aerodynamic conductance (3.7)
        self.fh = ln(813. / max(self.hveg, 0.25) - 5.45)
        self.ku1 = 0.305 / (self.fh * (self.fh + 2.3))

        self.logger.info("Starting Dynamic run...")

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
            exec("self." + s + " = cover(1.0)")


  def default_summarymaps(self):
      """
      *Optional*
      Return a default list of variables to report as summary maps in the outsum dir.
      """
        return []

    def parameters(self):
        """
        Define all model parameters here that the framework should handle for the model
        See wf_updateparameters and the parameters section of the ini file
        If you use this make sure to all wf_updateparameters at the start of the dynamic section
        and at the start/end of the initial section
        :returns modelparameters: list of model parameters
        """
        modelparameters = []

        # Static model parameters e.g.
        # modelparameters.append(self.ParamType(name="RunoffGeneratingGWPerc",stack="intbl/RunoffGeneratingGWPerc.tbl",type="static",default=0.1))
        # 3: Input time series ###################################################
        # self.P_mapstack = self.Dir + configget(self.config, "inputmapstacks", "Precipitation",
        #                                       "/inmaps/P")  # timeseries for rainfall
        # self.PET_mapstack = self.Dir + configget(self.config, "inputmapstacks", "EvapoTranspiration",
        #                                         "/inmaps/PET")  # timeseries for rainfall"/inmaps/PET"          # potential evapotranspiration
        # self.TEMP_mapstack = self.Dir + configget(self.config, "inputmapstacks", "Temperature",
        #                                          "/inmaps/TEMP")  # timeseries for rainfall "/inmaps/TEMP"          # global radiation
        # self.Inflow_mapstack = self.Dir + configget(self.config, "inputmapstacks", "Inflow",
        #                                            "/inmaps/IF")  # timeseries for rainfall "/inmaps/IF" # in/outflow locations (abstractions)

        # Meteo and other forcing
        # modelparameters.append(self.ParamType(name="Precipitation",stack=self.P_mapstack,type="timeseries",default=0.0,verbose=True,lookupmaps=[]))
        # modelparameters.append(self.ParamType(name="PotenEvap",stack=self.PET_mapstack,type="timeseries",default=0.0,verbose=True,lookupmaps=[]))
        # modelparameters.append(self.ParamType(name="Temperature",stack=self.TEMP_mapstack,type="timeseries",default=10.0,verbose=True,lookupmaps=[]))
        # modelparameters.append(self.ParamType(name="Inflow",stack=self.Inflow_mapstack,type="timeseries",default=0.0,verbose=False,lookupmaps=[]))
        return modelparameters

    def dynamic(self):
        """
          *Required*
          This is where all the time dependent functions are executed. Time dependent
          output should also be saved here.
        """
        # print 'useETPdata' , self.UseETPdata
        # Put the W3RA here. Stuff from W3RA_timestep_model.m
        # read meteo from file
        self.logger.debug("Running for: " + str(self.currentdatetime))
        self.PRECIP = cover(
            self.wf_readmap(self.PRECIP_mapstack, 0.0), scalar(0.0)
        )  # mm

        if self.UseETPdata == 1:
            self.TDAY = cover(
                self.wf_readmap(self.TDAY_mapstack, 10.0), scalar(10.0)
            )  # T in degC
            self.EPOT = cover(
                self.wf_readmap(self.EPOT_mapstack, 0.0), scalar(0.0)
            )  # mm
            # print "Using climatology for wind, air pressure and albedo."
        elif self.UseETPdata == 0:
            self.TMIN = cover(
                self.wf_readmap(self.TMIN_mapstack, 10.0), scalar(10.0)
            )  # T in degC
            self.TMAX = cover(
                self.wf_readmap(self.TMAX_mapstack, 10.0), scalar(10.0)
            )  # T in degC
            self.RAD = cover(
                self.wf_readmap(self.RAD_mapstack, 10.0), scalar(10.0)
            )  # W m-2 s-1
            self.WINDSPEED = cover(
                self.wf_readmap(self.WINDSPEED_mapstack, 10.0), scalar(10.0)
            )  # ms-1
            self.AIRPRESS = cover(
                self.wf_readmap(self.AIRPRESS_mapstack, 10.0), scalar(10.0)
            )  # Pa
            self.ALBEDO = cover(
                self.wf_readmapClimatology(self.ALBEDO_mapstack, default=0.1),
                scalar(0.1),
            )

        self.wf_multparameters()
        doy = self.currentdatetime.timetuple().tm_yday

        # conversion daylength
        setglobaloption("radians")
        m = scalar(1) - tan((self.latitude * scalar(pi) / scalar(180))) * tan(
            (
                (scalar(23.439) * scalar(pi) / scalar(180))
                * cos(scalar(2) * scalar(pi) * (doy + scalar(9)) / scalar(365.25))
            )
        )
        self.fday = min(
            max(
                scalar(0.02),
                scalar(acos(scalar(1) - min(max(scalar(0), m), scalar(2))))
                / scalar(pi),
            ),
            scalar(1),
        )  # fraction daylength

        # Assign forcing and estimate effective meteorological variables

        Pg = self.PRECIP  # mm

        if self.UseETPdata == 1:
            Ta = self.TDAY  # T in degC
            T24 = self.TDAY  # T in degC
        elif self.UseETPdata == 0:
            Rg = max(
                self.RAD, scalar(0.0001)
            )  # already in W m-2 s-1; set minimum of 0.01 to avoid numerical problems
            Ta = self.TMIN + scalar(0.75) * (self.TMAX - self.TMIN)  # T in degC
            T24 = self.TMIN + scalar(0.5) * (self.TMAX - self.TMIN)  # T in degC
            pex = min(
                scalar(17.27) * (self.TMIN) / (scalar(237.3) + self.TMIN), scalar(10)
            )  # T in degC
            pe = min(
                scalar(610.8) * (exp(pex)), scalar(10000.0)
            )  # Mean actual vapour pressure, from dewpoint temperature
        # windspeed is at 1m
        # u2 = scalar(WindFactor)*self.WINDSPEED*(scalar(1)-(scalar(1)-self.fday)*scalar(0.25))/self.fday
        self.u1 = (
            self.WINDSPEED
            * (scalar(1) - (scalar(1) - self.fday) * scalar(0.25))
            / self.fday
        )
        pair = self.AIRPRESS  # already in Pa

        # diagnostic equations

        self.LAI = self.SLA * self.Mleaf  # (5.3)
        fveg = max(1 - exp(-self.LAI / self.LAIref), 0.000001)  # (5.3)

        # Vc = max(0,EVI-0.07)/fveg
        fsoil = 1 - fveg
        w0 = self.S0 / self.S0max  # (2.)
        ws = self.Ss / self.Ssmax  # (2.1)
        wd = self.Sd / self.Sdmax  # (2.1)

        TotSnow = self.FreeWater + self.DrySnow
        # wSnow = self.FreeWater/(TotSnow+1e-5)

        # Spatialise catchment fractions
        # Sgfree =   max(self.Sg,0.0)
        # JS: Not sure if this is translated properly....
        # for i=1:par.Nhru
        ChannelSurface = min(0, (0.007 * self.Sr ** 0.75))
        OpenWaterFrac = max(ChannelSurface, self.OpenWaterFrac)
        
        # !! HANDometric functions go here !!

        # !! HANDometric functions go here !!

        # fsat =   min(1.0,max(min(0.005,0.007*self.Sr**0.75),Sgfree/self.Sgref))
        # Sghru =   self.Sg

        # CALCULATION OF PET
        # Conversions and coefficients (3.1)
        pesx = min((scalar(17.27) * Ta / (scalar(237.3) + Ta)), scalar(10))
        pes = min(
            scalar((scalar(610.8)) * exp(pesx)), scalar(10000)
        )  # saturated vapour pressure
        # fRH = pe/pes  # relative air humidity                                  -------------- check
        cRE = 0.03449 + 4.27e-5 * Ta
        # Caero = self.fday*0.176*(1+Ta/209.1)*(pair-0.417*pe)*(1-fRH)         -------------- check
        # keps = 1.4e-3*((Ta/187)**2+Ta/107+1)*(6.36*pair+pe)/pes
        ga = max(0.001, self.ku1 * self.u1)

        if self.UseETPdata == 1:
            self.E0 = max(self.EPOT, 0)
            keps = (
                0.655E-3 * pair / pes
            )  # See Appendix A3 (http://www.clw.csiro.au/publications/waterforahealthycountry/2010/wfhc-aus-water-resources-assessment-system.pdf) --------------------------------   check!

        elif self.UseETPdata == 0:
            # Aerodynamic conductance (3.7)

            ns_alb = self.ALBEDO  # available as parameter climatology
            Rgeff = Rg / self.fday
            # shortwave radiation balance (3.2)
            # alb_veg = 0.452*Vc
            # alb_soil = alb_wet+(alb_dry-alb_wet)*exp(-w0/w0ref_alb)
            # new equations for snow albedo
            alb_snow = self.SNOWALBEDO  # available as parameter climatology
            fsnow = min(
                1.0, 0.05 * TotSnow
            )  # assumed; ideally some lit research needed
            # alb = fveg*alb_veg+(fsoil-fsnow)*alb_soil +fsnow*alb_snow
            # alb = albedo
            alb = (1 - fsnow) * ns_alb + fsnow * alb_snow
            RSn = (1 - alb) * Rgeff
            # long wave radiation balance (3.3 to 3.5)
            StefBolz = 5.67e-8
            Tkelv = Ta + 273.16
            #            self.RLin = (0.65*(pe/Tkelv)**0.14)*StefBolz*Tkelv**4     # (3.3)
            self.RLin - self.LWDOWN  # available from meteo forcing
            RLout = StefBolz * Tkelv ** 4.0  # (3.4)
            self.RLn = self.RLin - RLout

            self.fGR = self.Gfrac_max * (1 - exp(-fsoil / self.fvegref_G))
            self.Rneff = max(1, (RSn + self.RLn) * (1 - self.fGR))

            fRH = pe / pes  # relative air humidity
            Caero = (
                self.fday * 0.176 * (1 + Ta / 209.1) * (pair - 0.417 * pe) * (1 - fRH)
            )  # -------------- check
            keps = 1.4e-3 * ((Ta / 187) ** 2 + Ta / 107 + 1) * (6.36 * pair + pe) / pes

            #  Potential evaporation
            kalpha = 1 + Caero * ga / self.Rneff
            self.E0 = cRE * (1 / (1 + keps)) * kalpha * self.Rneff * self.fday
            self.E0 = max(self.E0, 0)
            self.Eeq = cRE * (1 / (1 + keps)) * 1 * self.Rneff * self.fday
            self.Ept = cRE * (1 / (1 + keps)) * 1.26 * self.Rneff * self.fday

        # CALCULATION OF ET FLUXES AND ROOT WATER UPTAKE
        # Root water uptake constraint (4.4)
        U0max = scalar(0)
        Usmax = max(
            0, self.Us0 * min(1, ws / self.wslimU)
        )  ##0-waarden omdat ws1 bevat 0-waarden (zie regel 116)
        Udmax = max(
            0, self.Ud0 * min(1, wd / self.wdlimU)
        )  ##0-waarden omdat wd1 bevat 0-waarden (zie regel 118)
        Ugmax = max(0, self.Ug0 * max(0, fUg - fsat))
        Umax = max(Usmax, max(Udmax, Ugmax))

        # Maximum transpiration (4.3)
        Gsmax = self.Gs_scalar*self.cGsmax*self.Vc
        VPD = max(0,pes-pe)
        fD = self.Cg/(1+VPD/self.D50)
        gs = fveg*fD*Gsmax
        ft = 1/(1+(keps/(1+keps))*ga/gs)
        Etmax = ft*self.E0

        # Actual transpiration (4.1)
        Et = min(Umax, Etmax)

        # # Root water uptake distribution (2.3)
        U0 = max(min((U0max / (U0max + Usmax + Udmax + Ugmax)) * Et, self.S0 - 1e-2), 0)
        Us = max(min((Usmax / (U0max + Usmax + Udmax + Ugmax)) * Et, self.Ss - 1e-2), 0)
        Ud = max(min((Udmax / (U0max + Usmax + Udmax + Ugmax)) * Et, self.Sd - 1e-2), 0)
        Ug = max(min((Ugmax / (U0max + Usmax + Udmax + Ugmax)) * Et, self.Sd - 1e-2), 0)
        Et = U0 + Us + Ud + Ug  # to ensure mass balance

        # Soil evaporation (4.5)
        w0x = max(0, (self.S0 - U0) / self.S0max)  # (2.1)
        fsoilE = self.FsoilEmax * min(1, w0x / self.w0limE)
        Es0 = (1 - fsat) * fsoilE * (max(0, self.Eeq - Et))
        # Groundwater evaporation (4.6)
        Eg0 = max(0, fsat - fwater) * self.FsoilEmax * max(0, self.Eeq - Et)
        Es = Es0 + Eg0
        # Open water evaporation (4.7)
        Erl = fw_local*self.FwaterE*self.Ept
        Err = (fwater-fw_local)*self.FwaterE*self.Ept
        Er = Erl + Err

        # Rainfall interception evaporation (4.2)
        Sveg = self.S_sls * self.LAI
        fER = fveg * self.ER_coeff * max(0.05, self.hveg) ** self.ER_exp
        Pwet = max(
            0,
            (
                scalar(Sveg > 0 & fER > 0 & (fER / fveg) < 1)
                * -ln(1 - fER / fveg)
                * Sveg
                / fER
            ),
        )
        Ei = scalar(self.T24 > 0) * (
            scalar(Pg < Pwet) * fveg * Pg
            + scalar(Pg >= Pwet) * (fveg * Pwet + fER * (Pg - Pwet))
        )

        Edry = Et + Es + Er
        self.EACT = Edry + Ei  # for output only

        # HBV snow routine
        # Matlab: function [FreeWater,DrySnow,InSoil]=snow_submodel(Precipitation,Temperature,FreeWater,DrySnow)
        # derived from HBV-96 shared by Jaap Schellekens (Deltares) in May 2011
        # original in PCraster, adapted to Matlab by Albert van Dijk
        # HBV snow routine
        Pn = Pg - Ei

        # Snow routine parameters
        # parameters

        # Partitioning into fractions rain and snow
        Temperature = T24  # Dimmie, let op: tijdelijke regel!!
        RainFrac = max(0, min((Temperature - (TT - TTI / 2)) / TTI, 1))
        SnowFrac = 1 - RainFrac  # fraction of precipitation which falls as snow

        # Snowfall/melt calculations
        SnowFall = SnowFrac * Pn  # snowfall depth
        RainFall = RainFrac * Pn  # rainfall depth
        PotSnowMelt = Cfmax * max(
            0, Temperature - TT
        )  # Potential snow melt, based on temperature
        PotRefreezing = (
            Cfmax * CFR * max(TT - Temperature, 0)
        )  # Potential refreezing, based on temperature
        Refreezing = min(PotRefreezing, self.FreeWater)  # actual refreezing
        SnowMelt = min(PotSnowMelt, self.DrySnow)  # actual snow melt
        self.DrySnow = (
            self.DrySnow + SnowFall + Refreezing - SnowMelt
        )  # dry snow content
        self.FreeWater = self.FreeWater - Refreezing  # free water content in snow
        MaxFreeWater = self.DrySnow * WHC
        self.FreeWater = self.FreeWater + SnowMelt + RainFall
        InSoil = max(
            self.FreeWater - MaxFreeWater, 0
        )  # abundant water in snow pack which goes into soil
        self.FreeWater = self.FreeWater - InSoil
        # End of Snow Module
        Rmelt = scalar(Temperature < 0) * InSoil  # runs off if soil still frozen
        Ps = scalar(Temperature >= 0) * InSoil

        # CALCULATION OF WATER BALANCES
        # surface water fluxes (2.2)
        Rsof = fsat * Ps
        Pi = max(0, Ps - self.InitLoss)
        Rhof_soil = max(0, 1 - fsat - self.fImp) * (
            Pi - self.Pref * pcr_tanh(Pi / self.Pref)
        )  # CHECK IF THIS GOES OK IN PYTHON
        Rhof_imp = self.fImp * (
            Pi - self.Pref_imp * pcr_tanh(Pi / self.Pref_imp)
        )  # CHECK IF THIS GOES OK IN PYTHON
        Rhof = Rhof_soil + Rhof_imp
        QR = Rhof + Rsof + Rmelt
        I = Ps - Rhof - Rsof

        # SOIL WATER BALANCES (2.1 & 2.4)
        # Topsoil water balance (S0)

        # top soil layer hydrology
        Kr_0s = self.K0sat / self.Kssat
        Rh_0s = pcr_tanh(self.slope_coeff * self.slope * w0) * pcr_tanh(
            self.Kr_coeff * (Kr_0s - 1) * w0
        )
        # general default case
        Km = (self.K0sat * self.Kssat) ** 0.5
        A = Km / (self.S0max ** 2)
        B = 1
        C = -(self.S0 + I + Es)
        S0 = (-B + ((B ** 2 - 4 * A * C) ** 0.5)) / (2 * A)
        D0 = (1 - Rh_0s) * Km * ((self.S0 / self.S0max) ** 2)
        IF0 = Rh_0s * Km * ((self.S0 / self.S0max) ** 2)
        # depletion case
        imap = (self.S0 + I) <= Es
        Es = ifthenelse(imap, (self.S0 + I), Es)
        S0 = ifthenelse(imap, 0, S0)
        D0 = ifthenelse(imap, 0, D0)
        IF0 = ifthenelse(imap, 0, IF0)
        # saturation case
        imap = (self.S0max - self.S0 - self.K0sat) <= (I - Es)
        D0 = ifthenelse(imap, (1 - Rh_0s) * self.K0sat, D0)
        IF0 = ifthenelse(
            imap, Rh_0s * self.K0sat + (self.S0 - self.S0max - self.K0sat + I - Es), IF0
        )
        S0 = ifthenselse(imap, self.S0max, S0)
        # enforce mass balance (for numerical & rounding errors
        S0 = max(0, min(S0, self.S0max))
        massbal = self.S0 + I - Es - D0 - IF0 - S0
        D0 = D0 + (1 - Rh_0s) * massbal
        IF0 = IF0 + Rh_0s * massbal
        self.S0 = S0  # Update state

        # shallow soil layer hydrology
        Kr_sd = self.Kssat / self.Kdsat
        Rh_sd = pcr_tanh(self.slope_coeff * self.slope * ws) * pcr_tanh(
            self.Kr_coeff * (Kr_sd - 1) * ws
        )
        # general default case
        Km = (self.Kssat * self.Kdsat) ** 0.5
        A = Km / (self.Ssmax ** 2)
        B = 1
        C = -(self.Ss + D0 - Us)
        Ss = (-B + ((B ** 2 - 4 * A * C) ** 0.5)) / (2 * A)
        Ds = (1 - Rh_sd) * Km * ((Ss / self.Ssmax) ** 2)
        IFs = Rh_sd * Km * ((Ss / self.Ssmax) ** 2)
        # depletion case
        imap = (self.Ss + D0) <= Us
        Us = ifthenelse(imap, (self.Ss + D0), Us)
        Ss = ifthenelse(imap, 0, Ss)
        Ds = ifthenelse(imap, 0, Ds)
        IFs = ifthenelse(imap, 0, IFs)
        # saturation case
        imap = (self.Ssmax - self.Ss - self.Kssat) <= (D0 - Us)
        Ds = ifthenelse(imap, (1 - Rh_sd) * self.Kssat, Ds)
        IFs = ifthenelse(
            imap,
            Rh_sd * self.Kssat + (self.Ss - self.Ssmax - self.Kssat + D0 - Us),
            IFs,
        )
        Ss = ifthenselse(imap, self.Ssmax, Ss)
        # enforce mass balance (for numerical & rounding errors
        Ss = max(0, min(Ss, self.Ssmax))
        massbal = self.Ss + D0 - Us - Ds - IFs - Ss
        Ds = Ds + (1 - Rh_sd) * massbal
        IFs = IFs + Rh_sd * massbal
        self.Ss = Ss  # Update state

        # Deep soil layer hydrology
        A = self.Kdsat / (self.Sdmax ** 2)
        B = 1
        C = -(self.Sd + Ds - Ud)
        Sd = (-B + ((B ** 2 - 4 * A * C) ** 0.5)) / (2 * A)
        Dd = self.Kdsat * ((Sd / self.Sdmax) ** 2)
        IFd = 0 * Dd
        # depletion case
        imap = (self.Sd + Ds) <= Ud
        Ud = ifthenelse(imap, (self.Sd + Ds), Ud)
        Sd = ifthenelse(imap, 0, Sd)
        Dd = ifthenelse(imap, 0, Dd)
        IFd = ifthenelse(imap, 0, IFd)
        # saturation case
        imap   = (self.Sdmax-self.Sd-self.Kdsat)<=(Ds-Ud)
        Dd = ifthenelse(imap,self.Kdsat,Dd)
        IFd = ifthenelse(imap,self.Sd-self.Sdmax-self.Kdsat+Ds-Ud,IFd)
        Sd = ifthenselse(imap,self.Sdmax,Sd)
        # enforce mass balance (for numerical & rounding errors
        Sd = max(0, min(Sd, self.Sdmax))
        massbal = self.Sd + Ds - Ud - Dd - IFd - Sd
        Dd = Dd + massbal
        self.Sd = Sd  # Update state

        IFs = IFs + IFd
        # add up to interflow
        QR = QR + IF0 + IFs

        # CATCHMENT WATER BALANCE
        # Groundwater store water balance (Sg) (2.5)
        NetGf = Dd - Eg - Ug
        self.Sg = self.Sg + NetGf
        Sgfree = max(self.Sg, 0)
        Qg = min(Sgfree, (1 - exp(-self.K_gw)) * Sgfree)
        self.Sg = self.Sg - Qg

        # Surface water store water balance (Sr) (2.7)
        self.Sr = max(0, self.Sr + QR - Er + Qg)
        self.Qtot = max(0, min(self.Sr, (1 - exp(-self.K_rout)) * self.Sr))
        self.Sr = self.Sr - self.Qtot

        # VEGETATION ADJUSTMENT (5)
        fvmax = 1 - exp(-max(self.LAImax, 0.002778) / self.LAIref)
        fveq = (
            (1 / max((self.E0 / Umax) - 1, 1e-3))
            * (keps / (1 + keps))
            * (ga / (fd * Gsmax))
        )
        fveq = min(fveq, fvmax)
        dMleaf = -ln(1 - fveq) * self.LAIref / self.SLA - self.Mleaf

        # Mleafnet1 = dMleaf1 * (dMleaf1/self.Tgrow1) + dMleaf1 * dMleaf1/self.Tsenc1
        # Mleafnet2 = dMleaf2 * (dMleaf1/self.Tgrow2) + dMleaf2 * dMleaf2/self.Tsenc2
        Mleafnet = (
            scalar(dMleaf > 0) * (dMleaf / self.Tgrow)
            + scalar(dMleaf < 0) * dMleaf / self.Tsenc
        )

        self.Mleaf = self.Mleaf + Mleafnet
        self.LAI = self.SLA * self.Mleaf  # (5.3)

        fveg = 1 - exp(-self.LAI1 / self.LAIref1)  # (5.3)
        # in case this is desired as output:
        self.w0 = self.S0 / self.S0max
        self.TotSnow = self.DrySnow + self.FreeWater


# The main function is used to run the program from the command line


def main(argv=None):
    """
    *Optional*
    
    Perform command line execution of the model. This example uses the getopt
    module to parse the command line options.
    
    The user can set the caseName, the runDir, the timestep and the configfile.
    """
    global multpars
    caseName = (
        "../openstreams_w3ra"
    )  # "D:/trambaue/_Projects/GLOFFIS/201501/GLOFFIS_SA/Modules/openstreams_w3ra/"
    runId = "run_default"
    configfile = "wflow_W3RA.ini"
    _lastTimeStep = 15
    _firstTimeStep = 0
    timestepsecs = 86400
    fewsrun = False
    wflow_cloneMap = "wflow_subcatch.map"
    runinfoFile = "runinfo.xml"
    _NoOverWrite = False
    loglevel = logging.DEBUG
    LogFileName = "wflow.log"

    # This allows us to use the model both on the command line and to call
    # the model usinge main function from another python script.

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    opts, args = getopt.getopt(argv, "C:S:T:c:s:R:F:")

    for o, a in opts:
        if o == "-F":
            runinfoFile = a
            fewsrun = True
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

    if len(opts) <= 1:
        usage()

    if fewsrun:
        ts = getTimeStepsfromRuninfo(runinfoFile, timestepsecs)
        starttime = getStartTimefromRuninfo(runinfoFile)
        if ts:
            _lastTimeStep = ts
            _firstTimeStep = 1
        else:
            print("Failed to get timesteps from runinfo file: " + runinfoFile)
            exit(2)
    else:
        starttime = dt.datetime(1990,0o1,0o1)

    if _lastTimeStep < _firstTimeStep:
        print("The starttimestep (" + str(_firstTimeStep) + ") is smaller than the last timestep (" + str(
            _lastTimeStep) + ")")
        usage()

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(
        myModel, _lastTimeStep, firstTimestep=_firstTimeStep, datetimestart=starttime
    )
    dynModelFw.createRunId(
        NoOverWrite=_NoOverWrite,
        level=loglevel,
        logfname=LogFileName,
        model="wflow_W3RA",
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
            configset(myModel.config, "model", "reinit", "1", overwrite=True)
        if o == "-i":
            configset(myModel.config, "model", "intbl", a, overwrite=True)
        if o == "-s":
            configset(myModel.config, "model", "timestepsecs", a, overwrite=True)

    dynModelFw.setupFramework()

    dynModelFw._runInitial()
    dynModelFw._runResume()
    dynModelFw._runDynamic(0, 0)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
