#!/usr/bin/python

"""
Definition of the wflow_W3RA model.
---------------------------------------
The model is modified from the Australian Water Resources Assessment
Landscape (AWRA-L) model version 0.5
W3RA is documented in
van Dijk et al. (2013), Water Resour. Res., 49, 2729-2746, doi:10.1002/wrcr.20251
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

    -T: Set end time of the run: yyyy-mm-dd hh:mm:ss

    -S: Set start time of the run: yyyy-mm-dd hh:mm:ss
    
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
    for msg in args:
        print msg
    print __doc__
    sys.exit(0)


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
            "S01",
            "Ss1",
            "Sd1",
            "Mleaf1",
            "FreeWater1",
            "DrySnow1",
            "LAI1",
            "EVI1",
            "Sg",
            "Sr",
            "S02",
            "Ss2",
            "Sd2",
            "Mleaf2",
            "FreeWater2",
            "DrySnow2",
            "LAI2",
            "EVI2",
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

        # SET GLBOAL PARAMETER VALUES (however not used in original script)
        # Nhru=2
        # K_gw_scale=0.0146
        # K_gw_shape=0.0709
        # K_rout_scale=0.1943
        # K_rout_int=0.0589
        # FdrainFC_scale=0.2909
        # FdrainFC_shape=0.5154
        # Sgref_scale=3.2220
        # Sgref_shape=3.2860
        # fday=0.5000
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

        self.latitude = ycoordinate(boolean(self.Altitude))

        # Add reading of parameters here

        self.K_gw = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/k_gw.map"), 0.0, fail=True
        )
        self.K_rout = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/k_rout.map"), 0.0, fail=True
        )
        self.Sgref = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/sgref.map"), 0.0, fail=True
        )
        self.alb_dry1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/alb_dry.map"), 0.0, fail=True
        )
        self.alb_wet1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/alb_wet.map"), 0.0, fail=True
        )
        self.beta1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/beta.map"), 0.0, fail=True
        )
        self.cGsmax1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/cgsmax.map"), 0.0, fail=True
        )
        self.ER_frac_ref1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/er_frac_ref.map"), 0.0, fail=True
        )
        self.FdrainFC1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fdrainfc.map"), 0.0, fail=True
        )
        self.Fgw_conn1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fgw_conn.map"), 0.0, fail=True
        )
        self.Fhru1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fhru.map"), 0.0, fail=True
        )
        self.SLA1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/sla.map"), 0.0, fail=True
        )
        self.LAIref1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/lairef.map"), 0.0, fail=True
        )
        self.FsoilEmax1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fsoilemax.map"), 0.0, fail=True
        )
        self.fvegref_G1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fvegref_g.map"), 0.0, fail=True
        )
        self.FwaterE1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fwatere.map"), 0.0, fail=True
        )
        self.Gfrac_max1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/gfrac_max.map"), 0.0, fail=True
        )
        self.hveg1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hveg.map"), 0.0, fail=True
        )
        self.InitLoss1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/initloss.map"), 0.0, fail=True
        )
        self.LAImax1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/laimax.map"), 0.0, fail=True
        )
        self.PrefR1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/prefr.map"), 0.0, fail=True
        )
        self.S_sls1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/s_sls.map"), 0.0, fail=True
        )
        self.S0FC1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/s0fc.map"), 0.0, fail=True
        )
        self.SsFC1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/ssfc.map"), 0.0, fail=True
        )
        self.SdFC1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/sdfc.map"), 0.0, fail=True
        )
        self.Vc1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/vc.map"), 0.0, fail=True
        )
        self.w0ref_alb1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/w0ref_alb.map"), 0.0, fail=True
        )
        self.Us01 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/us0.map"), 0.0, fail=True
        )
        self.Ud01 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/ud0.map"), 0.0, fail=True
        )
        self.wslimU1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/wslimu.map"), 0.0, fail=True
        )
        self.wdlimU1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/wdlimu.map"), 0.0, fail=True
        )
        self.w0limE1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/w0lime.map"), 0.0, fail=True
        )
        self.Tgrow1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/tgrow.map"), 0.0, fail=True
        )
        self.Tsenc1 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/tsenc.map"), 0.0, fail=True
        )

        self.alb_dry2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/alb_dry2.map"), 0.0, fail=True
        )
        self.alb_wet2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/alb_wet2.map"), 0.0, fail=True
        )
        self.beta2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/beta2.map"), 0.0, fail=True
        )
        self.cGsmax2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/cgsmax2.map"), 0.0, fail=True
        )
        self.ER_frac_ref2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/er_frac_ref2.map"), 0.0, fail=True
        )
        self.FdrainFC2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fdrainfc2.map"), 0.0, fail=True
        )
        self.Fgw_conn2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fgw_conn2.map"), 0.0, fail=True
        )
        self.Fhru2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fhru2.map"), 0.0, fail=True
        )
        self.SLA2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/sla2.map"), 0.0, fail=True
        )
        self.LAIref2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/lairef2.map"), 0.0, fail=True
        )
        self.FsoilEmax2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fsoilemax2.map"), 0.0, fail=True
        )
        self.fvegref_G2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fvegref_g2.map"), 0.0, fail=True
        )
        self.FwaterE2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/fwatere2.map"), 0.0, fail=True
        )
        self.Gfrac_max2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/gfrac_max2.map"), 0.0, fail=True
        )
        self.hveg2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/hveg2.map"), 0.0, fail=True
        )
        self.InitLoss2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/initloss2.map"), 0.0, fail=True
        )
        self.LAImax2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/laimax2.map"), 0.0, fail=True
        )
        self.PrefR2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/prefr2.map"), 0.0, fail=True
        )
        self.S_sls2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/s_sls2.map"), 0.0, fail=True
        )
        self.S0FC2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/s0fc2.map"), 0.0, fail=True
        )
        self.SsFC2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/ssfc2.map"), 0.0, fail=True
        )
        self.SdFC2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/sdfc2.map"), 0.0, fail=True
        )
        self.Vc2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/vc2.map"), 0.0, fail=True
        )
        self.w0ref_alb2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/w0ref_alb2.map"), 0.0, fail=True
        )
        self.Us02 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/us02.map"), 0.0, fail=True
        )
        self.Ud02 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/ud02.map"), 0.0, fail=True
        )
        self.wslimU2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/wslimu2.map"), 0.0, fail=True
        )
        self.wdlimU2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/wdlimu2.map"), 0.0, fail=True
        )
        self.w0limE2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/w0lime2.map"), 0.0, fail=True
        )
        self.Tgrow2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/tgrow2.map"), 0.0, fail=True
        )
        self.Tsenc2 = self.wf_readmap(
            os.path.join(self.Dir, "staticmaps/tsenc2.map"), 0.0, fail=True
        )

        self.wf_multparameters()
        # Static, for the computation of Aerodynamic conductance (3.7)
        self.fh1 = ln(813. / self.hveg1 - 5.45)
        self.fh2 = ln(813. / self.hveg2 - 5.45)
        self.ku2_1 = 0.305 / (self.fh1 * (self.fh1 + 2.3))
        self.ku2_2 = 0.305 / (self.fh2 * (self.fh2 + 2.3))

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
                exec "self." + s + " = cover(1.0)"

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
            self.WINDSPEED = cover(
                self.wf_readmap(self.WINDSPEED_mapstack, default=1.0), scalar(1.0)
            )
            self.AIRPRESS = cover(
                self.wf_readmap(self.AIRPRESS_mapstack, default=980.0), scalar(980.0)
            )
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
        # rescale factor because windspeed climatology is at 2m
        WindFactor = 1.0
        # u2 = scalar(WindFactor)*self.WINDSPEED*(scalar(1)-(scalar(1)-self.fday)*scalar(0.25))/self.fday
        self.u2 = (
            scalar(WindFactor)
            * self.WINDSPEED
            * (scalar(1) - (scalar(1) - self.fday) * scalar(0.25))
            / self.fday
        )
        pair = self.AIRPRESS  # already in Pa

        # diagnostic equations

        self.LAI1 = self.SLA1 * self.Mleaf1  # (5.3)
        self.LAI2 = self.SLA2 * self.Mleaf2  # (5.3)
        fveg1 = max(1 - exp(-self.LAI1 / self.LAIref1), 0.000001)  # (5.3)
        fveg2 = max(1 - exp(-self.LAI2 / self.LAIref2), 0.000001)

        # Vc = max(0,EVI-0.07)/fveg
        fsoil1 = 1 - fveg1
        fsoil2 = 1 - fveg2
        w01 = self.S01 / self.S0FC1  # (2.1)
        w02 = self.S02 / self.S0FC2
        ws1 = self.Ss1 / self.SsFC1  # (2.1)
        ws2 = self.Ss2 / self.SsFC2
        wd1 = self.Sd1 / self.SdFC1  # (2.1)
        wd2 = self.Sd2 / self.SdFC2  # (2.1)

        TotSnow1 = self.FreeWater1 + self.DrySnow1
        TotSnow2 = self.FreeWater2 + self.DrySnow2
        wSnow1 = self.FreeWater1 / (TotSnow1 + 1e-5)
        wSnow2 = self.FreeWater2 / (TotSnow2 + 1e-5)

        # Spatialise catchment fractions
        Sgfree = max(self.Sg, 0.0)
        # JS: Not sure if this is translated properly....
        # for i=1:par.Nhru
        fwater1 = min(0.005, (0.007 * self.Sr ** 0.75))
        fwater2 = min(0.005, (0.007 * self.Sr ** 0.75))
        fsat1 = min(1.0, max(min(0.005, 0.007 * self.Sr ** 0.75), Sgfree / self.Sgref))
        fsat2 = min(1.0, max(min(0.005, 0.007 * self.Sr ** 0.75), Sgfree / self.Sgref))
        Sghru1 = self.Sg
        Sghru2 = self.Sg

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
        ga1 = self.ku2_1 * self.u2
        ga2 = self.ku2_2 * self.u2

        if self.UseETPdata == 1:
            self.E01 = max(self.EPOT, 0)
            self.E02 = max(self.EPOT, 0)
            keps = (
                0.655E-3 * pair / pes
            )  # See Appendix A3 (http://www.clw.csiro.au/publications/waterforahealthycountry/2010/wfhc-aus-water-resources-assessment-system.pdf) --------------------------------   check!

        elif self.UseETPdata == 0:
            # Aerodynamic conductance (3.7)

            ns_alb = self.ALBEDO
            Rgeff = Rg / self.fday
            # shortwave radiation balance (3.2)
            # alb_veg = 0.452*Vc
            # alb_soil = alb_wet+(alb_dry-alb_wet)*exp(-w0/w0ref_alb)
            # new equations for snow albedo
            alb_snow1 = 0.65 - 0.2 * wSnow1  # assumed; ideally some lit research needed
            alb_snow2 = 0.65 - 0.2 * wSnow2
            fsnow1 = min(
                1.0, 0.05 * TotSnow1
            )  # assumed; ideally some lit research needed
            fsnow2 = min(1.0, 0.05 * TotSnow2)
            # alb = fveg*alb_veg+(fsoil-fsnow)*alb_soil +fsnow*alb_snow
            # alb = albedo
            alb1 = (1 - fsnow1) * ns_alb + fsnow1 * alb_snow1
            alb2 = (1 - fsnow2) * ns_alb + fsnow2 * alb_snow2
            RSn1 = (1 - alb1) * Rgeff
            RSn2 = (1 - alb2) * Rgeff
            # long wave radiation balance (3.3 to 3.5)
            StefBolz = 5.67e-8
            Tkelv = Ta + 273.16
            self.RLin = (0.65 * (pe / Tkelv) ** 0.14) * StefBolz * Tkelv ** 4  # (3.3)
            RLout = StefBolz * Tkelv ** 4.0  # (3.4)
            self.RLn = self.RLin - RLout

            self.fGR1 = self.Gfrac_max1 * (1 - exp(-fsoil1 / self.fvegref_G1))
            self.fGR2 = self.Gfrac_max2 * (1 - exp(-fsoil2 / self.fvegref_G2))  # (3.5)
            self.Rneff1 = (RSn1 + self.RLn) * (1 - self.fGR1)
            self.Rneff2 = (RSn2 + self.RLn) * (1 - self.fGR2)

            fRH = pe / pes  # relative air humidity
            Caero = (
                self.fday * 0.176 * (1 + Ta / 209.1) * (pair - 0.417 * pe) * (1 - fRH)
            )  # -------------- check
            keps = 1.4e-3 * ((Ta / 187) ** 2 + Ta / 107 + 1) * (6.36 * pair + pe) / pes

            #  Potential evaporation
            kalpha1 = 1 + Caero * ga1 / self.Rneff1
            kalpha2 = 1 + Caero * ga2 / self.Rneff2
            self.E01 = cRE * (1 / (1 + keps)) * kalpha1 * self.Rneff1 * self.fday
            self.E02 = cRE * (1 / (1 + keps)) * kalpha2 * self.Rneff2 * self.fday
            self.E01 = max(self.E01, 0)
            self.E02 = max(self.E02, 0)

        # CALCULATION OF ET FLUXES AND ROOT WATER UPTAKE
        # Root water uptake constraint (4.4)
        Usmax1 = max(
            0, self.Us01 * min(1, ws1 / self.wslimU1)
        )  ##0-waarden omdat ws1 bevat 0-waarden (zie regel 116)
        Usmax2 = max(
            0, self.Us02 * min(1, ws2 / self.wslimU2)
        )  ##0-waarden omdat ws2 bevat 0-waarden (zie regel 117)
        Udmax1 = max(
            0, self.Ud01 * min(1, wd1 / self.wdlimU1)
        )  ##0-waarden omdat wd1 bevat 0-waarden (zie regel 118)
        Udmax2 = max(
            0, self.Ud02 * min(1, wd2 / self.wdlimU2)
        )  ##0-waarden omdat wd2 bevat 0-waarden (zie regel 119)
        # U0max = max(0, Us0*min(1,w0/wslimU))
        U0max1 = scalar(0)
        U0max2 = scalar(0)
        Utot1 = max(Usmax1, max(Udmax1, U0max1))
        Utot2 = max(Usmax2, max(Udmax2, U0max2))

        # Maximum transpiration (4.3)
        Gsmax1 = self.cGsmax1 * self.Vc1
        gs1 = fveg1 * Gsmax1
        ft1 = 1 / (1 + (keps / (1 + keps)) * ga1 / gs1)
        Etmax1 = ft1 * self.E01
        Gsmax2 = self.cGsmax2 * self.Vc2
        gs2 = fveg2 * Gsmax2
        ft2 = 1 / (1 + (keps / (1 + keps)) * ga2 / gs2)
        Etmax2 = ft2 * self.E02

        # Actual transpiration (4.1)
        Et1 = min(Utot1, Etmax1)
        Et2 = min(Utot2, Etmax2)

        # # Root water uptake distribution (2.3)
        U01 = max(min((U0max1 / (U0max1 + Usmax1 + Udmax1)) * Et1, self.S01 - 1e-2), 0)
        Us1 = max(min((Usmax1 / (U0max1 + Usmax1 + Udmax1)) * Et1, self.Ss1 - 1e-2), 0)
        Ud1 = max(min((Udmax1 / (U0max1 + Usmax1 + Udmax1)) * Et1, self.Sd1 - 1e-2), 0)
        Et1 = U01 + Us1 + Ud1  # to ensure mass balance

        U02 = max(min((U0max2 / (U0max2 + Usmax2 + Udmax2)) * Et2, self.S02 - 1e-2), 0)
        Us2 = max(min((Usmax2 / (U0max2 + Usmax2 + Udmax2)) * Et2, self.Ss2 - 1e-2), 0)
        Ud2 = max(min((Udmax2 / (U0max2 + Usmax2 + Udmax2)) * Et2, self.Sd2 - 1e-2), 0)
        Et2 = U02 + Us2 + Ud2

        # Soil evaporation (4.5)
        self.S01 = max(0, self.S01 - U01)
        self.S02 = max(0, self.S02 - U02)
        w01 = self.S01 / self.S0FC1  # (2.1)
        w02 = self.S02 / self.S0FC2  # (2.1)
        fsoilE1 = self.FsoilEmax1 * min(1, w01 / self.w0limE1)
        fsoilE2 = self.FsoilEmax2 * min(1, w02 / self.w0limE2)
        Es1 = max(0, min(((1 - fsat1) * fsoilE1 * (self.E01 - Et1)), self.S01 - 1e-2))
        Es2 = max(0, min(((1 - fsat2) * fsoilE2 * (self.E02 - Et2)), self.S02 - 1e-2))
        # Groundwater evaporation (4.6)
        Eg1 = min((fsat1 - fwater1) * self.FsoilEmax1 * (self.E01 - Et1), Sghru1)
        Eg2 = min((fsat2 - fwater2) * self.FsoilEmax2 * (self.E02 - Et2), Sghru2)
        # Open water evaporation (4.7)
        Er1 = min(fwater1 * self.FwaterE1 * max(0, self.E01 - Et1), self.Sr)
        Er2 = min(fwater2 * self.FwaterE2 * max(0, self.E02 - Et2), self.Sr)
        # Rainfall interception evaporation (4.2)
        Sveg1 = self.S_sls1 * self.LAI1
        fER1 = self.ER_frac_ref1 * fveg1
        Pwet1 = -ln(1 - fER1 / fveg1) * Sveg1 / fER1
        Ei1 = scalar(Pg < Pwet1) * fveg1 * Pg + scalar(Pg >= Pwet1) * (
            fveg1 * Pwet1 + fER1 * (Pg - Pwet1)
        )

        Sveg2 = self.S_sls2 * self.LAI2
        fER2 = self.ER_frac_ref2 * fveg2
        Pwet2 = -ln(1 - fER2 / fveg2) * Sveg2 / fER2
        Ei2 = scalar(Pg < Pwet2) * fveg2 * Pg + scalar(Pg >= Pwet2) * (
            fveg2 * Pwet2 + fER2 * (Pg - Pwet2)
        )

        self.EACT1 = (Et1 + Es1 + Eg1 + Er1 + Ei1) * self.Fhru1
        self.EACT2 = (Et2 + Es2 + Eg2 + Er2 + Ei2) * self.Fhru2
        self.EACT = self.EACT1 + self.EACT2

        # HBV snow routine
        # Matlab: function [FreeWater,DrySnow,InSoil]=snow_submodel(Precipitation,Temperature,FreeWater,DrySnow)
        # derived from HBV-96 shared by Jaap Schellekens (Deltares) in May 2011
        # original in PCraster, adapted to Matlab by Albert van Dijk
        # HBV snow routine
        Pn1 = Pg - Ei1
        Pn2 = Pg - Ei2
        Precipitation1 = Pn1
        Precipitation2 = Pn2

        # Snow routine parameters
        # parameters
        # TODO: Check this, not sure if this works.......
        x = scalar(Pg)
        Cfmax1 = 0.6 * 3.75653 * scalar(x >= 0)
        Cfmax2 = 3.75653 * scalar(x >= 0)
        TT1 = -1.41934 * scalar(
            x >= 0
        )  # critical temperature for snowmelt and refreezing
        TT2 = -1.41934 * scalar(x >= 0)
        TTI1 = 1.00000 * scalar(
            x >= 0
        )  # defines interval in which precipitation falls as rainfall and snowfall
        TTI2 = 1.00000 * scalar(x >= 0)
        CFR1 = 0.05000 * scalar(
            x >= 0
        )  # refreezing efficiency constant in refreezing of freewater in snow
        CFR2 = 0.05000 * scalar(x >= 0)
        WHC1 = 0.10000 * scalar(x >= 0)
        WHC2 = 0.10000 * scalar(x >= 0)

        # Partitioning into fractions rain and snow
        Temperature = T24  # Dimmie, let op: tijdelijke regel!!
        RainFrac1 = max(0, min((Temperature - (TT1 - TTI1 / 2)) / TTI1, 1))
        RainFrac2 = max(0, min((Temperature - (TT2 - TTI2 / 2)) / TTI2, 1))
        SnowFrac1 = 1 - RainFrac1  # fraction of precipitation which falls as snow
        SnowFrac2 = 1 - RainFrac2

        # Snowfall/melt calculations
        SnowFall1 = SnowFrac1 * Precipitation1  # snowfall depth
        SnowFall2 = SnowFrac2 * Precipitation2
        RainFall1 = RainFrac1 * Precipitation1  # rainfall depth
        RainFall2 = RainFrac2 * Precipitation2
        PotSnowMelt1 = Cfmax1 * max(
            0, Temperature - TT1
        )  # Potential snow melt, based on temperature
        PotSnowMelt2 = Cfmax2 * max(0, Temperature - TT2)
        PotRefreezing1 = (
            Cfmax1 * CFR1 * max(TT1 - Temperature, 0)
        )  # Potential refreezing, based on temperature
        PotRefreezing2 = Cfmax2 * CFR2 * max(TT2 - Temperature, 0)
        Refreezing1 = min(PotRefreezing1, self.FreeWater1)  # actual refreezing
        Refreezing2 = min(PotRefreezing2, self.FreeWater2)
        SnowMelt1 = min(PotSnowMelt1, self.DrySnow1)  # actual snow melt
        SnowMelt2 = min(PotSnowMelt2, self.DrySnow2)
        self.DrySnow1 = (
            self.DrySnow1 + SnowFall1 + Refreezing1 - SnowMelt1
        )  # dry snow content
        self.DrySnow2 = self.DrySnow2 + SnowFall2 + Refreezing2 - SnowMelt2
        self.FreeWater1 = self.FreeWater1 - Refreezing1  # free water content in snow
        self.FreeWater2 = self.FreeWater2 - Refreezing2
        MaxFreeWater1 = self.DrySnow1 * WHC1
        MaxFreeWater2 = self.DrySnow2 * WHC2
        self.FreeWater1 = self.FreeWater1 + SnowMelt1 + RainFall1
        self.FreeWater2 = self.FreeWater2 + SnowMelt2 + RainFall2
        InSoil1 = max(
            self.FreeWater1 - MaxFreeWater1, 0
        )  # abundant water in snow pack which goes into soil
        InSoil2 = max(self.FreeWater2 - MaxFreeWater2, 0)
        self.FreeWater1 = self.FreeWater1 - InSoil1
        self.FreeWater2 = self.FreeWater2 - InSoil2
        # End of Snow Module

        # CALCULATION OF WATER BALANCES
        # surface water fluxes (2.2)
        NetInSoil1 = max(0, (InSoil1 - self.InitLoss1))
        NetInSoil2 = max(0, (InSoil2 - self.InitLoss2))
        Rhof1 = (1 - fsat1) * (NetInSoil1 / (NetInSoil1 + self.PrefR1)) * NetInSoil1
        Rhof2 = (1 - fsat2) * (NetInSoil2 / (NetInSoil2 + self.PrefR2)) * NetInSoil2
        Rsof1 = fsat1 * NetInSoil1
        Rsof2 = fsat2 * NetInSoil2
        QR1 = Rhof1 + Rsof1
        QR2 = Rhof2 + Rsof2
        I1 = InSoil1 - QR1
        I2 = InSoil2 - QR2
        # SOIL WATER BALANCES (2.1 & 2.4)
        # Topsoil water balance (S0)
        self.S01 = self.S01 + I1 - Es1 - U01
        self.S02 = self.S02 + I2 - Es2 - U02
        SzFC1 = self.S0FC1
        SzFC2 = self.S0FC2
        Sz1 = self.S01
        Sz2 = self.S02
        wz1 = max(1e-2, Sz1) / SzFC1
        wz2 = max(1e-2, Sz2) / SzFC2
        self.TMP = SzFC1

        # TODO: Check if this works
        fD1 = scalar(wz1 > 1) * max(self.FdrainFC1, 1 - 1 / wz1) + scalar(
            wz1 <= 1
        ) * self.FdrainFC1 * exp(self.beta1 * scalar(wz1 - 1))
        fD2 = scalar(wz2 > 1) * max(self.FdrainFC2, 1 - 1 / wz2) + scalar(
            wz2 <= 1
        ) * self.FdrainFC2 * exp(self.beta2 * scalar(wz2 - 1))
        Dz1 = max(0, min(fD1 * Sz1, Sz1 - 1e-2))
        Dz2 = max(0, min(fD2 * Sz2, Sz2 - 1e-2))
        D01 = Dz1
        D02 = Dz2
        self.S01 = self.S01 - D01
        self.S02 = self.S02 - D02
        # Shallow root zone water balance (Ss)
        self.Ss1 = self.Ss1 + D01 - Us1
        self.Ss2 = self.Ss2 + D02 - Us2
        SzFC1 = self.SsFC1
        SzFC2 = self.SsFC2
        Sz1 = self.Ss1
        Sz2 = self.Ss2
        wz1 = max(1e-2, Sz1) / SzFC1
        wz2 = max(1e-2, Sz2) / SzFC2
        fD1 = scalar(wz1 > 1) * max(self.FdrainFC1, 1 - 1 / wz1) + scalar(
            wz1 <= 1
        ) * self.FdrainFC1 * exp(self.beta1 * scalar(wz1 - 1))
        fD2 = scalar(wz2 > 1) * max(self.FdrainFC2, 1 - 1 / wz2) + scalar(
            wz2 <= 1
        ) * self.FdrainFC2 * exp(self.beta2 * scalar(wz2 - 1))
        Dz1 = max(0, min(fD1 * Sz1, Sz1 - 1e-2))
        Dz2 = max(0, min(fD2 * Sz2, Sz2 - 1e-2))
        Ds1 = Dz1
        Ds2 = Dz2
        self.Ss1 = self.Ss1 - Ds1
        self.Ss2 = self.Ss2 - Ds2
        # Deep root zone water balance (Sd) (2.6)
        self.Sd1 = self.Sd1 + Ds1 - Ud1
        self.Sd2 = self.Sd2 + Ds2 - Ud2
        SzFC1 = self.SdFC1
        SzFC2 = self.SdFC2
        Sz1 = self.Sd1
        Sz2 = self.Sd2
        wz1 = max(1e-2, Sz1) / SzFC1
        wz2 = max(1e-2, Sz2) / SzFC2
        fD1 = scalar(wz1 > 1) * max(self.FdrainFC1, 1 - 1 / wz1) + scalar(
            wz1 <= 1
        ) * self.FdrainFC1 * exp(self.beta1 * scalar(wz1 - 1))
        fD2 = scalar(wz2 > 1) * max(self.FdrainFC2, 1 - 1 / wz2) + scalar(
            wz2 <= 1
        ) * self.FdrainFC2 * exp(self.beta2 * scalar(wz2 - 1))
        Dz1 = max(0, min(fD1 * Sz1, Sz1 - 1e-2))
        Dz2 = max(0, min(fD2 * Sz2, Sz2 - 1e-2))
        Dd1 = Dz1
        Dd2 = Dz2
        self.Sd1 = self.Sd1 - Dd1
        self.Sd2 = self.Sd2 - Dd2
        Y1 = min(
            self.Fgw_conn1 * max(0, self.wdlimU1 * self.SdFC1 - self.Sd1), Sghru1 - Eg1
        )
        Y2 = min(
            self.Fgw_conn2 * max(0, self.wdlimU2 * self.SdFC2 - self.Sd2), Sghru2 - Eg2
        )
        # Y = Fgw_conn.*max(0,wdlimU.*SdFC-Sd); #nog matlab script
        self.Sd1 = self.Sd1 + Y1
        self.Sd2 = self.Sd2 + Y2

        # CATCHMENT WATER BALANCE
        # Groundwater store water balance (Sg) (2.5)
        NetGf = (self.Fhru1 * (Dd1 - Eg1 - Y1)) + (self.Fhru2 * (Dd2 - Eg2 - Y2))
        self.Sg = self.Sg + NetGf
        Sgfree = max(self.Sg, 0)
        Qg = min(Sgfree, (1 - exp(-self.K_gw)) * Sgfree)
        self.Sg = self.Sg - Qg

        # Surface water store water balance (Sr) (2.7)
        self.Sr = self.Sr + (self.Fhru1 * (QR1 - Er1)) + (self.Fhru2 * (QR2 - Er2)) + Qg
        self.Qtot = min(self.Sr, (1 - exp(-self.K_rout)) * self.Sr)
        self.Sr = self.Sr - self.Qtot

        # VEGETATION ADJUSTMENT (5)

        fveq1 = (
            (1 / max((self.E01 / Utot1) - 1, 1e-3))
            * (keps / (1 + keps))
            * (ga1 / Gsmax1)
        )
        fveq2 = (
            (1 / max((self.E02 / Utot2) - 1, 1e-3))
            * (keps / (1 + keps))
            * (ga2 / Gsmax2)
        )
        fvmax1 = 1 - exp(-self.LAImax1 / self.LAIref1)
        fvmax2 = 1 - exp(-self.LAImax2 / self.LAIref2)
        fveq1 = min(fveq1, fvmax1)
        fveq2 = min(fveq2, fvmax2)
        dMleaf1 = -ln(1 - fveq1) * self.LAIref1 / self.SLA1 - self.Mleaf1
        dMleaf2 = -ln(1 - fveq2) * self.LAIref2 / self.SLA2 - self.Mleaf2

        # Mleafnet1 = dMleaf1 * (dMleaf1/self.Tgrow1) + dMleaf1 * dMleaf1/self.Tsenc1
        # Mleafnet2 = dMleaf2 * (dMleaf1/self.Tgrow2) + dMleaf2 * dMleaf2/self.Tsenc2
        Mleafnet1 = (
            scalar(dMleaf1 > 0) * (dMleaf1 / self.Tgrow1)
            + scalar(dMleaf1 < 0) * dMleaf1 / self.Tsenc1
        )
        Mleafnet2 = (
            scalar(dMleaf2 > 0) * (dMleaf2 / self.Tgrow2)
            + scalar(dMleaf2 < 0) * dMleaf2 / self.Tsenc2
        )

        self.Mleaf1 = self.Mleaf1 + Mleafnet1
        self.Mleaf2 = self.Mleaf2 + Mleafnet2
        self.LAI1 = self.SLA1 * self.Mleaf1  # (5.3)
        self.LAI2 = self.SLA2 * self.Mleaf2

        # Updating diagnostics
        self.LAI1 = self.SLA1 * self.Mleaf1  # (5.3)
        self.LAI2 = self.SLA2 * self.Mleaf2
        fveg1 = 1 - exp(-self.LAI1 / self.LAIref1)  # (5.3)
        fveg2 = 1 - exp(-self.LAI2 / self.LAIref2)
        fsoil1 = 1 - fveg1
        fsoil2 = 1 - fveg2
        w01 = self.S01 / self.S0FC1  # (2.1)
        w02 = self.S02 / self.S0FC2
        ws1 = self.Ss1 / self.SsFC1  # (2.1)
        ws2 = self.Ss2 / self.SsFC2
        wd1 = self.Sd1 / self.SdFC1  # (2.1)
        wd2 = self.Sd2 / self.SdFC2


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
        "default_w3ra"
    )  # "D:/trambaue/_Projects/GLOFFIS/201501/GLOFFIS_SA/Modules/openstreams_w3ra/"
    runId = "run_default"
    configfile = "wflow_W3RA.ini"
    _lastTimeStep = 0
    _firstTimeStep = 0
    timestepsecs = 86400

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

    opts, args = getopt.getopt(argv, "C:S:T:c:s:R:")

    for o, a in opts:
        if o == "-C":
            caseName = a
        if o == "-R":
            runId = a
        if o == "-c":
            configfile = a
        if o == "-s":
            timestepsecs = int(a)

    # if (len(opts) <=1):
    #    usage()

    starttime = dt.datetime(1990, 01, 01)

    if _lastTimeStep < _firstTimeStep:
        print "The starttimestep (" + str(
            _firstTimeStep
        ) + ") is smaller than the last timestep (" + str(_lastTimeStep) + ")"
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
        if o == "-T":
            configset(myModel.config, "run", "endtime", a, overwrite=True)
        if o == "-S":
            configset(myModel.config, "run", "starttime", a, overwrite=True)

    dynModelFw.setupFramework()

    dynModelFw._runInitial()
    dynModelFw._runResume()
    # dynModelFw._runDynamic(0,0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
