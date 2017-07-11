#!/usr/bin/python

#
# Wflow is Free software, see below:
# 
# Copyright (c) J. Schellekens 2005-2011
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
Definition of the wflow_wave model.
-----------------------------------

Runs the pcraster dynamic wave based on the output from the
kinematic wave.

The wflow\_sbm|hbv model must have saved mapstacks for
water level and discharge for each timestep (run*****.*** and lev*****.***).
If the name of you Q and/or H maps are different specify these in the
[inputmapstacks] section, e.g:
    
::

    [inputmapstacks]
    Q = run
    H = lev
    

The Settings for the dynamic wave itself should be put in the [dynamicwave] section
of the ini file.

::
        
    [dynamicwave]
    # Switch on dynamic wave for main rivers
        
    # Number of timeslices per dynamic waven substep
    TsliceDyn=900
    
    # number of substeps for the dynamic wave with respect to the model timesteps
    dynsubsteps=24
    
    # map with level boundary points
    #wflow_hboun = staticmaps/wflow_outlet.map
    
    # Optional river map for the dynamic wave that must be the same size or smaller as that of the
    # kinematic wave
    wflow_dynriver = staticmaps/wflow_dynriver.map
    
    # a fixed water level for each non-zero point in the wflow_hboun map 
    # level > 0.0 use that level
    # level == 0.0 use timeseries
    # level < 0.0 use upstream water level
    fixedLevel = -8.0
    
    # If this is set to one the program will try to keep the volume at the pits at a constant level
    lowerflowbound=1
    
    # instead of a fixed level a tss file with levels for each timesteps and each 
    # non-zero value in the wflow_hboun map
    #levelTss=intss/Hboun.tss
    #AdaptiveTimeStepping=1
    
    


Usage:
wflow_wave  -C case -R Runid -c inifile -h

    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -I: generate initial conditions from scratchs
    
    -c name of the config file (in the case directory)

    
    -h displays help information
    
    -l: loglevel (most be one of DEBUG, WARNING, ERROR)


    
$Author: schelle $
$Id: wflow_wave.py 913 2014-02-04 13:10:51Z schelle $
$Rev: 913 $

"""
import numpy
import os
import os.path
import shutil, glob
import getopt

from wflow.wf_DynamicFramework import *

#import scipy

from wflow.wflow_adapt  import *


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
      Initialize the object
      
      """
      self.thestep=0
      DynamicModel.__init__(self)

      self.caseName = os.path.abspath(Dir)
      self.clonemappath = os.path.join(os.path.abspath(Dir),"staticmaps",cloneMap)
      setclone(self.clonemappath)
      self.runId = RunDir
      self.Dir = os.path.abspath(Dir)
      self.configfile = configfile
      self.SaveDir = os.path.join(self.Dir,self.runId)


  def runDynamicWave(self):
        """
        Runs the dynamic wave for the main river
        Beware: Experimental, *very* slow and unstable
        """ 
 
        # Determine all the inflow points into the main river
        setglobaloption('manning')
        comb = ordinal(cover(self.DynRiver,0))
        dst = downstream(self.Ldd,comb)
        inf = ifthen(boolean(self.DynRiver),ordinal(0))
        inf = cover(inf,dst)
        self.Qin = ifthenelse(inf > 0,self.SurfaceRunoff * self.timestepsecs,scalar(0.0))
        self.Qin = upstream(self.Ldd,self.Qin)
        
        self.Qin = self.Qin/self.dynsubsteps
        
        # level boundary, fixed or TSS
        if self.fixed_h == 0.0:
            levelBoun = timeinputscalar(self.caseName + self.fixed_h_tss,ordinal(self.dynHBoundary))
        else:
            levelBoun = self.fixed_h
            
        self.oldTsliceDynDyn =  self.TsliceDyn 
        self.WaterLevelDyn = ifthen(self.DynRiver,self.WaterLevelDyn)
        for step in range(self.dynsubsteps):

            self.logger.debug("Dynamic wave substep: " + str(step))
            ChannelSurface = (self.ChannelBottomWidth + (self.ChannelForm * self.WaterLevelDyn * 2.0) + self.ChannelBottomWidth)/2.0
            self.AChannel = min(self.WaterLevelDyn,self.ChannelDepth) * ChannelSurface
            self.AFloodplain = max((self.WaterLevelDyn - self.ChannelDepth) * self.FloodplainWidth,0.0)
                        
            
            self.A = max(self.AChannel + self.AFloodplain,0.0001)
            self.velocity = self.SurfaceRunoffDyn/self.A
            self.crt = abs((self.timestepsecs/(self.TsliceDyn * self.dynsubsteps) * self.velocity)/self.ChannelLength)
            Vol = self.A * self.ChannelLength
            if  self.lowerflowbound:
                Qout = ifthen(boolean(pit(self.Ldd)),Vol * 0.9) 
                self.Qin = cover(Qout,self.Qin)
            
            crt = numpy.max(pcr2numpy(self.crt,0.0)) * 0.5
            if self.AdaptiveTimeStepping:
                self.TsliceDynDyn = numpy.max([self.TsliceDyn,numpy.min([numpy.max([1.0,self.TsliceDyn * crt]),self.timestepsecs/self.dynsubsteps/self.mintimestep])])
                
                self.logger.debug("Estimated timestep: " + str(self.timestepsecs/self.dynsubsteps/self.TsliceDynDyn))
            else:
                self.TsliceDynDyn = self.TsliceDyn            
            #self.TsliceDynDyn = 3600
                

            self.oldTsliceDynDyn = self.TsliceDynDyn
            self.crtsum=self.crtsum + self.crt
            
            self.EffectiveRoughness = ifthenelse(self.WaterLevelDyn>self.ChannelDepth,self.FloodplainRoughness,self.ChannelRoughness)
    
            self.SurfaceRunoffDyn=dynamicwaveq(self.LddIn,self.Qin,self.WaterLevelDyn,self.ChannelBottomLevel,
                                   self.EffectiveRoughness,self.ChannelLength,self.ChannelBottomWidth,
                                   self.ChannelDepth,self.ChannelForm,self.FloodplainWidth,
                                   self.timestepsecs/self.dynsubsteps,self.TsliceDynDyn,self.Structures,
                                   self.StructureA,self.StructureB,self.StructureCrestLevel)/ self.timestepsecs * self.dynsubsteps
            self.WaterLevelDyn=dynamicwaveh(self.LddIn,self.Qin,self.WaterLevelDyn,self.ChannelBottomLevel,
                                   self.EffectiveRoughness,self.ChannelLength,self.ChannelBottomWidth,
                                   self.ChannelDepth,self.ChannelForm,self.FloodplainWidth,
                                   self.timestepsecs/self.dynsubsteps,self.TsliceDynDyn,self.Structures,
                                   self.StructureA,self.StructureB,self.StructureCrestLevel)
           
            if self.fixed_h < 0.0:
                upstr1 = upstream(self.LddIn,self.WaterLevelDyn)
                upstr2 =  upstream(self.LddIn,upstr1)
                upstr3 =  upstream(self.LddIn,upstr2)
                upstr = (upstr1 + upstr2 + upstr3)/3.0
                levelBoun = upstr
                
            self.FloodPlainVol=self.AFloodplain * self.ChannelLength
            self.ChannelVol=self.AChannel * self.ChannelLength
            fxboun = ifthen(self.dynHBoundary > 0,scalar(levelBoun))                
            self.WaterLevelDyn = cover(ifthen(fxboun>0,fxboun),self.WaterLevelDyn)
     


  def stateVariables(self):
      """ 
      *Required*
      
      Returns a list of state variables that are essential to the model. 
      This list is essential for the resume and suspend functions to work.
      
      This function is specific for each model and **must** be present. This is
      where you specify the state variables of you model. If your model is stateless
      this function must return and empty array (states = [])
      

      
      :var SurfaceRunoffDyn.map: Discharge in m^3/sec
      :var WaterLevelDyn.map: Discharge in m
      
      """
      states = ['SurfaceRunoffDyn','WaterLevelDyn']
      
      return states
      
      
  def supplyCurrentTime(self):
      """
      *Optional*
      
      Supplies the current time in seconds after the start of the run
      This function is optional. If it is not set the framework assumes
      the model runs with daily timesteps.
      
      Ouput:
      
          - time in seconds since the start of the model run
          
      """
      
      return self.currentTimeStep() * int(configget(self.config,'model','timestepsecs','86400'))
  
  def suspend(self):
    """      
      Suspends the model to disk. All variables needed to restart the model
      are saved to disk as pcraster maps. Use resume() to re-read them
      
    """
        
    #self.logger.info("Saving initial conditions...")
    #: It is advised to use the wf_suspend() function 
    #: here which will suspend the variables that are given by stateVariables 
    #: function.
    self.logger.info("Saving initial conditions...")
    self.wf_suspend(os.path.join(self.SaveDir,"outstate"))


      
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

    #: Note the use of the configget functione below. This way you sepcify a default
    #: for a parameter but it can be overwritten by the uses in the ini file.
    self.timestepsecs = int(configget(self.config,'model','timestepsecs','86400'))
    self.reinit = int(configget(self.config,"run","reinit","0"))
    
    Qname=configget(self.config,"inputmapstacks","Q","run")
    Hname=configget(self.config,"inputmapstacks","H","lev")
    self.TsliceDyn = int(configget(self.config,"dynamicwave","TsliceDyn","900"))
    self.dynsubsteps = int(configget(self.config,"dynamicwave","dynsubsteps","24"))
    sizeinmetres = int(configget(self.config,"layout","sizeinmetres","0"))
    wflow_dynhboun = configget(self.config,"dynamicwave","wflow_hboun","staticmaps/wflow_hboun.map")
    wflow_dynriver = configget(self.config,"dynamicwave","wflow_dynriver","not_set")
    self.mintimestep=float(configget(self.config,'dynamicwave','mintimestep','1.0'))
    self.fixed_h = float(configget(self.config,"dynamicwave","fixedLevel","8.0"))
    self.dynHBoundary = pcrut.readmapSave(os.path.join(self.Dir, wflow_dynhboun),0.0)
    self.lowerflowbound=configget(self.config,"dynamicwave","lowerflowbound","0")
    self.fixed_h_tss=configget(self.config,"dynamicwave","levelTss","intss/Hboun.tss")
    self.logger.info("Dynamic wave timestep is: " + str(self.timestepsecs/self.dynsubsteps/self.TsliceDyn) + " seconds")
    self.logger.info("Lower boundary file: " + self.fixed_h_tss)
    self.AdaptiveTimeStepping = int(configget(self.config,"dynamicwave","AdaptiveTimeStepping","0"))

    
    self.basetimestep=86400
    self.SaveMapDir = os.path.join(self.Dir,self.runId,"outmaps")
    self.WL_mapstack=self.Dir + "/" + self.runId + "/outmaps/"  + Hname
    self.Q_mapstack=self.Dir + "/" + self.runId + "/outmaps/" + Qname
    self.Altitude=readmap(self.Dir + "/staticmaps/wflow_dem")
    self.River=readmap(self.Dir + "/staticmaps/wflow_river")
    self.Ldd=readmap(self.Dir + "/staticmaps/wflow_ldd")
    self.RiverWidth=readmap(os.path.join(self.Dir,self.runId,"outsum","RiverWidth.map"))
    self.DCL=readmap(os.path.join(self.Dir,self.runId,"outsum","DCL.map"))
    self.ZeroMap=0.0*scalar(self.Altitude)
    self.OutputLoc=readmap(self.Dir + "/staticmaps/wflow_gauges.map")
    self.OutputId=readmap(self.Dir + "/staticmaps/wflow_subcatch.map")

    if wflow_dynriver == "not_set":
        self.DynRiver = boolean(self.River)
    else:
        self.DynRiver = boolean(readmap(os.path.join(self.Dir, wflow_dynriver)))
        
    self.ChannelDepth=pcrut.readmapSave(self.Dir + "/staticmaps/ChannelDepth.map",8.0)
    self.ChannelDepth = self.ChannelDepth * scalar(boolean(self.DynRiver))
    self.ChannelBottomLevel = self.Altitude * scalar(boolean(self.DynRiver)) - self.ChannelDepth 
    self.FloodplainRoughness = pcrut.readmapSave(self.Dir + "/staticmaps/FloodplainRoughness.map",0.4)    
    self.FloodplainRoughness = self.FloodplainRoughness * scalar(boolean(self.DynRiver))
    self.ChannelRoughness = pcrut.readmapSave(self.Dir + "/staticmaps/ChannelRoughness.map",0.03)
    self.ChannelRoughness = self.ChannelRoughness * scalar(boolean(self.DynRiver))
    # COnvert to chezy        
    #self.ChannelRoughness = 1.49/self.ChannelRoughness
    self.ChannelLength = self.DCL * scalar(boolean(self.DynRiver))
    self.ChannelBottomWidth = self.RiverWidth
    self.ChannelForm = pcrut.readmapSave(self.Dir + "/staticmaps/ChannelForm.map",1.0)
    self.ChannelForm = self.ChannelForm * scalar(boolean(self.DynRiver))        
    self.FloodplainWidth =  pcrut.readmapSave(self.Dir + "/staticmaps/FloodplainWidth.map",300.0)  * scalar(boolean(self.DynRiver)) 
    self.FloodplainWidth = max(self.FloodplainWidth,(self.ChannelBottomWidth + (self.ChannelDepth*self.ChannelForm * 2.0))  * scalar(boolean(self.DynRiver)))
    
             
    self.Structures = boolean(self.ZeroMap * scalar(boolean(self.DynRiver))) 
    self.StructureA = self.ZeroMap * scalar(boolean(self.DynRiver))
    self.StructureB = self.ZeroMap * scalar(boolean(self.DynRiver))
    self.StructureCrestLevel = self.ZeroMap * scalar(boolean(self.DynRiver))
    
    report(self.FloodplainWidth,self.Dir + "/" + self.runId + "/outsum/FloodplainWidth.map")
    report(self.ChannelLength,self.Dir + "/" + self.runId + "/outsum/ChannelLength.map")
    report(self.ChannelDepth,self.Dir + "/" + self.runId + "/outsum/ChannelDepth.map")
    report(self.ChannelBottomLevel,self.Dir + "/" + self.runId + "/outsum/ChannelBottomLevel.map")
    report(self.ChannelRoughness,self.Dir + "/" + self.runId + "/outsum/ChannelRoughness.map")
    report(self.FloodplainRoughness,self.Dir + "/" + self.runId + "/outsum/FloodplainRoughness.map")
    report(self.ChannelForm,self.Dir + "/" + self.runId + "/outsum/ChannelForm.map")
    report(self.ChannelBottomWidth,self.Dir + "/" + self.runId + "/outsum/ChannelBottomWidth.map")

 
    # Make seperate LDD for Dynamic Wave
    self.LddIn= lddrepair(ifthen(boolean(self.DynRiver),self.Ldd))
    self.crtsum = self.ZeroMap
 

    self.logger.info("End of initial section.")


  def resume(self):
    """ 

    reads the initial conditions:
    
    :var self.WaterLevelDyn: Dynamic wave waterlevel [m]
    :var self.SurfaceRunoffDyn: Dynamic wave surface runoff [m^3/s]
    """
    #self.logger.info("Reading initial conditions...")
    #: It is advised to use the wf_resume() function 
    #: here which pick upt the variable save by a call to wf_suspend()
    if self.reinit == 1:
        self.logger.info("Setting initial conditions to default (zero!)")
        self.WaterLevelDyn=(self.ZeroMap + 0.1) * scalar(boolean(self.River)) 
        self.SurfaceRunoffDyn=self.ZeroMap * scalar(boolean(self.River))        

    else:
        self.wf_resume(os.path.join(self.Dir, "instate"))


    
  def dynamic(self):
      """
      *Required*
      
      This is where all the time dependent functions are executed. Time dependent
      output should also be saved here.
      
      :var self.FloodPlainVol: Volume of water in the floodplain [m^3]
      :var self.ChannelVol: Volume of water in the channel [m^3]
      :var self.WaterLevelDyn: Water Level [m]
      :var self.SurfaceRunoffDyn: Discharge [m^3/s]
      """

      self.logger.debug("Step: "+str(int(self.thestep + self._d_firstTimeStep))+"/"+str(int(self._d_nrTimeSteps)))
      self.thestep = self.thestep + 1
      self.SurfaceRunoff = self.wf_readmap(self.Q_mapstack,0.0)  
      self.WaterLevel = self.wf_readmap(self.WL_mapstack,0.0)
      
      self.runDynamicWave()


# The main function is used to run the program from the command line

def main(argv=None):  
    """
    *Optional*
    
    Perform command line execution of the model. This example uses the getopt
    module to parse the command line options.
    
    The user can set the caseName, the runDir, the timestep and the configfile.
    """      
    global multpars
    caseName = "default"
    runId = "run_default"
    configfile="wflow_wave.ini"
    _lastTimeStep = 0
    _firstTimeStep = 0
    timestepsecs=86400
    wflow_cloneMap = 'wflow_subcatch.map'

    runinfoFile="runinfo.xml"
    loglevel = logging.DEBUG
    
    # This allows us to use the model both on the command line and to call 
    # the model usinge main function from another python script.
    
    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return     

    opts, args = getopt.getopt(argv, 'C:S:T:c:s:R:fIs:hl:')
    
    for o, a in opts:
        if o == '-C': caseName = a
        if o == '-R': runId = a
        if o == '-c': configfile = a
        if o == '-s': timestepsecs = int(a)
        if o == '-T': _lastTimeStep=int(a)
        if o == '-S': _firstTimeStep=int(a)
        if o == '-l': exec "loglevel = logging." + a            
        if o == '-h': 
            usage()
            return
                    


        
    if (len(opts) <=1):
        usage()

    if _lastTimeStep < _firstTimeStep:
        print "The starttimestep (" + str(_firstTimeStep) +") is smaller than the last timestep (" + str(_lastTimeStep) + ")"
        usage()

        
    myModel = WflowModel(wflow_cloneMap, caseName,runId,configfile)
    dynModelFw = wf_DynamicFramework(myModel, _lastTimeStep,firstTimestep=_firstTimeStep)
    dynModelFw.createRunId(NoOverWrite=False,level=loglevel) 
    for o, a in opts:
        if o == '-I': configset(myModel.config,'model','reinit','1',overwrite=True) 
        if o == '-s': configset(myModel.config,'model','timestepsecs',a,overwrite=True)
 
    
    dynModelFw._runInitial()
    dynModelFw._runResume()
    dynModelFw._runDynamic(_firstTimeStep,_lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()
    

if __name__ == "__main__":
    main()
