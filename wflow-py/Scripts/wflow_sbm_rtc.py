
# coding: utf-8

# <h1>Basic test of the wflow BMI interface + reservoir

# In[10]:
import numpy
import os
import os.path
import shutil, glob
import getopt
import sys
import getopt

from wflow.wf_DynamicFramework import *

import wflow.wflow_bmi as bmi
import logging
import wflow.wflow_adapt as adapter
import datetime, calendar
import os
import numpy as np
import ConfigParser

#reload(bmi)

''' Todo: include error handling:
        - When WFlow_ids in the mapping section of config are not in the map
        - When RTC doesn't have data on the WFlow simulation period
        - When a RTC id in the id mapping section is not part of the model
        - ...
    Todo: add logging
'''

def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


def gettimestepfname(name,path,timestep):
    """
    Get the pcraster filename fro this step
    :param name:
    :param path:
    :param timestep:
    :return:
    """

    below_thousand = timestep % 1000
    above_thousand = timestep / 1000
    fname  = str(name + '%0' + str(8-len(name)) + '.f.%03.f') % (above_thousand, below_thousand)
    fname = os.path.join(path,fname)

    return fname

########################################################################
## Process command-line options                                        #
########################################################################

argv = sys.argv

try:
   opts, args = getopt.getopt(argv[1:], 'c:w:I:')
   print opts
except getopt.error, msg:
   print 'cannot parse commandline'
   sys.exit

for o, a in opts:
   if o == '-c': configfile = a
   if o == '-w' : cur_dir = os.path.abspath(a)
   if o == '-I' : IniFile = a

Config = ConfigParser.ConfigParser()
#inifile = Config.read('c:\FEWS\SI-WAMI\SI-WAMI\Modules\RTC\wflow_rtctools.ini')
inifile = Config.read(configfile)
########################################################################
## Parse ini-file                                                      #
########################################################################

#cur_dir = os.getcwd()
Config.sections()
os.chdir(cur_dir)
dir_rtc = os.path.join(cur_dir,(ConfigSectionMap("model")['dir_rtc_model']))
dir_wflow = os.path.join(cur_dir,(ConfigSectionMap("model")['dir_wflow_model']))
Bin_RTC = os.path.join(cur_dir,(ConfigSectionMap("RTC wrapper engine")['bin_rtc']))
inflow_map=os.path.join(dir_wflow,(ConfigSectionMap("id_maps")['samplemap_in']))
outflow_map=os.path.join(dir_wflow,(ConfigSectionMap("id_maps")['samplemap_out']))
ldd_map=os.path.join(dir_wflow,(ConfigSectionMap("ldd")['ldd_res']))

# id's wflow and RTC reservoir inflow points
id_in_rtc=[]
id_in_wflow=list(ConfigSectionMap("id_mapping_inflow"))
for index in range(len(id_in_wflow)):
    id_in_rtc.append(ConfigSectionMap("id_mapping_inflow")[id_in_wflow[index]])

# id's wflow and RTC reservoir outflow points
id_out_rtc=[]
id_out_wflow=list(ConfigSectionMap("id_mapping_outflow"))
for index in range(len(id_out_wflow)):
    id_out_rtc.append(ConfigSectionMap("id_mapping_outflow")[id_out_wflow[index]])

########################################################################
## Initialize models                                                   #
########################################################################

# In[]: Initialize the RTC-Tools model
os.chdir(Bin_RTC)

from wflow.wrappers.rtc.wrapperExtended import BMIWrapperExtended
#RTC_model = BMIWrapperExtended(engine=os.path.join(Bin_RTC,"RTCTools_BMI"))
RTC_model = BMIWrapperExtended(engine=os.path.join(Bin_RTC,"RTCTools_BMI"))
print 'RTCmodel', Bin_RTC,RTC_model
RTC_model.initialize('..')


# In[]: Initialize the WFlow model
os.chdir(dir_wflow)
LA_model = bmi.wflowbmi_csdms()
LA_model.initialize((IniFile), loglevel=logging.WARN)

# now get the forcings that wflow expects
# The bmi is such that you can get the input variables and the output variables. However, the
# input variable list also contains the in/out variables. So to
# get the input only we subtract the two lists.
invars = LA_model.get_input_var_names()
outvars = LA_model.get_output_var_names()
inputmstacks =  list(set(invars) - set(outvars))



# In[]: Investigate start time, end time and time step of both models

print 'WFlow:'
LA_dt = LA_model.get_time_step()

#LA_start = LA_model.get_start_time()
timeutc = adapter.getStartTimefromRuninfo('inmaps/runinfo.xml')
print timeutc
LA_start = calendar.timegm(timeutc.timetuple())
timeutc = adapter.getEndTimefromRuninfo('inmaps/runinfo.xml')
LA_end = calendar.timegm(timeutc.timetuple())
#LA_end = LA_model.get_end_time()
print LA_dt
print timeutc
print LA_start
print LA_end

print 'RTC-Tools'
RTC_dt = RTC_model.get_time_step()
RTC_start = RTC_model.get_start_time()
RTC_end = RTC_model.get_end_time()

print RTC_dt
print RTC_start
print RTC_end

if LA_start != RTC_start:
    print 'Error: start time of both model is not identical !!!'

if LA_dt != RTC_dt:
    print 'Error: time step of both models is not identical !!!'


# In[]:  Read and map reservoir inflow and outflow locations
Reservoir_inflow = pcr2numpy(scalar(pcraster.readmap(os.path.abspath(inflow_map))),np.NaN)
Reservoir_outflow = pcr2numpy(scalar(pcraster.readmap(os.path.abspath(outflow_map))),np.NaN)

inflow_list = list(np.unique(Reservoir_inflow[~np.isnan(Reservoir_inflow)]))
outflow_list = list(np.unique(Reservoir_outflow[~np.isnan(Reservoir_outflow)]))


# In[]:  Overwrite TopoLdd with modified version
ldd = pcraster.pcr2numpy(pcraster.readmap(os.path.join(dir_wflow,ldd_map)), np.NaN).astype(np.float32)
LA_model.set_value("TopoLdd",flipud(ldd).copy())


########################################################################
## Run models                                                          #
########################################################################

t = LA_start
timecounter = 0

while t < min(LA_end, RTC_end):
    #print "timestep = " + str(t)
    # run the WFlow model

    # first read forcing mapstacks (set in API section) and give to the model
    for thisstack in inputmstacks:
        toread =  gettimestepfname(thisstack,os.path.join(dir_wflow,'inmaps'),timecounter+1)
        nptoset = flipud(pcr2numpy(scalar(pcraster.readmap(os.path.abspath(toread))),-999.0)).copy()
        LA_model.set_value(thisstack,nptoset)

    LA_model.update()
    print "calculation timestep = " + str(timecounter)

    # Get the inflow from the wflow model runoff map and map
    inflowQ = flipud(LA_model.get_value("SurfaceRunoff")).copy()

    # Map the sum of WFlow Inflow to RTC
    for idx, wflow_id in enumerate(id_in_wflow):
        value = np.ndarray(shape=(1,1), dtype=float, order='F')
        value[0][0] = np.sum(inflowQ[np.where(Reservoir_inflow==int(wflow_id))])
        rtc_id = id_in_rtc[id_in_wflow.index(str(wflow_id))]
        print rtc_id + ' = ' + str(value[0][0])
        RTC_model.set_value(rtc_id, value)

    # run the RTC-Tools model
    RTC_model.update(-1.0)

    # Extract RTC outflow and supply on WFlow 'inflowfield'
    inflowfield = zeros_like(inflowQ).copy()
    for idx, wflow_id in enumerate(id_out_wflow):
        rtc_id = id_out_rtc[id_out_wflow.index(str(wflow_id))]
        Qout = RTC_model.get_var(rtc_id)
        if isfinite(Qout): # no nan's into wflow
            inflowfield[Reservoir_outflow==int(wflow_id)] = Qout

    LA_model.set_value("IF",flipud(inflowfield).copy())
    # This not not bmi but needed to update the kinematic wave reservoir
    LA_model.myModel.updateRunOff()

    #t = LA_model.get_current_time()
    t += LA_dt
    timecounter += 1


# In[]: Finalize....

LA_model.finalize()
RTC_model.finalize()


