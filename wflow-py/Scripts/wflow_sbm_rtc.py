
# coding: utf-8

# <h1>Basic test of the wflow BMI interface + reservoir

# In[10]:
import numpy
import os
import os.path
import shutil, glob
import getopt

from wflow.wf_DynamicFramework import *

import wflow.wflow_bmi as bmi
import logging
import wflow.wflow_adapt as adapter
import datetime, calendar
import os
import numpy as np
import ConfigParser

reload(bmi)

''' Todo: include error handling:
        - When WFlow_ids in the mapping section of config are not in the map
        - When RTC doesn't have data on the WFlow simulation period
        - When a RTC id in the id mapping section is not part of the model
        - ...
    Todo: add logging
'''

# In[]: Parse config-file
Config = ConfigParser.ConfigParser()

cur_dir = os.getcwd()
inifile = Config.read(os.path.join(cur_dir,'wflow_rtctools.ini'))
print inifile
Config.sections()
print Config.sections()

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
    
dir_rtc = os.path.join(cur_dir,(ConfigSectionMap("model")['dir_rtc_model']))
dir_wflow = os.path.join(cur_dir,(ConfigSectionMap("model")['dir_wflow_model']))  
Bin_RTC = os.path.join(cur_dir,(ConfigSectionMap("RTC wrapper engine")['bin_rtc']))
inflow_map=ConfigSectionMap("id_maps")['samplemap_in'] 
outflow_map=ConfigSectionMap("id_maps")['samplemap_out']  
ldd_map = ConfigSectionMap("LDD")['ldd_res']

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

# In[]: Initialize the RTC-Tools model
os.chdir(Bin_RTC)
from wrapperExtended import BMIWrapperExtended
RTC_model = BMIWrapperExtended(engine=os.path.join(Bin_RTC,"RTCTools_BMI"))
RTC_model.initialize('..')


# In[]: Initialize the WFlow model
os.chdir(dir_wflow)
LA_model = bmi.wflowbmi_csdms()
LA_model.initialize(('wflow_sbm.ini'), loglevel=logging.ERROR)


# In[]: Investigate start time, end time and time step of both models

print 'WFlow:'
LA_dt = LA_model.get_value("timestepsecs")
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
Reservoir_inflow = pcr2numpy(scalar(pcraster.readmap(os.path.join(dir_wflow,inflow_map))),np.NaN)
Reservoir_outflow = pcr2numpy(scalar(pcraster.readmap(os.path.join(dir_wflow,outflow_map))),np.NaN)

inflow_list = list(np.unique(Reservoir_inflow[~np.isnan(Reservoir_inflow)]))
outflow_list = list(np.unique(Reservoir_outflow[~np.isnan(Reservoir_outflow)]))


# In[]:  Overwrite TopoLdd with modified version
ldd_map = 'staticmaps/wflow_ldd_rtc.map'
ldd = pcraster.pcr2numpy(pcraster.readmap(os.path.join(dir_wflow,ldd_map)), np.NaN)

LA_model.set_value("TopoLdd",flipud(ldd))


# In[]: Now run the models

t = LA_start
timecounter = 0

while t < min(LA_end, RTC_end):

    # run the WFlow model
    LA_model.update()
    
    # Get the inflow from the wflow model runoff map and map
    inflowQ = flipud(LA_model.get_value("SurfaceRunoff"))
    
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
    inflowfield = zeros_like(inflowQ)
    for idx, wflow_id in enumerate(id_out_wflow):
        rtc_id = id_out_rtc[id_out_wflow.index(str(wflow_id))]
        Qout = RTC_model.get_var(rtc_id)
        print rtc_id + ' = ' + str(Qout)
        inflowfield[np.where(Reservoir_outflow==int(wflow_id))] = Qout    

    LA_model.set_value("IF",flipud(inflowfield))
    # This not not bmi but needed to update the kinematic wave reservoit
    LA_model.myModel.updateRunOff()
    
    #t = LA_model.get_current_time()
    print t
    t += LA_dt


# In[]: Finalize....

LA_model.finalize()
RTC_model.finalize()

