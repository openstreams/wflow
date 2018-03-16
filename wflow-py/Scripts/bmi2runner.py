"""

bmi2runner - runs multiple linked bmi models

Usage:

        bmi2runner -c configfile -l loglevel

        -l: loglevel (most be one of DEBUG, WARNING, ERROR)

Example ini file:

::

    [models]
    wflow_sbm=wflow_sbm/wflow_sbm_comb.ini
    wflow_routing=wflow_routing/wflow_routing_comb.ini

    [exchanges]
    # From_model/var -> To_model/var
    wflow_sbm@InwaterMM=wflow_routing@IW

"""
import wflow
import wflow.wflow_bmi_combined as wfbmi
import wflow.pcrut as pcrut
import os
import time



"""
Perform command line execution of the model.
"""

configfile = "bmi2runner.ini"

loglevel = 'INFO'
combilogger = pcrut.setlogger('bmi2runner.log','bmi2runner_logging',thelevel=loglevel)

# Construct object and initilize the models
combilogger.info('Starting combined bmi object')
bmiobj = wfbmi.wflowbmi_csdms()

#this line is needed when running from batch script
#os.sys.path.append(os.getcwd() +'\\rtc\\rtc_brantas\\bin\\')

bmiobj.initialize_config(configfile,loglevel=loglevel)
bmiobj.initialize_model()

#Get and set start and end times
start = bmiobj.get_start_time()
end = bmiobj.get_end_time()
bmiobj.set_start_time(start)
bmiobj.set_end_time(end)

#Update models (if necessary) to start time
bmiobj.update_to_start_time(start)

#Number of steps to run models
ts = bmiobj.get_time_step()
steps = int((end - start)/ts + 1)

print 'start = ', start#time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(start))
print 'start time rtc =', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(bmiobj.bmimodels['RTC-Tools'].get_start_time()))
print 'start time wflow =', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(bmiobj.bmimodels['wflow_sbm'].get_start_time()))
print 'start time lintul =', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(bmiobj.bmimodels['wflow_lintul'].get_start_time()))


cts = bmiobj.currenttimestep
# Loop over the time duration    
while cts < steps:

    bmiobj.update()
    cts = bmiobj.currenttimestep


#else: 
#    bmiobj.bmimodels['RTC-Tools'].update()
#    bmiobj.bmimodels['wflow_sbm'].update()
#    bmiobj.bmimodels['wflow_lintul'].update()
#    cts = bmiobj.currenttimestep
        
bmiobj.bmimodels['RTC-Tools'].finalize()
bmiobj.bmimodels['wflow_sbm'].finalize()
bmiobj.bmimodels['wflow_lintul'].finalize()
combilogger.info('Finishing run')






