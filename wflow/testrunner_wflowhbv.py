"""
Controlling the model via the API
---------------------------------
Example of how to use the framework to control the wflow\_\* models
from other software, in this case a python script.

In this example case the controlling piece of software (this script)
provides the forcing data using the API. It also interrogates the model
to get results and display them on screen.

Each wflow\_\* model should provide a def supplyVariableNamesAndRoles function
that returns a list of variables and their roles. This function can than
be used by the controlling program to interrogate the model. 

Some important things to consider:
    
+ when initialzing the framework you must do so for the maximum number of \
timesteps you want to run the model for. This is needed for the pcraster\
timeseries output functions.

+ the framework internally saves the state variable for the last timesteps. As\
 such it is possible to redo the previous timesteps after calling the \
quickresume function. If you want to go back more than one timesteps you \
will need to call the normal resume function wich will restore the states \
to the last time the suspend function was call. Administering this is left \
to the controlling application.

+ Inputs (forcing variables) must be set at the start of a timesteps, results\
 should be read after each timestep.
 
"""


from .wflow_hbv import *

npmap0 = []
ltt = []


def main():

    global npmap0
    global ltt
    # define start and stop time of the run
    startTime = 1
    stopTime = 5
    currentTime = 1

    # set runid, cl;onemap and casename. Also define the ini file
    runId = "memtest"
    configfile = "wflow_hbv_mem.ini"
    wflow_cloneMap = "wflow_subcatch.map"
    caseName = "../../examples/wflow_rhine_hbv"
    # Mske a usermodel object
    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    # initialise the framework
    dynModelFw = wf_DynamicFramework(myModel, stopTime, startTime)

    # Load model config from files and check directory structure
    dynModelFw.createRunId(NoOverWrite=False)
    # Run the initial part of the model (reads parameters and sets initial values)
    dynModelFw._runInitial()  # Runs initial part

    dynModelFw._runResume()  # gets the state variables

    # Get list of variables supplied by the model

    dd = dynModelFw.wf_supplyVariableNamesAndRoles()
    print(dd)
    dynModelFw.wf_setValueLdd("TopoLdd", 5.0, 6.46823, 51.6821)
    npmap0 = dynModelFw.wf_supplyMapAsNumpy("TopoLdd")
    ltt = dynModelFw.wf_supplyMapAsList("SurfaceRunoff")

    for ts in range(startTime, stopTime):

        # Get value at pit

        inflowQ = dynModelFw.wf_supplyScalar("SurfaceRunoff", 6.46823, 51.6821)
        outflowQ = dynModelFw.wf_supplyScalar("SurfaceRunoff", 6.43643, 51.7226)

        # Ass inflow to outflow
        # dynModelFw.wf_setValue("ForecQ_qmec", -1.0 * inflowQ  ,6.46823,51.6821)
        Resoutflow = inflowQ
        dynModelFw.wf_setValue("ForecQ_qmec", Resoutflow, 6.43643, 51.7226)
        dynModelFw.wf_setValues("P", scalar(ts) * 0.1)
        # dynModelFw.wf_setValue("ForecQ_qmec",inflowQ * 1000 ,6.47592,51.7288)
        # update runoff ONLY NEEDED IF YOU FIDDLE WITH THE KIN_WAVE RESERVOIR
        myModel.updateRunOff()
        dynModelFw._runDynamic(ts, ts)  # runs for all timesteps
        # dynModelFw.wf_setValue("SurfaceRunoff",0.0,6.46823,51.6821)
        # dynModelFw.wf_setValue("SurfaceRunoff",0.0,6.11535,51.8425)
        npmap0 = dynModelFw.wf_supplyMapAsNumpy("ForecQ_qmec")
        npmap1 = dynModelFw.wf_supplyMapAsNumpy("P")
        dynModelFw.wf_setValuesAsNumpy("xx", npmap1)
        npmap2 = dynModelFw.wf_supplyMapAsNumpy("DezeBestaatNiet")
        # myModel.updateRunOff()

    dynModelFw._runSuspend()  # saves the state variables
    os.chdir("../../")


if __name__ == "__main__":
    main()
