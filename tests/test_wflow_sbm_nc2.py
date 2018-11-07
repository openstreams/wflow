__author__ = "schelle"

import unittest
import wflow.wflow_sbm as wf

"""
Run sceleton for 10 steps and checks if the outcome is approx that of the reference run
"""


class MyTest(unittest.TestCase):
    def testapirun_netcfd2(self):
        startTime = 1
        stopTime = 30
        currentTime = 1

        # set runid, clonemap and casename. Also define the ini file
        runId = "unittest"
        configfile = "wflow_sbm_NC2.ini"
        wflow_cloneMap = "wflow_catchment.map"
        caseName = "../examples/wflow_rhine_sbm_nc"

        myModel = wf.WflowModel(wflow_cloneMap, caseName, runId, configfile)
        # initialise the framework
        dynModelFw = wf.wf_DynamicFramework(myModel, stopTime, startTime)

        # Load model config from files and check directory structure
        dynModelFw.createRunId(NoOverWrite=False, level=wf.logging.DEBUG)
        # Run the initial part of the model (reads parameters and sets initial values)
        dynModelFw._runInitial()  # Runs initial part

        dynModelFw._runResume()  # gets the state variables
        sump = 0.0
        for ts in range(startTime, stopTime + 1):
            dynModelFw._runDynamic(ts, ts)  # runs for all timesteps

        dynModelFw._runSuspend()  # saves the state variables
        dynModelFw._wf_shutdown()


if __name__ == "__main__":
    unittest.main()
