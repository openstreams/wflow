__author__ = 'schelle'

import unittest
import wflow.wflow_sceleton as wf
import os
"""
Run sceleton for 10 steps and checks if the outcome is approx that of the reference run
"""

class MyTest(unittest.TestCase):

    def testapirun(self):
        startTime = 1
        stopTime = 10
        currentTime = 1

          # set runid, clonemap and casename. Also define the ini file
        runId = "unittest"
        configfile="wflow_sceleton.ini"
        wflow_cloneMap = 'wflow_catchment.map'
        caseName="wflow_sceleton"

        myModel = wf.WflowModel(wflow_cloneMap, caseName,runId,configfile)
         # initialise the framework
        dynModelFw = wf.wf_DynamicFramework(myModel, stopTime,startTime)

          # Load model config from files and check directory structure
        dynModelFw.createRunId(NoOverWrite=False,level=wf.logging.ERROR)
        # Run the initial part of the model (reads parameters and sets initial values)
        dynModelFw._runInitial() # Runs initial part

        dynModelFw._runResume() # gets the state variables

        for ts in range(startTime,stopTime):
            dynModelFw._runDynamic(ts,ts) # runs for all timesteps
            dynModelFw.logger.info("Doing step: " + str(ts))
        dynModelFw._runSuspend() # saves the state variables
        dynModelFw._wf_shutdown()

        # Now read the csv results and check of they match the first run
        # Sum should be approx c 4.569673676
        my_data = wf.genfromtxt(os.path.join(caseName,runId,"tes.csv"), delimiter=',')

        self.assertAlmostEquals(4.569673676,my_data[:,2].sum())


if __name__ == '__main__':
    unittest.main()
