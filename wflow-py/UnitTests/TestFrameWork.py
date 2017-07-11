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
        stopTime = 20
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
        dynModelFw.createRunId(NoOverWrite=False,level=wf.logging.DEBUG)
        # Run the initial part of the model (reads parameters and sets initial values)
        dynModelFw._runInitial() # Runs initial part

        dynModelFw._runResume() # gets the state variables

        for ts in range(startTime,stopTime):
            dynModelFw._runDynamic(ts,ts) # runs for all timesteps
            dynModelFw.logger.info("Doing step: " + str(ts))
        dynModelFw._runSuspend() # saves the state variables
        dynModelFw._wf_shutdown()

        my_data = wf.genfromtxt(os.path.join(caseName,runId,"tes.csv"), delimiter=',')

        self.assertAlmostEquals(14.885358393192291,my_data[:,2].sum())
        my_data_mean = wf.genfromtxt(os.path.join(caseName, runId, "tes_mean_5.csv"), delimiter=',')
        self.assertAlmostEquals( 20.727288454771042, my_data_mean[:, 2].sum())


if __name__ == '__main__':
    unittest.main()
