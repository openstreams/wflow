__author__ = "schelle"

import os
import unittest

import numpy as np
import wflow.wflow_sbm as wf

"""
Run wflow_sbm for 10 steps and checks if the outcome is approx that of the reference run
"""


class MyTest(unittest.TestCase):
    def testapirunexample(self):
        startTime = 1
        stopTime = 30
        currentTime = 1

        orgdir = os.getcwd()
        os.chdir("../examples/")
        # set runid, clonemap and casename. Also define the ini file
        runId = "unittestsbm"
        configfile = "wflow_sbm.ini"
        wflow_cloneMap = "wflow_subcatch.map"
        caseName = "wflow_rhine_sbm/"

        myModel = wf.WflowModel(wflow_cloneMap, caseName, runId, configfile)
        # initialise the framework
        dynModelFw = wf.wf_DynamicFramework(myModel, stopTime, startTime)

        # Load model config from files and check directory structure
        dynModelFw.createRunId(NoOverWrite=False, level=wf.logging.WARN)
        # Run the initial part of the model (reads parameters and sets initial values)
        dynModelFw._runInitial()  # Runs initial part

        dynModelFw._runResume()  # gets the state variables
        sump = 0.0
        for ts in range(startTime, stopTime):
            dynModelFw._runDynamic(ts, ts)  # runs for all timesteps
            dynModelFw.logger.info("Doing step: " + str(ts))
        dynModelFw._runSuspend()  # saves the state variables
        dynModelFw._wf_shutdown()

        # now read the csv results acn check of they match the first run
        # Sum should be approx c 4.569673676
        my_data = np.genfromtxt(
            os.path.join(caseName, runId, "specrun.csv"), delimiter=","
        )

        os.chdir(orgdir)

        print("Checking specific runoff ....")
        self.assertAlmostEqual(0.05023503426991738, my_data[:, 2].sum(), places=4)


if __name__ == "__main__":
    unittest.main()
