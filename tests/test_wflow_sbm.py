__author__ = "schelle"

import os
import unittest

import numpy as np
import wflow.wflow_sbm as wf

"""
Run wflow_sbm for 10 steps and checks if the outcome is approx that of the reference run
"""


class MyTest(unittest.TestCase):
    def testapirun(self):
        startTime = 1
        stopTime = 30
        currentTime = 1

        # set runid, clonemap and casename. Also define the ini file
        runId = "unittestsbm"
        configfile = "wflow_sbm.ini"
        wflow_cloneMap = "wflow_catchment.map"
        caseName = "wflow_sbm"

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
            if ts < 10:
                dynModelFw.wf_setValues("P", 0.0)
            elif ts <= 15:
                dynModelFw.wf_setValues("P", 5.0)
                sump = sump + 5.0
            else:
                dynModelFw.wf_setValues("P", 0.0)

            dynModelFw.wf_setValues("PET", 3.0)
            dynModelFw.wf_setValues("TEMP", 10.0)
            dynModelFw._runDynamic(ts, ts)  # runs for all timesteps
            dynModelFw.logger.info("Doing step: " + str(ts))
        dynModelFw._runSuspend()  # saves the state variables
        dynModelFw._wf_shutdown()

        # now read the csv results acn check of they match the first run
        # Sum should be approx c 4.569673676
        my_data = np.genfromtxt(
            os.path.join(caseName, runId, "wbsurf.csv"), delimiter=","
        )

        print("Checking surface water budget ....")
        self.assertAlmostEqual(-1.003901559215592e-06, my_data[:, 2].sum(), places=4)
        my_data = np.genfromtxt(
            os.path.join(caseName, runId, "wbsoil.csv"), delimiter=","
        )
        print("Checking soil water budget ....")
        self.assertAlmostEqual( 0.0007127678763936274, my_data[:, 2].sum(), places=4)
        print("Checking precip sum ....")
        my_data = np.genfromtxt(os.path.join(caseName, runId, "P.csv"), delimiter=",")
        self.assertAlmostEqual(sump, my_data[:, 2].sum())


if __name__ == "__main__":
    unittest.main()
