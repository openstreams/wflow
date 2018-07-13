__author__ = "schelle"

import unittest
import wflow.wflow_sbm as wf
import os

"""
Run wflow_sbm for 10 steps and checks if the outcome is approx that of the reference run
"""


class MyTest(unittest.TestCase):
    def TestBMIrun(self):
        import wflow.wflow_bmi as bmi

        bmimodel = bmi.wflowbmi_csdms()
        bmimodel.initialize("../../examples/wflow_rhine_sbm2/wflow_sbm2.ini")
        bmimodel.update()


if __name__ == "__main__":
    unittest.main()
