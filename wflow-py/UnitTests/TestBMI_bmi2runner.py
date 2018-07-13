__author__ = "schelle"

import unittest
import logging
import sys

sys.path = ["../wflow"] + ["../"] + sys.path
sys.path = ["../Scripts"] + sys.path
import bmi2runner as bmirun
import time
import os
import wflow_bmi_combined as bmi

"""
Simple test for wflow bmi framework
"""


class MyTest(unittest.TestCase):
    def testbmirunner(self):

        bmirun.main(["-c", "combined/bmirunner.ini"])
        self.assertTrue(os.path.exists("combined/run_default/wflow.log"))


if __name__ == "__main__":
    unittest.main()
