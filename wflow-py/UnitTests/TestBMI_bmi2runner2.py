__author__ = "schelle"

import unittest
import sys

sys.path = ["../wflow"] + ["../"] + sys.path
sys.path = ["../Scripts"] + sys.path
import bmi2runner as bmirun
import wflow_bmi_combined as bmi

"""
Simple test for wflow bmi framework
"""


class MyTest(unittest.TestCase):
    def testbmirunner_set(self):

        configfile = "combined/bmirunner.ini"
        bmiobj = bmi.wflowbmi_csdms()
        bmiobj.initialize_config(configfile)
        bmiobj.initialize_model()
        start = bmiobj.get_start_time()
        end = bmiobj.get_end_time()
        bmiobj.set_start_time(start)
        # bmiobj.set_end_time(end)
        # Get time for the loop

        ts = bmiobj.get_time_step()
        curtime = bmiobj.get_current_time()
        # Loop over the time duration

        while curtime < end:
            bmiobj.update_until(curtime + ts)
            curtime = bmiobj.get_current_time()

        bmiobj.finalize()


if __name__ == "__main__":
    unittest.main()
