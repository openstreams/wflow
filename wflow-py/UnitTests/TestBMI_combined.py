__author__ = "schelle"

import unittest
import logging
import sys

sys.path = ["../wflow"] + ["../"] + sys.path
import wflow_bmi_combined as bmi
import time
import os

"""
Simple test for wflow bmi framework
"""


class MyTest(unittest.TestCase):
    def testbmifuncs(self):

        bmiobj = bmi.wflowbmi_csdms()
        bmiobj.initialize("bmirunner.ini", loglevel=logging.INFO)

        print bmiobj.get_component_name().split(",")
        print bmiobj.get_input_var_names()
        print bmiobj.get_output_var_names()
        print bmiobj.get_start_time()
        endtime = bmiobj.get_end_time()
        curtime = bmiobj.get_current_time()
        print endtime
        print curtime
        print bmiobj.get_time_step()
        print bmiobj.get_attribute_names()
        steps = 0
        print steps
        while curtime < endtime:
            bmiobj.update()
            steps = steps + 1
            curtime = bmiobj.get_current_time()

        atn = bmiobj.get_attribute_names()
        print atn[0]
        print bmiobj.get_attribute_value(atn[0])
        bmiobj.finalize()
        self.assertEquals(steps, 29)
        self.assertEquals(curtime, bmiobj.get_current_time())


if __name__ == "__main__":
    unittest.main()
