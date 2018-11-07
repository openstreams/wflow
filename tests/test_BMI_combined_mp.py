__author__ = "schelle"

import unittest
import sys

sys.path = ["../wflow"] + ["../"] + sys.path
import wflow_bmi_combined_mp as bmi

"""
Simple test for wflow bmi framework
"""


class MyTest(unittest.TestCase):
    def testbmifuncs(self):

        bmiobj = bmi.wflowbmi_csdms()
        bmiobj.initialize("bmirunner.ini")

        print(bmiobj.get_component_name().split(","))
        print(bmiobj.get_input_var_names())
        print(bmiobj.get_output_var_names())
        print(bmiobj.get_start_time())
        print(bmiobj.get_end_time())
        print(bmiobj.get_current_time())
        print(bmiobj.get_time_step())
        print(bmiobj.get_attribute_names())
        steps = (
            bmiobj.get_end_time() - bmiobj.get_start_time()
        ) / bmiobj.get_time_step() + 1
        print(steps)
        for a in range(0, int(steps)):
            bmiobj.update()

        atn = bmiobj.get_attribute_names()
        print(atn[0])
        print(bmiobj.get_attribute_value(atn[0]))
        bmiobj.finalize()


if __name__ == "__main__":
    unittest.main()
