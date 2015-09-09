__author__ = 'schelle'

import unittest
import logging
import wflow.wflow_bmi as bmi
import time
import os
import datetime
"""
Simple test for wflow bmi framework
"""

class MyTest(unittest.TestCase):

    def testbmiinitfuncs(self):
        bmiobj = bmi.wflowbmi_csdms()

        bmiobj.initialize_config('wflow_sceleton/wflow_sceleton.ini',loglevel=logging.ERROR)
        print("-------------- Time units: ")
        print(bmiobj.get_time_units())
        print("-------------- Default current time: ")
        print(datetime.datetime.utcfromtimestamp(bmiobj.get_current_time()))

        print("-------------- Set start time to 5: ")
        bmiobj.set_start_time(5 * 86400)
        print("-------------- Get current time: ")
        print(datetime.datetime.utcfromtimestamp(bmiobj.get_current_time()))

        print("-------------- Get timestep: ")
        print(bmiobj.get_time_step())
        print("-------------- get endtime: ")
        print(datetime.datetime.utcfromtimestamp(bmiobj.get_end_time()))
        print("-------------- set endtime: ")
        bmiobj.set_end_time(10.0 * 86400)
        print("-------------- get endtime: ")
        print(datetime.datetime.utcfromtimestamp(bmiobj.get_end_time()))

        print("-------------- get_var_units: ")
        print(bmiobj.get_var_units("SS"))
        bmiobj.finalize()


if __name__ == '__main__':
    unittest.main()
