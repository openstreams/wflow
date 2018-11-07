__author__ = "schelle"

import unittest
import logging
import wflow.wflow_bmi as bmi
import datetime

"""
Simple test for wflow bmi framework
"""


class MyTest(unittest.TestCase):
    def testbmiinitfuncs(self):
        bmiobj = bmi.wflowbmi_csdms()

        bmiobj.initialize_config(
            "wflow_sceleton/wflow_sceleton_notime.ini", loglevel=logging.INFO
        )
        print("-------------- Time units: ")
        print((bmiobj.get_time_units()))
        print("-------------- Default current time: ")
        print((datetime.datetime.utcfromtimestamp(bmiobj.get_current_time())))
        print("-------------- Default start time: ")
        print((datetime.datetime.utcfromtimestamp(bmiobj.get_start_time())))
        print("-------------- Default end time: ")
        print((datetime.datetime.utcfromtimestamp(bmiobj.get_end_time())))
        print("-------------- Set start time to 5: ")
        bmiobj.set_start_time(5 * 86400)
        print("-------------- Set start time to 10: ")
        bmiobj.set_end_time(10 * 86400)

        print("-------------- Updated start time: ")
        stime = datetime.datetime.utcfromtimestamp(bmiobj.get_start_time())
        print("-------------- Updated end time: ")
        etime = datetime.datetime.utcfromtimestamp(bmiobj.get_end_time())

        print("Init model...")
        bmiobj.initialize_model()
        print((bmiobj.dynModel._userModel().config.get("run", "starttime")))

        nstime = datetime.datetime.utcfromtimestamp(bmiobj.get_start_time())

        netime = datetime.datetime.utcfromtimestamp(bmiobj.get_end_time())

        self.assertEqual(stime, nstime)
        self.assertEqual(etime, netime)
        print("-------------- Set start time to 5: ")
        bmiobj.set_start_time(5 * 86400)
        print("-------------- Get current time: ")
        print((datetime.datetime.utcfromtimestamp(bmiobj.get_current_time())))

        print("-------------- Get timestep: ")
        print((bmiobj.get_time_step()))
        print("-------------- get endtime: ")
        print((datetime.datetime.utcfromtimestamp(bmiobj.get_end_time())))
        print("-------------- set endtime: ")
        bmiobj.set_end_time(10.0 * 86400)
        print("-------------- get endtime: ")
        print((datetime.datetime.utcfromtimestamp(bmiobj.get_end_time())))

        print("-------------- get_var_units: ")
        print((bmiobj.get_var_units("SS")))
        bmiobj.finalize()


if __name__ == "__main__":
    unittest.main()
