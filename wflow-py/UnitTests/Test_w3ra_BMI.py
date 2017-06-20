__author__ = 'schelle'

import unittest
import logging
import sys
import datetime
sys.path = ['../'] + sys.path
import wflow.wflow_bmi as bmi
import time
from dateutil import parser
import calendar
import os
import wflow.wflow_sceleton as wf
import numpy as np

"""
Simple test for wflow bmi framework
"""


class MyTest(unittest.TestCase):




    def testbmirunnetcdfw3ra(self):
        bmiobj = bmi.wflowbmi_csdms()
        bmiobj.initialize_config('../../examples/openstreams_w3ra_usa/wflow_w3ra.ini',loglevel=logging.DEBUG)
        bmiobj.set_attribute_value('run:runlengthdetermination','intervals')

        stime= calendar.timegm(parser.parse("2014-05-13 00:00:00").utctimetuple())
        etime = calendar.timegm(parser.parse("2014-05-21 00:00:00").utctimetuple())
        bmiobj.set_start_time(stime)
        bmiobj.set_end_time(etime)
        st = bmiobj.get_start_time()
        #print st
        ett = bmiobj.get_end_time()
        ts = bmiobj.get_time_step()

        bmiobj.initialize_model()
        curtime = bmiobj.get_current_time()
        cnt = 0
        lastcurtime = bmiobj.get_current_time()
        while curtime < ett:
            avar = bmiobj.get_value('LAI1')
            bmiobj.set_value('PRECIP',avar)
            cnt = cnt + 1
            bmiobj.update_until(curtime + ts)
            #print (curtime + ts)/ts
            curtime = bmiobj.get_current_time()
            #print bmiobj.get_current_time() - lastcurtime
            lastcurtime = bmiobj.get_current_time()


        bmiobj.finalize()
        # Check the values in a state file as a refrence. This is what the baselien model gives
        x, y, data, FillVal = wf.readMap('../../examples/openstreams_w3ra_usa/run_default/outstate/Sd2.map','PCRaster')
        tmean = np.ma.masked_invalid(data.astype(np.float64)).mean()
        tmax = np.ma.masked_invalid(data.astype(np.float64)).max()
        tmin = np.ma.masked_invalid(data.astype(np.float64)).min()
        self.assertAlmostEquals(266.18075561523438, tmax)
        self.assertAlmostEquals(-7.8522729383405979e+37, tmean)
        self.assertAlmostEquals(-3.4028234663852886e+38, tmin)



if __name__ == '__main__':
    unittest.main()
