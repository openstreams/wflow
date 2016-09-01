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
        print st
        ett = bmiobj.get_end_time()
        ts = bmiobj.get_time_step()

        bmiobj.initialize_model()
        curtime = st
        cnt = 0
        lastcurtime = bmiobj.get_current_time()
        while curtime < ett:
            avar = bmiobj.get_value('LAI1')
            bmiobj.set_value('PRECIP',avar)
            cnt = cnt + 1
            bmiobj.update_until(curtime + ts)
            print (curtime + ts)/ts
            curtime = bmiobj.get_current_time()
            print bmiobj.get_current_time() - lastcurtime
            lastcurtime = bmiobj.get_current_time()


        bmiobj.finalize()


if __name__ == '__main__':
    unittest.main()
