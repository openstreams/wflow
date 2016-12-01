__author__ = 'schelle'

import unittest
import logging
import sys
sys.path = ['../'] + sys.path
import wflow.wflow_bmi as bmi
import time
import os

"""
Simple test for wflow bmi framework
"""


class MyTest(unittest.TestCase):


    def testbmifuncs(self):

        bmiobj = bmi.wflowbmi_csdms()
        bmiobj.initialize('wflow_sceleton/wflow_sceleton.ini',loglevel=logging.ERROR)

        print("-------------- Grid origin: ")
        gorigin = bmiobj.get_grid_origin('Altitude')
        #print(gorigin)
        self.assertAlmostEquals(sum([45.875934703275561, 5.2088299822062254]), sum(gorigin),places=4)

        print("-------------- Grid shape: ")
        print(bmiobj.get_grid_shape('Altitude'))
        self.assertAlmostEquals(sum([169L, 187L]), sum(bmiobj.get_grid_shape('Altitude')),places=4)

        print("-------------- Grid spacing: ")
        print(bmiobj.get_grid_spacing('Altitude'))

        print("-------------- Grid X: ")
        print(bmiobj.get_grid_x('Altitude'))

        print("-------------- Grid Y: ")
        print(bmiobj.get_grid_y('Altitude'))

        print("-------------- Grid Z: ")
        print(bmiobj.get_grid_z('Altitude'))

        print("-------------- Name: ")
        print(bmiobj.get_component_name())

        print("-------------- Input var names: ")
        print(bmiobj.get_input_var_names())

        print("-------------- UNit of var TEMP: ")
        print(bmiobj.get_var_units('TEMP'))

        print("-------------- UNit of var P: ")
        print(bmiobj.get_var_units('P'))

        print("-------------- Output var names: ")
        print(bmiobj.get_output_var_names())

        print("-------------- Time units: ")
        print(bmiobj.get_time_units())

        print("-------------- Time step: ")
        print(bmiobj.get_time_step())

        print("-------------- Start time: ")
        print(bmiobj.get_start_time())

        print("-------------- Current time: ")
        print(bmiobj.get_current_time())
        a = bmiobj.get_current_time()
        # print(time.localtime(bmiobj.get_current_time()))

        os.environ['TZ'] = 'Europe/London'

        print("-------------- Current time (set to london): ")
        print(bmiobj.get_current_time())
        b = bmiobj.get_current_time()

        self.assertAlmostEquals(a,b)

        print("-------------- update: ")
        bmiobj.update()

        print("-------------- Current time after update: ")
        print(bmiobj.get_current_time())
        print(time.localtime(bmiobj.get_current_time()))

        print("-------------- Start time: ")
        print(bmiobj.get_start_time())
        print(time.localtime(bmiobj.get_start_time()))

        print("-------------- End time: ")
        print(bmiobj.get_end_time())
        print(time.localtime(bmiobj.get_end_time()))

        print("-------------- Grid type: ")
        print(bmiobj.get_grid_type('Altitude'))

        print("-------------- Var type: ")
        print(bmiobj.get_var_type('Altitude'))

        print("-------------- Var rank: ")
        print(bmiobj.get_var_rank('Altitude'))

        print("-------------- Var size: ")
        print(bmiobj.get_var_size('Altitude'))

        print("-------------- Var nbytes: ")
        print(bmiobj.get_var_nbytes('Altitude'))

        print("-------------- Getvalue: ")
        print(bmiobj.get_value('Altitude'))

        print("-------------- Getvalue: ")
        print(bmiobj.get_value('timestepsecs'))

        print("-------------- get_attribute_names: ")
        names = bmiobj.get_attribute_names()
        print names

        print("-------------- get_attribute_value: ")
        print names[0]
        print(bmiobj.get_attribute_value(names[0]))

        print("-------------- set_attribute_value: ")
        print names[0]
        bmiobj.set_attribute_value(names[0],"SET By TEST")
        print(bmiobj.get_attribute_value(names[0]))
        self.assertEquals("SET By TEST",bmiobj.get_attribute_value(names[0]))

        print("-------------- set_start_time: ")
        bmiobj.set_start_time(0)
        print(bmiobj.get_attribute_value("run:starttime"))

        print("-------------- save the state:")
        bmiobj.save_state(".")
        self.assertTrue(os.path.exists("TSoil.map"))
        os.remove("TSoil.map")

        bmiobj.finalize()

    def testbmirun(self):
        bmiobj = bmi.wflowbmi_csdms()
        bmiobj.initialize('wflow_sceleton/wflow_sceleton.ini',loglevel=logging.DEBUG)
        print(bmiobj.get_var_type("IF"))
        et = bmiobj.get_end_time()
        st = bmiobj.get_start_time()
        ts = 86400
        # Do timesteps and chak

        bmiobj.update_until(et)
        bmiobj.get_current_time()
        bmiobj.finalize()
        print et - bmiobj.get_current_time()
        self.assertEquals(et,bmiobj.get_current_time())



    def testbmirun_l(self):
        print 'Run with update(-1)'
        bmiobj = bmi.wflowbmi_light()
        bmiobj.initialize('wflow_sceleton/wflow_sceleton.ini',loglevel=logging.ERROR)
        et = bmiobj.get_end_time()
        st = bmiobj.get_start_time()
        bmiobj.update(et - st)
        bmiobj.finalize()


    def testbmirun_space_in_name(self):
        print 'Run with update(-1)'
        bmiobj = bmi.wflowbmi_light()
        bmiobj.initialize('wflow sceleton/wflow sceleton.ini',loglevel=logging.ERROR)
        et = bmiobj.get_end_time()
        st = bmiobj.get_start_time()
        bmiobj.update(et - st)
        bmiobj.finalize()


    def testbmirunnetcdf(self):
        bmiobj = bmi.wflowbmi_csdms()
        bmiobj.initialize_config('wflow_sbm/wflow_sbm_nc.ini',loglevel=logging.ERROR)
        #bmiobj.set_attribute_value('run:starttime','2014-05-01 01:00:00')
        #bmiobj.set_attribute_value('run:endtime','2014-05-17 00:00:00')
        #bmiobj.set_attribute_value('run:runlengthdetermination','steps')


        bmiobj.set_start_time(1399597200)
        bmiobj.set_end_time(1399597200 + (4 * 3600))

        st = bmiobj.get_start_time()
        ett = bmiobj.get_end_time()
        ts = bmiobj.get_time_step()

        bmiobj.initialize_model()
        tt = bmiobj.get_value('timestepsecs')
        curtime = st
        cnt = 0
        lastcurtime = bmiobj.get_current_time()
        while curtime < ett:
            avar = bmiobj.get_value('PET')
            bmiobj.set_value('PET',avar + 10.0)
            cnt = cnt + 1
            bmiobj.update_until(curtime + ts)
            print (curtime + ts)/ts
            curtime = bmiobj.get_current_time()
            print bmiobj.get_current_time() - lastcurtime
            lastcurtime = bmiobj.get_current_time()


        bmiobj.finalize()


if __name__ == '__main__':
    unittest.main()
