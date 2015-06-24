__author__ = 'schelle'

import unittest
import logging
import wflow.wflow_bmi as bmi
import os
"""
Simple test for wflow bmi framework
"""

class MyTest(unittest.TestCase):

    def testbmirun(self):
        bmiobj = bmi.wflowbmi_csdms()
        bmiobj.initialize('wflow_sceleton/wflow_sceleton.ini',loglevel=logging.ERROR)

        et = bmiobj.get_end_time()
        bmiobj.update_until(et)
        bmiobj.finalize()

    def testbmifuncs(self):

        bmiobj = bmi.wflowbmi_csdms()
        bmiobj.initialize('wflow_sceleton/wflow_sceleton.ini',loglevel=logging.ERROR)

        print("-------------- Grid origin: ")
        gorigin = bmiobj.get_grid_origin('Altitude')
        print(gorigin)

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

        print("-------------- Output var names: ")
        print(bmiobj.get_output_var_names())

        print("-------------- Time units: ")
        print(bmiobj.get_time_units())

        print("-------------- Time step: ")
        print(bmiobj.get_time_step())

        print("-------------- Current time: ")
        print(bmiobj.get_current_time())

        print("-------------- Grid shape: ")
        print(bmiobj.get_grid_shape('Altitude'))

        print("-------------- End time: ")
        print(bmiobj.get_end_time())

        bmiobj.finalize()




if __name__ == '__main__':
    unittest.main()
