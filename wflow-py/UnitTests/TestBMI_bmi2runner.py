__author__ = 'schelle'

import unittest
import logging
import sys
sys.path = ['../wflow'] + ['../'] + sys.path
sys.path = ['../Scripts'] + sys.path
import bmi2runner as bmirun
import time
import os

"""
Simple test for wflow bmi framework
"""


class MyTest(unittest.TestCase):

    def testbmifuncs(self):

        bmirun.main(['-c','combined/bmirunner.ini'])




if __name__ == '__main__':
    unittest.main()
