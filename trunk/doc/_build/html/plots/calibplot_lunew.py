# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 14:42:47 2012

@author: schelle
"""


import numpy
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import datetime
import pandas
import pcrut
from stats import *

Case = "rhine_sbm"
#Case = "rhineNew"
warmup = 500 # timesteps to skip in warmup phase for calculations
cooldown = 50

location = 7
obs,head=pcrut.readtss(Case + "/calib_new.tss")
obs = obs[1:3601,:]


    
for location in range(0,14):
    sim,hd=pcrut.readtss(Case + "/newlu_test/run.tss")
    #shift one day!!!!
    pers = numpy.size(obs,axis=0)
    
    a = get_nash_sutcliffe(obs[warmup:len(obs)-cooldown,location],sim[warmup:len(sim)-cooldown,location],NoData=numpy.nan)
    
    #trange = pandas.DatetimeIndex(datetime.datetime(1985,1,1),periods=pers,offset=pandas.DateOffset())
    ts = pandas.Series(obs[:,location],index=pandas.date_range('1/1/1985',periods=pers))
    tssim = pandas.Series(sim[:,location],index=pandas.date_range('1/1/1985',periods=pers))
    
    plt.figure(location)
    plt.autoscale(enable=True)
    ts.plot(label='Observed',color='blue')
    plt.autoscale(enable=False)
    tssim.plot(label='Simulated',color='black')
    plt.title(head[location] + ": NS = " + str(a[0])) 
    plt.legend()