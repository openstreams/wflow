#!/usr/bin/env python
# we need the netcdf4 library. This allows us to read netcdf from OpeNDAP or files
import netCDF4

# Uncomment these if needed
from numpy import *
from scipy import *
from matplotlib import *
from pylab import *

# specify an url

baseurl =  'http://wci.earth2observe.eu/thredds/dodsC/ecmwf/met_forcing_v0/%d/Rainf_daily_E2OBS_%d%02d.nc'

months = arange(1,13,1)
years= arange(1979,2013,1)


tot= []

for year in years:
    for mon in months:
        rainurl = baseurl % (year, year, mon)
        print "processing: " + rainurl    
        # create a dataset object
        raindataset = netCDF4.Dataset(rainurl)
        lat = raindataset.variables['lat'][:]
        lon = raindataset.variables['lon'][:]
        precip = raindataset.variables['Rainf']
        
        # Select lat and long for our cathcment
        # Bounding box for our catchment
        BB = dict(
           lon=[ 143, 150],
           lat=[-37, -33]
           )
         
        (latidx,) = logical_and(lat >= BB['lat'][0], lat < BB['lat'][1]).nonzero()
        (lonidx,) = logical_and(lon >= BB['lon'][0], lon < BB['lon'][1]).nonzero()
        
        # get rid of the non used lat/lon now
        lat = lat[latidx]
        lon = lon[lonidx]
        # Now get the time for the x-axis
        time = raindataset.variables['time']
        timeObj = netCDF4.num2date(time[:], units=time.units, calendar=time.calendar)
        
        #Now determine area P for each timestep and display in a graph
        # first  the mean per area lat, next average those also
        # Multiply with timestep in seconds to get mm
        p_select = precip[:,latidx.min():latidx.max(),lonidx.min():lonidx.max()]
        #print p_select
     
        # PLot the sum over this month for the subcatchment
        
        Lon,Lat = meshgrid(lon, lat)
        #mesh = pcolormesh(Lon,Lat,p_select.sum(axis=0))
        #title("Cumulative precipitation") 
        p_mean = p_select.mean(axis=1).mean(axis=1)
        print len(timeObj)
        print p_mean
        raindataset.close()
        
        if len(tot) == 0:
            tot = p_mean.copy()
        else:
            tot = hstack((tot,p_mean))
    #f, ax = plt.subplots(figsize=(16, 6))
    #ax.plot(timeObj,p_mean)
    #title("Precipitation")
    
plot(tot * 86400)    
