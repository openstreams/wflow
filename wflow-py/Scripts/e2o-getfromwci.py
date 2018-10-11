#!/usr/bin/env python
# we need the netcdf4 library. This allows us to read netcdf from OpeNDAP or files
import netCDF4
import cftime
import osgeo.gdal as gdal
from osgeo.gdalconst import *

# Uncomment these if needed
from numpy import *
from scipy import *
from matplotlib import *
from pylab import *


# specify an url

"""
http://'http://wci.earth2observe.eu/thredds/dodsC/ecmwf/met_forcing_v0/1980/Tair_daily_E2OBS_198001.nc

"""

baseurl = "http://wci.earth2observe.eu/thredds/dodsC/ecmwf/met_forcing_v0/%d/Rainf_daily_E2OBS_%d%02d.nc"
# baseurl =  'http://wci.earth2observe.eu/thredds/dodsC/ecmwf/met_forcing_v0/%d/Tair_daily_E2OBS_%d%02d.nc'

months = arange(1, 13, 1)

years = arange(2010, 2013, 1)
clonemap = "wflow_dem.map"
lowresout = "lowres"
finalout = "inmaps"
mapstackname = "P"
# mapstackname = "TEMP"
ncdatafield = "Rainf"
# ncdatafield =  'Tair'


def writeMap(fileName, fileFormat, x, y, data, FillVal):
    """ Write geographical data into file"""

    verbose = False
    gdal.AllRegister()
    driver1 = gdal.GetDriverByName("GTiff")
    driver2 = gdal.GetDriverByName(fileFormat)

    # Processing
    if verbose:
        print "Writing to temporary file " + fileName + ".tif"
    # Create Output filename from (FEWS) product name and date and open for writing
    TempDataset = driver1.Create(
        fileName + ".tif", data.shape[1], data.shape[0], 1, gdal.GDT_Float32
    )
    # Give georeferences
    xul = x[0] - (x[1] - x[0]) / 2
    yul = y[0] + (y[0] - y[1]) / 2

    print xul
    print yul
    TempDataset.SetGeoTransform([xul, x[1] - x[0], 0, yul, 0, y[1] - y[0]])
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data, 0, 0)
    TempBand.FlushCache()
    TempBand.SetNoDataValue(FillVal)
    # Create data to write to correct format (supported by 'CreateCopy')
    if verbose:
        print "Writing to " + fileName + ".map"
    outDataset = driver2.CreateCopy(fileName, TempDataset, 0)
    TempDataset = None
    outDataset = None
    if verbose:
        print "Removing temporary file " + fileName + ".tif"
    os.remove(fileName + ".tif")

    if verbose:
        print "Writing to " + fileName + " is done!"


tot = []
cnt = 0
for year in years:
    for mon in months:
        rainurl = baseurl % (year, year, mon)
        print "processing: " + rainurl
        # create a dataset object
        ncdataset = netCDF4.Dataset(rainurl)
        lat = ncdataset.variables["lat"][:]
        lon = ncdataset.variables["lon"][:]
        ncdata = ncdataset.variables[ncdatafield]

        # Select lat and long for our cathcment
        # Bounding box for our catchment
        BB = dict(lon=[143, 150], lat=[-37, -33])

        (latidx,) = logical_and(lat >= BB["lat"][0], lat < BB["lat"][1]).nonzero()
        (lonidx,) = logical_and(lon >= BB["lon"][0], lon < BB["lon"][1]).nonzero()

        print lonidx
        print latidx
        print lat[latidx]
        print lon[lonidx]
        # get rid of the non used lat/lon now
        lat = lat[latidx]
        lon = lon[lonidx]
        # Now get the time for the x-axis
        time = ncdataset.variables["time"]
        timeObj = cftime.num2date(time[:], units=time.units, calendar=time.calendar)

        # Now determine area P for each timestep and display in a graph
        # first  the mean per area lat, next average those also
        # Multiply with timestep in seconds to get mm

        # unfortunateley Tair also has  heigh dimension and Precip not
        if mapstackname == "P":
            p_select = (
                ncdata[:, latidx.min() : latidx.max(), lonidx.min() : lonidx.max()]
                * 86400
            )
        if mapstackname == "TEMP":
            p_select = (
                ncdata[:, 0, latidx.min() : latidx.max(), lonidx.min() : lonidx.max()]
                - 273.15
            )
        # print p_select

        # PLot the sum over this month for the subcatchment

        Lon, Lat = meshgrid(lon, lat)
        # mesh = pcolormesh(Lon,Lat,p_select.sum(axis=0))
        # title("Cumulative precipitation")
        p_mean = p_select.mean(axis=1).mean(axis=1)
        print lon
        print lat
        ncdataset.close()

        if len(tot) == 0:
            tot = p_mean.copy()
        else:
            tot = hstack((tot, p_mean))

        arcnt = 0
        for a in timeObj:
            cnt = cnt + 1
            below_thousand = cnt % 1000
            above_thousand = cnt / 1000
            mapname = str(
                mapstackname + "%0" + str(8 - len(mapstackname)) + ".f.%03.f"
            ) % (above_thousand, below_thousand)
            print "saving map: " + os.path.join(lowresout, mapname)
            writeMap(
                os.path.join(lowresout, mapname),
                "PCRaster",
                lon,
                lat[::-1],
                flipud(p_select[arcnt, :, :]),
                -999.0,
            )
            arcnt = arcnt + 1
            execstr = (
                "resample --clone "
                + clonemap
                + " "
                + os.path.join(lowresout, mapname)
                + " "
                + os.path.join(finalout, mapname)
            )
            print "resampling map: " + execstr
            os.system(execstr)
