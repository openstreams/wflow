#!/usr/bin/python

# pcr2netcdf is Free software, see below:
#
# Copyright (c) J. Schellekens 2014
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
syntax:
    pcr2netcdf -S date -E date - N mapstackname -I mapstack_folder 
               -O netcdf_name [-b buffersize] [-c inifile]

    -S startdate in "%d-%m-%Y %H:%M:%S" e.g. 31-12-1990 00:00:00
    -E endDate in "%d-%m-%Y %H:%M:%S"
    -N Mapstack-name (prefix)
       You can sepcify multiple input mapstack  to merge them into one netcdf
       e.g. -M P -M TEMP -M PET
    -I input mapstack folder
    -O output netcdf file
    -b maxbuf - maximum number of timesteps to buffer before writing (default = 600)
    -t timestep - (set timestep in seconds, default = 86400) Only 86400 and 3600 supported
    -F FORMAT (default = NETCDF4, use NETCDF3_CLASSIC for OpenDA)
    -z switch on zlib compression (default=Off)
    -c inifile (contains netcdf meta data)

This utility is made to simplify running wflow models with OpenDA. The
OpenDA link needs the forcing timeseries to be in netcdf format. Use this to convert
all the input mapstacks to netcdf.


(c) J. Schellekens

Created on Tue May 13 07:37:04 2014

based on GLOFRIS_utils.py by H Winsemius

"""

    
import time
import datetime as dt
import getopt
import sys
from numpy import *
import netCDF4 as nc4
import osgeo.gdal as gdal
import os
import logging
import logging.handlers
import ConfigParser
import wflow.pcrut as _pcrut

def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)


def readMap(fileName, fileFormat,logger):
    """ 
    Read geographical file into memory
    """
    import osgeo.gdal as gdal 
    # Open file for binary-reading
    mapFormat = gdal.GetDriverByName(fileFormat)
    mapFormat.Register()
    ds = gdal.Open(fileName)
    if ds is None:
        logger.error('Could not open ' + fileName + '. Something went wrong!! Shutting down')
        sys.exit(1)
        # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX    = geotrans[1]
    resY    = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = linspace(originX+resX/2,originX+resX/2+resX*(cols-1),cols)
    y = linspace(originY+resY/2,originY+resY/2+resY*(rows-1),rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1) # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0,0,cols,rows)
    FillVal = RasterBand.GetNoDataValue()
    RasterBand = None
    ds = None
    return x, y, data, FillVal
    

def getnetcdfmetafromini(inifile):
    """
    Gets a netcdf mete data dictionary from an ini file

    :param inifile: inifile with a metadata section
    :return : dictionary with meta data
    """
    metadata = {}

    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    if os.path.exists(inifile):
        config.read(inifile)
    else:
        print ("Cannot open ini file: " +  inifile)
        exit(1)

    metadata = dict(config.items('metadata'))

    return metadata


def getvarmetadatafromini(inifile,var):
    """

    :param inifile: inifile
    :param var: Name of the mapstack
    :return: dictionary
    """
    metadata = {}

    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    if os.path.exists(inifile):
        config.read(inifile)
    else:
        print ("Cannot open ini file: " +  inifile)
        exit(1)

    metadata = dict(config.items(var))

    return metadata


def write_netcdf_timeseries(srcFolder, srcPrefix, trgFile, trgVar, trgUnits, trgName, timeList, metadata, logger,maxbuf=600,Format="NETCDF4",zlib=False,startidx=0):
    """
    Write pcraster mapstack to netcdf file. Taken from GLOFRIS_Utils.py
    
    - srcFolder - Folder with pcraster mapstack
    - srcPrefix - name of the mapstack
    - trgFile - target netcdf file
    - tgrVar - variable in nc file
    - trgUnits - units for the netcdf file
    - timeLists - list of times
    - metadata - dict with metedata for var

    Optional argumenrs
    - maxbuf = 600: number of timesteps to buffer before writing
    
    """

    # if necessary, make trgPrefix maximum of 8 characters    
    if len(srcPrefix) > 8:
        srcPrefix = srcPrefix[0:8]
    # Open target netCDF file
    nc_trg = nc4.Dataset(trgFile, 'a',format=Format,zlib=zlib)
    # read time axis and convert to time objects
    logger.debug("Creating time object..")
    time = nc_trg.variables['time']
    nc_trg.set_fill_off()
    timeObj = nc4.num2date(time[:], units=time.units, calendar=time.calendar)

    
    try:
        nc_var = nc_trg.variables[trgVar]
    except:
        # prepare the variable
        nc_var = nc_trg.createVariable(trgVar, 'f4', ('time', 'lat', 'lon',), fill_value=-9999., zlib=zlib)
        nc_var.units = trgUnits
        nc_var.standard_name = trgName
        print metadata
        for attr in metadata:
            #print metadata[attr]
            nc_var.setncattr(attr, metadata[attr])

    
    nc_Fill = nc_var._FillValue
    # Create a buffer of a number of timesteps to speed-up writing    
    bufsize = minimum(len(timeList),maxbuf) 
    timestepbuffer = zeros((bufsize,len(nc_trg.variables['lat']),len(nc_trg.variables['lon'])))

    # now loop over all time steps, check the date and write valid dates to a list, write time series to PCRaster maps
    for nn, curTime in enumerate(timeList):
        logger.debug("Adding time: " + str(curTime))
        idx = where(timeObj==curTime)[0]
        count = nn + startidx

        below_thousand = count % 1000
        above_thousand = count / 1000
        # read the file of interest
        pcraster_file  = str(srcPrefix + '%0' + str(8-len(srcPrefix)) + '.f.%03.f') % (above_thousand, below_thousand)
        pcraster_path = os.path.join(srcFolder, pcraster_file)
        # write grid to PCRaster file
        logger.debug("processing map: " + pcraster_file)
        x, y, data, FillVal = readMap(pcraster_path, 'PCRaster',logger)
        logger.debug("Setting fillval...")
        data[data==FillVal] = nc_Fill
        
        buffreset = (idx + 1) % maxbuf
        bufpos = (idx) % maxbuf
        #timestepbuffer[bufpos,:,:] =  flipud(data)
        # Weird, the flupud is no longer needed!!!!
        logger.debug("Adding data to array...")
        timestepbuffer[bufpos,:,:] =  data

        if buffreset == 0 or idx ==  bufsize -1 or nn + 1 == len(timeList):
            logger.info("Writing buffer to file at: " + str(curTime) + " " + str(int(bufpos) + 1) + " timesteps")
            nc_var[idx-bufsize+1:idx+1,:,:] = timestepbuffer
            nc_trg.sync()


    nc_trg.sync()
    nc_trg.close()

    
    
def prepare_nc(trgFile, timeList, x, y, metadata, logger, units='Days since 1900-01-01 00:00:00', calendar='gregorian',Format="NETCDF4",zlib=False):
    """
    This function prepares a NetCDF file with given metadata, for a certain year, daily basis data
    The function assumes a gregorian calendar and a time unit 'Days since 1900-01-01 00:00:00'
    """
    import datetime as dt
    
    logger.info('Setting up "' + trgFile + '"')
    startDayNr = nc4.date2num(timeList[0], units=units, calendar=calendar)
    endDayNr   = nc4.date2num(timeList[-1], units=units, calendar=calendar)
    time       = arange(startDayNr,endDayNr+1)
    nc_trg     = nc4.Dataset(trgFile,'w',format=Format,zlib=zlib)

    logger.info('Setting up dimensions and attributes. lat: ' + str(len(y))+ " lon: " + str(len(x)))
    nc_trg.createDimension('time', 0) #NrOfDays*8
    nc_trg.createDimension('lat', len(y))
    nc_trg.createDimension('lon', len(x))
    DateHour = nc_trg.createVariable('time','f8',('time',))
    DateHour.units = units
    DateHour.calendar = calendar
    DateHour.standard_name = 'time'
    DateHour.long_name = 'time'
    DateHour.axis = 'T'
    DateHour[:] = time
    y_var = nc_trg.createVariable('lat','f4',('lat',))
    y_var.standard_name = 'latitude'
    y_var.long_name = 'latitude'
    y_var.units = 'degrees_north'
    y_var.axis = 'Y'
    x_var = nc_trg.createVariable('lon','f4',('lon',))
    x_var.standard_name = 'longitude'
    x_var.long_name = 'longitude'
    x_var.units = 'degrees_east'
    x_var.axis = 'X'
    y_var[:] = y
    x_var[:] = x
    projection= nc_trg.createVariable('projection','c')
    projection.long_name = 'wgs84'
    projection.EPSG_code = 'EPSG:4326'
    projection.proj4_params = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    projection.grid_mapping_name = 'latitude_longitude'

    # now add all attributes from user-defined metadata
    for attr in metadata:
        nc_trg.setncattr(attr, metadata[attr])
    nc_trg.sync()
    nc_trg.close()


def setlogger(logfilename,loggername, thelevel=logging.INFO):
    """
    Set-up the logging system and return a logger object. Exit if this fails
    """

    try:
        #create logger
        logger = logging.getLogger(loggername)
        if not isinstance(thelevel, int):
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(thelevel)
        ch = logging.FileHandler(logfilename,mode='w')
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        #create formatter
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        #add formatter to ch
        ch.setFormatter(formatter)
        console.setFormatter(formatter)
        #add ch to logger
        logger.addHandler(ch)
        logger.addHandler(console)
        logger.debug("File logging to " + logfilename)
        return logger
    except IOError:
        print "ERROR: Failed to initialize logger with logfile: " + logfilename
        sys.exit(2)

def date_range(start, end, tdelta="days"):



    if tdelta == "days":
        r = (end+dt.timedelta(days=1)-start).days
        return [start+dt.timedelta(days=i) for i in range(r)]
    else:
        r = (end+dt.timedelta(days=1)-start).days * 24
        return [start+dt.timedelta(hours=i) for i in range(r)]

def date_range_peryear(start, end, tdelta="days"):

    ret = []
    for yrs in range(start.year,end.year+1):
        if tdelta == "days":
            r = (dt.datetime(yrs+1,1,1)-dt.datetime(yrs,1,1)).days
            ret.append([dt.datetime(yrs,1,1)+dt.timedelta(days=i) for i in range(r)])
        else:
            r = (dt.datetime(yrs+1,1,1)-dt.datetime(yrs,1,1)).days * 24
            ret.append([dt.datetime(yrs,1,1)+dt.timedelta(hours=i) for i in range(r)])

    return ret

def main(argv=None):
    """
    Perform command line execution of the model.
    """
    # initiate metadata entries
    metadata = {}
    metadata['title'] = 'wflow input mapstack'
    metadata['institution'] = 'Deltares'
    metadata['source'] = 'pcr2netcdf'
    metadata['history'] = time.ctime()
    metadata['references'] = 'http://wflow.googlecode.com'
    metadata['Conventions'] = 'CF-1.4'
     
    ncoutfile = "inmaps.nc"
    mapstackfolder="inmaps"
    inifile = None
    mapstackname=[]
    var=[]
    varname=[]
    unit="mm"
    startstr="1-1-1990 00:00:00"
    endstr="2-2-1990 00:00:00"
    mbuf=600
    timestepsecs = 86400
    
    clonemap=None
    Format="NETCDF4"
    zlib=False
    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return    

    ## Main model starts here
    ########################################################################
    try:
        opts, args = getopt.getopt(argv, 'c:S:E:N:I:O:b:t:F:z')
    except getopt.error, msg:
        usage(msg)

    for o, a in opts:
        if o == '-S': startstr = a
        if o == '-E': endstr = a
        if o == '-O': ncoutfile = a
        if o == '-c': inifile = a
        if o == '-I': mapstackfolder = a
        if o == '-b': mbuf = int(a)
        if o == '-z': zlib=True
        if o == '-F': Format=a
        if o == '-t': 
            timestepsecs = int(a)
        if o == '-N': 
            mapstackname.append(a)
            var.append(a)
            varname.append(a)


    # Use first timestep as clone-map
    logger = setlogger('pcr2netcdf.log','pcr2netcdf', thelevel = logging.DEBUG)

    count = 1
    below_thousand = count % 1000
    above_thousand = count / 1000
    clonemapname  = str(mapstackname[0] + '%0' + str(8-len(mapstackname[0])) + '.f.%03.f') % (above_thousand, below_thousand)
    clonemap = os.path.join(mapstackfolder, clonemapname)
    _pcrut.setclone(clonemap)
   
    x = _pcrut.pcr2numpy(_pcrut.xcoordinate(_pcrut.boolean(_pcrut.cover(1.0))),NaN)[0,:]
    y = _pcrut.pcr2numpy(_pcrut.ycoordinate(_pcrut.boolean(_pcrut.cover(1.0))),NaN)[:,0]
    

    start=dt.datetime.strptime(startstr,"%d-%m-%Y %H:%M:%S")
    end=dt.datetime.strptime(endstr,"%d-%m-%Y %H:%M:%S")
    if timestepsecs == 86400:
        timeList = date_range_peryear(start, end, tdelta="days")
    else:   
        timeList = date_range_peryear(start, end, tdelta="hours")

    if inifile is not None:
        inimetadata = getnetcdfmetafromini(inifile)
        metadata.update(inimetadata)


    # break up into separate years

    startmapstack = 1
    for yr_timelist in timeList:
        ncoutfile_yr = os.path.splitext(ncoutfile)[0] + "_" + str(yr_timelist[0].year) + os.path.splitext(ncoutfile)[1]
        prepare_nc(ncoutfile_yr, yr_timelist, x, y, metadata, logger,Format=Format,zlib=zlib)

        idx = 0
        for mname in mapstackname:
            logger.info("Converting mapstack: " + mname + " to " + ncoutfile)
            # get variable attributes from ini file here
            varmeta = getvarmetadatafromini(inifile,var[idx])

            write_netcdf_timeseries(mapstackfolder, mname, ncoutfile_yr, var[idx], unit, varname[idx], yr_timelist, varmeta, logger,maxbuf=mbuf,Format=Format,zlib=zlib,startidx=startmapstack)
            startmapstack = startmapstack + len(yr_timelist)
            idx = idx + 1


if __name__ == "__main__":
    main()
