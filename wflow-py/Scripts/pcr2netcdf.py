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
    For mapstacks:

    pcr2netcdf -S date -E date -N mapstackname -I mapstack_folder
               -O netcdf_name [-b buffersize] [-c inifile][-s start][-d digit][-Y][-P EPSG]
               [-i inputformat][-F netcdfformat]

    For single maps (no series)
    pcr2netcdf -M -S date  -N mapname
               -O netcdf_name [-b buffersize] [-c inifile][-s start][-d digit][-Y][-P EPSG][-C clone]

    -M single map mode
    -S startdate in "%d-%m-%Y %H:%M:%S" e.g. 31-12-1990 00:00:00
    -E endDate in "%d-%m-%Y %H:%M:%S"
    -s startstep (in the mapstack, default = 1)
    -Y Make seperate files per year
    -i GDAL Input file format string (default: PCRaster)
    -P set the EPSG string. default: "EPSG:4326"
    -N Mapstack-name (prefix)
       You can specify multiple input mapstack  to merge them into one netcdf
       e.g. -N P -N TEMP -N PET
       In combination with the -M option you can use wildcards

    -I input mapstack folder
    -O output netcdf file
    -b maxbuf - maximum number of timesteps to buffer before writing (default = 600)
    -t timestep - (set timestep in seconds, default = 86400) Only 86400 and 3600 supported
    -F FORMAT (default = NETCDF4, use NETCDF3_CLASSIC for OpenDA)
    -z switch on zlib compression (default=On)
    -d least_significant_digit (set precision for more compression)
    -c inifile (contains netcdf meta data)
    -C clonemap.map use maps as mask and clone

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
import osgeo.gdal as gdal
import wflow.wf_netcdfio as ncdf
import glob

def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)

def writeMap(fileName, fileFormat, x, y, data, FillVal):
    """
    Write geographical data into file. Also replace NaN by FillVall

    :param fileName:
    :param fileFormat:
    :param x:
    :param y:
    :param data:
    :param FillVal:
    :return:
    """


    verbose = False
    gdal.AllRegister()
    driver1 = gdal.GetDriverByName('GTiff')
    driver2 = gdal.GetDriverByName(fileFormat)

    data[isnan(data)] = FillVal
    # Processing
    if verbose:
        print 'Writing to temporary file ' + fileName + '.tif'
        print "Output format: " + fileFormat
    # Create Output filename from (FEWS) product name and date and open for writing
    TempDataset = driver1.Create(fileName + '.tif',data.shape[1],data.shape[0],1,gdal.GDT_Float32)
    # Give georeferences
    xul = x[0]-(x[1]-x[0])/2
    yul = y[0]+(y[0]-y[1])/2

    TempDataset.SetGeoTransform( [ xul, x[1]-x[0], 0, yul, 0, y[1]-y[0] ] )
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data.astype(float32),0,0)
    TempBand.FlushCache()
    TempBand.SetNoDataValue(FillVal)

    # Create data to write to correct format (supported by 'CreateCopy')
    if verbose:
        print 'Writing to ' + fileName + '.map'
    if fileFormat == 'GTiff':
        outDataset = driver2.CreateCopy(fileName, TempDataset, 0 ,options = ['COMPRESS=LZW'])
    else:
        outDataset = driver2.CreateCopy(fileName, TempDataset, 0)
    TempDataset = None
    outDataset = None
    if verbose:
        print 'Removing temporary file ' + fileName + '.tif'
    os.remove(fileName + '.tif');

    if verbose:
        print 'Writing to ' + fileName + ' is done!'

def readMap(fileName, fileFormat,logger,unzipcmd='pigz -d -k'):
    """ 
    Read PCRaster geographical file into memory
    """
    unzipped = 0
    if not os.path.exists(fileName):
        # try and unzip
        if os.path.exists(fileName + ".gz"):
            os.system(unzipcmd + ' ' + fileName + ".gz")
            logger.info("unzipping: " + fileName + ".gz")
            unzipped = 1

    pcrdata =  _pcrut.readmap(fileName)
    x = _pcrut.pcr2numpy(_pcrut.xcoordinate(_pcrut.boolean(_pcrut.cover(1.0))),NaN)[0,:]
    y = _pcrut.pcr2numpy(_pcrut.ycoordinate(_pcrut.boolean(_pcrut.cover(1.0))),NaN)[:,0]

    FillVal =  float(1E31)
    data = _pcrut.pcr2numpy(pcrdata,FillVal)
    if unzipped:
        #Delete uncompressed file if compressed exsists
        if os.path.exists(fileName + ".gz"):
            logger.info("Removing: " + fileName)
            os.remove(fileName)


    return x, y, data, FillVal
    

def _readMap(fileName, fileFormat,logger):
    """
    Read geographical file into memory
    """

    #Open file for binary-reading

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
    #if resY < 0.0:
    #    y = linspace(originY+abs(resY)/2,originY+abs(resY)/2+abs(resY)*(rows-1),rows)[::-1]
    #    y = linspace(originY + resY / 2, originY + resY / 2 + resY * (rows - 1), rows)
    #else:
    y = linspace(originY + resY / 2, originY + resY / 2 + resY * (rows - 1), rows)
    RasterBand = ds.GetRasterBand(1) # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0,0,cols,rows)
    FillVal = RasterBand.GetNoDataValue()
    RasterBand = None
    del ds
    # #ds = None

    return x, y, data.copy(), FillVal






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


def write_netcdf_timeseries(srcFolder, srcPrefix, trgFile, trgVar, trgUnits, trgName, timeList, metadata,
                            logger,clone,maxbuf=600,Format="NETCDF4",zlib=True,least_significant_digit=None,startidx=0,EPSG="EPSG:4326",
                            FillVal=1E31):
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
    complevel=9

    # if necessary, make trgPrefix maximum of 8 characters    
    if len(srcPrefix) > 8:
        srcPrefix = srcPrefix[0:8]
    # Open target netCDF file
    nc_trg = nc4.Dataset(trgFile, 'a',format=Format)
    # read time axis and convert to time objects
    logger.debug("Creating time object..")
    time = nc_trg.variables['time']
    nc_trg.set_fill_off()
    timeObj = nc4.num2date(time[:], units=time.units, calendar=time.calendar)

    
    try:
        nc_var = nc_trg.variables[trgVar]
    except:
        # prepare the variable

        if EPSG.lower() == "epsg:4326":
            nc_var = nc_trg.createVariable(trgVar, 'f4', ('time', 'lat', 'lon',), fill_value=FillVal, zlib=zlib,
                                                        complevel=complevel, least_significant_digit=least_significant_digit)
            nc_var.coordinates = "lat lon"
        else:
            nc_var = nc_trg.createVariable(trgVar, 'f4', ('time', 'y', 'x',), fill_value=FillVal, zlib=zlib,
                                                        complevel=complevel, least_significant_digit=least_significant_digit)
            nc_var.coordinates = "lat lon"
            nc_var.grid_mapping = "crs"

        nc_var.units = trgUnits
        nc_var.standard_name = trgName
        #print metadata
        for attr in metadata:
            #print metadata[attr]
            nc_var.setncattr(attr, metadata[attr])

    
    nc_Fill = nc_var._FillValue
    # Create a buffer of a number of timesteps to speed-up writing    
    bufsize = minimum(len(timeList),maxbuf)

    if len(shape(nc_trg.variables['lat'])) == 2:
        latlen = shape(nc_trg.variables['lat'])[0]
        lonlen = shape(nc_trg.variables['lon'])[1]
    else:
        latlen = len(nc_trg.variables['lat'])
        lonlen = len(nc_trg.variables['lon'])

    timestepbuffer = zeros((bufsize,latlen,lonlen))

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
        #x, y, data, FillVal = readMap(pcraster_path, 'PCRaster',logger)
        x, y, data, FFillVal = _readMap(pcraster_path, 'PCRaster',logger)
        logger.debug("Setting fillval...")


        data[data==FFillVal] = float(nc_Fill)
        data[isinf(data)] = float(nc_Fill)
        data[isnan(data)] = float(nc_Fill)
        data[clone <= -999] = float(nc_Fill)
        data[clone == FFillVal] = float(nc_Fill)

        buffreset = (idx + 1) % maxbuf
        bufpos = (idx) % maxbuf
        logger.debug("Adding data to array...")
        timestepbuffer[bufpos,:,:] =  data

        logger.debug("index: " + str(idx-bufpos) + " index: " + str(idx) + "bufpos: " + str(bufpos) + "idx: " + str(idx))
        if buffreset == 0 or idx ==  bufsize -1 or nn + 1 == len(timeList):
            logger.info("Writing buffer to file at: " + str(curTime) + " " + str(int(bufpos) + 1) + " timesteps")
            nc_var[idx-bufpos:idx+1,:,:] = timestepbuffer[0:bufpos+1,:,:]
            nc_trg.sync()


    #nc_trg.sync()
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

# def date_range(start, end, tdelta="days"):
#
#
#
#     if tdelta == "days":
#         r = (end+dt.timedelta(days=1)-start).days
#         return [start+dt.timedelta(days=i) for i in range(r)]
#     else:
#         r = (end+dt.timedelta(days=1)-start).days * 24
#         return [start+dt.timedelta(hours=i) for i in range(r)]

def date_range_peryear(start, end, tdelta="days"):

    ret = []
    for yrs in range(start.year,end.year):
        ed = min(dt.datetime(yrs+1,1,1), end)

        if tdelta == "days":
            r = (ed-dt.datetime(yrs,1,1)).days
            ret.append([dt.datetime(yrs,1,1)+dt.timedelta(days=i) for i in range(r)])
        else:
            r = ((ed-dt.datetime(yrs,1,1)).days) * 24
            ret.append([dt.datetime(yrs,1,1)+dt.timedelta(hours=i) for i in range(r)])


    return ret

def date_range(start, end, timestepsecs):
        r = int((end + dt.timedelta(seconds=timestepsecs) - start).total_seconds()/timestepsecs)
        return [start + dt.timedelta(seconds=(timestepsecs * i)) for i in range(r)]


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
    inifile = "not set"
    mapstackname=[]
    var=[]
    varname=[]
    unit="mm"
    startstr="1-1-1990 00:00:00"
    endstr="2-2-1990 00:00:00"
    mbuf=600
    timestepsecs = 86400

    outputFillVal = 1E31
    clonemap=None
    OFormat="NETCDF4"
    IFormat = 'PCRaster'
    EPSG="EPSG:4326"
    Singlemap = False
    zlib=True
    least_significant_digit=None
    clonemapname = 'None'
    startstep = 1
    perYear = False
    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return    


    ## Main model starts here
    ########################################################################
    try:
        opts, args = getopt.getopt(argv, 'c:S:E:N:I:O:b:t:F:zs:d:YP:Mi:C:')
    except getopt.error, msg:
        usage(msg)

    for o, a in opts:
        if o == '-S': startstr = a
        if o == '-s':
            startstep = int(a)
        if o == '-E': endstr = a
        if o == '-i': IFormat = a
        if o == '-O': ncoutfile = a
        if o == '-c': inifile = a
        if o == '-I': mapstackfolder = a
        if o == '-b': mbuf = int(a)
        if o == '-Y': perYear = True
        if o == '-z': zlib=True
        if o == '-P': EPSG = a
        if o == '-M': Singlemap = True
        if o == '-F': OFormat=a
        if o == '-d': least_significant_digit = int(a)
        if o == '-C': clonemapname = a
        if o == '-t': 
            timestepsecs = int(a)
        if o == '-N':
            flst = glob.glob(a)
            if len(flst) == 0:
                mapstackname.append(a)
                var.append(a)
                varname.append(a)
            else:
                mapstackname = flst
                var =flst
                varname =flst

    # Use first timestep as clone-map
    logger = setlogger('pcr2netcdf.log','pcr2netcdf', thelevel = logging.DEBUG)

    count = 1
    below_thousand = count % 1000
    above_thousand = count / 1000

    if clonemapname == 'None':
        clonemapname = str(mapstackname[0] + '%0' + str(8 - len(mapstackname[0])) + '.f.%03.f') % (above_thousand, below_thousand)

    clonemap = os.path.join(mapstackfolder, clonemapname)


    if Singlemap:
        clonemap = mapstackname[0]


    if IFormat == 'PCRaster':
        _pcrut.setclone(clonemap)

    x, y, clone, FillVal = _readMap(clonemap, IFormat, logger)



    start=dt.datetime.strptime(startstr,"%d-%m-%Y %H:%M:%S")

    if Singlemap:
        end = start
    else:
        end=dt.datetime.strptime(endstr,"%d-%m-%Y %H:%M:%S")

    if timestepsecs == 86400:
        if perYear:
            timeList = date_range_peryear(start, end, tdelta="days")
        else:
            timeList = date_range(start, end, timestepsecs)
    else:
        if perYear:
            timeList = date_range_peryear(start, end, tdelta="hours")
        else:
            timeList = date_range(start, end,timestepsecs)

    if os.path.exists(inifile):
        inimetadata = getnetcdfmetafromini(inifile)
        metadata.update(inimetadata)


    # break up into separate years


    if not Singlemap:
        varmeta = {}

        startmapstack = startstep

        if perYear:
            for yr_timelist in timeList:

                ncoutfile_yr = os.path.splitext(ncoutfile)[0] + "_" + str(yr_timelist[0].year) + os.path.splitext(ncoutfile)[1]

                if os.path.exists(ncoutfile_yr):
                    logger.info("Skipping file: " + ncoutfile_yr)
                else:
                    ncdf.prepare_nc(ncoutfile_yr, yr_timelist, x, y, metadata, logger,Format=OFormat,zlib=zlib,
                                    EPSG=EPSG,FillValue=outputFillVal)

                    idx = 0
                    for mname in mapstackname:
                        logger.info("Converting mapstack: " + mname + " to " + ncoutfile)
                        # get variable attributes from ini file here
                        if os.path.exists(inifile):
                            varmeta = getvarmetadatafromini(inifile,var[idx])

                        write_netcdf_timeseries(mapstackfolder, mname, ncoutfile_yr, var[idx], unit, varname[idx], \
                                                yr_timelist, varmeta, logger, clone,maxbuf=mbuf,Format=OFormat,
                                                zlib=zlib,least_significant_digit=least_significant_digit,
                                                startidx=startmapstack,EPSG=EPSG,FillVal=outputFillVal)
                        idx = idx + 1

                logger.info("Old stack: " + str(startmapstack) + " new startpoint " + str(startmapstack + len(yr_timelist) -1))
                startmapstack = startmapstack + len(yr_timelist)

        else:
             #ncoutfile_yr = os.path.splitext(ncoutfile)[0] + "_" + str(yr_timelist[0].year) + os.path.splitext(ncoutfile)[1]
             ncdf.prepare_nc(ncoutfile, timeList, x, y, metadata, logger,Format=OFormat,zlib=zlib,EPSG=EPSG)
             idx = 0
             for mname in mapstackname:
                logger.info("Converting mapstack: " + mname + " to " + ncoutfile)
                # get variable attributes from ini file here
                if os.path.exists(inifile):
                    varmeta = getvarmetadatafromini(inifile,var[idx])

                write_netcdf_timeseries(mapstackfolder, mname, ncoutfile, var[idx], unit, varname[idx], timeList, varmeta,\
                                        logger,clone,maxbuf=mbuf,Format=OFormat,zlib=zlib,least_significant_digit=least_significant_digit,\
                                        startidx=startmapstack,EPSG=EPSG,FillVal=outputFillVal)
                idx = idx + 1
    else:
        NcOutput = ncdf.netcdfoutputstatic(ncoutfile, logger, timeList[0],1,timestepsecs=timestepsecs,
                                                     maxbuf=1, metadata=metadata, EPSG=EPSG,Format=OFormat,
                                                     zlib=zlib)

        for file in mapstackname:
            pcrdata = _pcrut.readmap(file)
            thevar = os.path.basename(file)
            NcOutput.savetimestep(1, pcrdata, unit="mm", var=thevar, name=file)







if __name__ == "__main__":
    main()
