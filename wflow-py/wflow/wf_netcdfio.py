"""
wf_netcdfio
-----------

netcdf reading and writing for wflow

$Author: schelle $
$Id: wf_DynamicFramework.py 915 2014-02-10 07:33:56Z schelle $
$Rev: 915 $
"""


import netCDF4
# the two below are needed fpr bbfreeze
import netCDF4.utils
import netcdftime

from pcraster import *
from numpy import *

import time
import datetime as dt


import wflow.wflow_lib as wflow_lib
import wflow.pcrut as _pcrut

globmetadata = {}
globmetadata['title'] = 'wflow output mapstack'
globmetadata['institution'] = 'Deltares'
globmetadata['source'] = 'wflow'
globmetadata['history'] = time.ctime()
globmetadata['references'] = 'http://wflow.googlecode.com'
globmetadata['Conventions'] = 'CF-1.4'
# 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', or 'NETCDF3_CLASSIC'
netcdfformat = "NETCDF4"

def prepare_nc(trgFile, timeList, x, y, metadata, logger, units='Days since 1900-01-01 00:00:00', calendar='gregorian',Format="NETCDF4",complevel=1,zlib=True):
    """
    This function prepares a NetCDF file with given metadata, for a certain year, daily basis data
    The function assumes a gregorian calendar and a time unit 'Days since 1900-01-01 00:00:00'
    """
    import datetime as dt

    logger.info('Setting up netcdf output: ' + trgFile)
    print timeList[0]
    startDayNr = netCDF4.date2num(timeList[0].replace(tzinfo=None), units=units, calendar=calendar)
    endDayNr   = netCDF4.date2num(timeList[-1].replace(tzinfo=None), units=units, calendar=calendar)
    time       = arange(startDayNr,endDayNr+1)
    nc_trg     = netCDF4.Dataset(trgFile,'w',format=Format,zlib=zlib)

    logger.info('Setting up dimensions and attributes.Steps: ' + str(len(timeList)) + ' lat: ' + str(len(y))+ " lon: " + str(len(x)))
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



class netcdfoutput():

    def __init__(self,netcdffile,logger,starttime,timesteps,timestepsecs=86400,metadata={},maxbuf=25):
        """
        Under construction
        """

        def date_range(start, end, tdelta="days"):
            if tdelta == "days":
                r = (end+dt.timedelta(days=1)-start).days
                return [start+dt.timedelta(days=i) for i in range(r)]
            else:
                r = (end+dt.timedelta(days=1)-start).days * 24
                return [start+dt.timedelta(hours=i) for i in range(r)]

        self.logger = logger
        # Do not allow a max buffer larger than the number of timesteps
        self.maxbuf =  maxbuf if timesteps >= maxbuf else timesteps
        self.ncfile = netcdffile
        self.timesteps = timesteps
        rows = pcraster._pcraster.clone().nrRows()
        cols = pcraster._pcraster.clone().nrCols()
        cellsize = pcraster._pcraster.clone().cellSize()
        yupper = pcraster._pcraster.clone().north()
        xupper = pcraster._pcraster.clone().west()
        x = _pcrut.pcr2numpy(_pcrut.xcoordinate(_pcrut.boolean(_pcrut.cover(1.0))),NaN)[0,:]
        y = _pcrut.pcr2numpy(_pcrut.ycoordinate(_pcrut.boolean(_pcrut.cover(1.0))),NaN)[:,0]

        end=starttime + dt.timedelta(seconds=timestepsecs * self.timesteps-1)
        if timestepsecs == 86400:
            timeList = date_range(starttime, end, tdelta="days")
        else:
            timeList = date_range(starttime, end, tdelta="hours")

        self.timestepbuffer = zeros((self.maxbuf,len(y),len(x)))
        self.bufflst={}

        globmetadata.update(metadata)

        prepare_nc(self.ncfile,timeList,x,y,globmetadata,logger,Format=netcdfformat)


    def savetimestep(self,timestep,pcrdata,unit="mm",var='P',name="Precipitation"):
        """
        save a single timestep for a variable

        input:
            - timestep - current timestep
            - pcrdata - pcraster map to save
            - unit - unit string
            - var - variable string
            - name - name of the variable
        """
        # Open target netCDF file
        var = os.path.basename(var)
        self.nc_trg = netCDF4.Dataset(self.ncfile, 'a',format=netcdfformat,zlib=True,complevel=1)
        self.nc_trg.set_fill_off()
        # read time axis and convert to time objects
        time = self.nc_trg.variables['time']
        timeObj = netCDF4.num2date(time[:], units=time.units, calendar=time.calendar)

        idx = timestep -1

        buffreset = (idx + 1) % self.maxbuf
        bufpos = (idx) % self.maxbuf

        try:
             nc_var = self.nc_trg.variables[var]
        except:
            self.logger.debug("Creating variable " + var + " in netcdf file. Format: " + netcdfformat)
            nc_var = self.nc_trg.createVariable(var, 'f4', ('time', 'lat', 'lon',), fill_value=-9999.0, zlib=True,complevel=1)
            nc_var.units = unit
            nc_var.standard_name = name
            self.nc_trg.sync()

        miss =  float(nc_var._FillValue)
        data = pcr2numpy(pcrdata,miss)

        if self.bufflst.has_key(var):
            self.bufflst[var][bufpos,:,:] =  data
        else:
            self.bufflst[var] = self.timestepbuffer.copy()
            self.bufflst[var][bufpos,:,:] =  data

        # Write out timestep buffer.....
        if buffreset == 0 or idx ==  self.maxbuf -1 or self.timesteps <= timestep:
            spos = idx-bufpos
            self.logger.debug("Writing buffer for " + var + " to file at: " + str(spos) + " " + str(int(bufpos) + 1) + " timesteps")
            nc_var[spos:idx+1,:,:] = self.bufflst[var][0:bufpos+1,:,:]
            self.nc_trg.sync()


    def finish(self):
        """
        Flushes and closes the netcdf file

        :return: Nothing
        """
        if hasattr(self,"nc_trg"):
            self.nc_trg.sync()
            self.nc_trg.close()


class netcdfinput():

    def __init__(self,netcdffile,logging,vars=[]):
        """
        First try to setup a class read netcdf files
        (converted with pcr2netcdf.py)

        netcdffile: file to read the forcing data from
        logging: python logging object
        vars: list of variables to get from file
        """

        if os.path.exists(netcdffile):
            self.dataset = netCDF4.Dataset(netcdffile,mode='r')
        else:
            logging.error(os.path.abspath(netcdffile) + " not found!")
            exit(ValueError)

        logging.info("Reading input from netCDF file: " + netcdffile + ": " + str(self.dataset).replace('\n',' '))
        self.alldat ={}
        a = pcr2numpy(cover(0.0),0.0).flatten()
        # Determine steps to load in mem based on estimated memory usage
        floatspermb = 1048576/4
        maxmb = 4000
        self.maxsteps = maxmb * len(a)/floatspermb + 1
        self.fstep = 0
        self.lstep = self.fstep + self.maxsteps

        for var in vars:
            try:
                self.alldat[var] = self.dataset.variables[var][self.fstep:self.maxsteps]
            except:
                self.alldat.pop(var, None)
                logging.warn("Variable " + var + " not found in netcdf file: " + netcdffile)


    def gettimestep(self,timestep,logging,var='P'):
        """
        Gets a map for a single timestep. reads data in blocks assuming sequential access

        timestep: framework timestep (1-based)
        logging: python logging object
        var: variable to get from the file
        """
        ncindex = timestep -1
        if self.alldat.has_key(var):
            if ncindex == self.lstep: # Read new block of data in mem
                logging.debug("reading new netcdf data block starting at: " + str(ncindex))
                for vars in self.alldat:
                    self.alldat[vars] = self.dataset.variables[vars][ncindex:ncindex + self.maxsteps]
                self.fstep = ncindex
                self.lstep = ncindex + self.maxsteps
            np_step = self.alldat[var][ncindex-self.fstep,:,:]
            return numpy2pcr(Scalar, np_step, 1E31)
        else:
            logging.debug("Var (" + var + ") not found returning 0")
            return cover(scalar(0.0))
