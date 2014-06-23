__author__ = 'schelle'


import netCDF4
from pcraster import *
from numpy import *





def prepare_nc(trgFile, timeList, x, y, metadata, logger, units='Days since 1900-01-01 00:00:00', calendar='gregorian',Format="NETCDF4",zlib=False):
    """
    This function prepares a NetCDF file with given metadata, for a certain year, daily basis data
    The function assumes a gregorian calendar and a time unit 'Days since 1900-01-01 00:00:00'
    """
    import datetime as dt

    logger.info('Setting up "' + trgFile + '"')
    startDayNr = netCDF4.date2num(timeList[0], units=units, calendar=calendar)
    endDayNr   = netCDF4.date2num(timeList[-1], units=units, calendar=calendar)
    time       = arange(startDayNr,endDayNr+1)
    nc_trg     = netCDF4.Dataset(trgFile,'w',format=Format,zlib=zlib)

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



class netcdfoutput():

    def __init__(self,netcdffile,vars=[]):
        """
        Under construction
        """


    def savetimestep(self,timestep,var='P'):
        """
        ss
        """



class netcdfinput():
    def __init__(self,netcdffile,logging,vars=[]):
        """
        First try to setup a class read netcdf files
        (converted with pcr2netcdf.py)

        netcdffile: file to read the forcing data from
        logging: python logging object
        vars: list of variables to get from file
        """
        self.dataset = netCDF4.Dataset(netcdffile,mode='r')
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
                print ncindex - self.fstep
            np_step = self.alldat[var][ncindex-self.fstep,:,:]
            return numpy2pcr(Scalar, np_step, 1E31)
        else:
            logging.debug("Var (" + var + ") not found returning 0")
            return cover(scalar(0.0))
