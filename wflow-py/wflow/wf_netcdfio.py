__author__ = 'schelle'


import netCDF4
from pcraster import *


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
