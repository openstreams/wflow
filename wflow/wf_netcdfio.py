"""
wf_netcdfio
-----------

netcdf reading and writing for wflow

$Author: schelle $
$Id: wf_DynamicFramework.py 915 2014-02-10 07:33:56Z schelle $
$Rev: 915 $
"""

import datetime as dt
import os
import sys

import cftime
import netCDF4
import osgeo
import osgeo.ogr
import pyproj
import numpy as np
import pcraster as pcr

globmetadata = {}
globmetadata["title"] = "wflow output mapstack"
globmetadata["institution"] = "Deltares"
globmetadata["source"] = "wflow"
globmetadata["history"] = dt.datetime.now().isoformat()
globmetadata["references"] = "https://github.com/openstreams/wflow"
globmetadata["Conventions"] = "CF-1.4"


def convertCoord(proj_src, proj_trg, x, y):
    """
    Convert a list of x,y pairs in a certain projection to another projection
    input:
        proj_src:   string, EPSG or proj4 string referring to projection of source coordinates
        proj_trg:   string, EPSG or proj4 string referring to projection of target coordinates
        x:          NumPy array, vector or 2D array of x-coordinates (source)
        y:          NumPy array, vector or 2D array of y-coordinates (source)
    output:
        X:          NumPy array, vector or 2D array of x-coordinates (target)
        Y:          NumPy array, vector or 2D array of y-coordinates (target)
    """
    srs1 = pyproj.Proj(proj_src)  # OPT['proj4_params'])
    srs2 = pyproj.Proj(proj_trg)  # wgs84
    X, Y = pyproj.transform(srs1, srs2, x, y)  # Do add 0. to avoid trunc issues.
    return X, Y


def prepare_nc(
    trgFile,
    timeList,
    x,
    y,
    metadata,
    logger,
    EPSG="EPSG:4326",
    units=None,
    calendar="gregorian",
    Format="NETCDF4",
    complevel=9,
    zlib=True,
    least_significant_digit=None,
    FillValue=1e31,
):
    """
    This function prepares a NetCDF file with given metadata, for a certain year, daily basis data
    The function assumes a gregorian calendar and a time unit 'Days since 1900-01-01 00:00:00'
    """

    logger.info("Setting up netcdf output: " + trgFile)

    if units == None:  # Use start of the run
        epoch = timeList[0]
        units = "seconds since %04d-%02d-%02d %02d:%02d:%02d.0 00:00" % (
            epoch.year,
            epoch.month,
            epoch.day,
            epoch.hour,
            epoch.minute,
            epoch.second,
        )

    startDayNr = cftime.date2num(
        timeList[0].replace(tzinfo=None), units=units, calendar=calendar
    )
    endDayNr = cftime.date2num(
        timeList[-1].replace(tzinfo=None), units=units, calendar=calendar
    )

    timeAR = np.linspace(startDayNr, endDayNr, num=len(timeList))

    if os.path.exists(trgFile):
        os.remove(trgFile)

    nc_trg = netCDF4.Dataset(
        trgFile, "w", format=Format, zlib=zlib, complevel=complevel
    )

    logger.info(
        "Setting up dimensions and attributes. Steps: "
        + str(len(timeList))
        + " lat: "
        + str(len(y))
        + " lon: "
        + str(len(x))
    )
    if len(timeAR) == 1:
        nc_trg.createDimension("time", 1)
    else:
        nc_trg.createDimension("time", 0)  # NrOfDays*8

    DateHour = nc_trg.createVariable(
        "time", "f8", ("time",), fill_value=FillValue, zlib=zlib, complevel=complevel
    )
    DateHour.units = units
    DateHour.calendar = calendar
    DateHour.standard_name = "time"
    DateHour.long_name = "time"
    DateHour.axis = "T"
    DateHour[:] = timeAR

    # make a proj4 string
    srs = osgeo.osr.SpatialReference()
    res = srs.ImportFromEPSG(int(EPSG[5:]))
    if res != 0:
        logger.error(
            "EPGS not converted correctly: "
            + EPSG
            + ". Is the GDAL_DATA environment variable set correctly?"
        )
        sys.exit(1)

    projStr = srs.ExportToProj4()
    proj_src = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

    if srs.IsProjected() == 0:  # ONly lat lon needed
        nc_trg.createDimension("lat", len(y))
        nc_trg.createDimension("lon", len(x))
        y_var = nc_trg.createVariable(
            "lat", "f4", ("lat",), fill_value=FillValue, zlib=zlib, complevel=complevel
        )
        y_var.standard_name = "latitude"
        y_var.long_name = "latitude"
        y_var.units = "degrees_north"
        y_var.axis = "Y"
        x_var = nc_trg.createVariable(
            "lon", "f4", ("lon",), fill_value=FillValue, zlib=zlib, complevel=complevel
        )
        x_var.standard_name = "longitude"
        x_var.long_name = "longitude"
        x_var.units = "degrees_east"
        x_var.axis = "X"
        y_var[:] = y
        x_var[:] = x
        crs = nc_trg.createVariable("crs", "c")
        crs.long_name = "wgs84"
        crs.proj4_params = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
        crs.grid_mapping_name = "latitude_longitude"
    else:  # Assume regular grid in m
        nc_trg.createDimension("y", len(y))
        nc_trg.createDimension("x", len(x))
        y_var = nc_trg.createVariable(
            "y", "f4", ("y",), fill_value=FillValue, zlib=zlib, complevel=complevel
        )
        y_var.standard_name = "projection_y_coordinate"
        y_var.long_name = "y-coordinate in Cartesian system"
        y_var.units = "m"
        y_var.axis = "Y"
        x_var = nc_trg.createVariable(
            "x", "f4", ("x",), fill_value=FillValue, zlib=zlib, complevel=complevel
        )
        x_var.standard_name = "projection_x_coordinate"
        x_var.long_name = "x-coordinate in Cartesian system"
        x_var.units = "m"
        x_var.axis = "X"
        y_var[:] = y
        x_var[:] = x
        crs = nc_trg.createVariable("crs", "c")
        crs.long_name = EPSG
        crs.grid_mapping_name = "universal_transverse_mercator"
        crs.utm_zone_number = srs.GetUTMZone()
        crs.semi_major_axis = srs.GetSemiMajor()
        crs.inverse_flattening = srs.GetInvFlattening()
        crs._CoordinateTransformType = "Projection"
        crs._CoordinateAxisTypes = "y x"
        crs.proj4_params = projStr
        # Also write lat lon fields
        XI, YI = np.meshgrid(x, y)
        lon_vals, lat_vals = convertCoord(projStr, proj_src, XI, YI)
        # Need to create lat-lon fields
        lat = nc_trg.createVariable("lat", "f4", ("y", "x"))
        lat.standard_name = "latitude"
        lat.long_name = "latitude coordinate"
        lat.units = "degrees_north"
        lat.coordinates = "lat lon"
        lat.grid_mapping = "wgs84"
        # lat._CoordinateAxisType = "Lat"
        lat[:, :] = lat_vals
        lon = nc_trg.createVariable("lon", "f4", ("y", "x"))
        lon.standard_name = "longitude"
        lon.long_name = "longitude coordinate"
        lon.units = "degrees_east"
        lon.coordinates = "lat lon"
        lon.grid_mapping = "wgs84"
        # lon._CoordinateAxisType = "Lon"
        lon[:, :] = lon_vals

    crs.EPSG_code = EPSG

    # now add all attributes from user-defined metadata
    for attr in metadata:
        nc_trg.setncattr(attr, metadata[attr])
    nc_trg.sync()
    nc_trg.close()


class netcdfoutput:
    def __init__(
        self,
        netcdffile,
        logger,
        starttime,
        timesteps,
        EPSG="EPSG:4326",
        timestepsecs=86400,
        metadata={},
        zlib=True,
        Format="NETCDF4",
        maxbuf=25,
        least_significant_digit=None,
    ):
        """
        Under construction
        """

        self.EPSG = EPSG
        self.zlib = zlib
        self.Format = Format
        self.least_significant_digit = least_significant_digit

        def date_range(start, end, timestepsecs):
            r = int(
                (end + dt.timedelta(seconds=timestepsecs) - start).total_seconds()
                / timestepsecs
            )
            return [start + dt.timedelta(seconds=(timestepsecs * i)) for i in range(r)]

        self.logger = logger
        # Do not allow a max buffer larger than the number of timesteps
        self.maxbuf = maxbuf if timesteps >= maxbuf else timesteps
        self.ncfile = netcdffile
        self.timesteps = timesteps
        rows = pcr.clone().nrRows()
        cols = pcr.clone().nrCols()
        cellsize = pcr.clone().cellSize()
        yupper = pcr.clone().north()
        xupper = pcr.clone().west()
        x = pcr.pcr2numpy(pcr.xcoordinate(pcr.boolean(pcr.cover(1.0))), np.nan)[0, :]
        y = pcr.pcr2numpy(pcr.ycoordinate(pcr.boolean(pcr.cover(1.0))), np.nan)[:, 0]

        # Shift one timestep as we output at the end
        # starttime = starttime + dt.timedelta(seconds=timestepsecs)
        end = starttime + dt.timedelta(seconds=timestepsecs * (self.timesteps - 1))

        timeList = date_range(starttime, end, timestepsecs)
        self.timestepbuffer = np.zeros((int(self.maxbuf), len(y), len(x)))
        self.bufferdirty = True
        self.bufflst = {}

        globmetadata.update(metadata)

        prepare_nc(
            self.ncfile,
            timeList,
            x,
            y,
            globmetadata,
            logger,
            Format=self.Format,
            EPSG=EPSG,
            zlib=self.zlib,
            least_significant_digit=self.least_significant_digit,
        )

        self.nc_trg = None

    def savetimestep(
        self,
        timestep,
        pcrdata,
        unit="mm",
        var="P",
        name="Precipitation",
        flushonly=False,
    ):
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
        if not self.nc_trg:
            self.nc_trg = netCDF4.Dataset(
                self.ncfile, "a", format=self.Format, zlib=self.zlib, complevel=9
            )
            self.nc_trg.set_fill_off()
        # read time axis and convert to time objects
        # TODO: use this to append time
        # time = self.nc_trg.variables['time']
        # timeObj = cftime.num2date(time[:], units=time.units, calendar=time.calendar)

        idx = timestep - 1

        buffreset = int((idx + 1) % self.maxbuf)
        bufpos = int((idx) % self.maxbuf)

        try:
            nc_var = self.nc_trg.variables[var]
        except:
            self.logger.debug(
                "Creating variable " + var + " in netcdf file. Format: " + self.Format
            )
            if self.EPSG.lower() == "epsg:4326":
                nc_var = self.nc_trg.createVariable(
                    var,
                    "f4",
                    ("time", "lat", "lon"),
                    fill_value=-9999.0,
                    zlib=self.zlib,
                    complevel=9,
                    least_significant_digit=self.least_significant_digit,
                )
                nc_var.coordinates = "lat lon"
            else:
                nc_var = self.nc_trg.createVariable(
                    var,
                    "f4",
                    ("time", "y", "x"),
                    fill_value=-9999.0,
                    zlib=self.zlib,
                    complevel=9,
                    least_significant_digit=self.least_significant_digit,
                )
                nc_var.coordinates = "lat lon"
                nc_var.grid_mapping = "crs"

            nc_var.units = unit
            nc_var.standard_name = name
            self.nc_trg.sync()

        miss = float(nc_var._FillValue)
        data = pcr.pcr2numpy(pcr.scalar(pcrdata), miss)

        if var in self.bufflst:
            self.bufflst[var][bufpos, :, :] = data
        else:
            self.bufflst[var] = self.timestepbuffer.copy()
            self.bufflst[var][bufpos, :, :] = data

        # Write out timestep buffer.....
        self.bufferdirty = True

        if buffreset == 0 or idx == self.maxbuf - 1 or self.timesteps <= timestep:
            spos = idx - bufpos
            self.logger.debug(
                "Writing buffer for "
                + var
                + " to file at: "
                + str(spos)
                + " "
                + str(int(bufpos) + 1)
                + " timesteps"
            )
            nc_var[spos : idx + 1, :, :] = self.bufflst[var][0 : bufpos + 1, :, :]
            self.bufferdirty = False
            self.nc_trg.sync()

    def finish(self):
        """
        Flushes and closes the netcdf file

        :return: Nothing
        """
        if hasattr(self, "nc_trg"):
            if self.bufferdirty:
                self.logger.warning(
                    "Finishing before expected run-length exceeded. Buffer not flushed"
                )
            self.nc_trg.sync()
            self.nc_trg.close()


class netcdfoutputstatic:
    def __init__(
        self,
        netcdffile,
        logger,
        starttime,
        timesteps,
        EPSG="EPSG:4326",
        timestepsecs=86400,
        metadata={},
        zlib=True,
        Format="NETCDF4",
        maxbuf=25,
        least_significant_digit=None,
    ):
        """
        Under construction
        """

        self.EPSG = EPSG
        self.zlib = zlib
        self.Format = Format
        self.least_significant_digit = least_significant_digit

        def date_range(start, end, timestepsecs):
            r = int(
                (end + dt.timedelta(seconds=timestepsecs) - start).total_seconds()
                / timestepsecs
            )
            return [start + dt.timedelta(seconds=(timestepsecs * i)) for i in range(r)]

        self.logger = logger
        # Do not allow a max buffer larger than the number of timesteps
        self.maxbuf = maxbuf if timesteps >= maxbuf else timesteps
        self.ncfile = netcdffile
        self.timesteps = timesteps
        rows = pcr.clone().nrRows()
        cols = pcr.clone().nrCols()
        cellsize = pcr.clone().cellSize()
        yupper = pcr.clone().north()
        xupper = pcr.clone().west()
        x = pcr.pcr2numpy(pcr.xcoordinate(pcr.boolean(pcr.cover(1.0))), np.nan)[0, :]
        y = pcr.pcr2numpy(pcr.ycoordinate(pcr.boolean(pcr.cover(1.0))), np.nan)[:, 0]

        # Shift one timestep as we output at the end
        # starttime = starttime + dt.timedelta(seconds=timestepsecs)
        end = starttime + dt.timedelta(seconds=timestepsecs * (self.timesteps - 1))

        timeList = date_range(starttime, end, timestepsecs)
        self.timestepbuffer = np.zeros((self.maxbuf, len(y), len(x)))
        self.bufflst = {}
        self.buffdirty = False

        globmetadata.update(metadata)

        prepare_nc(
            self.ncfile,
            timeList,
            x,
            y,
            globmetadata,
            logger,
            Format=self.Format,
            EPSG=EPSG,
            zlib=self.zlib,
            least_significant_digit=self.least_significant_digit,
        )

    def savetimestep(self, timestep, pcrdata, unit="mm", var="P", name="Precipitation"):
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
        self.nc_trg = netCDF4.Dataset(
            self.ncfile, "a", format=self.Format, zlib=self.zlib, complevel=9
        )
        self.nc_trg.set_fill_off()
        # read time axis and convert to time objects
        # TODO: use this to append time
        # time = self.nc_trg.variables['time']
        # timeObj = cftime.num2date(time[:], units=time.units, calendar=time.calendar)

        idx = timestep - 1

        buffreset = (idx + 1) % self.maxbuf
        bufpos = (idx) % self.maxbuf

        try:
            nc_var = self.nc_trg.variables[var]
        except:
            self.logger.debug(
                "Creating variable " + var + " in netcdf file. Format: " + self.Format
            )
            if self.EPSG.lower() == "epsg:4326":
                nc_var = self.nc_trg.createVariable(
                    var,
                    "f4",
                    ("time", "lat", "lon"),
                    fill_value=-9999.0,
                    zlib=self.zlib,
                    complevel=9,
                    least_significant_digit=self.least_significant_digit,
                )
                nc_var.coordinates = "lat lon"
            else:
                nc_var = self.nc_trg.createVariable(
                    var,
                    "f4",
                    ("time", "y", "x"),
                    fill_value=-9999.0,
                    zlib=self.zlib,
                    complevel=9,
                    least_significant_digit=self.least_significant_digit,
                )
                nc_var.coordinates = "lat lon"
                nc_var.grid_mapping = "crs"

            nc_var.units = unit
            nc_var.standard_name = name
            self.nc_trg.sync()

        miss = float(nc_var._FillValue)
        data = pcr.pcr2numpy(pcr.scalar(pcrdata), miss)

        if var in self.bufflst:
            self.bufflst[var][bufpos, :, :] = data
            self.buffdirty = True
        else:
            self.bufflst[var] = self.timestepbuffer.copy()
            self.bufflst[var][bufpos, :, :] = data
            self.buffdirty = True

        # Write out timestep buffer.....

        if buffreset == 0 or idx == self.maxbuf - 1 or self.timesteps <= timestep:
            spos = idx - bufpos
            self.logger.debug(
                "Writing buffer for "
                + var
                + " to file at: "
                + str(spos)
                + " "
                + str(int(bufpos) + 1)
                + " timesteps"
            )
            nc_var[spos : idx + 1, :, :] = self.bufflst[var][0 : bufpos + 1, :, :]
            self.nc_trg.sync()
            self.buffdirty = False

    def finish(self):
        """
        Flushes and closes the netcdf file

        :return: Nothing
        """
        if hasattr(self, "nc_trg"):
            self.nc_trg.sync()
            self.nc_trg.close()
            if self.buffdirty:
                self.logger.error("Finishing with dirty netcdf write buffer...!")


class netcdfinput:
    def __init__(self, netcdffile, logging, vars=[]):
        """
        First try to setup a class read netcdf files
        (converted with pcr2netcdf.py)

        netcdffile: file to read the forcing data from
        logging: python logging object
        vars: list of variables to get from file
        """

        if os.path.exists(netcdffile):
            self.dataset = netCDF4.Dataset(netcdffile, mode="r")
        else:
            msg = os.path.abspath(netcdffile) + " not found!"
            logging.error(msg)
            raise ValueError(msg)

        logging.info("Reading input from netCDF file: " + netcdffile)
        self.alldat = {}
        a = pcr.pcr2numpy(pcr.cover(0.0), 0.0).flatten()
        # Determine steps to load in mem based on estimated memory usage
        floatspermb = 1048576 / 4
        maxmb = 40

        self.maxlentime = len(self.dataset.variables["time"])
        self.maxsteps = np.minimum(
            maxmb * len(a) / floatspermb + 1, self.maxlentime - 1
        )
        self.fstep = 0
        self.lstep = self.fstep + self.maxsteps
        self.offset = 0
        self.datetime = self.dataset.variables["time"][:]
        if hasattr(self.dataset.variables["time"], "units"):
            self.timeunits = self.dataset.variables["time"].units
        else:
            self.timeunits = "Seconds since 1970-01-01 00:00:00"
        if hasattr(self.dataset.variables["time"], "calendar"):
            self.calendar = self.dataset.variables["time"].calendar
        else:
            self.calendar = "gregorian"
        self.datetimelist = cftime.num2date(
            self.datetime, self.timeunits, calendar=self.calendar
        )

        try:
            self.x = self.dataset.variables["x"][:]
        except:
            self.x = self.dataset.variables["lon"][:]

        # Now check Y values to see if we must flip the data
        try:
            self.y = self.dataset.variables["y"][:]
        except:
            self.y = self.dataset.variables["lat"][:]

        # test if 1D or 2D array
        if len(self.y.shape) == 1:
            if self.y[0] > self.y[-1]:
                self.flip = False
            else:
                self.flip = True
        else:  # not sure if this works
            self.y = self.y[:][0]
            if self.y[0] > self.y[-1]:
                self.flip = False
            else:
                self.flip = True

        x = pcr.pcr2numpy(pcr.xcoordinate(pcr.boolean(pcr.cover(1.0))), np.nan)[0, :]
        y = pcr.pcr2numpy(pcr.ycoordinate(pcr.boolean(pcr.cover(1.0))), np.nan)[:, 0]

        # Get average cell size
        acc = (
            np.diff(x).mean() * 0.25
        )  # non-exact match needed becuase of possible rounding problems
        if self.flip:
            (self.latidx,) = np.logical_and(
                self.y[::-1] + acc >= y.min(), self.y[::-1] <= y.max() + acc
            ).nonzero()
            (self.lonidx,) = np.logical_and(
                self.x + acc >= x.min(), self.x <= x.max() + acc
            ).nonzero()
        else:
            (self.latidx,) = np.logical_and(
                self.y + acc >= y.min(), self.y <= y.max() + acc
            ).nonzero()
            (self.lonidx,) = np.logical_and(
                self.x + acc >= x.min(), self.x <= x.max() + acc
            ).nonzero()

        if len(self.lonidx) != len(x):
            logging.error("error in determining X coordinates in netcdf...")
            logging.error("model expects: " + str(x.min()) + " to " + str(x.max()))
            logging.error(
                "got coordinates  netcdf: "
                + str(self.x.min())
                + " to "
                + str(self.x.max())
            )
            logging.error(
                "got len from  netcdf x: "
                + str(len(x))
                + " expected "
                + str(len(self.lonidx))
            )
            raise ValueError("X coordinates in netcdf do not match model")

        if len(self.latidx) != len(y):
            logging.error("error in determining Y coordinates in netcdf...")
            logging.error("model expects: " + str(y.min()) + " to " + str(y.max()))
            logging.error(
                "got from  netcdf: " + str(self.y.min()) + " to " + str(self.y.max())
            )
            logging.error(
                "got len from  netcdf y: "
                + str(len(y))
                + " expected "
                + str(len(self.latidx))
            )
            raise ValueError("Y coordinates in netcdf do not match model")

        for var in vars:
            try:
                self.alldat[var] = self.dataset.variables[var]
            except:
                self.alldat.pop(var, None)
                logging.warning(
                    "Variable " + var + " not found in netcdf file: " + netcdffile
                )

    def gettimestep(self, timestep, logging, tsdatetime=None, var="P", shifttime=False):
        """
        Gets a map for a single timestep. reads data in blocks assuming sequential access

        :var timestep: framework timestep (1-based)
        :var logging: python logging object
        :var var: variable to get from the file
        :var shifttime: is True start at 1 in the NC file (instead of 0)
        :var tsdatetime: Assumed date/time of this timestep

            without ensembles (dims = 3)
            window = data[dpos,latidx.min():latidx.max()+1,lonidx.min():lonidx.max()+1]
            
            with ensembles from Delft-FEWS (dims = 4):
            window = data[dpos,realization,latidx.min():latidx.max()+1,lonidx.min():lonidx.max()+1]
        """
        if shifttime:
            ncindex = timestep
        else:
            ncindex = timestep - 1

        ncindex = ncindex + self.offset

        if self.datetimelist.size < ncindex + 1:
            ncindex = self.datetimelist.size - 1
        
        if tsdatetime != None:
            if tsdatetime.replace(tzinfo=None) != self.datetimelist[ncindex]:
                logging.warning(
                    "Date/time does not match. Wanted "
                    + str(tsdatetime)
                    + " got "
                    + str(self.datetimelist[ncindex])
                )
                import bisect

                pos = bisect.bisect_left(
                    self.datetimelist, tsdatetime.replace(tzinfo=None)
                )
                if pos >= self.datetimelist.size:
                    pos = self.datetimelist.size - 1
                    logging.warning(
                        "No matching date/time found using last date/time again..."
                    )
                self.offset = pos - ncindex
                logging.warning(
                    "Adjusting to the date/time at index and setting offset: "
                    + str(pos)
                    + ":"
                    + str(self.offset)
                    + ":"
                    + str(self.datetimelist[pos])
                )
                ncindex = pos

        if var in self.alldat:
            # if ncindex == self.lstep:  # Read new block of data in mem
            #    logging.debug("reading new netcdf data block starting at: " + str(ncindex))
            #    for vars in self.alldat:
            #        self.alldat[vars] = self.dataset.variables[vars][ncindex:ncindex + self.maxsteps]
            #
            # self.fstep = ncindex
            # self.lstep = ncindex + self.maxsteps

            if len(self.alldat[var].dimensions) == 3:
                np_step = self.alldat[var][
                    ncindex - self.fstep,
                    self.latidx.min() : self.latidx.max() + 1,
                    self.lonidx.min() : self.lonidx.max() + 1,
                ]
            if len(self.alldat[var].dimensions) == 4:
                np_step = self.alldat[var][
                    ncindex - self.fstep,
                    0,
                    self.latidx.min() : self.latidx.max() + 1,
                    self.lonidx.min() : self.lonidx.max() + 1,
                ]

            miss = float(self.dataset.variables[var]._FillValue)
            if self.flip:
                return pcr.numpy2pcr(pcr.Scalar, np.flipud(np_step).copy(), miss), True
            else:
                return pcr.numpy2pcr(pcr.Scalar, np_step, miss), True
        else:
            # logging.debug("Var (" + var + ") not found returning 0")
            return pcr.cover(pcr.scalar(0.0)), False


class netcdfinputstates:
    def __init__(self, netcdffile, logging, vars=[]):
        """
        First try to setup a class read netcdf files
        (converted with pcr2netcdf.py)

        netcdffile: file to read the forcing data from
        logging: python logging object
        vars: list of variables to get from file
        """

        self.fname = netcdffile
        if os.path.exists(netcdffile):
            self.dataset = netCDF4.Dataset(netcdffile, mode="r")
        else:
            msg = os.path.abspath(netcdffile) + " not found!"
            logging.error(msg)
            raise ValueError(msg)

        logging.info("Reading state input from netCDF file: " + netcdffile)
        self.alldat = {}
        a = pcr.pcr2numpy(pcr.cover(0.0), 0.0).flatten()
        # Determine steps to load in mem based on estimated memory usage
        floatspermb = 1048576 / 4
        maxmb = 40
        self.maxsteps = maxmb * len(a) / floatspermb + 1
        self.fstep = 0
        self.lstep = self.fstep + self.maxsteps

        self.datetime = self.dataset.variables["time"][:]
        if hasattr(self.dataset.variables["time"], "units"):
            self.timeunits = self.dataset.variables["time"].units
        else:
            self.timeunits = "Seconds since 1970-01-01 00:00:00"
        if hasattr(self.dataset.variables["time"], "calendar"):
            self.calendar = self.dataset.variables["time"].calendar
        else:
            self.calendar = "gregorian"
        self.datetimelist = cftime.num2date(
            self.datetime, self.timeunits, calendar=self.calendar
        )

        try:
            self.x = self.dataset.variables["x"][:]
        except:
            self.x = self.dataset.variables["lon"][:]

        # Now check Y values to see if we must flip the data
        try:
            self.y = self.dataset.variables["y"][:]
        except:
            self.y = self.dataset.variables["lat"][:]

        # test if 1D or 2D array
        if len(self.y.shape) == 1:
            if self.y[0] > self.y[-1]:
                self.flip = False
            else:
                self.flip = True
        else:  # not sure if this works
            self.y = self.y[:][0]
            if self.y[0] > self.y[-1]:
                self.flip = False
            else:
                self.flip = True

        x = pcr.pcr2numpy(pcr.xcoordinate(pcr.boolean(pcr.cover(1.0))), np.nan)[0, :]
        y = pcr.pcr2numpy(pcr.ycoordinate(pcr.boolean(pcr.cover(1.0))), np.nan)[:, 0]

        # Get average cell size
        acc = (
            np.diff(x).mean() * 0.25
        )  # non-exact match needed becuase of possible rounding problems
        if self.flip:
            (self.latidx,) = np.logical_and(
                self.y[::-1] + acc >= y.min(), self.y[::-1] <= y.max() + acc
            ).nonzero()
            (self.lonidx,) = np.logical_and(
                self.x + acc >= x.min(), self.x <= x.max() + acc
            ).nonzero()
        else:
            (self.latidx,) = np.logical_and(
                self.y + acc >= y.min(), self.y <= y.max() + acc
            ).nonzero()
            (self.lonidx,) = np.logical_and(
                self.x + acc >= x.min(), self.x <= x.max() + acc
            ).nonzero()

        if len(self.lonidx) != len(x):
            logging.error("error in determining X coordinates in netcdf...")
            logging.error("model expects: " + str(x.min()) + " to " + str(x.max()))
            logging.error(
                "got coordinates  netcdf: "
                + str(self.x.min())
                + " to "
                + str(self.x.max())
            )
            logging.error(
                "got len from  netcdf x: "
                + str(len(x))
                + " expected "
                + str(len(self.lonidx))
            )
            raise ValueError("X coordinates in netcdf do not match model")

        if len(self.latidx) != len(y):
            logging.error("error in determining Y coordinates in netcdf...")
            logging.error("model expects: " + str(y.min()) + " to " + str(y.max()))
            logging.error(
                "got from  netcdf: " + str(self.y.min()) + " to " + str(self.y.max())
            )
            logging.error(
                "got len from  netcdf y: "
                + str(len(y))
                + " expected "
                + str(len(self.latidx))
            )
            raise ValueError("Y coordinates in netcdf do not match model")

        for var in vars:
            try:
                self.alldat[var] = self.dataset.variables[var][
                    self.fstep : self.maxsteps
                ]
            except:
                self.alldat.pop(var, None)
                logging.warning(
                    "Variable " + var + " not found in netcdf file: " + netcdffile
                )

    def gettimestep(self, timestep, logging, var="P", tsdatetime=None):
        """
        Gets a map for a single timestep. reads data in blocks assuming sequential access

        timestep: framework timestep (1-based)
        logging: python logging object
        var: variable to get from the file
        """
        ncindex = timestep - 1

        if var in self.dataset.variables:
            if tsdatetime != None:
                if tsdatetime.replace(tzinfo=None) != self.datetimelist[
                    ncindex]:
                    logging.warning(
                        "Date/time of state ("
                        + var
                        + " in "
                        + self.fname
                        + ")does not match. Wanted "
                        + str(tsdatetime)
                        + " got "
                        + str(self.datetimelist[ncindex])
                    )

            np_step = self.dataset.variables[var][
                ncindex,
                self.latidx.min() : self.latidx.max() + 1,
                self.lonidx.min() : self.lonidx.max() + 1,
            ]

            miss = float(self.dataset.variables[var]._FillValue)
            return pcr.numpy2pcr(pcr.Scalar, np_step, miss), True
        else:
            # logging.debug("Var (" + var + ") not found returning map with 0.0")
            return pcr.cover(pcr.scalar(0.0)), False


class netcdfinputstatic:
    def __init__(self, netcdffile, logging):
        """
        First try to setup a class read netcdf files
        (converted with pcr2netcdf.py)

        netcdffile: file to read the forcing data from
        logging: python logging object
        vars: list of variables to get from file
        """

        if os.path.exists(netcdffile):
            self.dataset = netCDF4.Dataset(netcdffile, mode="r")
        else:
            msg = os.path.abspath(netcdffile) + " not found!"
            logging.error(msg)
            raise ValueError(msg)

        try:
            self.x = self.dataset.variables["x"][:]
        except:
            self.x = self.dataset.variables["lon"][:]
        # Now check Y values to see if we must flip the data
        try:
            self.y = self.dataset.variables["y"][:]
        except:
            self.y = self.dataset.variables["lat"][:]

        x = pcr.pcr2numpy(pcr.xcoordinate(pcr.boolean(pcr.cover(1.0))), np.nan)[0, :]
        y = pcr.pcr2numpy(pcr.ycoordinate(pcr.boolean(pcr.cover(1.0))), np.nan)[:, 0]

        (self.latidx,) = np.logical_and(self.x >= x.min(), self.x < x.max()).nonzero()
        (self.lonidx,) = np.logical_and(self.y >= x.min(), self.y < y.max()).nonzero()

        logging.info("Reading static input from netCDF file: " + netcdffile)

    def gettimestep(self, timestep, logging, var="P"):
        """
        Gets a map for a single timestep. reads data in blocks assuming sequential access

        timestep: framework timestep (1-based)
        logging: python logging object
        var: variable to get from the file
        """

        if var in self.dataset.variables:
            np_step = self.alldat[var][
                timestep - 1,
                self.latidx.min() : self.latidx.max() + 1,
                self.lonidx.min() : self.lonidx.max() + 1,
            ]
            miss = float(self.dataset.variables[var]._FillValue)
            return pcr.numpy2pcr(pcr.Scalar, np_step, miss), True
        else:
            logging.debug("Var (" + var + ") not found returning map with 0.0")
            return pcr.cover(pcr.scalar(0.0)), False
