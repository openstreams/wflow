
import wflow.bmi as bmi
import wflow
import logging
import os
from wflow.pcrut import setlogger



class runnerbmi(bmi.Bmi):
    """
    csdms BMI implementation for pcraster/python models

    .. todo::

        implement translation of units

    .. todo::

        implement translation of long_var_names
    """

    def __init__(self):
        """
        Initialises the object

        :return: nothing
        """
        self.loggingmode=logging.DEBUG
        self.bmi_routing = None
        self.bmi_hydro = None
        self.bmilogger = setlogger('bmi_runner.log','wflow_bmi_logging',thelevel=self.loggingmode)
        self.bmilogger.info("__init__: wflow_bmi object initialised.")

    def initialize_config(self, filename):
        """
        *Extended functionality*, see https://github.com/eWaterCycle/bmi/blob/master/src/main/python/bmi.py

        see initialize

        :param filename:
        :param loglevel:
        :return: nothing
        """

        self.bmimodels = {}

        self.datadir = os.path.dirname(filename)
        inifile = os.path.basename(filename)
        hydroconfigfilename = os.path.join(self.datadir,'wflow_sbm.ini')
        routingconfigfilename = os.path.join(self.datadir,'wflow_routing.ini')

        self.bmimodels[hydroconfigfilename] = wflow.wflow_bmi.wflowbmi_csdms()
        self.bmimodels[routingconfigfilename] = wflow.wflow_bmi.wflowbmi_csdms()

        # now get the to config files from the config file of overall bmi
        #self.bmimodels[hydroconfigfilename].initialize_config(hydroconfigfilename)
        #self.bmimodels[routingconfigfilename].initialize_config(routingconfigfilename)

        for key, value in self.bmimodels.iteritems():
            value.initialize_config(key)




    def initialize_model(self):
        """
        *Extended functionality*, see https://github.com/eWaterCycle/bmi/blob/master/src/main/python/bmi.py

        see initialize

        :param self:
        :return: nothing
        """

        for key, value in self.bmimodels.iteritems():
            value.initialize_model(key)


    def set_start_time(self, start_time):
        """

        :param start_time: time in units (seconds) since the epoch
        :return: nothing
        """
        for key, value in self.bmimodels.iteritems():
            self.bmimodels[key].set_start_time(start_time)

    def set_end_time(self, end_time):
        """
        :param end_time: time in units (seconds) since the epoch
        :return:
        """
        for key, value in self.bmimodels.iteritems():
            self.bmimodels[key].set_send_time(end_time)



    def get_attribute_names(self):
        """
        Get the attributes of the model return in the form of section_name:attribute_name

        :return: list of attributes
        """
        #TODO: Prefix with foreward slash of  model name!!!!
        attr = []

        for key, value in self.bmimodels.iteritems():
            attr = attr + self.bmimodels[key].get_attribute_names()

        return attr


    def get_attribute_value(self, attribute_name):
        """
        :param attribute_name:
        :return: attribute value
        """
        raise NotImplementedError

    def set_attribute_value(self, attribute_name, attribute_value):
        """
        :param attribute_name: name using the section:option notation
        :param attribute_value: string value of the option
        :return:
        """
        raise NotImplementedError



    def initialize(self, filename,loglevel=logging.DEBUG):
        """
        Initialise the model. Should be called before any other method.

        :var filename: full path to the wflow ini file
        :var loglevel: optional loglevel (default == DEBUG)

        Assumptions for now:

            - the configfile wih be a full path
            - we define the case from the basedir of the configfile


        """

        self.bmilogger.info("initialize: Initialising wflow bmi with ini: " + filename)
        self.initialize_config(filename,loglevel=loglevel)
        self.initialize_model()


    def update(self):
        """
        Propagate the model to the next model timestep
        """
        for key, value in self.bmimodels.iteritems():
            self.bmimodels[key].update()

        self.currenttimestep = self.currenttimestep + 1

    def update_until(self, time):
        """
        Update the model until and including the given time.

        - one or more timesteps foreward
        - max one timestep backward

        :var  double time: time in the units and epoch returned by the function get_time_units.
        """
        curtime = self.get_current_time()

        if abs(time - curtime)% self.get_time_step() != 0:
            self.bmilogger.error('update_until: timespan not dividable by timestep: ' + str(abs(time - curtime)) +
                                 ' and ' + str(self.get_time_step()))
            raise ValueError("Update in time not a multiple of timestep")

        if curtime > time:
            timespan = curtime - time
            nrstepsback = int(timespan/self.get_time_step())
            self.bmilogger.debug('update_until: update timesteps back ' + str(nrstepsback) + ' to ' + str(curtime + timespan))
            if nrstepsback > 1:
                raise ValueError("Time more than one timestep before current time.")
            for key, value in self.bmimodels.iteritems():
                self.bmimodels[key].dynModel.wf_QuickResume()

        else:
            timespan = time - curtime
            nrsteps = int(timespan/self.get_time_step())
            self.bmilogger.debug('update_until: update ' + str(nrsteps) + ' timesteps forward from ' + str(curtime) + ' to ' + str(curtime + timespan))
            self.bmilogger.debug('update_until: step ' + str(self.currenttimestep) + ' to ' + str(self.currenttimestep + nrsteps -1))

            #self.dynModel._runDynamic(self.currenttimestep, self.currenttimestep + nrsteps -1)
            for key, value in self.bmimodels.iteritems():
                self.bmimodels[key].update_until(time)

            self.currenttimestep = self.currenttimestep + nrsteps

    def update_frac(self, time_frac):
        """
        Not implemented. Raises a NotImplementedError
        """
        raise NotImplementedError

    def save_state(self, destination_directory):
        """
        Ask the model to write its complete internal current state to one or more state files in the given directory.
        Afterwards the given directory will only contain the state files and nothing else.
        Sates are save in the models' native format.

        :var destination_directory: the directory in which the state files should be written.
        """
        for key, value in self.bmimodels.iteritems():
            self.bmimodels[key].save_state(destination_directory)


    def load_state(self, source_directory):
        """
        Ask the model to load its complete internal current state from one or more
        state files in the given directory.

        :var  source_directory: the directory from which the state files should be
        read.
        """
        for key, value in self.bmimodels.iteritems():
            self.bmimodels[key].save_state(source_directory)

    def finalize(self):
        """
        Shutdown the library and clean up the model.
        Uses the default (model configured) state location to also save states.
        """
        # First check if the seconf initilize_states has run
        for key, value in self.bmimodels.iteritems():
            self.bmimodels[key].finalize()

    def get_component_name(self):
        """
        :return:  identifier of the model based on the name of the ini file
        """
        names = []
        for key, value in self.bmimodels.iteritems():
            names.append(self.bmimodels[key].get_component_name())

        return ",".join(names)


    def get_input_var_names(self):
        """

        :return: List of String objects: identifiers of all input variables of the model:
        """
        names = []
        for key, value in self.bmimodels.iteritems():
            names.append(self.bmimodels[key].get_input_var_names())

        return inames

    def get_output_var_names(self):
        """
        Returns the list of model output variables

        :return: List of String objects: identifiers of all output variables of the model:
        """
        namesroles = self.dynModel.wf_supplyVariableNamesAndRoles()
        inames = []

        for varrol in namesroles:
            if varrol[1] == 1 or varrol[1] == 2:
                inames.append(varrol[0])

        self.bmilogger.debug("get_output_var_names: " + str(inames))
        return inames

    def get_var_type(self, long_var_name):
        """
        Gets the variable type as a numpy type string

        :return: variable type string, compatible with numpy:
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(long_var_name)

        if hasattr(npmap,'dtype'):
            self.bmilogger.debug("get_var_type: " + str(npmap.dtype))
            return str(npmap.dtype)
        else:
            self.bmilogger.debug("get_var_type: " + str(None))
            return None

    def get_var_rank(self, long_var_name):
        """
        Gets the number of dimensions for a variable

        :var  String long_var_name: identifier of a variable in the model:
        :return: array rank or 0 for scalar (number of dimensions):
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(long_var_name)

        self.bmilogger.debug("get_var_rank: (" + long_var_name + ") " + str(len(npmap.shape)))
        return len(npmap.shape)

    def get_var_size(self, long_var_name):
        """
        Gets the number of elements in a variable (rows * cols)

        :var  String long_var_name: identifier of a variable in the model:
        :return: total number of values contained in the given variable (number of elements in map)
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(long_var_name)

        if hasattr(npmap,'size'):
            self.bmilogger.debug("get_var_size: " + str(npmap.size))
            return npmap.size
        else:
            self.bmilogger.debug("get_var_size: " + str(1))
            return 1

    def get_var_nbytes(self, long_var_name):
        """
        Gets the number of bytes occupied in memory for a given variable.

        :var  String long_var_name: identifier of a variable in the model:
        :return: total number of bytes contained in the given variable (number of elements * bytes per element)
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(long_var_name)

        self.bmilogger.debug("get_var_nbytes: " + str(npmap.size * npmap.itemsize))
        return npmap.size * npmap.itemsize

    def get_start_time(self):
        """
        Gets the start time of the model.

        :return: start time in the units and epoch returned by the function get_time_units
        """
        st = self.dynModel.wf_supplyStartTime()
        self.bmilogger.debug("get_start_time: " + str(st))
        return st

    def get_current_time(self):
        """
        Get the current time since the epoch of the model

        :return: current time of simulation n the units and epoch returned by the function get_time_units
        """

        st = self.dynModel.wf_supplyCurrentTime()
        self.bmilogger.debug("get_current_time: " + str(st) + " " + str(self.dynModel.DT.currentDateTime.strftime("%Y-%m-%d %H:%M:%S")))
        return st

    def get_end_time(self):
        """
        Get the end time of the model run in units since the epoch

        :return: end time of simulation n the units and epoch returned by the function get_time_units
        """
        et = self.dynModel.wf_supplyEndTime()
        self.bmilogger.debug("get_end_time: " + str(et))
        return et

    def get_time_step(self):
        """
        Get the model time steps in units since the epoch

        :return: duration of one time step of the model in the units returned by the function get_time_units
        """
        ts = self.dynModel.DT.timeStepSecs
        self.bmilogger.debug("get_time_step: " + str(ts))
        return ts

    def get_time_units(self):
        """
        Return the time units of the model as a string

        :return: Return a string formatted using the UDUNITS standard from Unidata.
        (http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/cf-conventions.html#time-coordinate)
        """
        tu = self.dynModel.wf_supplyEpoch()
        self.bmilogger.debug("get_time_units: " + str(tu))

        return tu

    def get_value(self, long_var_name):
        """
        Get the value(s) of a variable as a numpy array

        :var long_var_name: name of the variable
        :return: a np array of long_var_name
        """
        if long_var_name in self.inputoutputvars:
            ret = self.dynModel.wf_supplyMapAsNumpy(long_var_name)
            self.bmilogger.debug("get_value: " + long_var_name)
            if self.wrtodisk:
                fname = str(self.currenttimestep) + "_get_" + long_var_name + ".map"
                arpcr = self.dynModel.numpy2pcr(self.dynModel.Scalar, ret, -999)
                self.bmilogger.debug("Writing to disk: " + fname)
                self.dynModel.report(arpcr,fname)

            try:
                fret = np.flipud(ret)
                return fret
            except:
                return ret
        else:
            self.bmilogger.error("get_value: " + long_var_name + ' not in list of output values ' + str(self.inputoutputvars))
            return None

    def get_value_at_indices(self, long_var_name, inds):
        """
        Get a numpy array of the values at the given indices

        :var long_var_name: identifier of a variable in the model:
        :var inds: List of list each tuple contains one index for each dimension of the given variable, i.e. each tuple indicates one element in the multi-dimensional variable array:

        :return: numpy array of values in the data type returned by the function get_var_type.
        """

        if long_var_name in self.inputoutputvars:
            self.bmilogger.debug("get_value_at_indices: " + long_var_name + ' at ' + str(inds))
            npmap = np.flipud(self.dynModel.wf_supplyMapAsNumpy(long_var_name))
            return npmap[inds]
        else:
            self.bmilogger.error("get_value_at_indices: " + long_var_name + ' not in list of output values ' + str(self.inputoutputvars))
            return None


    def set_value_at_indices(self, long_var_name, inds, src):
        """
        Set the values in a variable using a numpy array of the values given indices

        :var long_var_name: identifier of a variable in the model:
        :var inds: List of Lists of integers inds each nested List contains one index for each dimension of the given variable,
                                        i.e. each nested List indicates one element in the multi-dimensional variable array,
                                        e.g. [[0, 0], [0, 1], [15, 19], [15, 20], [15, 21]] indicates 5 elements in a 2D grid.:
        :var src: Numpy array of values. one value to set for each of the indicated elements:
        """

        if long_var_name in self.outputonlyvars:
            self.bmilogger.error("set_value_at_indices: " + long_var_name + " is listed as an output only variable, cannot set. " + str(self.outputonlyvars))
            raise ValueError("set_value_at_indices: " + long_var_name + " is listed as an output only variable, cannot set. " + str(self.outputonlyvars))
        else:
            self.bmilogger.debug("set_value_at_indices: " + long_var_name + ' at ' + str(inds))
            npmap = np.flipud(self.dynModel.wf_supplyMapAsNumpy(long_var_name))
            npmap[inds] = src
            self.dynModel.wf_setValuesAsNumpy(long_var_name,np.flipud(npmap))

    def get_grid_type(self, long_var_name):
        """
        Get the grid type according to the enumeration in BmiGridType

        :var String long_var_name: identifier of a variable in the model.

        :return: BmiGridType type of the grid geometry of the given variable.
        """
        ret=BmiGridType()

        self.bmilogger.debug("get_grid_type: " + long_var_name + ' result: ' + str(ret.UNIFORM))

        return ret.UNIFORM

    def get_grid_shape(self, long_var_name):
        """
        Return the shape of the grid. Only return something for variables with a uniform, rectilinear or structured grid. Otherwise raise ValueError.

        :var long_var_name: identifier of a variable in the model.

        :return: List of integers: the sizes of the dimensions of the given variable, e.g. [500, 400] for a 2D grid with 500x400 grid cells.
        """
        dim =  self.dynModel.wf_supplyGridDim()
        #[ Xll, Yll, xsize, ysize, rows, cols]

        self.bmilogger.debug("get_grid_shape: " + long_var_name + ' result: ' + str([dim[4], dim[5]]))

        return [dim[4], dim[5]]

    def get_grid_spacing(self, long_var_name):
        """
        Only return something for variables with a uniform grid. Otherwise raise ValueError.

        :var long_var_name: identifier of a variable in the model.

        :return: The size of a grid cell for each of the dimensions of the given variable, e.g. [width, height]: for a 2D grid cell.
        """
        dims = self.dynModel.wf_supplyGridDim()[2:4]
        x = dims[0]
        y = dims[1]
        self.bmilogger.debug("get_grid_spacing: " + long_var_name + ' result: ' + str([y, x]))
        return [y, x]

    def get_grid_origin(self, long_var_name):
        """
        gets the origin of the model grid.

        :var String long_var_name: identifier of a variable in the model.

        :return: X, Y: ,the lower left corner of the grid.
        """
        dims = self.dynModel.wf_supplyGridDim()
        x = dims[0]
        y = dims[1]
        self.bmilogger.debug("get_grid_origin: " + long_var_name + ' result: ' + str([y, x]))
        return [y, x]

    def get_grid_x(self, long_var_name):
        """
        Give X coordinates of point in the model grid

        :var String long_var_name: identifier of a variable in the model.

        :return: Numpy array of doubles: x coordinate of grid cell center for each grid cell, in the same order as the
        values returned by function get_value.
        """
        self.bmilogger.debug("get_grid_x: " + long_var_name)
        return self.dynModel.wf_supplyMapXAsNumpy()

    def get_grid_y(self, long_var_name):
        """
        Give Y coordinates of point in the model grid

        :var String long_var_name: identifier of a variable in the model.

        :return: Numpy array of doubles: y coordinate of grid cell center for each grid cell, in the same order as the
        values returned by function get_value.

        """
        self.bmilogger.debug("get_grid_y: " + long_var_name)
        return np.flipud(self.dynModel.wf_supplyMapYAsNumpy())

    def get_grid_z(self, long_var_name):
        """
        Give Z coordinates of point in the model grid

        :var String long_var_name: identifier of a variable in the model.

        :return: Numpy array of doubles: z coordinate of grid cell center for each grid cell, in the same order as the values returned by function get_value.
        """
        self.bmilogger.debug("get_grid_z: " + long_var_name)
        return np.flipud(self.dynModel.wf_supplyMapZAsNumpy())

    def get_var_units(self, long_var_name):
        """
        Supply units as defined in the API section of the ini file

        :var long_var_name: identifier of a variable in the model.

        :return:   String: unit of the values of the given variable. Return a string formatted
        using the UDUNITS standard from Unidata. (only if set properly in the ini file)
        """

        nru = self.dynModel.wf_supplyVariableNamesAndRoles()

        unit ='mm'

        for it in nru:
            if long_var_name == it[0]:
                unit = it[2]

        self.bmilogger.debug("get_var_units: " + long_var_name + ' result: ' + unit)
        return unit

    def set_value(self, long_var_name, src):
        """
        Set the values(s) in a map using a numpy array as source

        :var long_var_name: identifier of a variable in the model.
        :var src: all values to set for the given variable. If only one value
                  is present a uniform map will be set in the wflow model.
        """

        if self.wrtodisk:
            fname = str(self.currenttimestep) + "_set_" + long_var_name + ".map"
            arpcr = self.dynModel.numpy2pcr(self.dynModel.Scalar, src, -999)
            self.bmilogger.debug("Writing to disk: " + fname)
            self.dynModel.report(arpcr,fname)

        if long_var_name in self.outputonlyvars:
            self.bmilogger.error("set_value: " + long_var_name + " is listed as an output only variable, cannot set. " + str(self.outputonlyvars))
            raise ValueError("set_value: " + long_var_name + " is listed as an output only variable, cannot set. " + str(self.outputonlyvars))
        else:
            if len(src) == 1:
                self.bmilogger.debug("set_value: (uniform value) " + long_var_name + '(' +str(src) + ')')
                self.dynModel.wf_setValues(long_var_name,float(src))
            else:
                self.bmilogger.debug("set_value: (grid) " + long_var_name)
                self.dynModel.wf_setValuesAsNumpy(long_var_name, np.flipud(src))

    def get_grid_connectivity(self, long_var_name):
        """
        Not applicable, raises NotImplementedError
        Should return the ldd if present!!
        """
        raise NotImplementedError

    def get_grid_offset(self, long_var_name):
        """
        Not applicable raises NotImplementedError
        """
        raise NotImplementedError
