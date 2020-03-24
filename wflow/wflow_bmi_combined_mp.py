import configparser
import logging
import os

import numpy as np
import wflow.bmi as bmi
import wflow.wflow_bmi as wfbmi
import pcraster as pcr
from wflow.pcrut import setlogger


def iniFileSetUp(configfile):
    """
    Reads .ini file and returns a config object.


    """
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(configfile)
    return config


def configsection(config, section):
    """
    gets the list of lesy in a section

    Input:
        - config
        - section

    Output:
        - list of keys in the section
    """
    try:
        ret = config.options(section)
    except:
        ret = []

    return ret


class wflowbmi_csdms(bmi.Bmi):
    """
    csdms BMI implementation runner for combined pcraster/python models


    + all variables are identified by: component_name@variable_name
    + this version is only tested for a one-way link
    + get_component_name returns a comma separated list of components

    """

    def __init__(self):
        """
        Initialises the object

        :return: nothing
        """
        from collections import OrderedDict

        self.bmimodels = OrderedDict()
        self.currenttimestep = 0
        self.exchanges = []
        self.comp_sep = "@"
        self.wrtodisk = False
        if os.getenv("wflow_bmi_combined_writetodisk", "False") in "True":
            self.wrtodisk = True

        self.loggingmode = logging.ERROR
        logstr = os.getenv("wflow_bmi_loglevel", "ERROR")
        if logstr in "ERROR":
            self.loggingmode = logging.ERROR
        if logstr in "WARNING":
            self.loggingmode = logging.WARNING
        if logstr in "INFO":
            self.loggingmode = logging.INFO
        if logstr in "DEBUG":
            self.loggingmode = logging.DEBUG

        self.bmilogger = setlogger(
            "wflow_bmi_combined.log",
            "wflow_bmi_combined_logging",
            thelevel=self.loggingmode,
        )
        self.bmilogger.info("__init__: wflow_bmi_combined object initialised.")
        if self.wrtodisk:
            self.bmilogger.warning("Will write all bmi set- and get- grids to disk!...")

    def __getmodulenamefromvar__(self, long_var_name):
        """

        :param long_var_name:
        :return: name of the module
        """
        return long_var_name.split(self.comp_sep)[0]

    def initialize_config(self, filename, loglevel=logging.DEBUG):
        """
        *Extended functionality*, see https://github.com/eWaterCycle/bmi/blob/master/src/main/python/bmi.py

        Read the ini file for the comnined bmi model and initializes all the bmi models
        listed in the config file.

        :param filename:
        :return: nothing
        """

        self.currenttimestep = 1

        fullpathname = os.path.abspath(filename)
        self.config = iniFileSetUp(fullpathname)
        self.datadir = os.path.dirname(fullpathname)
        inifile = os.path.basename(filename)

        self.models = configsection(self.config, "models")
        self.exchanges = configsection(self.config, "exchanges")

        for mod in self.models:
            self.bmimodels[mod] = wfbmi.wflowbmi_csdms()

        # Initialize all bmi model objects
        for key, value in self.bmimodels.items():
            modconf = os.path.join(self.datadir, self.config.get("models", key))
            self.bmimodels[key].initialize_config(modconf, loglevel=loglevel)

    def initialize_model(self):
        """
        *Extended functionality*, see https://github.com/eWaterCycle/bmi/blob/master/src/main/python/bmi.py

        initalises all bmi models listed in the config file. Als does the first (timestep 0) data exchange

        :param self:
        :return: nothing
        """

        for key, value in self.bmimodels.items():
            self.bmimodels[key].initialize_model()

        # Copy and set the variables to be exchanged for step 0
        for key, value in self.bmimodels.items():
            # step one update first model
            curmodel = self.bmimodels[key].get_component_name()
            for item in self.exchanges:
                supplymodel = self.__getmodulenamefromvar__(item)
                if curmodel == supplymodel:
                    outofmodel = self.get_value(item).copy()
                    tomodel = self.config.get("exchanges", item)
                    self.set_value(tomodel, outofmodel)

        self.bmilogger.info(self.bmimodels)

    def set_start_time(self, start_time):
        """
        Sets the start time for all bmi models

        :param start_time: time in units (seconds) since the epoch
        :return: nothing
        """
        for key, value in self.bmimodels.items():
            self.bmimodels[key].set_start_time(start_time)

    def set_end_time(self, end_time):
        """
        sets the end time for all bmi models

        :param end_time: time in units (seconds) since the epoch
        :return:
        """
        for key, value in self.bmimodels.items():
            self.bmimodels[key].set_end_time(end_time)

    def get_attribute_names(self):
        """
        Get the attributes of the model return in the form of model_name@section_name:attribute_name

        :return: list of attributes
        """
        names = []
        for key, value in self.bmimodels.items():
            names.append(self.bmimodels[key].get_attribute_names())
            names[-1] = [
                self.bmimodels[key].get_component_name() + self.comp_sep + s
                for s in names[-1]
            ]

        ret = [item for sublist in names for item in sublist]
        return ret

    def get_attribute_value(self, attribute_name):
        """
        gets the attribute value for the name. The name should adhere to the following convention::
            model_name@section_name:attribute_name

        :param attribute_name:
        :return: attribute value
        """
        cname = attribute_name.split(self.comp_sep)
        return self.bmimodels[cname[0]].get_attribute_value(cname[1])

    def set_attribute_value(self, attribute_name, attribute_value):
        """
        :param attribute_name: name using the model_name@section:option notation
        :param attribute_value: string value of the option
        :return:
        """
        cname = attribute_name.split(self.comp_sep)
        self.bmimodels[cname[0]].set_attribute_value(cname[1], attribute_value)

    def initialize(self, filename, loglevel=logging.DEBUG):
        """
        Initialise the model. Should be called before any other method.

        :var filename: full path to the combined model ini file
        """

        self.initialize_config(filename, loglevel=loglevel)
        self.initialize_model()

    def update(self):
        """
        Propagate the model to the next model timestep

        The function iterates over all models
        """
        for key, value in self.bmimodels.items():
            # step one update model
            self.bmimodels[key].update()
            # do all exchanges
            curmodel = self.bmimodels[key].get_component_name()
            for item in self.exchanges:
                supplymodel = self.__getmodulenamefromvar__(item)
                if curmodel == supplymodel:
                    outofmodel = self.get_value(item)
                    tomodel = self.config.get("exchanges", item)
                    self.set_value(tomodel, outofmodel)

        self.currenttimestep = self.currenttimestep + 1

    def update_until(self, time):
        """
        Update the model until and including the given time.

        - one or more timesteps foreward
        - max one timestep backward

        :var  double time: time in the units and epoch returned by the function get_time_units.
        """
        curtime = self.get_current_time()

        if abs(time - curtime) % self.get_time_step() != 0:
            raise ValueError("Update in time not a multiple of timestep")

        if curtime > time:
            timespan = curtime - time
            nrstepsback = int(timespan / self.get_time_step())
            if nrstepsback > 1:
                raise ValueError("Time more than one timestep before current time.")
            for key, value in self.bmimodels.items():
                self.bmimodels[key].dynModel.wf_QuickResume()

        else:
            timespan = time - curtime
            nrsteps = int(timespan / self.get_time_step())

            # self.dynModel._runDynamic(self.currenttimestep, self.currenttimestep + nrsteps -1)
            for st in range(0, nrsteps):
                # for key, value in self.bmimodels.iteritems():
                self.update()

            # self.currenttimestep = self.currenttimestep + nrsteps

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
        for key, value in self.bmimodels.items():
            self.bmimodels[key].save_state(destination_directory)

    def load_state(self, source_directory):
        """
        Ask the model to load its complete internal current state from one or more
        state files in the given directory.

        :var  source_directory: the directory from which the state files should be
        read.
        """
        for key, value in self.bmimodels.items():
            self.bmimodels[key].save_state(source_directory)

    def finalize(self):
        """
        Shutdown the library and clean up the model.
        Uses the default (model configured) state location to also save states.
        """
        for key, value in self.bmimodels.items():
            self.bmimodels[key].finalize()

        self.bmilogger.info("finalize.")

    def get_component_name(self):
        """
        :return:  identifier of the models separated by a comma (,)
        """
        names = []
        for key, value in self.bmimodels.items():
            names.append(self.bmimodels[key].get_component_name())

        return ",".join(names)

    def get_input_var_names(self):
        """

        :return: List of String objects: identifiers of all input variables of the model:
        """
        names = []
        for key, value in self.bmimodels.items():
            names.append(self.bmimodels[key].get_input_var_names())
            names[-1] = [
                self.bmimodels[key].get_component_name() + self.comp_sep + s
                for s in names[-1]
            ]

        ret = [item for sublist in names for item in sublist]
        return ret

    def get_output_var_names(self):
        """
        Returns the list of model output variables

        :return: List of String objects: identifiers of all output variables of the model:
        """
        names = []
        for key, value in self.bmimodels.items():
            names.append(self.bmimodels[key].get_output_var_names())
            names[-1] = [
                self.bmimodels[key].get_component_name() + self.comp_sep + s
                for s in names[-1]
            ]

        ret = [item for sublist in names for item in sublist]

        return ret

    def get_var_type(self, long_var_name):
        """
        Gets the variable type as a numpy type string

        :return: variable type string, compatible with numpy:
        """

        # first part should be the component name
        cname = long_var_name.split(self.comp_sep)
        for key, value in self.bmimodels.items():
            nn = self.bmimodels[key].get_component_name()
            if nn == cname[0]:
                ret = self.bmimodels[key].get_var_type(cname[1])

        return ret

    def get_var_rank(self, long_var_name):
        """
        Gets the number of dimensions for a variable

        :var  String long_var_name: identifier of a variable in the model:
        :return: array rank or 0 for scalar (number of dimensions):
        """
        # first part should be the component name
        cname = long_var_name.split(self.comp_sep)
        for key, value in self.bmimodels.items():
            nn = self.bmimodels[key].get_component_name()
            if nn == cname[0]:
                ret = self.bmimodels[key].get_var_rank(cname[1])

        return ret

    def get_var_size(self, long_var_name):
        """
        Gets the number of elements in a variable (rows * cols)

        :var  String long_var_name: identifier of a variable in the model:
        :return: total number of values contained in the given variable (number of elements in map)
        """
        # first part should be the component name
        cname = long_var_name.split(self.comp_sep)
        for key, value in self.bmimodels.items():
            nn = self.bmimodels[key].get_component_name()
            if nn == cname[0]:
                ret = self.bmimodels[key].get_var_size(cname[1])

        return ret

    def get_var_nbytes(self, long_var_name):
        """
        Gets the number of bytes occupied in memory for a given variable.

        :var  String long_var_name: identifier of a variable in the model:
        :return: total number of bytes contained in the given variable (number of elements * bytes per element)
        """
        # first part should be the component name
        cname = long_var_name.split(self.comp_sep)
        for key, value in self.bmimodels.items():
            nn = self.bmimodels[key].get_component_name()
            if nn == cname[0]:
                ret = self.bmimodels[key].get_var_nbytes(cname[1])

        return ret

    def get_start_time(self):
        """
        Gets the start time of the model.

        :return: start time of last model in the list. Tiem sare assumed to be identical
        """
        st = []
        for key, value in self.bmimodels.items():
            st.append(self.bmimodels[key].get_start_time())

        return np.array(st).max()

    def get_current_time(self):
        """
        Get the current time since the epoch of the model

        :return: current time of simulation n the units and epoch returned by the function get_time_units
        """

        st = []
        for key, value in self.bmimodels.items():
            st.append(self.bmimodels[key].get_current_time())

        return st[-1]

    def get_end_time(self):
        """
        Get the end time of the model run in units since the epoch

        :return: end time of simulation n the units and epoch returned by the function get_time_units
        """
        st = []
        for key, value in self.bmimodels.items():
            st.append(self.bmimodels[key].get_end_time())

        return np.array(st).min()

    def get_time_step(self):
        """
        Get the model time steps in units since the epoch

        :return: duration of one time step of the model with the largest! timestep.
        """
        st = []
        for key, value in self.bmimodels.items():
            st.append(self.bmimodels[key].get_time_step())

        return max(st)

    def get_time_units(self):
        """
        Return the time units of the model as a string

        :return: Return a string formatted using the UDUNITS standard from Unidata.
        (http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/cf-conventions.html#time-coordinate)
        """
        st = []
        for key, value in self.bmimodels.items():
            st.append(self.bmimodels[key].get_time_units())

        return st[-1]

    def get_value(self, long_var_name):
        """
        Get the value(s) of a variable as a numpy array

        :var long_var_name: name of the variable
        :return: a np array of long_var_name
        """
        # first part should be the component name
        self.bmilogger.debug("get_value: " + long_var_name)
        cname = long_var_name.split(self.comp_sep)
        if cname[0] in self.bmimodels:
            tmp = self.bmimodels[cname[0]].get_value(cname[1])
            if self.wrtodisk:
                pcr.report(
                    pcr.numpy2pcr(pcr.Scalar, tmp, -999),
                    long_var_name + "_get_" + str(self.get_current_time()) + ".map",
                )
            return tmp
        else:
            self.bmilogger.error("get_value: " + long_var_name + " returning None!!!!")
            return None

    def get_value_at_indices(self, long_var_name, inds):
        """
        Get a numpy array of the values at the given indices

        :var long_var_name: identifier of a variable in the model:
        :var inds: List of list each tuple contains one index for each dimension of the given variable, i.e. each tuple indicates one element in the multi-dimensional variable array:

        :return: numpy array of values in the data type returned by the function get_var_type.
        """
        cname = long_var_name.split(self.comp_sep)

        if cname[0] in self.bmimodels:
            return self.bmimodels[cname[0]].get_value_at_indices(cname[1], inds)
        else:
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
        cname = long_var_name.split(self.comp_sep)
        if cname[0] in self.bmimodels:
            self.bmimodels[cname[0]].set_value_at_indices(cname[1], inds, src)

    def get_grid_type(self, long_var_name):
        """
        Get the grid type according to the enumeration in BmiGridType

        :var String long_var_name: identifier of a variable in the model.

        :return: BmiGridType type of the grid geometry of the given variable.
        """
        cname = long_var_name.split(self.comp_sep)

        if cname[0] in self.bmimodels:
            return self.bmimodels[cname[0]].get_grid_type(cname[1])
        else:
            return None

    def get_grid_shape(self, long_var_name):
        """
        Return the shape of the grid. Only return something for variables with a uniform, rectilinear or structured grid. Otherwise raise ValueError.

        :var long_var_name: identifier of a variable in the model.

        :return: List of integers: the sizes of the dimensions of the given variable, e.g. [500, 400] for a 2D grid with 500x400 grid cells.
        """
        cname = long_var_name.split(self.comp_sep)

        if cname[0] in self.bmimodels:
            return self.bmimodels[cname[0]].get_grid_shape(cname[1])
        else:
            return None

    def get_grid_spacing(self, long_var_name):
        """
        Only return something for variables with a uniform grid. Otherwise raise ValueError.

        :var long_var_name: identifier of a variable in the model.

        :return: The size of a grid cell for each of the dimensions of the given variable, e.g. [width, height]: for a 2D grid cell.
        """
        cname = long_var_name.split(self.comp_sep)

        if cname[0] in self.bmimodels:
            return self.bmimodels[cname[0]].get_grid_spacing(cname[1])
        else:
            return None

    def get_grid_origin(self, long_var_name):
        """
        gets the origin of the model grid.

        :var String long_var_name: identifier of a variable in the model.

        :return: X, Y: ,the lower left corner of the grid.
        """
        cname = long_var_name.split(self.comp_sep)

        if cname[0] in self.bmimodels:
            return self.bmimodels[cname[0]].get_grid_origin(cname[1])
        else:
            return None

    def get_grid_x(self, long_var_name):
        """
        Give X coordinates of point in the model grid

        :var String long_var_name: identifier of a variable in the model.

        :return: Numpy array of doubles: x coordinate of grid cell center for each grid cell, in the same order as the
        values returned by function get_value.
        """
        cname = long_var_name.split(self.comp_sep)

        if cname[0] in self.bmimodels:
            return self.bmimodels[cname[0]].get_grid_x(cname[1])
        else:
            return None

    def get_grid_y(self, long_var_name):
        """
        Give Y coordinates of point in the model grid

        :var String long_var_name: identifier of a variable in the model.

        :return: Numpy array of doubles: y coordinate of grid cell center for each grid cell, in the same order as the
        values returned by function get_value.

        """
        cname = long_var_name.split(self.comp_sep)

        if cname[0] in self.bmimodels:
            return self.bmimodels[cname[0]].get_grid_y(cname[1])
        else:
            return None

    def get_grid_z(self, long_var_name):
        """
        Give Z coordinates of point in the model grid

        :var String long_var_name: identifier of a variable in the model.

        :return: Numpy array of doubles: z coordinate of grid cell center for each grid cell, in the same order as the values returned by function get_value.
        """
        cname = long_var_name.split(self.comp_sep)

        if cname[0] in self.bmimodels:
            return self.bmimodels[cname[0]].get_grid_z(cname[1])
        else:
            return None

    def get_var_units(self, long_var_name):
        """
        Supply units as defined in the API section of the ini file

        :var long_var_name: identifier of a variable in the model.

        :return:   String: unit of the values of the given variable. Return a string formatted
        using the UDUNITS standard from Unidata. (only if set properly in the ini file)
        """

        cname = long_var_name.split(self.comp_sep)

        if cname[0] in self.bmimodels:
            return self.bmimodels[cname[0]].get_var_units(cname[1])
        else:
            return None

    def set_value(self, long_var_name, src):
        """
        Set the values(s) in a map using a numpy array as source

        :var long_var_name: identifier of a variable in the model.
        :var src: all values to set for the given variable. If only one value
                  is present a uniform map will be set in the wflow model.
        """
        # first part should be the component name
        self.bmilogger.debug("set_value: " + long_var_name)
        cname = long_var_name.split(self.comp_sep)
        if cname[0] in self.bmimodels:
            self.bmimodels[cname[0]].set_value(cname[1], src)
            if self.wrtodisk:
                pcr.report(
                    pcr.numpy2pcr(pcr.Scalar, src, -999),
                    long_var_name + "_set_" + str(self.get_current_time()) + ".map",
                )

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
        
    def get_grid_rank(self, grid_id):
        """
        Not applicable raises NotImplementedError
        """
        raise NotImplementedError

    def get_grid_size(self, grid_id):
        """
        Not applicable raises NotImplementedError
        """
        raise NotImplementedError

    def get_var_grid(self, long_var_name):
        """
        Not applicable raises NotImplementedError
        """
        raise NotImplementedError
        
    def get_var_itemsize(self, long_var_name):
        """
        Not applicable raises NotImplementedError
        """
        raise NotImplementedError
