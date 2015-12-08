__author__ = 'schelle'

import os
import logging
import datetime
import parser

import wflow.bmi as bmi
import numpy as np
from wflow.pcrut import setlogger

# TODO: Set log level also ini to be able to make quiet or non-quiet runs
# TODO: set re-init in the ini file to be able to make cold start runs
# TODO: Rework framework to get rid of max timesteps shit

class wflowbmi_ligth(object):
    """
    Deltares specific light version of the BMI. Used for internal model linkage
    """

    def __init__(self):
        """

        :return:
        """
        self.loggingmode = logging.ERROR
        logstr = os.getenv('wflow_bmi_loglevel', 'ERROR')

        if logstr in 'ERROR':
            self.loggingmode = logging.ERROR
        if logstr in 'WARNING':
            self.loggingmode = logging.WARNING
        if logstr in 'INFO':
            self.loggingmode = logging.INFO
        if logstr in 'DEBUG':
            self.loggingmode = logging.DEBUG

        self.bmilogger = setlogger('wflow_bmi.log','wflow_bmi_logging',thelevel=self.loggingmode)
        self.bmilogger.info("__init__: wflow_bmi object initialised.")

    def initialize(self, configfile=None,loglevel=logging.DEBUG):
        """
        Assumptions for now:
        - the configfile wih be a full path
        - we define the case from the basedir of the configfile
        """
        retval = 0
        self.currenttimestep = 1
        wflow_cloneMap = 'wflow_subcatch.map'
        datadir = os.path.dirname(configfile)
        inifile = os.path.basename(configfile)
        runid = "run_default"
        # The current pcraster framework needs a max number of timesteps :-(
        # This cannot be avoided at the moment (needs a rework)
        # set to 10000 for now
        # .. todo::
        #       Get name of module from ini file name

        maxNrSteps = 10000
        if "wflow_sbm.ini" in configfile:
            import wflow_sbm as wf
            self.name = "wflow_sbm"
        elif "wflow_hbv.ini" in configfile:
            import wflow_hbv as wf
            self.name = "wflow_hbv"
        elif "wflow_routing.ini" in configfile:
            import wflow_routing as wf
            self.name = "wflow_routing"
        else:
            modname = os.path.splitext(os.path.basename(configfile))[0]
            exec "import wflow." + modname + " as wf"
            self.name = modname

        self.bmilogger.info("initialize: Initialising wflow bmi with ini: " + configfile)
        myModel = wf.WflowModel(wflow_cloneMap, datadir, runid, inifile)

        self.dynModel = wf.wf_DynamicFramework(myModel, maxNrSteps, firstTimestep = 1)
        self.dynModel.createRunId(NoOverWrite=0,level=loglevel,model=os.path.basename(configfile))
        self.dynModel._runInitial()
        self.dynModel._runResume()

        return retval

    def finalize(self):
        """
        Shutdown the library and clean up the model.

        """
        self.dynModel._runSuspend()
        self.dynModel._wf_shutdown()
        self.bmilogger.debug("finalize: shutting done bmi")

    def update(self, dt):
        """
        Return type string, compatible with numpy.
        Propagate the model dt timesteps or ( if dt == -1) to the end of the model run
        """
        # TODO: fix dt = -1 problem, what do we want here?
        #curstep = self.dynModel.wf_
        if dt != -1:
            self.dynModel._runDynamic(self.currenttimestep, self.currenttimestep + dt)
            self.currenttimestep = self.currenttimestep + dt
        else:
            nrsteps = int(dt/self.dynModel.timestepsecs)
            self.dynModel._runDynamic(self.currenttimestep, self.currenttimestep + nrsteps)
            self.currenttimestep = self.currenttimestep + nrsteps

    def get_time_units(self):
        """

        :return: time units as a CF convention string
        (http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/cf-conventions.html#time-coordinate)
        """

        return self.dynModel.wf_supplyEpoch()

    def get_var_count(self):
        """
        Return number of variables
        """
        return self.dynModel.wf_supplyVariableCount()

    def get_var_name(self, i):
        """
        Return variable name
        """

        names = self.dynModel.wf_supplyVariableNames()
        return names[i]

    def get_var_type(self, name):
        """
        Return type string, compatible with numpy.
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(name)

        if hasattr(npmap,'dtype'):
            return str(npmap.dtype)
        else:
            return None

    def get_var_rank(self, name):
        """
        Return array rank or 0 for scalar.
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(name)

        return len(npmap.shape)

    def get_var_shape(self, name):
        """
        Return shape of the array.
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(name)

        return npmap.shape

    def get_start_time(self):
        """
        returns start time
        """
        return self.dynModel.wf_supplyStartTime()

    def get_end_time(self):
        """
        returns end time of simulation
        """
        return self.dynModel.wf_supplyEndTime()

    def get_current_time(self):
        """
        returns current time of simulation
        """
        return self.dynModel.wf_supplyCurrentTime()

    def get_var(self, name):
        """
        Return an nd array from model library
        """
        return np.flipud(self.dynModel.wf_supplyMapAsNumpy(name))

    def set_var(self, name, var):
        """
        Set the variable name with the values of var
        Assume var is a numpy array
        """
        #TODO: check the numpy type
        self.dynModel.wf_setValuesAsNumpy(name, np.flipud(var))


    def set_var_slice(self, name, start, count, var):
        """
        Overwrite the values in variable name with data
        from var, in the range (start:start+count).
        Start, count can be integers for rank 1, and can be
        tuples of integers for higher ranks.
        For some implementations it can be equivalent and more efficient to do:
        `get_var(name)[start[0]:start[0]+count[0], ..., start[n]:start[n]+count[n]] = var`
        """
        tmp = np.flipud(self.get_var(name).copy())
        try:
            # if we have start and count as a number we can do this
            tmp[start:(start+count)] = var
        except:
            # otherwise we have to loop over all dimensions
            slices = [np.s_[i:(i+n)] for i,n in zip(start, count)]
            tmp[slices]

        self.set_var(name, name, np.flipud(tmp))

    def set_var_index(self, name, index, var):
        """
        Overwrite the values in variable "name" with data
        from var, at the flattened (C-contiguous style)
        indices. Indices is a vector of 0-based
        integers, of the same length as the vector var.
        For some implementations it can be equivalent
        and more efficient to do:
        `get_var(name).flat[index] = var`
        """
        tmp = np.flipud(self.get_var(name).copy())
        tmp.flat[index] = var
        self.set_var(name, name, np.flipud(tmp))

    def inq_compound(self, name):
        """
        Return the number of fields of a compound type.
        """
        pass

    def inq_compound_field(self, name, index):
        """
        Lookup the type,rank and shape of a compound field
        """
        pass



class LookupNames():
    """

    """
    def __init__(self, filename):
        """
        :param filename: filename with the translation table, format: long_var_name:model_var_name
        :return: nothing
        """
        raise NotImplementedError


class wflowbmi_csdms(bmi.Bmi):
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
        self.currenttimestep = 0
        self.name = "undefined"
        self.myModel = None
        self.dynModel = None

        self.loggingmode = logging.ERROR
        logstr = os.getenv('wflow_bmi_loglevel', 'ERROR')

        if logstr in 'ERROR':
            self.loggingmode = logging.ERROR
            print logstr
        if logstr in 'WARNING':
            self.loggingmode = logging.WARNING
            print logstr
        if logstr in 'INFO':
            self.loggingmode = logging.INFO
            print logstr
        if logstr in 'DEBUG':
            self.loggingmode = logging.DEBUG
            print logstr

        self.bmilogger = setlogger('wflow_bmi.log','wflow_bmi_logging',thelevel=self.loggingmode)
        self.bmilogger.info("__init__: wflow_bmi object initialised.")

    def initialize_config(self, filename, loglevel=logging.DEBUG):
        """
        *Extended functionality*, see https://github.com/eWaterCycle/bmi/blob/master/src/main/python/bmi.py

        see initialize

        :param filename:
        :param loglevel:
        :return: nothing
        """
        def get_output_only_var_names(self):
            """
            Returns the list of model output only variables. This is not BMI but needed to check

            :return: List of String objects: identifiers of all output variables of the model:
            """
            namesroles = self.dynModel.wf_supplyVariableNamesAndRoles()
            inames = []

            for varrol in namesroles:
                if varrol[1] == 1:
                    inames.append(varrol[0])

            return inames


        self.currenttimestep = 1
        wflow_cloneMap = 'wflow_subcatch.map'
        self.datadir = os.path.dirname(filename)
        inifile = os.path.basename(filename)
        runid = "run_default"
        # The current pcraster framework needs a max number of timesteps :-(
        # This cannot be avoided at the moment (needs a rework)
        # set to 10000 for now
        #
        maxNrSteps = 10000
        self.bmilogger.info("initialize_config: Initialising wflow bmi with ini: " + filename)

        if "wflow_sbm" in filename:
            import wflow.wflow_sbm as wf
            self.name = "wflow_sbm"
        elif "wflow_hbv" in filename:
            import wflow.wflow_hbv as wf
            self.name = "wflow_hbv"
        elif "wflow_routing" in filename:
            import wflow.wflow_routing as wf
            self.name = "wflow_routing"
        else:
            modname = os.path.splitext(os.path.basename(filename))[0]
            exec "import wflow." + modname + " as wf"
            self.name = modname

        self.myModel = wf.WflowModel(wflow_cloneMap, self.datadir, runid, inifile)

        self.dynModel = wf.wf_DynamicFramework(self.myModel, maxNrSteps, firstTimestep = 1)
        self.dynModel.createRunId(doSetupFramework=False,NoOverWrite=0,level=loglevel,model=os.path.basename(filename))

        self.outputonlyvars = get_output_only_var_names(self)
        self.inputoutputvars = self.get_output_var_names()

    def initialize_model(self):
        """
        *Extended functionality*, see https://github.com/eWaterCycle/bmi/blob/master/src/main/python/bmi.py

        see initialize

        :param self:
        :return: nothing
        """
        self.bmilogger.info("initialize_model: Initialising wflow bmi with ini, loading initial state")
        self.dynModel.setupFramework()
        self.dynModel._runInitial()
        self.dynModel._runResume()


    def set_start_time(self, start_time):
        """

        :param start_time: time in units (seconds) since the epoch
        :return: nothing
        """

        dateobj = datetime.datetime.utcfromtimestamp(start_time)
        datestrimestr = dateobj.strftime("%Y-%m-%d %H:%M:%S")

        self.dynModel._userModel().config.set("run",'starttime',datestrimestr)
        self.dynModel._userModel().datetime_firststep=dateobj
        self.dynModel.datetime_firststep=dateobj
        self.dynModel._userModel().currentdatetime = self.dynModel._userModel().datetime_firststep
        self.bmilogger.debug("set_start_time: " + str(start_time) + " " + str(datestrimestr))

    def set_end_time(self, end_time):
        """
        :param end_time: time in units (seconds) since the epoch
        :return:
        """

        dateobj = datetime.datetime.utcfromtimestamp(end_time)
        datestrimestr = dateobj.strftime("%Y-%m-%d %H:%M:%S")
        self.dynModel._userModel().config.set("run",'endtime',datestrimestr)
        self.dynModel._userModel().datetime_laststep=dateobj
        self.dynModel.datetime_laststep=dateobj
        self.bmilogger.debug("set_end_time: " + str(end_time) + " " + str(datestrimestr))



    def get_attribute_names(self):
        """
        Get the attributes of the model return in the form of section_name:attribute_name

        :return: list of attributes
        """
        self.bmilogger.debug("get_attribute_names")
        attr = []
        for sect in self.dynModel._userModel().config.sections():
            for opt in self.dynModel._userModel().config.options(sect):
                attr.append(sect + ":" + opt)
        return attr

    def get_attribute_value(self, attribute_name):
        """
        :param attribute_name:
        :return: attribute value
        """
        self.bmilogger.debug("get_attribute_value: " + attribute_name)
        attrpath = attribute_name.split(":")

        if len(attrpath) == 2:
            return self.dynModel._userModel().config.get(attrpath[0],attrpath[1])
        else:
            raise Warning("attributes should follow the name:option  convention")

    def set_attribute_value(self, attribute_name, attribute_value):
        """
        :param attribute_name: name using the section:option notation
        :param attribute_value: string value of the option
        :return:
        """
        self.bmilogger.debug("set_attribute_value: " + attribute_value)
        attrpath = attribute_name.split(":")
        if len(attrpath) == 2:
            self.dynModel._userModel().config.set(attrpath[0],attrpath[1],attribute_value)
        else:
            self.bmilogger.warn("Attributes should follow the name:option  convention")
            raise Warning("attributes should follow the name:option  convention")



    def initialize(self, filename,loglevel=logging.DEBUG):
        """
        Initialise the model. Should be called before any other method.

        :var filename: full path to the wflow ini file
        :var loglevel: optional loglevel (default == DEBUG)

        Assumptions for now:

            - the configfile wih be a full path
            - we define the case from the basedir of the configfile

        .. todo::

            Get name of module from ini file name
        """

        self.bmilogger.info("initialize: Initialising wflow bmi with ini: " + filename)
        self.initialize_config(filename,loglevel=loglevel)
        self.initialize_model()
        self.dynModel.wf_resume(os.path.join(self.datadir,'instate'))

    def update(self):
        """
        Propagate the model to the next model timestep
        """
        self.bmilogger.debug('update: update one timestep: ' + str(self.currenttimestep) + ' to ' + str(self.currenttimestep + 1))
        self.dynModel._runDynamic(self.currenttimestep, self.currenttimestep)
        self.currenttimestep = self.currenttimestep + 1

    def update_until(self, time):
        """
        Update the model until and including the given time.

        - one or more timesteps foreward
        - max one timestep backward

        :var  double time: time in the units and epoch returned by the function get_time_units.
        """
        curtime = self.get_current_time()

        if abs(time - curtime)% self.dynModel.timestepsecs != 0:
            self.bmilogger.error('update_until: timespan not dividable by timestep: ' + str(abs(time - curtime)) +
                                 ' and ' + str(self.dynModel.timestepsecs))
            raise ValueError("Update in time not a multiple of timestep")

        if curtime > time:
            timespan = curtime - time
            nrstepsback = int(timespan/self.dynModel.timestepsecs)
            self.bmilogger.debug('update_until: update timesteps back ' + str(nrstepsback) + ' to ' + str(curtime))
            if nrstepsback > 1:
                raise ValueError("Time more than one timestep before current time.")
            self.dynModel.wf_QuickResume()
        else:
            timespan = time - curtime
            nrsteps = int(timespan/self.dynModel.timestepsecs)
            self.bmilogger.debug('update_until: update timesteps foreward ' + str(nrsteps) + ' to ' + str(curtime))
            self.dynModel._runDynamic(self.currenttimestep, self.currenttimestep + nrsteps)
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
        self.bmilogger.debug("save_state: " + destination_directory)
        self.dynModel.wf_suspend(destination_directory)

    def load_state(self, source_directory):
        """
        Ask the model to load its complete internal current state from one or more
        state files in the given directory.

        :var  source_directory: the directory from which the state files should be
        read.
        """

        if os.path.isabs(source_directory):
            new_source_directory = source_directory
        else:
            new_source_directory = os.path.join(self.datadir,source_directory)

        self.bmilogger.debug("load_state: " + new_source_directory)
        self.dynModel.wf_resume(new_source_directory)

    def finalize(self):
        """
        Shutdown the library and clean up the model.
        Uses the default (model configured) state location to also save states.
        """
        # First check if the seconf initilize_states has run
        self.bmilogger.info("finalize.")
        if hasattr(self.dynModel,"framework_setup"):
            self.dynModel._runSuspend()
            self.dynModel._wf_shutdown()


    def get_component_name(self):
        """
        :return:  identifier of the model based on the name of the ini file
        """
        self.bmilogger.debug("get_component_name: " + str(self.name))
        return self.name

    def get_input_var_names(self):
        """

        :return: List of String objects: identifiers of all input variables of the model:
        """
        namesroles = self.dynModel.wf_supplyVariableNamesAndRoles()

        # input variable are all forcing variables and all state variables
        # 0 = input (to the model)
        # 1 = is output (from the model)
        # 2 = input/output (state information)
        # 3 = model parameter
        inames = []

        for varrol in namesroles:
            if varrol[1] == 0 or varrol[1] == 2:
                inames.append(varrol[0])

        self.bmilogger.debug("get_input_var_names: " + str(inames))
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

        self.bmilogger.debug("get_var_rank: " + str(len(npmap.shape)))
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
        return self.dynModel.wf_supplyCurrentTime()

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
        ts = self.dynModel.timestepsecs
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

class BmiGridType(object):
    UNKNOWN = 0
    UNIFORM = 1
    RECTILINEAR = 2
    STRUCTURED = 3
    UNSTRUCTURED = 4

