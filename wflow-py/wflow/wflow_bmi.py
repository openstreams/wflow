__author__ = 'schelle'

import os
import logging

import wflow.bmi as bmi

#TODO: Set log level also ini to be able to make quiet or non-quiet runs
#TODO: set re-init in the ini file to be able to make cold start runs
#TODO: Rework framework to get rid of max timesteps shit

class wflowbmi_ligth(object):

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
        elif "wflow_hbv.ini" in configfile:
            import wflow_sbm as wf
        elif "wflow_routing.ini" in configfile:
            import wflow_routing as wf
        else:
            raise NotImplementedError

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



    def update(self, dt):
        """
        Return type string, compatible with numpy.
        Propagate the model one timestep?
        """
        #curstep = self.dynModel.wf_
        if dt == -1:
            self.dynModel._runDynamic(self.currenttimestep, self.currenttimestep)
            self.currenttimestep = self.currenttimestep + 1
        else:
            nrsteps = int(dt/self.dynModel.timestepsecs)
            self.dynModel._runDynamic(self.currenttimestep, self.currenttimestep + nrsteps -1)
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

        return npmap.dtype


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
        return self.dynModel.wf_supplyMapAsNumpy(name)


    def set_var(self, name, var):
        """
        Set the variable name with the values of var
        Assume var is a numpy array
        """
        #TODO: check the numpy type
        self.dynModel.wf_setValuesAsNumpy(name, var)



    def set_var_slice(self, name, start, count, var):
        """
        Overwrite the values in variable name with data
        from var, in the range (start:start+count).
        Start, count can be integers for rank 1, and can be
        tuples of integers for higher ranks.
        For some implementations it can be equivalent and more efficient to do:
        `get_var(name)[start[0]:start[0]+count[0], ..., start[n]:start[n]+count[n]] = var`
        """
        tmp = self.get_var(name).copy()
        try:
            # if we have start and count as a number we can do this
            tmp[start:(start+count)] = var
        except:
            # otherwise we have to loop over all dimensions
            slices = [np.s_[i:(i+n)] for i,n in zip(start, count)]
            tmp[slices]
        self.set_var(name, name, tmp)


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
        tmp = self.get_var(name).copy()
        tmp.flat[index] = var
        self.set_var(name, name, tmp)




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



class wflowbmi_csdms(bmi.Bmi):

    def initialize(self, filename,loglevel=logging.DEBUG):
        """
        Assumptions for now:
        - the configfile wih be a full path
        - we define the case from the basedir of the configfile
        """
        retval = 0
        self.currenttimestep = 1
        wflow_cloneMap = 'wflow_subcatch.map'
        datadir = os.path.dirname(filename)
        inifile = os.path.basename(filename)
        runid = "run_default"
        # The current pcraster framework needs a max number of timesteps :-(
        # This cannot be avoided at the moment (needs a rework)
        # set to 10000 for now
        # .. todo::
        #       Get name of module from ini file name

        maxNrSteps = 10000
        if "wflow_sbm.ini" in filename:
            import wflow_sbm as wf
        elif "wflow_hbv.ini" in filename:
            import wflow_sbm as wf
        elif "wflow_routing.ini" in filename:
            import wflow_routing as wf
        else:
            raise NotImplementedError

        self.myModel = wf.WflowModel(wflow_cloneMap, datadir, runid, inifile)

        self.dynModel = wf.wf_DynamicFramework(self.myModel, maxNrSteps, firstTimestep = 1)
        self.dynModel.createRunId(NoOverWrite=0,level=loglevel,model=os.path.basename(filename))
        self.dynModel._runInitial()
        self.dynModel._runResume()
        return retval

    def update(self):
        """
        Propagate the model to the next timestep
        """

        self.dynModel._runDynamic(self.currenttimestep, self.currenttimestep)
        self.currenttimestep = self.currenttimestep + 1



    def update_until(self, time):
        """
        Update the model until the given time.

        :var  double time: time in the units and epoch returned by the function get_time_units.
        """

        timespan = time - self.get_current_time()
        nrsteps = int(timespan/self.dynModel.timestepsecs)
        self.dynModel._runDynamic(self.currenttimestep, self.currenttimestep + nrsteps -1)
        self.currenttimestep = self.currenttimestep + nrsteps

    def update_frac(self, time_frac):
        """
        No idea what to do with this one.....

        Input parameters:
        double time_frac: ???
        """
        raise NotImplementedError


    def save_state(self, destination_directory):
        """
        Ask the model to write its complete internal current state to one or more state files in the given directory.
        Afterwards the given directory should only contain the state files and nothing else.

        :var destination_directory: the directory in which the state files should be written.
        """

        self.dynModel.wf_suspend(destination_directory)

    def finalize(self):
        """
        Shutdown the library and clean up the model.
        Uses the default (model configured) state location to also save states.

        """
        self.dynModel._runSuspend()
        self.dynModel._wf_shutdown()


    def get_component_name(self):
        """
        :return  identifier of the model:
        """

        return self.dynModel.name


    def get_input_var_names(self):
        """
        :return List of String objects: identifiers of all input variables of the model:
        """
        namesroles = self.dynModel.wf_supplyVariableNamesAndRoles()

        # input variable are all forcing variables and all state variables
        #0 = input (to the model)
        #1 = is output (from the model)
        #2 = input/output (state information)
        #3 = model parameter
        inames = []

        for varrol in namesroles:
            if varrol[1] == 0 or varrol[1] == 2:
                inames.append(varrol[0])
        return inames

    def get_output_var_names(self):
        """
        :return List of String objects: identifiers of all output variables of the model:
        """
        namesroles = self.dynModel.wf_supplyVariableNamesAndRoles()

        # input variable are all forcing variables and all state variables
        #0 = input (to the model)
        #1 = is output (from the model)
        #2 = input/output (state information)
        #3 = model parameter
        inames = []

        for varrol in namesroles:
            if varrol[1] == 1 or varrol[1] == 2:
                inames.append(varrol[0])
        return inames


    def get_var_type(self, long_var_name):
        """
        :return type string, compatible with numpy:
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(long_var_name)

        return npmap.dtype


    def get_var_rank(self, long_var_name):
        """
        :var  String long_var_name: identifier of a variable in the model:

        :return array rank or 0 for scalar (number of dimensions):
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(long_var_name)

        return len(npmap.shape)

    def get_var_size(self, long_var_name):
        """
        :var  String long_var_name: identifier of a variable in the model:

        :return total number of values contained in the given variable (number of elements in map):
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(long_var_name)

        return npmap.size

    def get_var_nbytes(self, long_var_name):
        """
        :var  String long_var_name: identifier of a variable in the model:

        :return total number of bytes contained in the given variable (number of elements * bytes per element):
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(long_var_name)

        return npmap.size * npmap.itemsize


    def get_start_time(self):
        """
        :return start time in the units and epoch returned by the function get_time_units:
        """
        return self.dynModel.wf_supplyStartTime()

    def get_current_time(self):
        """
        :return current time of simulation n the units and epoch returned by the function get_time_units:
        """
        return self.dynModel.wf_supplyCurrentTime()


    def get_end_time(self):
        """
        :return end time of simulation n the units and epoch returned by the function get_time_units:
        """
        return self.dynModel.wf_supplyEndTime()

    def get_time_step(self):
        """
        :return duration of one time step of the model in the units returned by the function get_time_units:
        """
        return self.dynModel.timestepsecs

    def get_time_units(self):
        """
        :return: time units as a CF convention string
        (http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/cf-conventions.html#time-coordinate)
        """

        return self.dynModel.wf_supplyEpoch()

    def get_value(self, long_var_name):
        """
        :var long_var_name name of the variable
        :return an np array of long_var_name
        """
        return self.dynModel.wf_supplyMapAsNumpy(long_var_name)


    def get_value_at_indices(self, long_var_name, inds):
        """
        :var long_var_name: identifier of a variable in the model:
        :var List of list each tuple contains one index for each dimension of the given variable, i.e. each tuple indicates one element in the multi-dimensional variable array:

        :return  numpy array of values in the data type returned by the function get_var_type.
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(long_var_name)

        return npmap[inds]

    def set_value_at_indices(self, long_var_name, inds, src):
        """
        :var long_var_name: identifier of a variable in the model:
        :var inds: List of Lists of integers inds each nested List contains one index for each dimension of the given variable,
                                        i.e. each nested List indicates one element in the multi-dimensional variable array,
                                        e.g. [[0, 0], [0, 1], [15, 19], [15, 20], [15, 21]] indicates 5 elements in a 2D grid.:
        :var src: Numpy array of values. one value to set for each of the indicated elements:
        """

        npmap = self.dynModel.wf_supplyMapAsNumpy(long_var_name)
        npmap[inds] = src
        self.dynModel.wf_setValuesAsNumpy(long_var_name,npmap)


    def get_grid_type(self, long_var_name):
        """
        :var String long_var_name: identifier of a variable in the model.

        :return BmiGridType type of the grid geometry of the given variable:
        """

        ret=BmiGridType()

        return ret.UNIFORM


    def get_grid_shape(self, long_var_name):
        """
        Only return something for variables with a uniform, rectilinear or structured grid. Otherwise raise ValueError.

        :var long_var_name: identifier of a variable in the model.

        :return List of integers: the sizes of the dimensions of the given variable, e.g. [500, 400] for a 2D grid with 500x400 grid cells.
        """

        dim =  self.dynModel.wf_supplyGridDim()
        #[ Xul, Yul, xsize, ysize, rows, cols]

        return dim[4,5]

   def get_grid_spacing(self, long_var_name):
        """
        Only return something for variables with a uniform grid. Otherwise raise ValueError.

        :var long_var_name: identifier of a variable in the model.

        :return the size of a grid cell for each of the dimensions of the given variable, e.g. [width, height]: for a 2D grid cell.
        """
        raise NotImplementedError


    def get_grid_origin(self, long_var_name):
        """
        Only return something for variables with a uniform grid. Otherwise raise ValueError.

        Input parameters:
        String long_var_name: identifier of a variable in the model.

        Return value:
        List of doubles: the coordinate of the grid origin for each of the dimensions of the given variable. For a 2D grid this must be the lower left corner of the grid.
        """
        raise NotImplementedError

    def get_grid_x(self, long_var_name):
        """
        Only return something for variables with a rectilinear, structured or unstructured grid. Otherwise raise ValueError.

        Input parameters:
        String long_var_name: identifier of a variable in the model.

        Return value:
        Numpy array of doubles: x coordinate of grid cell center for each grid cell, in the same order as the values returned by function get_value.
                         For a rectilinear grid: x coordinate of column center for each column.
        """
        return self.myModel.pcr2numpy(self.myModel.xcoordinate(1),0.0)

    def get_grid_y(self, long_var_name):
        """
        Only return something for variables with a rectilinear, structured or unstructured grid. Otherwise raise ValueError.

        Input parameters:
        String long_var_name: identifier of a variable in the model.

        Return value:
        Numpy array of doubles: x coordinate of grid cell center for each grid cell, in the same order as the values returned by function get_value.
                         For a rectilinear grid: x coordinate of column center for each column.
        """
        return self.myModel.pcr2numpy(self.myModel.ycoordinate(1),0.0)


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


    def get_var_shape(self, name):
        """
        Return shape of the array.
        """
        npmap = self.dynModel.wf_supplyMapAsNumpy(name)

        return npmap.shape





    def set_var(self, name, var):
        """
        Set the variable name with the values of var
        Assume var is a numpy array
        """
        #TODO: check the numpy type
        self.dynModel.wf_setValuesAsNumpy(name, var)


    def set_var_slice(self, name, start, count, var):
        """
        Overwrite the values in variable name with data
        from var, in the range (start:start+count).
        Start, count can be integers for rank 1, and can be
        tuples of integers for higher ranks.
        For some implementations it can be equivalent and more efficient to do:
        `get_var(name)[start[0]:start[0]+count[0], ..., start[n]:start[n]+count[n]] = var`
        """
        tmp = self.get_var(name).copy()
        try:
            # if we have start and count as a number we can do this
            tmp[start:(start+count)] = var
        except:
            # otherwise we have to loop over all dimensions
            slices = [np.s_[i:(i+n)] for i,n in zip(start, count)]
            tmp[slices]
        self.set_var(name, name, tmp)


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
        tmp = self.get_var(name).copy()
        tmp.flat[index] = var
        self.set_var(name, name, tmp)



class BmiGridType(object):
    UNKNOWN = 0
    UNIFORM = 1
    RECTILINEAR = 2
    STRUCTURED = 3
    UNSTRUCTURED = 4

