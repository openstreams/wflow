__author__ = 'schelle'

import os
import logging

#TODO: Set log level also ini to be able to make quiet runs
#TODO: set re-init in the ini file to be able to make cold start runs
#TODO: Rework framework to get rid of max timesteps shit

class wflowbmi(object):

    def initialize(self, configfile=None):
        """
        Assumptions for now:
        - the configfile wih be a full path
        - we define the case from the basedir of the configfile
        """
        retval = 0;
        wflow_cloneMap = 'wflow_subcatch.map'
        datadir = os.path.dirname(configfile)
        inifile = os.path.basename(configfile)
        runid = "run_default"
        # The current pcraster framework needs a max number of timesteps :-(
        # This cannot be avoided at the moment (needs a rework)
        # set to 10000 for now
        maxNrSteps = 10000
        if "_sbm" in configfile:
            import wflow_sbm as wf
        elif "hbv" in configfile:
            import wflow_sbm as wf

        myModel = wf.WflowModel(wflow_cloneMap, datadir, runid, inifile)

        self.dynModel = wf.wf_DynamicFramework(myModel, maxNrSteps, firstTimestep = 1)
        self.dynModel.createRunId(NoOverWrite=0,level=logging.ERROR)
        self.dynModel._runInitial()
        self.dynModel._runResume()


        return retval

    def finalize(self):
        """
        Shutdown the library and clean up the model.
        Note that the Fortran library's cleanup code is not up to snuff yet,
        so the cleanup is not perfect. Note also that the working directory is
        changed back to the original one.
        """
        self.dynModel._runSuspend()
        self.dynModel._wf_shutdown()



    def update(self, dt):
        """
        Return type string, compatible with numpy.
        """
        pass


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
