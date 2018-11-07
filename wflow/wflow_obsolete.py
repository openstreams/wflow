from .wflow_hbv import *


class wflow_model:
    def initialize(self, timeStepsToRun, timeStepInSeconds, dataDirectory, iniFile):
        self.currentTime = 1
        runId = "run_default"
        configfile = iniFile
        wflow_cloneMap = "wflow_subcatch.map"
        # Initializes the model class
        myModel = WflowModel(wflow_cloneMap, dataDirectory, runId, configfile)
        myModel.timestepsecs = timeStepInSeconds
        self.dynModelFw = wf_DynamicFramework(myModel, timeStepsToRun, firstTimestep=1)
        self.dynModelFw.createRunId(NoOverWrite=0)
        self.dynModelFw._runInitial()
        self.dynModelFw._runResume()

    def run_last_time_step(self):
        """
        Redo the last timestep
        """
        self.dynModelFw.wf_QuickResume()
        self.currentTime = self.currentTime - 1
        self.dynModelFw._runDynamic(self.currentTime, self.currentTime)
        self.currentTime = self.currentTime + 1

    def get_grid_parameters(self):
        """
        return the dimension of the current model grid as list
      
        [Xul, Yul, xsize, ysize, rows, cols]
        """
        return self.dynModelFw.wf_supplyGridDim()

    def get_variable_values(self, name):
        return self.dynModelFw.wf_supplyMapAsList(name)

    def get_variable_names(self):
        return self.dynModelFw.wf_supplyVariableNames()

    def get_variable_roles(self):
        return self.dynModelFw.wf_supplyVariableRoles()

    def get_variable_units(self):
        return self.dynModelFw.wf_supplyVariableUnits()

    def set_variable_values(self, name, values):
        self.dynModelFw.wf_setValues(name, values)

    def run_time_step(self):
        self.dynModelFw._runDynamic(self.currentTime, self.currentTime)
        self.currentTime = self.currentTime + 1

    def finalize(self):
        self.dynModelFw._runSuspend()
