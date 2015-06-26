# -*- coding: utf-8 -*-
import csv
import math
import os
import random
import shutil
import sys
import numpy
import frameworkBase
import mcFramework



class ParticleFilterFramework(frameworkBase.FrameworkBase):
  """
  Framework class for the particle filter method.

  `userModel`
    Instance that models the :ref:`Particle Filter Model Concept <particleFilterModelConcept>`.
  """

  def __init__(self,
    userModel):
    frameworkBase.FrameworkBase.__init__(self)
    self._d_model = userModel
    self._testRequirements()
    self._d_totalTimesteps = self._userModel().nrTimeSteps()
    self._d_trackCloned = {}

    self._resetSampleWeights()

    self._addMethodToClass(self._setParticleWeight)
    self._addMethodToClass(self._runPremcloop)
    self._addMethodToClass(self._runPostmcloop)


    # \todo !!!test if filter timesteps are in interval of model timesteps...

  def _userModel(self):
    return self._d_model._userModel()

  ## \brief Returns the framework provided by the user
  def _userFramework(self):
    return self._d_model

  def _testRequirements(self):
    """
    Test whether the user model models the
    :ref:`Particle Filter Model Concept <particleFilterModelConcept>`.

    .. todo::

       The implementation should just perform the concept check and be done
       with it. Don't use framework code.
    """
    if hasattr(self._userModel(), "_userModel"):
      msg = "The _userModel method is deprecated and obsolete"
      self.showWarning(msg)
    #\todo test to dynamic framework model
    if not isinstance(self._userFramework(), mcFramework.MonteCarloFramework):
      msg = "Model must be instance of MonteCarloFramework"
      raise frameworkBase.FrameworkError(msg)

    if not hasattr(self._d_model, 'run'):
      msg = "No 'run' section defined."
      raise frameworkBase.FrameworkError(msg)

    if not hasattr(self._userModel(), 'updateWeight'):
      msg = "Cannot run particle filter framework: Implement 'updateWeight' method"
      raise frameworkBase.FrameworkError(msg)

    if not hasattr(self._userModel(), 'suspend'):
      msg = "Cannot run particle filter framework: Implement 'suspend' method"
      raise frameworkBase.FrameworkError(msg)

    if not hasattr(self._userModel(), 'resume'):
      msg = "Cannot run particle filter framework: Implement 'resume' method"
      raise frameworkBase.FrameworkError(msg)

  def _particleWeights(self):
    return self._userModel()._d_particleWeights

  ## \brief Creates the subdirectories for state variables
  # \todo test if mc dirs are there...
  def _initialiseSampleDirectories(self):
    sample = self._userModel()._firstSampleNumber()
    while sample <= self._userModel()._lastSampleNumber():
      cwd = os.getcwd()

      dirname = "%d" % (sample)
      varName = "stateVar"
      os.chdir(dirname)

      if not os.path.isdir(varName):
        # Create sample directory.
        os.mkdir(varName)
      else :
        #if not os.path.isdir(varName):
        #  # Remove existing file with name of sample directory.
        os.remove(varName)
        os.mkdir(varName)

      os.chdir(cwd)
      assert os.path.exists(os.path.join(dirname,"stateVar")) and \
        os.path.isdir(os.path.join(dirname,"stateVar"))
      sample += 1

  def setFilterTimesteps(self, filterTimesteps):
    """
    Set the filter moments.
    """
    assert isinstance(filterTimesteps, (list, numpy.ndarray))
    # \todo assert some more
    for filtertimestep in filterTimesteps:
      assert filtertimestep <= self._userModel().nrTimeSteps()
      assert filtertimestep > 0

    self._userModel()._d_filterTimesteps = filterTimesteps

  def filterTimesteps(self):
    """
    Return a list of filter moments.
    """
    return self._userModel()._d_filterTimesteps

  def run(self):
    """
    Run the user model in the filter mode.

    Re-implemented from ShellScript.
    """
    if(hasattr(self._userModel(), 'run')):
      self._userModel().run()
    else:
      self._atStartOfScript()
      self._initialiseSampleDirectories()

      lastPeriod = len(self._userModel()._d_filterTimesteps)

      if lastPeriod == 0:
        self.showError("No filter timesteps specified")
        sys.exit()



      # set the proposal/initial weight distribution by user
      if hasattr(self._userModel(), 'setInitialParticleWeights'):
        self._userModel()._d_particleWeights = self._userModel().setInitialParticleWeights()

        # check initial weights
        assert type(self._particleWeights()) == list
        assert len(self._particleWeights()) == self._userModel().nrSamples()

        for i in range(0, len(self._particleWeights())):
          assert type(self._particleWeights()[i]) == float

      # run the premc loop
      self._userModel()._runPremcloop()

      # looping over the filter periods
      for currentPeriod in range(0, len(self._userModel()._d_filterTimesteps) + 1):

        # \todo replace with a better solution...
        sumW = sum(self._particleWeights())
        assert abs(sumW - 1.0) < 0.00001

        # here we never execute premc and postmc
        self._runMonteCarlo(currentPeriod, lastPeriod)




        if not currentPeriod == lastPeriod:
          # update the weights
          # calling the "objective fuction" for each sample
          for sample in range(1, self._userModel().nrSamples() + 1):
            self._userModel()._setCurrentSample(sample)
            self._userModel()._d_inUpdateWeight = True
            fitnessValue = self._userModel().updateWeight()
            self._userModel()._d_inUpdateWeight = False
            assert type(fitnessValue) == float or type(fitnessValue) == int
            self._userModel()._d_particleWeights[sample - 1] = fitnessValue


          # determine samples to clone
          samplesToClone = self._samplesToClone(self._particleWeights())
          assert sum(samplesToClone) == self._userModel().nrSamples()


          # clone the data
          self._cloneData(samplesToClone)


          # reset the sample weights
          self._resetSampleWeights()

          currentPeriod += 1
          self._userModel()._d_filterPeriod += 1

    # run the postmc loop
    self._userModel()._setFirstTimeStep(self._userModel().firstTimeStep()) #1)
    self._userModel()._runPostmcloop()
    self._createGraph()
    return 0

  def _createGraph(self):
    oFile = open("samples.dot","w")
    oFile.write("digraph G {\n")
    oFile.write("nodesep=.002;\n")
    oFile.write("ranksep=1.00;\n")
    oFile.write("rankdir=LR;\n")
    oFile.write("node [shape=plaintext];\n")
    for i in range(0, len(self._userModel()._d_filterTimesteps) + 1):
      oFile.write("{range=same;")
      for j in range(1, self._userModel().nrSamples() + 1):
        oFile.write("\"%d-%d\";" % (i, j))
      oFile.write("}\n")

    for key in self._d_trackCloned.keys():
      if self._d_trackCloned[key] == 0:
        #oFile.write("\"%s\";\n" % (key))
        pass
      elif not type(self._d_trackCloned[key]) == list:
        oFile.write("\"%s\" -> \"%s\";\n" % (key, self._d_trackCloned[key]))
      else:
        for i in self._d_trackCloned[key]:
          oFile.write("\"%s\" -> \"%s\";\n" % (key, i)) #self._d_trackCloned[key]))
    oFile.write("}")
    oFile.close()

  def _cloneData(self,
    samplesToClone):
    # determine 'dead' samples
    notClonedSamples = []
    for i in range(1, len(samplesToClone) + 1):
      if samplesToClone[i - 1] == 0:
        notClonedSamples.append(i)
        # record the samples which are deleted
        el1 = str(self._userModel()._d_filterPeriod) + "-"+ str(i)
        el2 = 0
        self._d_trackCloned[el1] = el2
      elif samplesToClone[i - 1] == 1:
        # record the samples which are continued, but not cloned
        el1 = str(self._userModel()._d_filterPeriod) + "-"+ str(i)
        el2 = str(self._userModel()._d_filterPeriod + 1) + "-"+ str(i)
        self._d_trackCloned[el1] = el2

    # rm state in 'dead' samples
    for i in range(0, len(notClonedSamples)):
      dirName = os.path.join("%d" % notClonedSamples[i], "stateVar")
      shutil.rmtree(dirName)

    clonedSamples = []
    for i in range(1, len(samplesToClone) + 1):
      if samplesToClone[i - 1] > 1:
        clonedSamples.append(i)

    # each cloned sample
    for i in range(0, len(clonedSamples)):
      # record which sample is cloned
      clonedFrom = str(self._userModel()._d_filterPeriod) + "-"+ str(clonedSamples[i])
      clonedTo = []
      # how often current sample cloned? Do not clone yourself
      for j in range(1, int(samplesToClone[clonedSamples[i] - 1])):
        # \todo dir or subdir
        source = os.path.join(str(clonedSamples[i]), "stateVar")

        sample = str(notClonedSamples.pop())
        dest = os.path.join(sample, "stateVar")
        shutil.copytree(source, dest)

        if j == 1:
          # record the cloned sample as 'continuing itself'
          clonedTo.append(str(self._userModel()._d_filterPeriod + 1) + "-" + str(clonedSamples[i]))
          # record the clone
          clonedTo.append(str(self._userModel()._d_filterPeriod + 1) + "-"+ sample)
        else:
          # record the clone
          clonedTo.append(str(self._userModel()._d_filterPeriod + 1) + "-"+ sample)
      self._d_trackCloned[clonedFrom] = clonedTo

  def _normaliseWeights(self,
    weights):
    assert weights
    sumWeights = sum(weights)
    if sumWeights == 0:
      norm = [1.0 / self._userModel().nrSamples()] * self._userModel().nrSamples()
    else:
      norm = [0.0] * len(weights)
      for i in range(0, len(weights)):
        norm[i] =  weights[i] / sumWeights

    return norm

  def _resetSampleWeights(self):
    assert self._userModel().nrSamples() > 0
    self._userModel()._d_particleWeights = [1.0 / self._userModel().nrSamples()] * self._userModel().nrSamples()

  def _cumulativeWeights(self,
    weights):
    cumulative = [0.0] * self._userModel().nrSamples()
    value = 0.0
    for i in range(len(weights)):
      value += weights[i]
      cumulative[i] = value
    return cumulative

  def _startEndOfPeriod(self,
    currentPeriod,
    lastPeriod):
    # determine start end end timestep of current period
    if currentPeriod == 0:
      startTimestep = 1
      endTimestep = self._userModel()._d_filterTimesteps[currentPeriod]
    elif currentPeriod == lastPeriod:
      startTimestep = self._userModel()._d_filterTimesteps[currentPeriod -1] + 1
      endTimestep = self._d_totalTimesteps
    else:
      startTimestep = self._userModel()._d_filterTimesteps[currentPeriod - 1] + 1
      endTimestep = self._userModel()._d_filterTimesteps[currentPeriod]

    assert startTimestep <= endTimestep
    return startTimestep, endTimestep

  def _executePrePostMc(self,
    currentPeriod,
    lastPeriod):
    if currentPeriod == 0:
      # execute premc
      premc = True
      postmc = False
    elif currentPeriod == lastPeriod:
      # execute postmc
      premc = False
      postmc = True
    else:
      # without pre/postmc
      premc = False
      postmc = False
    # \todo assert something
    return premc, postmc

  def _runMonteCarlo(self,
    currentPeriod,
    lastPeriod):
    #  get user model and (re)set start and end time
    startTimestep, endTimestep = self._startEndOfPeriod(currentPeriod, lastPeriod)

    self._userModel()._setNrTimeSteps(endTimestep)
    self._userModel()._setFirstTimeStep(startTimestep)
    # \todo test this one, why is this endTimestep?!
    # otherwise timestep == 0 in updateWeight?!
    self._userModel()._setCurrentTimeStep(endTimestep)

    # determine MC execution with(out) pre/postloop
    #premc, postmc = self._executePrePostMc(currentPeriod, lastPeriod)

    # run the model in mc mode for current filter period
    self._incrementIndentLevel()
    self._atStartOfFilterPeriod(currentPeriod)
    self._d_model.run(False, False)
    self._atEndOfFilterPeriod()
    self._decrementIndentLevel()

  def particleWeight(self,
    sample):
    """
    Return the weight of a particle.
    """
    assert self._d_particleWeights
    assert sample >= self._firstSampleNumber()
    assert sample <= self._lastSampleNumber()
    return self._d_particleWeights[sample - 1]

  def _setParticleWeight(self,
    sample,
    weight):
    assert self._d_particleWeights
    assert sample >= self._firstSampleNumber()
    assert sample <= self._lastSampleNumber()
    self._d_particleWeights[sample - 1] = weight

  def readmap(self, name):
    """
    Read sample data from disk.

    Returns the map of the current time step from the current sample directory.
    """
    return self._readmapNew(name)

  def readDeterministic(self, name):
    """
    Read deterministic data from disk.

    Returns the map of the current time step from the current working directory.
    """
    if self._userModel()._inPremc() or self._inPostmc() or self._inInitial():
      newName = name + ".map"
    else:
      newName = frameworkBase.generateNameT(name, self._userModel().currentTimeStep())
    import pcraster
    return pcraster.readmap(newName)



## \brief Sequential importance resampling algorithm
class SequentialImportanceResamplingFramework(ParticleFilterFramework):
  ## \brief Constructor
  def __init__(self, userModel):
    ParticleFilterFramework.__init__(self, userModel)
    self._addMethodToClass(self.optimalSampleNumber)


  def _samplesToClone(self, weights):
    # normalise weights
    normalisedWeights = self._normaliseWeights(self._particleWeights())

    cumulativeWeights = self._cumulativeWeights(normalisedWeights)

    samplesToClone = [0] * self._userModel().nrSamples()
    for i in range(0, self._userModel().nrSamples()):
      lower = 0.0
      uniformReal = random.uniform(0.0, 1.0)
      for j in range(0, len(cumulativeWeights)):
        upper = cumulativeWeights[j]
        if uniformReal > lower and uniformReal <= upper:
          samplesToClone[j] += 1
        lower = upper

    self._writeFilterStatistics(normalisedWeights, cumulativeWeights, samplesToClone)
    return samplesToClone


  def optimalSampleNumber(self, filterTimestep):
    filename = "filter%s.csv" % (filterTimestep)

    handle = open(filename, 'r')
    col = numpy.zeros(self._userModel().nrSamples())
    row = 0
    for line in handle.xreadlines():
      col1, col2, col3, col4 = line.split(';')
      if row > 0:
        col[row - 1] = float(col2)
      row += 1

    return ((numpy.nonzero(col == max(col))[0])[0]) + 1



  def _writeFilterStatistics(self, normalisedWeights, cumulativeWeights, samplesToClone):
    filename = "filterSIR_%s.csv" % (self._userModel().currentTimeStep())
    csvFile = csv.writer(open(filename, "w"), delimiter=";",quoting=csv.QUOTE_NONNUMERIC)

    csvFile.writerow(["sample", "normalised weight", "cumulative weight", "resampled particles"])
    for i in range(0, self._userModel().nrSamples()):
      csvFile.writerow([(i+1), normalisedWeights[i], cumulativeWeights[i], samplesToClone[i]])








## \brief Residual resampling algorithm
class ResidualResamplingFramework(ParticleFilterFramework):
  ## \brief Constructor
  def __init__(self, userModel):
    ParticleFilterFramework.__init__(self, userModel)


  def _samplesToClone(self, sampleWeights):
    # normalise weights
    weights = self._normaliseWeights(sampleWeights)

    nrSamples = self._userModel().nrSamples()
    samplesToClone = [int(0)] * nrSamples
    cdfResidualWeights = numpy.zeros(nrSamples)

    # resampling factor
    for i in range(0, self._userModel().nrSamples()):
      samplesToClone[i] = int(math.floor(nrSamples * weights[i]))
    # bookkeeping for reporting statistics
    resamplingFactor = samplesToClone[:]

    # obtain residual number of samples to clone
    nrNewSamples = sum(samplesToClone)
    if nrNewSamples < nrSamples:
      nrMissingSamples = int(nrSamples - nrNewSamples)
      residualWeights = [0.0] * nrSamples

      for i in range(0, len(residualWeights)):
        residualWeights[i] = (nrSamples * weights[i] - math.floor(nrSamples * weights[i])) / nrMissingSamples

      #\todo replace this
      normalisedResidualWeights = [0.0] * nrSamples
      sumWeights = sum(residualWeights)
      for i in range(0, len(residualWeights)):
        normalisedResidualWeights[i] = residualWeights[i] / sumWeights

      assert normalisedResidualWeights
      cdfResidualWeights = self._cumulativeWeights(normalisedResidualWeights)

      for i in range(0, nrMissingSamples):
        lower = 0.0
        uniformReal = random.uniform(0.0, 1.0)
        for j in range(0, len(cdfResidualWeights)):
          upper = cdfResidualWeights[j]
          if uniformReal > lower and uniformReal <= upper:
            samplesToClone[j] += 1
          lower = upper

    assert sum(samplesToClone) == nrSamples
    self._writeFilterStatistics(weights, resamplingFactor, cdfResidualWeights, samplesToClone)
    return samplesToClone


  def _writeFilterStatistics(self, normalisedWeights, resamplingFactor, cdfResidualWeights, samplesToClone):
    filename = "filterRR_%s.csv" % (self._userModel().currentTimeStep())
    csvFile = csv.writer(open(filename, "w"), delimiter=";",quoting=csv.QUOTE_NONNUMERIC)

    csvFile.writerow(["sample", "normalised weight", "resampling factor", "cdf residual weights", "resampled particles"])
    for i in range(0, self._userModel().nrSamples()):
      csvFile.writerow([(i+1), normalisedWeights[i], resamplingFactor[i], cdfResidualWeights[i], samplesToClone[i]])

