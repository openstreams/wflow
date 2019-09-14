# -*- coding: utf-8 -*-
import math
import os
import shutil
import string
import numpy
import numpy.ma
from pcraster import *
from .frameworkBase import generateNameS, generateNameT, generateNameST
from . import generalfunctions, regression



def _percentile(
  array,
  level):
  k = (len(array)-1) * level
  f = math.floor(k)
  c = math.ceil(k)
  if f == c:
    return array[int(k)]
  else:
    d0 = array[int(f)] * (c - k)
    d1 = array[int(c)] * (k - f)
    return d0 + d1



def selectSArray(name, sampleNumbers, row, col):
  """Selects values at row, col from raster name in Monte Carlo samples.

  name -- Name of raster.
  sampleNumber -- Numbers of MC samples to use.
  row -- Row index of cell to read.
  col -- Col index of cell to read.
  The returned array does not contain missing values so the size is maximimal
  sampleNumbers but possibly smaller.

  Returned array has elements of type numpy.float32"""
  mask = numpy.zeros(len(sampleNumbers)).astype(numpy.bool_)
  array = numpy.zeros(len(sampleNumbers)).astype(numpy.float32)
  i = 0
  while i < len(sampleNumbers):
    filename = generateNameS(name, sampleNumbers[i])
    array[i], mask[i] = readFieldCell(filename, row, col)
    i += 1
  array = numpy.compress(numpy.logical_not(mask), array)
  return array



def selectSArrays(
  name,
  sampleNumbers):
  masks = []
  arrays = []
  nrCells = clone().nrRows() * clone().nrCols()

  # For each cell.
  # Create empty arrays for cell values and 'missing valueness'.
  c = 0
  while c < nrCells:
    masks.append(numpy.zeros(len(sampleNumbers)).astype(numpy.bool_))
    arrays.append(numpy.zeros(len(sampleNumbers)).astype(numpy.float32))
    c += 1

  # For each sample.
  # Read raster and assign cell values to arrays.
  s = 0
  while s < len(sampleNumbers):
    raster = readmap(generateNameS(name, sampleNumbers[s]))

    # For each cell.
    c = 0
    while c < nrCells:
      arrays[c][s], masks[c][s] = cellvalue(raster, c + 1)
      c += 1
    s += 1

  # For each cell.
  # Compress each array with cell values.
  c = 0
  while c < nrCells:
    arrays[c] = numpy.compress(masks[c], arrays[c])
    c += 1

  return arrays



# Helper function for aggregation of simulation results. This function
# creates an array for each cell location in the rasters for the variables
# in the names variable. With this array the function is called which can
# calculate statistics or whatever.
def aggregateSPerCell(name, sampleNumbers, calculator):
  for row in range(clone().nrRows()):
    for col in range(clone().nrCols()):
      calculator.run(row, col, selectSArray(name, sampleNumbers, row, col))
  return calculator.result()



def aggregateS(
  name,
  sampleNumbers,
  calculator):
  calculator.run(selectSArrays(name, sampleNumbers))
  return calculator.result()



# Helper class for aggeration of simulation results.
class PercentileCalculator:

  def __init__(self, percentiles):
    if isinstance(percentiles, int):
      self.d_percentiles = [percentiles]
    else:
      self.d_percentiles = percentiles
    assert isinstance(self.d_percentiles, list)
    self.d_fields = []
    for percentile in self.d_percentiles:
      self.d_fields.append(newScalarField())

  def run(self,
    arrays):
    for c in range(clone().nrRows() * clone().nrCols()):
      if len(arrays[c]) > 0:
        arrays[c].sort()
        for p in range(len(self.d_percentiles)):
          value = _percentile(arrays[c], self.d_percentiles[p])
          self.d_fields[p].setCell(numpy.float64(value), c)

  def result(self):
    assert len(self.d_fields) >= 1
    if len(self.d_fields) == 1:
      return self.d_fields[0]
    else:
      return self.d_fields



def probability(name, sampleNumbers):
  """
  Calculates the probability that a cell is TRUE.

  name
    Name of the (boolean) raster for which each sample has a realization.

  sampleNumbers 
    List of numbers of samples to aggregate.

  Returns a raster with probabilities.
  """
  present = scalar(0)
  count = scalar(0)
  for sample in sampleNumbers:
    filename = generateNameS(name, sample)
    raster   = readmap(filename)
    present  = ifthenelse(raster, present + 1, present)
    count    = ifthen(defined(raster), count + 1)
  return present / count



def average(name, sampleNumbers):
  """
  Calculates the average value of each cell.

  name
    Name of the scalar raster for which each sample has a realization.

  sampleNumbers
    List of numbers of samples to aggregate.

  Returns a raster with average values.
  """
  sum = scalar(0)
  count = scalar(0)
  for sample in sampleNumbers:
    filename = generateNameS(name, sample)
    raster   = readmap(filename)
    sum      = sum + raster
    count    = ifthen(defined(raster), count + 1)
  return sum / count



def variance(name, sampleNumbers):
  """
  Calculates the variance of each cell.

  name
    Name of the scalar raster for which each sample has a realization.

  sampleNumbers
    List of numbers of samples to aggregate.

  Returns a raster with variances.
  """
  sumOfSquaredValues, sumOfValues, count = scalar(0), scalar(0), scalar(0)
  for sample in sampleNumbers:
    filename           = generateNameS(name, sample)
    raster             = readmap(filename)
    sumOfSquaredValues = sumOfSquaredValues + raster ** 2
    sumOfValues        = sumOfValues + raster
    count              = ifthen(defined(raster), count + 1)
  return (count * sumOfSquaredValues - sumOfValues ** 2) / (count * (count - 1))



def stddev(name, sampleNumbers):
  """
  Calculates the standard deviation of each cell.

  name
    Name of the scalar raster for which each sample has a realization.

  sampleNumbers
    List of numbers of samples to aggregate.

  Returns a raster with standard deviations.
  """
  return sqrt(variance(name, sampleNumbers))



def percentile(
  name,
  sampleNumbers,
  percentiles):
  """
  Calculates a percentile for each cell.

  name
    Name of the scalar raster for which each sample has a realization.

  sampleNumbers
    List of numbers of samples to aggregate.

  percentiles
    Percentile or list of percentiles to calculate. Percentiles range from
    [0.0, 1.0].

  Returns a raster or a list of rasters with percentiles.
  """
  return aggregateS(name, sampleNumbers, PercentileCalculator(percentiles))



# row >= 1, <= nrRows, col >= 1, <= nrCols
def timeseries(name, timeSteps, row, col):
  mask = numpy.zeros(len(timeSteps)).astype(numpy.bool_)
  steps = numpy.zeros(len(timeSteps)).astype(numpy.int32)
  array = numpy.zeros(len(timeSteps)).astype(numpy.float32)
  i = 0
  while i < len(timeSteps):
    filename = generateNameT(name, timeSteps[i])
    steps[i] = timeSteps[i]
    array[i], mask[i] = readFieldCell(filename, row, col)
    i += 1
  steps = numpy.compress(numpy.logical_not(mask), steps)
  array = numpy.compress(numpy.logical_not(mask), array)
  return steps, array



#def varmean(names,sampleNumbers, timeSteps):
#  nrSamples=scalar(len(sampleNumbers))
#  for name in names:
#    for step in timeSteps:
#      sumSquared=scalar(0.0)
#      sum=scalar(0.0)
#      for sample in sampleNumbers:
#        realization=scalar(generateNameST(name,sample,step))
#        sumSquared=sumSquared+realization*realization
#        sum=sum+realization
#      std=(scalar(1.0)/nrSamples)*sqrt(nrSamples*sumSquared-sum*sum)
#      var=std*std
#      mean=sum/nrSamples
#      report(var, generateNameT(name + '-var', step))
#      report(mean, generateNameT(name + '-ave', step))



def correlation(location, independentName, dependentName, locationName, sampleNumbers, timeSteps):
  location=boolean(location)
  name=independentName + '_' + dependentName + '_' + locationName
  tssFileIntercept = file("%s%s.tss" % (name,'_int'), "w")
  tssFileSlope = file("%s%s.tss" % (name,'_slope'), "w")
  tssFileRSquared = file("%s%s.tss" % (name,'_rSq'), "w")
  for step in timeSteps:
    values=[]
    for sample in sampleNumbers:
      smallValue=0.0000000000000000001
      fileNameOne=generateNameST(dependentName,sample,step)
      valueOne=generalfunctions.getCellValueAtBooleanLocation(location,scalar(fileNameOne))
      pairList=[valueOne + smallValue]
      fileNameTwo=generateNameST(independentName,sample,step)
      valueTwo=generalfunctions.getCellValueAtBooleanLocation(location,scalar(fileNameTwo))
      pairList.append(valueTwo + smallValue)
      values.append(pairList)
      #print valueOne + smallValue, valueTwo + smallValue
    reg=regression.linearRegression(values, 1)
    rSq=regression.linearRSquared(values,reg)
    #print step, reg, rSq
    tssFileIntercept.write("%d %g\n" % (step, reg[0]))
    tssFileSlope.write("%d %g\n" % (step, reg[1]))
    tssFileRSquared.write("%d %g\n" % (step, rSq))
  tssFileIntercept.close()
  tssFileSlope.close()
  tssFileRSquared.close()



def staticInput(timeSteps):
  return len(timeSteps) == 1 and timeSteps[0] == 0



def deterministicInput(sampleNumbers):
  return len(sampleNumbers) == 1 and sampleNumbers[0] == 0



def sampleMin(name, sampleNumbers):
  """
  Calculates the minimum value of each cell.

  name
    Name of the scalar raster for which each sample has a realization.

  sampleNumbers
    List of numbers of samples to aggregate.

  Returns a raster with minimum values.
  """
  minimum = scalar(1e31)
  for sample in sampleNumbers:
    filename = generateNameS(name, sample)
    raster   = scalar(readmap(filename))
    minimum      = ifthenelse(pcrlt(raster,minimum),raster,minimum)
  return minimum



def sampleMax(name, sampleNumbers):
  """
  Calculates the maximum value of each cell.

  name
    Name of the scalar raster for which each sample has a realization.

  sampleNumbers
    List of numbers of samples to aggregate.

  Returns a raster with maximum values.
  """
  maximum = scalar(-1e31)
  for sample in sampleNumbers:
    filename = generateNameS(name, sample)
    raster   = scalar(readmap(filename))
    maximum      = ifthenelse(pcrgt(raster,maximum),raster,maximum)
  return maximum



def uniquesamples(name, sampleNumbers):
  """
  Retrieves the unique samples.

  name
    Name of the raster for which each sample has a realization.

  sampleNumbers
    List of numbers of samples to aggregate.

  Returns a list with sets of corresponding loops.
  """
  uniqueSets=[]
  for sample in sampleNumbers:
    filename = generateNameS(name, sample)
    raster   = readmap(filename)
    setNumber=0
    sampleAddedToExistingSet=False
    for uniqueSet in uniqueSets:
      if generalfunctions.mapeq(uniqueSet[0],raster):
        uniqueSets[setNumber].append(sample)
        sampleAddedToExistingSet=True
        break
      setNumber += 1
    if sampleAddedToExistingSet==False:
      uniqueSets.append([raster,sample])
  firstLoopOfEachUniqueSet=[]
  for uniqueSet in uniqueSets:
    firstLoopOfEachUniqueSet.append(uniqueSet[1])
  return len(firstLoopOfEachUniqueSet),firstLoopOfEachUniqueSet



def mcaveragevariance(names,sampleNumbers, timeSteps):
  if staticInput(timeSteps):
    for name in names:
      mean=average(name + '.map', sampleNumbers)
      var=variance(name + '.map', sampleNumbers)
      minimum=sampleMin(name + '.map', sampleNumbers)
      maximum=sampleMax(name + '.map', sampleNumbers)
      #std=stddev(name + '.map', sampleNumbers)
      report(mean, name + '-ave.map')
      report(var, name + '-var.map')
      report(minimum, name + '-min.map')
      report(maximum, name + '-max.map')
      report(sqrt(var)/mean, name + '-err.map')
  else:
    nrSamples=scalar(len(sampleNumbers))
    for name in names:
      for step in timeSteps:
        var=variance(generateNameT(name,step),sampleNumbers)
        mean=average(generateNameT(name,step), sampleNumbers)
        report(mean, generateNameT(name + '-ave', step))
        report(var, generateNameT(name + '-var', step))
        report(sqrt(var)/mean, generateNameT(name + '-err', step))



def mcpercentiles(
  names,
  percentiles,
  sampleNumbers,
  timeSteps):
  if staticInput(timeSteps):
    for name in names:
      results = percentile(name + ".map", sampleNumbers, percentiles)
      for i in range(len(percentiles)):
        report(results[i], "%s_%s.map" % (name, percentiles[i]))
  else:
    for name in names:
      for step in timeSteps:
        results = percentile(generateNameT(name, step), sampleNumbers,
          percentiles)
        assert len(results) == len(percentiles)
        for i in range(len(percentiles)):
          report(results[i], "%s_%d_%s.map" % (name, step, percentiles[i]))



def createtimeseries(names, nameExtension, locations,sampleNumbers,timeSteps):
  if deterministicInput(sampleNumbers):
    for name in names:
      tssFile = file(name + nameExtension + '.tss', "w")
      tssFile.write("timeseries scalar\n")
      tssFile.write("2\n")
      tssFile.write("timestep\n")
      tssFile.write("%s\n" % (name))
      for step in timeSteps:
        timeseriesValue=mapmaximum(ifthen(locations,generateNameT(name,step)))
        value, valid = cellvalue(timeseriesValue, 1, 1); assert valid
        tssFile.write("%d %g\n" % (step, value))
      tssFile.close()
  else:
    for name in names:
      for sample in sampleNumbers:
        tssFile = file(generateNameS("%s%s.tss" % (name,nameExtension), sample), "w")
        tssFile.write("timeseries scalar\n")
        tssFile.write("2\n")
        tssFile.write("timestep\n")
        tssFile.write("%s\n" % (name))
        for step in timeSteps:
          timeseriesValue=mapmaximum(ifthen(locations,generateNameST(name,sample,step)))
          value, valid = cellvalue(timeseriesValue, 1, 1); assert valid
          tssFile.write("%d %g\n" % (step, value))
        tssFile.close()



def createtimeseriesnewfileformat(names,locations,sampleNumbers,timeSteps, quantiles):
  if deterministicInput(sampleNumbers):
    for quantile in quantiles:
      for name in names:
        tssFile = file(name + "_" + str(quantile) + '.tss', "w")
        tssFile.write("timeseries scalar\n")
        tssFile.write("2\n")
        tssFile.write("timestep\n")
        tssFile.write("%s\n" % (name))
        for step in timeSteps:
          filename = name + ("_%d_" % (step)) + str(quantile) + ".map"
          timeseriesValue=mapmaximum(ifthen(locations,filename))
          value, valid = cellvalue(timeseriesValue, 1, 1); assert valid
          tssFile.write("%d %g\n" % (step, value))
        tssFile.close()
  else:
    print('timeseries for monte carlo loops not yet available')



def createGstatRealizations(setOfRealizations, nameCommandFile, nameOutMapList):
  # number of realizations required
  nSim=len(setOfRealizations)
  # open template gstat script and replace
  #print nameCommandFile
  gstatTemplate=file(nameCommandFile + '.gst','r')
  gstatTemplateString=gstatTemplate.read()
  gstatTemplate.close()
  gstatString=gstatTemplateString.replace('NSIM',str(nSim))
  gstatFile=file('tmpGstat.gst','w')
  gstatFile.write(gstatString)
  gstatFile.close()
  # run gstat
  os.system('gstat tmpGstat.gst')
  os.remove('tmpGstat.gst')
  # rename files
  i=1
  for realization in setOfRealizations:
    #print realization
    item=0
    for name in nameOutMapList:
      gstatOutputFileName=generateNameT('g_' + name,i)
      #print gstatOutputFileName, realization[item]
      shutil.move(gstatOutputFileName,realization[item])
      item=item+1
    i = i + 1



def createAllGstatRealizations(nameCommandFile,nameOutMapList,nrRealPerGstatCall,sampleNumbers,timeSteps):
  # create list names with filenames required
  names=[]
  names.append([])
  i = 0
  j = 0
  for sample in sampleNumbers:
    for step in timeSteps:
      if j == nrRealPerGstatCall:
        j = 1
        i = i+1
        names.append([])
      else:
        j = j+1
      namesOneSampleOneTimeStep=[]
      for nameOutMap in nameOutMapList:
        if staticInput(timeSteps):
          fileName=generateNameS(nameOutMap,sample) + '.map'
        else:
          fileName=generateNameST(nameOutMap,sample,step)
        namesOneSampleOneTimeStep.append(fileName)
      names[i].append(namesOneSampleOneTimeStep)
  for setOfRealizations in names:
    createGstatRealizations(setOfRealizations,nameCommandFile, nameOutMapList)
    #print setOfRealizations

