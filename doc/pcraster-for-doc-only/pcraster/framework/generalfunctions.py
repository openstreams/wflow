#!/usr/bin/env python
# -*- coding: utf-8 -*-


from pcraster import *

# time in hours

def getCellValue(Map, Row, Column):
  Value, Valid=cellvalue(Map, Row, Column)
  if Valid:
    return Value
  else:
    print 'missing value in input of getCellValue'

def getCellValueAtBooleanLocation(location,map):
  # map can be any type, return value always float
  valueMap=mapmaximum(ifthen(location,scalar(map)))
  value=getCellValue(valueMap,1,1)
  return value

def printCellValue(self, mapVariable, variableNameToPrint, unit, row, column):
  cellValue=getCellValue(mapVariable,row,column)
  print variableNameToPrint + ' (' + unit + ') at row ' + str(row) + ', column: ' + str(column) + ' is: ' + str(cellValue)

def onePeriod(self, startTime, endTime, timeStepDuration, currentTimeStep):
  # this could be separated in two functions, one converting hours to
  # time steps, one creating the period
  time = float(currentTimeStep) * float(timeStepDuration)
  period = (time > startTime) & (time < endTime)
  return period

def mapeq(mapOne, mapTwo):
  mapOneScalar=scalar(mapOne)
  mapTwoScalar=scalar(mapTwo)
  difference=mapOneScalar-mapTwoScalar
  cellEqual=pcreq(difference,scalar(0))
  mapEqual=pcrgt(mapminimum(scalar(cellEqual)),scalar(0.5))
  return getCellValue(mapEqual,1,1)

def slopeToDownstreamNeighbour(dem, ldd):
  slopeToDownstreamNeighbour=(dem-downstream(ldd,dem))/downstreamdist(ldd)
  return slopeToDownstreamNeighbour

def slopeToDownstreamNeighbourNotFlat(dem,ldd,minSlope):
  slopeToDownstreamNeighbourMap=slopeToDownstreamNeighbour(dem,ldd)
  lddArea=defined(ldd)
  minSlopeCover=ifthen(lddArea,scalar(minSlope))
  slopeToDownstreamNeighbourNotFlat=cover(max(minSlopeCover,slopeToDownstreamNeighbourMap), minSlopeCover)
  return slopeToDownstreamNeighbourNotFlat

def distancetodownstreamcell(Ldd):
  distanceToDownstreamCell=max(downstreamdist(Ldd),celllength())
  return distanceToDownstreamCell

def createTimeSeriesList(timeSeriesFile):
  file = open(timeSeriesFile,'r')
  piet = file.readlines()
  newList=[]
  for line in piet:
    lineList=string.split(line)
    newList.append(lineList)
  file.close()
  return newList

def timeInputSparse(fileName):
  return os.path.exists(fileName)

def normalcorrelated(normalX, normalY, correlation):
  # returns realizations of two normal variables with
  # mean zero and var 1 having correlation of correlation
  # based on:
  # x=normal()
  # y=ax+b*normal()
  # correlation = a /  sqrt( sqr(a) + sqr(b) )
  x=scalar(normalX)
  y=(x+sqrt((1/sqr(correlation))-1)*scalar(normalY)) * scalar(correlation)
  return x,y
