#!/usr/bin/env python
# -*- coding: utf-8 -*-


from numpy import *
import numpy as np
import sys



#----------------------------------------------------------------------#
#  function:   linearRegression                                        #
#----------------------------------------------------------------------#
def linearRegression(rawData, equationOrder):

   # Purpose: Regress the coefficients for an equation of the form:
   #     Y = C0 + C1*X + C2*X^2 + ... + Cn*X^n
   # The regression outputs a list of the coefficients ([C0,C1,..,C[n])

   # rawData: A supplied list of Y and X values which serve as a basis
   # for computing the coefficients.  The general form of the list is:
   # [[Y0, X0],[Y1, X1],...,[Ym,Xm]].

   # equationOrder: The equation order is the power function for the X
   # variable where the number of coefficients returned is the equation + 1.
   # For example, an equation of 1 maps to Y = C0 + C1X.  An equation order
   # of 2 maps to Y = C0 + C1*X + C2*X*X.

   # require part of contract
   assert(type(rawData) == types.ListType), \
      "Raw data input must be a list"
   assert(equationOrder > 0), \
      "Equation order must be greater than 0th order"
   assert(len(rawData) >= equationOrder), \
      "Number of data points must be >= to number of coefficients be calculated"
   for each in (rawData):
      assert(type(each) == types.ListType), \
         "Raw data input must be a list of data values"
      assert(len(each) > 1), \
         "More than one data point is required for the raw data"
      assert(len(each) == len(rawData[0])), \
         "All data points in raw data must have the same number items"

   xEquationForm = [lambda rawItem, coefIndex: pow(rawItem[1], coefIndex)]
   return regression(rawData, xEquationForm * (equationOrder+1))

#----------------------------------------------------------------------#
#  function:   regression                                              #
#----------------------------------------------------------------------#
def regression(rawData, xEquationForm, yEquationForm = lambda rawItem: rawItem[0]):

   # Purpose: Regress the coefficients for an equation of a generalized form
   # as supplied by the xEquationForm and yEquationForm functions:
   #     Y = C0*X0 + C1*X1 + C2*X2 + ... + Cn*Xn
   # Where Y = yEquationForm() and Xn = xEquationForm[n]()
   # The regression outputs a list of the coefficients ([C0,C1,..,C[n])

   # rawData: A supplied list of Y and X values which serve as a basis
   # for computing the coefficients.  The general form of the list is:
   # [[Y0, X0],[Y1, X1],...,[Ym,Xm]].

   # equationOrder: The equation order is the power function for the X
   # variable where the number of coefficients returned is the equation + 1.
   # For example, an equation of 1 maps to Y = C0 + C1X.  An equation order
   # of 2 maps to Y = C0 + C1*X + C2*X*X.

   # require part of contract
   assert(type(rawData) == types.ListType), \
      "Raw data input must be a list"
   assert(type(xEquationForm) == types.ListType), \
      "X Equation form must be defined in a list"
   assert(type(yEquationForm == types.FunctionType)), \
      "Y Equation form must be a lambda function"
   assert(len(xEquationForm) > 0), \
      "X Equation form must not be an empty list"
   assert(len(rawData) >= len(xEquationForm)), \
      "Number of data points must be >= to number of coefficients be calculated"
   for each in (xEquationForm):
      assert(type(each) == types.FunctionType), \
         "X Equation form must be lambda functions"
   for each in (rawData):
      assert(type(each) == types.ListType), \
         "Raw data input must be a list of data values"
      assert(len(each) > 1), \
         "More than one data point is required for the raw data"
      assert(len(each) == len(rawData[0])), \
         "All data points in raw data must have the same number items"

   # number of coefficients being solved for
   numCoefficients = len(xEquationForm)

   # the value of each term for the equation
   term = [0] * numCoefficients

   # the matrices for the simultaneous equations
   B = [0] * numCoefficients
   A = []
   for i in range(numCoefficients):
      A.append([0] * numCoefficients)

   # loop through all the raw data samples
   for each in rawData:

      # sum the y values
      yCurrent = float(yEquationForm(each))
      B[0] = B[0] + yCurrent

      # sum the x values
      for i in range(numCoefficients):
         term[i] = float(xEquationForm[i](each, i))
         A[0][i] = A[0][i] + term[i]

      # now set up the rest of rows in the matrices
      # (multiplying each row by each term)
      for i in range(1, numCoefficients):
         B[i] = B[i] + yCurrent * term[i]
         for j in range(numCoefficients):
            A[i][j] = A[i][j] + term[i] * term[j]

   # solve for the coefficients
   return gauss(A, B)

#----------------------------------------------------------------------#
#  function:   linearRSquared                                          #
#----------------------------------------------------------------------#
def linearRSquared(rawData, coefficients):

   # Purpose: Compute the R-Squared statistic for the supplied coefficients

   # require part of contract
   assert(type(rawData) == types.ListType), \
      "Raw data input must be a list"
   assert(type(coefficients) == types.ListType), \
      "Computed coefficients must be a list"
   assert(len(coefficients) > 0), \
      "At least coefficient is required"
   assert(len(rawData) >= len(coefficients)), \
      "Number of data points must be >= to number of coefficients"
   for each in (rawData):
      assert(type(each) == types.ListType), \
         "Raw data input must be a list of data values"
      assert(len(each) > 1), \
         "More than one data point is required for the raw data"
      assert(len(each) == len(rawData[0])), \
         "All data points in raw data must have the same number items"

   xEquationForm = [lambda rawItem, coefIndex: pow(rawItem[1], coefIndex)]
   return solveRSquared(rawData, coefficients, xEquationForm * len(coefficients))

#----------------------------------------------------------------------#
#  function: solveRSquared                                             #
#----------------------------------------------------------------------#
def solveRSquared(rawData, coefficients, xEquationForm, \
   yEquationForm = lambda rawItem: rawItem[0]):

   # Purpose: Compute the R-Squared statistic for the supplied coefficients

   # require part of contract
   assert(type(rawData) == types.ListType), \
      "Raw data input must be a list"
   assert(type(xEquationForm) == types.ListType), \
      "X Equation form must be defined in a list"
   assert(type(yEquationForm == types.FunctionType)), \
      "Y Equation form must be a lambda function"
   assert(len(xEquationForm) > 0), \
      "X Equation form must not be an empty list"
   assert(len(rawData) >= len(xEquationForm)), \
      "Number of data points must be >= to number of coefficients be calculated"
   for each in (xEquationForm):
      assert(type(each) == types.FunctionType), \
         "X Equation form must be lambda functions"
   for each in (rawData):
      assert(len(each) == len(rawData[0])), \
         "All data points in raw data must have the same number items"
      assert(len(each) > 1), \
         "More than one data point is required for the raw data"

   # number of coefficients being solved for
   numCoefficients = len(xEquationForm)

   # sum of y*y
   ysquare = 0

   # number of raw data samples
   samples = len(rawData)

   # the value of each term for the equation
   term = [0] * numCoefficients

   # the matrices for the simultaneous equations
   B = [0] * numCoefficients
   A = []
   for i in range(numCoefficients):
      A.append([0] * numCoefficients)

   # use the first column as the default value for the dependent variable
   if yEquationForm == []:
      yEquationForm = lambda rawItem: rawItem[0]

   # loop through all the raw data samples
   for each in rawData:

      # sum the y values
      yCurrent = float(yEquationForm(each))
      B[0] = B[0] + yCurrent
      ysquare = ysquare + yCurrent * yCurrent

      # sum the x values
      for i in range(numCoefficients):
         term[i] = float(xEquationForm[i](each, i))
         A[0][i] = A[0][i] + term[i]

      # now set up the rest of rows in the matrices
      # (multiplying each row by each term)
      for i in range(1, numCoefficients):
         B[i] = B[i] + float(yCurrent * term[i])
         for j in range(numCoefficients):
            A[i][j] = A[i][j] + term[i] * term[j]

   # calculate the r-squared statistic
   sumsquare = 0
   yaverage = B[0] / samples
   for i in range(numCoefficients):
      xaverage = A[0][i] / samples
      sumsquare = sumsquare + coefficients[i] * (B[i] - (samples * xaverage * yaverage))
   return sumsquare / (ysquare - (samples * yaverage * yaverage))

#----------------------------------------------------------------------#
#  function:   gauss                                                   #
#----------------------------------------------------------------------#
def gauss(AMatrix, BMatrix):

   # solve linear equations of the form:
   #     |A00 A01 ... A0n|   |coefficient0|   |B0|
   #     |A10 A11 ... A1n| * |coefficient1| = |B1|
   #     |... ... ... ...|   |     ...    |   |..|
   #     |An0 An1 ... Ann|   |coefficientn|   |Bn|
   # where |A| and |B| are supplied and |coefficient| is the solution.

   # require part of contract
   assert(type(AMatrix) == types.ListType), \
      "Input matrix must be a list"
   assert(type(BMatrix) == types.ListType), \
      "Input matrix must be a list"
   assert(len(AMatrix) > 0), \
      "Input matrix must not be an empty list"
   assert(len(AMatrix) == len(BMatrix)), \
      "A and B matrix must have same number of rows"
   for each in (AMatrix):
      assert(len(each) == len(AMatrix)), \
         "A matrix must be of square (NxN) dimensions"

   # get the length of the matrix
   n = len(AMatrix)

   # copy the matrices to local variables - inplace substitution used
   A = map(lambda x: list(x), AMatrix)
   B = list(BMatrix)

   # initialize the output results
   coefficients = [0] * n

   # condition the matrix for the solution
   for i in range(n-1):
      pivot = A[i][i]
      assert(pivot != 0), \
         "Zero pivot value encountered"
      for j in range(i+1, n):
         multiplier = float(A[j][i]) / pivot
         for k in range(i+1, n):
            A[j][k] = A[j][k] - multiplier * A[i][k]
         B[j] = B[j] - multiplier * B[i]

   # solve for the coefficients
   assert(A[n-1][n-1] != 0), \
      "Divide by zero encountered in solution"
   coefficients[n-1] = float(B[n-1]) / A[n-1][n-1]
   for i in range(n-2, -1, -1):
      top = B[i]
      for j in range(i+1, n):
         top = top - (A[i][j]* coefficients[j])
      assert(A[i][i] != 0), \
         "Divide by zero encountered in solution"
      coefficients[i] = float(top) / A[i][i]

   return coefficients
