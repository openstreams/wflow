

.. index::
   single: timeoutput
.. _timeoutput:

**********
timeoutput
**********
.. topic:: timeoutput

   Expression value of an uniquely identified cell or cells written to a time series per timestep

::

  ResultTimeSeries = timeoutput(idexpression, expression)

idexpression
   spatial
   boolean, nominal, ordinal

expression
   spatial
   boolean, nominal, ordinal, scalar, directional, ldd

ResultTimeSeries
   
   time series

Operation
=========


This operation is used in the iterative sections (dynamic, storage and transport sections) of
a dynamic model script only. The keyword report precedes the operation. The
idexpression is an expression that contains one identified cell or several uniquely identified cells. For each timestep the value of expression at the identified cell or cells is written to ResultTimeSeries, which is an ascii formatted time series file.  After a model run, the time series contains for each identified cell a list of expression cell values per timestep. 



The idexpression must contain one or more uniquely identified cells, which are numbered with consecutive whole values, starting with 1. The remaining cells must have a value 0. For instance, if you want to write the expression values from three different cells to ResultTimeSeries, these cells must have the values 1, 2 and 3 on idexpression respectively, the remaining cells must have a value 0.  



The ResultTimeSeries is an ascii formatted time series with header. It has the following lay out:  



line 1: header, description: expression map name  



line 2: header, number of columns in the file 




line 3: header, time column description: model time 




line 4 up to and including line N+ 3: header, the numbers of the N uniquely
identified cells: 1,2,3,...N.





subsequent lines: data formatted in rows and columns. Each row represents one timestep
I at time t(I) in the model from which the time series is the result; the first
row contains data for timestep I = 1, the second row for timestep I =2, etc. The
first column contains the time t at the timesteps. At the first row which contains data
for the first time step (I = 1) it is always the starttime t(1).

The remaining columns (column number 2 up to and including :emphasis:`n` + 1, where :emphasis:`n` is the number of uniquely identified cells as above said) contain values that are taken from expression: column number :emphasis:`n` + 1 contains data of the cell that has a value :emphasis:`n`.   

Notes
=====


In principle, each unique identifier is represented by one cell on idexpression. This is the basic use: to sample certain locations. Alternatively, more than one cell on idexpression may have the same unique identifier value. In that case ResultTimeSeries contains for the id under consideration an aggregate expression value of the set of cells with that id. If expression is of data type scalar or directional the average expression value of cells with the id under consideration on idexpression is written to ResultTimeSeries. If the data type is directional cells without a direction (cellvalue -1) are discarded in the calculation of the average value. If all cells with the id under consideration have no direction, a -1 value is written to ResultTimeSeries. If expression is of data type boolean, nominal or ordinal the highest score (most ocurring cell value) of the cells with the id under consideration is written to the time series. If several values all have the highest score, the highest value is assigned.    

Group
=====
This operation belongs to the group of  Time operators 

Examples
========
