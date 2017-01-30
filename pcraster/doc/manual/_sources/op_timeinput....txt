

.. index::
   single: timeinput...
.. index::
   single: timeinputboolean
.. index::
   single: timeinputnominal
.. index::
   single: timeinputordinal
.. index::
   single: timeinputscalar
.. index::
   single: timeinputdirectional
.. index::
   single: timeinputldd
.. _timeinput...:

************
timeinput...
************
.. topic:: timeinput...

   Cell values per timestep read from a time series that is linked to a map with unique identifiers

::

  Result = timeinput...(type, TimeSeries, idexpression)

::

  Result = timeinputnominal(TimeSeries, idexpression)

::

  Result = timeinputboolean(TimeSeries, idexpression)

::

  Result = timeinputnominal(TimeSeries, idexpression)

::

  Result = timeinputordinal(TimeSeries, idexpression)

::

  Result = timeinputscalar(TimeSeries, idexpression)

::

  Result = timeinputdirectional(TimeSeries, idexpression)

::

  Result = timeinputldd(TimeSeries, idexpression)

TimeSeries
   
   ascii formatted time series

idexpression
   spatial, non spatial
   boolean, nominal, ordinal

Result
   spatial
   type is specified by the sort of command: timeinputboolean results in a boolean Result, timeinputnominal in a nominal Result etc.

Operation
=========
Example of a time series file with header
-----------------------------------------

  | rain (mm) per area for model with timer 1 6 1
  | 4
  | model time
  | Area station A
  | Area station B
  | Area station C
  | 1    2.0  0.0  0.0    
  | 2    3.0  2.0  1.0
  | 3    7.0  5.0  3.0
  | 4    9.0  12.0 6.0
  | 5    6.0  0.0  5.0
  | 6    0.0  1.0  2.0
  | 

This operation is used in the iterative sections (dynamic, storage and transport sections) of a dynamic model script only. TimeSeries is an ascii formatted time series that contains cell values formatted in rows and columns. During a run of a Dynamic Model TimeSeries is read from top to bottom: for each timestep, a row gives cell values that are assigned to Result at the timestep under consideration.  The cell values in a row are assigned to the cells of Result on basis of unique identifiers on idexpression: each column on TimeSeries gives cell values for an unique identifier. Each timestep, the cell value in a column is assigned to the cells on Result that have a unique identifier on idexpression that corresponds with the unique identifier of the column.  



The data type that is assigned to Result is specified by the sort of operator that is used.  



The contents and partly also the format (number of rows) of the
TimeSeries must match the dynamic model for which the TimeSeries is used, especially the time dimension of the model. For a description of the time dimension and the terms used, see section VI.2.2d. Two types of format for the TimeSeries can be used; the timeinput... operator detects the formats by itself:  



1) a time series file with a header 




line 1: header, description 




line 2: header, number of columns in the file 




line 3: header, time column description 




line 4 up to and including line :emphasis:`n` + 3: header, the names of the :emphasis:`n` unique identifiers.  




subsequent lines: data formatted in rows and columns. Each row represents
one timestep I at time t(I) in the model for which the time
series is used; the first row contains data for timestep I = 1, the
second row for timestep I =2, etc. The first column contains the time
t at the timesteps.
The remaining columns
(column number 2 and further) contain values that are assigned to
Result. These values must be in the domain of the data type that is  given to Result. Each column contains data for an unique identifier. Column number 2 is associated with an unique identifier 1, column number 3 with an unique identifier 2, etc. In general, starting with the second column, a column number :emphasis:`c` is related to an unique identifier :emphasis:`c`- 1. The columns must be separated by one or more whitespace characters (spaces, tabs), the number of characters does not matter.  



2) a plain time series file 



This is a file formatted like the time series file with header, but without
header lines.




During operation, the values in a column of TimeSeries are assigned to the cells of Result that have an unique identifier on idexpression that corresponds with the unique identifier value associated with the column of TimeSeries. So, idexpression must contain a set of whole unique identifier cell values, starting with a value one, that is related to the unique identifiers of the columns in TimeSeries: cells with a value 1 are assigned data from the second column in TimeSeries, cells with a value 2 are assigned data from the third column etc.   



An example of a time series file with header is given in the table explained above. It is meant for
a dynamic model with starttime 1, endtime 6 and timeslice 1. It gives
precipitation for three areas. For instance, values for the area of station B
are in the third column. These are assigned to cells of Result that have a value 2 on idexpression.   

Notes
=====


Cells with an unique identifier on idexpression that is not represented by a column associated with the same unique identifier on TimeSeries are assigned a missing value on Result. For instance: let TimeSeries contain three columns, the first column with the time and two columns with data associated with unique identifiers 1 and 2 respectively. All cells of Result that have a idexpression value different from 1 or 2 are assigned a missing value.  



A timeseries generated and stored in the database during a model run by
the report keyword (or the timeoutput operator) cannot be imported during the same model run with the timeinput... operator.  

Group
=====
This operation belongs to the group of  Time operators 

Examples
========
