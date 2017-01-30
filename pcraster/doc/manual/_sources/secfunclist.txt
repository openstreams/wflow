.. _secfunclist:

**********************************************
Functional list of applications and operations
**********************************************


.. _secpointop:

Point operators
===============


.. _groupbool:

Boolean operators
-----------------
:ref:`and`
  Performs a Boolean-and operation on two expressions, on a cell-by-cell basis.
:ref:`not`
 Performs a Boolean-not operation on two expressions, on a cell-by-cell basis.
:ref:`or`
 Performs a Boolean-or operation on two expressions, on a cell-by-cell basis.
:ref:`xor`
 Performs a Boolean-xor operation on two expressions, on a cell-by-cell basis.
 



.. _groupcomp:

Comparison operators
--------------------
:ref:`eq or == <eq>`
  Performs a relational-equal-to operation on two expressions, on a cell-by-cell basis.
:ref:`ge or >= <ge>`
 Performs a relational-greater-than-or-equal-to operation on two expressions, on a cell-by-cell basis.
:ref:`gt or > <gt>`
 Performs a relational-greater-than operation on two expressions, on a cell-by-cell basis.
:ref:`le or \<= <le>`
 Performs a relational-less-than-or-equal-to operation on two expressions, on a cell-by-cell basis.
:ref:`lt or \< <lt>`
 Performs a relational-less-than operation on two expressions, on a cell-by-cell basis.
:ref:`ne or != <ne>`
 Performs a relational-not-equal-to operation on two expressions, on a cell-by-cell basis.
 



.. _groupcond:

Conditional statements
----------------------
:ref:`if then else <ifthenelse>`
 For each cell a Boolean expression determines whether the value of the first expression or the value of a second expression is assigned to the result
:ref:`if then <ifthen>`
  Return missing values if condition is not met.
:ref:`repeat { ...} until <secrepeatuntil>` (pcrcalc)
  Conditional iteration in pcrcalc.
 



.. _groupmissvalue:

Missing value creation, detection, alteration
---------------------------------------------
:ref:`cover`
  Substitutes missing values on an expression for values selected from one or more different expression(s), on a cell-by-cell basis.
:ref:`defined`
 Assigns a Boolean TRUE for non missing values on the input expression and FALSE for missing values, on a cell-by-cell basis.
:ref:`lddmask`
 Cuts a local drain direction map resulting in a (smaller) sound local drain direction map.
:ref:`nodirection`
 For an expression of directional data type, returns TRUE for cells without a direction and FALSE otherwise for cells with a direction.
:ref:`if then <ifthen>`
  Return missing values if condition is not met.
:ref:`inversedistance`
 Interpolate values
 



.. _groupreltab:

Relations in tables
-------------------
:ref:`lookup... <lookup>`
 For each cell, compares the cell value(s) of one or more expression(s) with the search key in a table and assigns a new value linked to that record in the key which matches the value(s) of the input expression
:ref:`lookuplinear`
 Assigns table key values with possible interpolation between key values.
:ref:`table`
 Creates on basis of one or more maps a table with a score for each key in the table. The score is the total area of the cells that match the key in the table.
 



.. _grouporder:

Order
-----
:ref:`order`
  Returns ordinal numbers to cells in ascending order.
:ref:`markwhilesumge, markwhilesumle <markwhilesum>`
    Marks each cell in specified order until the sum reaches a specified limit.
:ref:`areaorder`
 Within each area ordinal numbers to cells in ascending order.
:ref:`argorder,argorderwithid <argorder>`
 Identify highest value by argument order
:ref:`argorderarealimited,argorderwithidarealimited <argorderarealimited>`
 Identify highest value by argument order with a limit per argument
:ref:`argorderaddarealimited argorderwithidaddarealimited <argorderaddarealimited>`
 Variation on argorder
:ref:`pred`
 For each cell returns an ordinal number which is the ordinal number of the next lower ordinal class (predecessor) on the expression.
:ref:`succ`
 For each cell returns an ordinal number which is the ordinal number of the next higher ordinal class (predecessor) on the expression.
 



.. _groupmaxmin:

Maximize, minimize
------------------
:ref:`max`
  For each cell, determines the maximum value of multiple expressions and assigns it to the corresponding cell for the result.
:ref:`min`
 For each cell, determines the minimum value of multiple expressions and assigns it to the corresponding cell for the result.
 



.. _groupmath:

Arithmic operators, trigonometric, exponential, logarithmic functions
---------------------------------------------------------------------
:ref:`* <ster>`
  Multiplies the values of two expressions and sends this product to the result, on a cel-by-cell basis.
:ref:`** <sterster>`
 Calculates the Nth power of the first expression, where N is the value on a second expression and sends it to the result, on a cell-by-cell basis.
:ref:`- <minus>`
 Subtracts the value of the second expression from the value of the first expression and assigns it to the result, on a cell-by-cell basis. 
:ref:`+ <plus>`
 Adds the values of two expressions and assigns this sum to the result, on a cell-by-cell basis.
:ref:`/ or div <slash>`
 Divides the value of a first expressions by the value of a second expression and assigns this quotient to the result, on a cell-by-cell basis.
:ref:`abs`
 Calculates the absolute value of an expression, on a cell-by-cell basis.
:ref:`acos`
 Calculates the inverse cosine value of an expression, on a cell-by-cell basis.
:ref:`asin`
 Calculates the inverse sine value of an expression, on a cell-by-cell basis.
:ref:`atan`
 Calculates the inverse tangent value of an expression, on a cell-by-cell basis.
:ref:`cos`
 Calculates the cosine of an expression, on a cell-by-cell basis.
:ref:`exp`
 Calculates the base\ :sub:`e` exponential of an expression, on a cell-by-cell basis.
:ref:`fac`
  Faculty or factorial of a natural positive number
:ref:`idiv`
 Divides (integer division) the values on a first expression by the values on a second expression and assigns this quotient to the result, on a cell-by-cell basis.
:ref:`ln`
 Calculates the natural logarithm (base\ :sub:`e`) exponential of an expression, on a cell-by-cell basis.
:ref:`log10`
 Calculates the (base\ :sub:`e`) logarithm of an expression, on a cell-by-cell basis.
:ref:`mod`
 Divides (integer division) the values on a first expression by the values on a second expression and assigns the remainder to the result, on a cell-by-cell basis.
:ref:`sin`
 Calculates the sine of an expression, on a cell-by-cell basis.
:ref:`sqr`
 Calculates the square of an expression, on a cell-by-cell basis.
:ref:`sqrt`
 Calculates the square root of an expression, on a cell-by-cell basis.
:ref:`tan`
 Calculates the tangent of an expression, on a cell-by-cell basis.
 



.. _groupround:

Rounding
--------
:ref:`roundup`
  For each cell, the value of an expression is rounded upwards. Values of the results will be whole numbers.
:ref:`rounddown`
 For each cell, the value of an expression is rounded downwards. Values of the results will be whole numbers.
:ref:`roundoff`
 For each cell, the value of an expression is rounded off. Values of the results will be whole numbers.
 



.. _grouptypeconvers:

Data types: Conversion and assignment
-------------------------------------
:ref:`boolean`
  Converts from nominal, ordinal, scalar, directional or ldd data type to a boolean data type or generates a map of boolean data type with one constant value.
:ref:`directional`
 Converts from boolean, nominal, ordinal, scalar or ldd data type to a directional data type or generates a map of directional data type with one constant value.
:ref:`ldd`
 Converts from boolean, nominal, ordinal, scalar or directional data type to a ldd data type or generates a map of ldd data type with one constant value.
:ref:`nominal`
 Converts from boolean, ordinal, scalar, directional or ldd data type to a nominal data type or generates a map of nominal data type with one constant value.
:ref:`ordinal`
 Converts from boolean, nominal, scalar, directional or ldd data type to a ordinal data type or generates a map of ordinal data type with one constant value.
:ref:`scalar`
 Converts from boolean, nominal, ordinal, directional or ldd data type to a scalar data type or generates a map of scalar data type with one constant value.
:ref:`spatial`
 Conversion of a non-spatial value to a spatial data type.
 



.. _grouppointfield:

Random number generation - cells
--------------------------------
:ref:`normal`
  For each cell that is TRUE on a Boolean expression, assigns a value taken from a normal distribution with mean 0 and standard deviation 1.
:ref:`uniform`
 For each cell that is TRUE on a Boolean expression, assigns a value taken from a uniform distribution between 0 and 1.
 



.. _groupcoord:

Coordinates, unique ID's
------------------------
:ref:`uniqueid`
  For each cell that is TRUE on a Boolean expression, assigns a unique whole value
:ref:`xcoordinate`
 For each cell that is TRUE on a Boolean expression, assigns the xcoordinate of the cell
:ref:`ycoordinate`
 For each cell that is TRUE on a Boolean expression, assigns the ycoordinate of the cell
 



.. _secneighop:

Neighbourhood operators
=======================


.. _groupwind:

Windows operations
------------------
:ref:`shift, shift0 <shift>`
 Shifts the value of each cell a specified number of cells in the assigned direction.
:ref:`window4total`
 Sum the values of the four surrounding cells.
:ref:`windowaverage`
 For each cell, finds the average of cell values within a specified square neighbourhood and assigns it to the corresponding cell for the result
:ref:`windowdiversity`
 For each cell, finds the number of unique values within a specified square neighbourhood and assigns it to the corresponding cell for the result
:ref:`windowhighpass`
 Increases spatial frequency within a specified square neighbourhood. For each cell, it calculates the sum of cell values of an expression in a specified surrounding neighbourhood; this is subtracted from the cell values itself multiplied by twice the number of cells in the surrounding neighbourhood. The result of this calculation is assigned to the corresponding cell for the result.
:ref:`windowmajority`
 For each cell, finds the most often occurring cell values within a specified square neighbourhood and assigns it to the corresponding cell for the result
:ref:`windowmaximum`
 For each cell, finds the maximum cell value within a specified square neighbourhood and assigns it to the corresponding cell for the result
:ref:`windowminimum`
 For each cell, finds the minimum cell value within a specified square neighbourhood and assigns it to the corresponding cell for the result
:ref:`windowtotal`
 For each cell, finds the sum of cell values within a specified square neighbourhood and assigns it to the corresponding cell for the result
 



.. _groupderelev:

Derivatives of elevation maps
-----------------------------
:ref:`aspect`
 For each cell, calculates the aspect using elevations from a digital elevation model.
:ref:`horizontan`
 Calculates the maximum tangent of the angles of neighbouring cells in the direction of the sun.
:ref:`lddcreate`
  Creates a local drain direction map expression using the 8 points pour algorithm with flow directions from each cell to its steepest downslope neighbour. Pits can be removed with pit removing threshold map(s).
:ref:`lddcreatedem`
 Creates a modified digital elevation model which fits the local drain direction map generated on the basis of the original digital elevation model (the elevation model is the input of the operation).
:ref:`plancurv`
 For each cell, calculates the planform curvature (i.e. curvature transverse to the slope) using elevations from a digital elevation model.
:ref:`profcurv`
 For each cell, calculates the profile curvature (i.e. curvature in the direction of the slope) using elevations from a digital elevation model.
:ref:`slope`
 For each cell, calculates the slope using elevations from a digital elevation model.
 



.. _groupspread:

Spread operations
-----------------
:ref:`influencesimplegauss`
  Simple unweighted gaussian decrease of influence from sources.
:ref:`spread`
  For each cell, calculates the friction-distance of the shortest material-distance path over a map with friction material from an identified source cell or cells to the cell under consideration.
:ref:`spreadmax`
 Like spread but with a maximum distance argument.
:ref:`spreadldd`
 For each cell, calculates the friction-distance of the shortest material-distance path over a map with friction material from an identified source cell or cells to the cell under consideration, where only paths are considered in downstream direction from the source cells.
:ref:`spreadlddzone`
 For each cell, determines the shortest friction-distance path over a map with friction from an identified source cell or cells to the cell under consideration, where only paths are considered in downstream direction from the source cells. The value of the source cell at the start of this shortest material-distance path is assigned to the cell under consideration.
:ref:`spreadzone`
 Determines for each cell the shortest friction-distance path over a map with friction material from an identified source cell or cells to the cell under consideration. The value of the source cell at the start of this shortest friction-distance is assigned to the cell under consideration.
 



.. _groupldd:

Operations with local drain direction maps
------------------------------------------

:ref:`ldddist`
 Calculates for each cell the material-distance of the path over a map with friction from the cell under consideration to the downstream nearest TRUE cell.
:ref:`slopelength`
 For each cell assigns the accumulative-material-distance of the longest accumulative-material-path upstream over the local drain direction network to one of the cells against the divide of its catchment.
 

:emphasis:`transport of material:`

:ref:`accuflux`
 For each cell calculates the accumulated amount of material that flows out of the cell into its neighbouring downstream cell. This accumulated amount is the amount of material in the cell itself plus the amount of material in upstream cells of the cell.
:ref:`accucapacityflux and accucapacitystate <accucapacity>`
 Transport input of material downstream over a local drain direction network; material is transported from one cell to its downstream cell when the transport capacity is exceeded, the remaining material is stored.
:ref:`accufractionflux and accufractionstate <accufraction>`
 Transport input of material downstream over a local drain direction network; a fraction of the material is transported to its downstream cell, the remaining material is stored.
:ref:`accuthresholdflux and accuthresholdstate <accuthreshold>`
 Transport material downstream over a local drain direction network when transport threshold is exceeded
:ref:`accutriggerflux and accutriggerstate <accutrigger>`
 Transport input of material downstream over a local drain direction network; transport from one cell to its downstream cell only takes place if the amount of the material input to the cell exceeds the transport trigger of the cell: if the trigger is exceeded, all material is transported, else all material is stored.
:ref:`accutraveltimestate, accutraveltimeflux <accutraveltime>`
 Transports material downstream over a distance dependent on a given velocity.
:ref:`accutraveltimefractionstate, accutraveltimefractionflux, accutraveltimefractionremoved <accutraveltimefraction>`
 Transports a fraction of material downstream over a distance dependent on a given velocity.
:ref:`dynamicwaveq, dynamicwaveh <dynamicwave>`
 Dynamic Wave equation
:ref:`dynwavestate, dynwaveflux, lookuppotential, lookupstate <dynwave>`
 Dynamic wave equation as state and flux.
:ref:`kinematic`
 Kinematic Wave equation
:ref:`kinwavestate, kinwaveflux <kinwave>`
 Kinematic Wave equation as state and flux
 

:emphasis:`miscellaneous operations:`

:ref:`catchment`
  Determines the catchment(s) (watershed, basin) of each one or more specified cells, subcatchments are not identified.
:ref:`catchmenttotal`
  Total catchment for the entire upstream area
:ref:`downstream`
 Returns the value of the neighbouring downstream cell.
:ref:`downstreamdist`
 Returns the distance to the first cell downstream.
:ref:`lddrepair`
 Repairs an unsound local drain direction map.
:ref:`path`
 Determines for each TRUE cell on a Boolean input expression the path over the local drain direction network downstream to its pit; on the result, each cell which is on a path is assigned a TRUE.
:ref:`pit`
 Assigns a Boolean TRUE to all pit cells on a local drain direction network.
:ref:`streamorder`
 Assigns the stream order index to all cells on a local drain direction network.
:ref:`subcatchment`
 Determines the (sub-)catchment(s) (watershed, basin) of each one or more specified cells, subcatchments are identified.
:ref:`transient`
 Simulates transient groundwater flow according to the implicit finite difference method.
:ref:`upstream`
 For each cell, assigns the sum of the cell values of its upstream cell(s).
 



.. _groupvisi:

Operations for visibility analysis
----------------------------------
:ref:`extentofview`
 Determines the total length of the lines in a number of directions from the cell under consideration to the first cell with a different value.
:ref:`view`
 Assigns a TRUE or FALSE value for each cell on the result according to the visibility of that cell from one or more viewpoint cells in a 3D landscape defined by a digital elevation model.
 



.. _secareaop:

Area operations
===============


.. _grouparea:

Operations over areas
---------------------
:ref:`areaarea`
  For each cell, assigns the area of the area to which the cell belongs. Areas are identified by cell values on a expression with classes.
:ref:`areaaverage`
 For each cell, assigns the average value of the cells that belong to the same area to the cell itself. Areas are identified by cell values on a expression with classes.
:ref:`areadiversity`
 For each cell, assigns the number of unique cell values that belong to the same area to the cell itself. Areas are identified by cell values on a expression with classes.
:ref:`areamajority`
 For each cell, assigns the most often occurring cell value of cells that belong to the same area to the cell itself. Areas are identified by cell values on a expression with classes.
:ref:`areamaximum`
 For each cell, assigns the maximum value of the cells that belong to the same area to the cell itself. Areas are identified by cell values on a expression with classes.
:ref:`areaminimum`
 For each cell, assigns the minimum value of the cells that belong to the same area to the cell itself. Areas are identified by cell values on a expression with classes.
:ref:`areatotal`
 For each cell, assigns the sum of cells of cells that belong to the same area to the cell itself. Areas are identified by cell values on a expression with classes.
:ref:`clump`
 Identifies all continuous groups of with the same value ('clumps'); cells belonging to one clump are assigned the same new unique value.  
 



.. _groupareafield:

Random number generation - areas
--------------------------------
:ref:`areanormal`
 Assigns to each area one value taken from a normal distribution with a mean 0 and a standard deviation 1.
:ref:`areauniform`
 Assigns to each area one value taken from a uniform distribution with a mean 0 and a standard deviation 1.
 



.. _secmapop:

Map operations
==============


.. _groupmap:

Operations over maps
--------------------
:ref:`mapmaximum`
 Determines the maximum cell value of all cells values.
:ref:`mapminimum`
 Determines the minimum cell value of all cells values.
:ref:`maptotal`
 Sums all cell values
:ref:`maparea`
 Calculates the total area represented by the non missing value cells.
:ref:`cellarea`
 Assigns the area of one cell.
:ref:`celllength`
 Assigns the length which is identical in vertical and horizontal direction of one cell.
 



.. _groupmapfield:

Random number generation - map
------------------------------
:ref:`mapnormal`
 Assigns to all cells one non spatial value taken from a normal distribution with a mean 0 and a standard deviation 1.
:ref:`mapuniform`
 Assigns to all cells one non spatial value taken from a uniform distribution with a mean 0 and a standard deviation 1.
 



.. _sectimeop:

Time operations
===============


.. _grouptime:

Time operations
---------------
:ref:`time`
  Assigns for each time step the time at that time step.
:ref:`timeinput`
 Assigns for each time step one of a set of maps in the database. Each time step, the map is taken with the extension that refers to the time at the time step.
:ref:`timeinput...`
 For each time step assigns cell values read from a time series that is linked to a map with unique identifiers. Per time step, the time series gives for each unique identifier a cell value that is assigned to cells on the map with a corresponding unique identifier.
:ref:`timeinputsparse`
   Returns a map for each timestep where map-timesteps combinations can be missing.
:ref:`timeinputmodulo`
   Returns a map for each timestep using a modulo index into a shorter map stack.
:ref:`lookupmapstack`
 Reads a variable assigned map number from a map stack.
:ref:`timeoutput`
 For each cell writes the expression value of an uniquely identified cell or cells to a time series. After a model run, the time series contains for each identified cell a list of expression cell values per time step.
:ref:`timeslice`
 Assigns the timeslice.
 



.. _secdatamanop:

Data management
===============


.. _groupcreaclone:

map creation, changing attributes
---------------------------------
:ref:`mapattr`
  Generates a new PCRaster map with attributes specified by the user or changes location attributes of an existing PCRaster map.
 



.. _groupconvdata:

Conversion of data
------------------
:ref:`asc2map`
 Converts from ascii file format (including ARC/INFO and GENAMAP ascii output) to PCRaster map format.
:ref:`col2map`
 Converts from column file format (including simplified Geo-EAS format used in the GSTAT module of PCRaster) to PCRaster map format.
:ref:`map2asc`
 Converts from PCRaster map format to ascii file format (including ascii input format for ARC/INFO.
:ref:`map2col`
 Converts from PCRaster map format to column file format (including simplified Geo-EAS format also used in the GSTAT module of PCRaster).
 



.. _groupcutjoin:

Cutting and joining together PCRaster maps
------------------------------------------
:ref:`resample`
 Cuts one PCRaster map or joins several PCRaster maps by resampling to the cells of the resulting PCRaster map.
 



.. _groupleg:

Generation of legends
---------------------
:ref:`legend`
 Attaches a legend to or changes the legend of one or more maps.
 



.. _groupscreen:

Screen output
-------------

Visualisation of maps and timeseries can be done with aguila
 
