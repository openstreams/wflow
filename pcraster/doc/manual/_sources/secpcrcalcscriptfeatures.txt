.. _secpcrcalcscriptfeatures:

***********************
pcrcalc script features
***********************

.. index::
   single: error messages

.. _pcrcalcerrormessage:

error messages
==============

All error message of pcrcalc have the format:

 | scriptFileName:lineNr:colNr:ERROR: A descriptive message

The message is a list of positions separated by the symbol :

scriptFileName
    Name of the script file or the ?-symbol if a statement is issued from the command line.
lineNr
    The line number within the script file.
colNr
    The column number where the error is detected.

In most cases, pcrcalc will point to the exact position where the error is present. There are a few occasions that the actual (syntax) error is before the position that pcrcalc reports. One of these occasions is the generic syntax error like:

::

   example.mod:5:36:ERROR: Syntax error at symbol ';'

In this example, the symbol ";" at line 5, character position 36 is not expected. This may be caused by earlier symbols in the script. 

.. index::
   single: repeat until
.. index::
   single: conditional iteration


.. _secrepeatuntil:

repeat until
============


Conditional iteration (looping) allows to perform model code a conditional
number of times. The syntax is:

::

  repeat {
  Statements
  } until BooleanCondition ;


The Statements are executed repeatedly until the BooleanCondition is true. For
example:

::

  nrLoopsExecuted=0;
  repeat {
   nrLoopsExecuted=nrLoopsExecuted+1;
  } until nrLoopsExecuted eq 4;

The statement nrLoopsExecuted=nrLoopsExecuted+1 is executed 4 times.

If the BooleanCondition is spatial (it evaluates to a map) then the Statements
are executed until no cells in BooleanCondition are false(0). In other words
execution of the loop is stopped if all cells are true(1) or missing value.
Note that the execution is not stopped at a per cell base, but when ALL cells
meet the specified condition. In theory, stopping at a per cell base might
decrease execution time, but the additonal bookkeeping of which cells needs one
iteration more at each loop will neutralise this gain.
A more realistic example is an iteration to fit a certain function, with a
series of succesive estimates that will stop if the next estimate only differs
a certain small epsilon value.

::

  prevEst= ...initial estimate ... ;
  repeat {
    nextEst= ... estimate ... ;
    difference = abs(nextEst-prevEst);
    prevEst=nextEst;
  } until difference < 0.00001;

If the result of nextEst is a map, then this loop will continue until each cell
in difference has a value smaller than 0.000001.

..
  ***** Current bug and possible work around *****
  In the oldcalc version, the repeat-until construction has a serious flaw: under
  certains conditions map data that is referenced within in the repeat code block
  for the last time in the script is prematurely released. This bug will
  terminate the model run with a message like: "RUN TIME ERROR prevEst has no
  value." The solution is to define a dummy operation on prevEst just after the
  repeat until construct:
  prevEst= ...initial estimate ... ;
  repeat {
    nextEst= ... estimate ... ;
    difference = abs(nextEst-prevEst);
    prevEst=nextEst;
  } until difference < 0.00001;
  # dummy operation
  prevEst = prevEst + 0;

.. index::
   pair: report; selective

.. _selectivereports:

selective reports
=================

In scripts one can define at which timesteps a
map must be written to the permanent database (in other words your harddisk).

This feature consists of the following extensions to the modelling language:

- Defining a sequence of report moments in the timer section of the script by a symbolic name reportName.
- Using that symbolic name in the report statement in the dynamic section: report(reportName).
- Definining a sequence of report moments in the report statement in the dynamic section: report(reportMoments).
- Introduction of the keyword endtime as the symbolic name for the last timestep.
- Introduction of the keyword reportdefault as the symbolic name for the reports without an explicit definition of the report moments.
- Note that selective reports only apply for map stacks, when writing to a timeseries file they are ignored, a timeseries file will always contains data for all timesteps except if the option :ref:`----noheader <noheader>` is specified.

To explain this feature, we use the following minimal script

::

  timer 1 1000 1;
  initial
  dynamic
   report stack1_ = input.map;
   report stack2_ = input.map;

The script simply copies the input.map to 2000 maps named stack1_0.001,
stack2_0.002 up to stack1_1.000 and stack2_0.001 up to stack2_1.000.

Now suppose we want to write to stack stack1\_ only at step 1, 10, 900, 1000 and
stack2\_ at every fifth step. We then define and use two report moment
definitions in the timer sections and apply them to the reports:

::

  timer 1 1000 1;
   rep1 = 1,10,900,endtime;
   rep2 = 5+5..endtime;
  initial
  dynamic
   report(rep1) stack1_ = input.map;
   report(rep2) stack2_ = input.map;

Different moments in a report moment definition are seperated by ',' as in
1,10,900,endtime that defines 4 moments. A range of moments can be given by the
syntax start+step..end. Step increases the moments up and including end. In the
example above stack2\_ will be reported at timestep 5,10,15,25, etc. until 1000.
If the + is omitted a step value of 1 is assumed. In other words 5..10 and
5+1..10 will both result in 5,6,7,8,9,10.

A report moment definition can also be placed within the report statements:

::

  timer 1 1000 1;
  initial
  dynamic
   report(endtime)   stack1_ = input.map;
   report(1,5+3..12) stack2_ = input.map;

In the example above stack1\_ is only reported at the last step and stack2\_ at
1,5,8,11.

A special report definition is reportdefault. Defining this one in the timer
section causes all stack reports without explicit moments to report only at
particular moments:

::

  timer 1 1000 1;
   reportdefault = 900+5..endtime;
  initial
  dynamic
   report stack1_ = input.map;
   report stack2_ = input.map;

Both stack1\_ and stack2\_ are reported at steps 900,905,910..1000.

.. index::
   single: argument substitution
.. index::
   single: $

.. _secsubstitution:

substitution of arguments
=========================

In both command line expressions and model scripts, parts of the model can be
substituted from shell-like arguments. For example when sum_2.mod contains:

::
  $1 = $2 + $3;

Then calling

::

  pcrcalc -f sum_2.mod sum_2.map add_1.map  add_2.map

Will yield execution of:

::

  sum_2.map = add_1.map + add_2.map;


On the command line, the model and the arguments need to be separated by ;; .
For example:

::

   pcrcalc "$1 = $2 + $3 ;;" sum_2.map add_1.map add_2.map


Basic substitution rules:

- everything, except comments, may substituted, even operators, functions, etc.
- if a substitution fails, blanks are returned: no error message is
  printed, the resulting string might be incorrect. (use -t to check
  substitution without running the model)
- A $-sign followed by a number refers to an argument on the command line, e.g. $1 $2
- If the $-sign is followed by a non-numerical string the argument is a
  shell or environment variable, for example $RESULT.
- Simple arguments may be enclosed in curly braces, e.g. $1 = $2 + $3
  equals ${1} = $2 + $3;
  This enables prefix-substitution of a variable such as ${1}.map = $2 + $3;

Advanced substitution rules:

 - A range of parameters can be given in the ${from,to} construct. For example: $1 = max(${2,3}) ;; max.map in1.map in2.map becomes max.map = max(in1.map,in2.map)
 - n in the first or second arguments denotes the number of arguments: $1 = max(${2,n}) ;; max.map in1.map in2.map in3.map becomes max.map = max(in1.map,in2.map,in3.map)
 - The ${from,to} construct prints an ',' between every argument. Another argument separator can be given explicitly: $1 = ${2,n,+} ;; sum.map in1.map in2.map in3.map becomes sum.map = in1.map+in2.map+in3.map
 - A wrapper around each arguments can be given: $1 = ${2,n,+,sqrt($)} ;; sumsqrt.map in1.map in2.map in3.map becomes sumsqrt.map = sqrt(in1.map)+sqrt(in2.map)+sqrt(in3.map)
 - In the 4th argument, the wrapper, a $ is given for where the argument should be inserted. That $-sign is optional, allowing the following construct: $1 = (${2,n,+}) / (${2,n,+,1});; av.map in1.map in2.map in3.map becomes av.map = (in1.map+in2.map+in3.map) / (1+1+1)
 -  white space in argument 3 and 4 is kept in the substitution: $1 = ${2,n,and,not $} ;; andnot.map in1.map in2.map in3.map becomes andnot.map = not in1.mapandnot in2.mapandnot in3.map while $1 = ${2,n, and ,not $} ;; andnot.map in1.map in2.map in3.map becomes the correct andnot.map = not in1.map and not in2.map and not in3.map

Additional notes:

 - comments (#) are not allowed within a substitution
 - Note that ${n} denotes the last argument, while $n means an environment variable named n.
 - Note that is an extra substitution level, so one needs to give $$1 in an UNIX shell, if the command line expression is enclosed in "-symbols. This is not neccessary if it is enclosed in '-symbols.
 - If you are trying $-constructs and you get the message: ERROR: parse error near line '?' near symbol '?' Then two errors are frequent:

   - you forgot the delimit the expression and the arguments by ;;
   - the substitution went wrong. First try -t, to see what the result of the substitution is. As said before, the substitution mechanism silently ignores most errors, and returns blanks. So be warned.


