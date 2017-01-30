.. _reference:

Reference
*********

.. toctree::
   :maxdepth: 2

   ClassIndex
   VariableCollection

Parameter files
---------------

Parameter values can be stored in external files. The syntax of the file contents is simliar to the one of oldcalc. The file must contain the following columns:

First column
  VarName:  the name of the collection used in the model script

One (or more) column(s):
  Columns holding the index type names. The number of columns must equal the order of Index variables passed to the collection

One column:
  holding the parameter value. In case of a string entry a map will be read from disk. Paths including spaces must be enclosed by quotation marks, e.g. "New Folder\input.map"

Columns need to be separated by space or tabulator characters.

