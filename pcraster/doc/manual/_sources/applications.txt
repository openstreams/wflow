.. _applications:

############
Applications
############
This chapter describes the usage of all applications in the PCRaster package.
To get a the usage on the command line type the command name without any options.
For example, typing:

::

   asc2map

Will give:

::

  asc2map version: Apr 23 2011 (linux/x86_64)
  USAGE:   asc2map AsciiInputFile MapOutputFile
   BLNOSDV  data type for map file.
   s c      separator character (default none or ',')
   m s      MV in AsciiInputFile (default 1E31)
   a        AsciiInputFile is Arc/Info gridascii output
   g        AsciiInputFile is genamap (report) format
   r #      nr. of lines to skip before each row
   h #      nr. of lines to skip for header
  --clone f clone file

.. toctree::
   :maxdepth: 1

   app_asc2map.rst
   app_col2map.rst
   app_legend.rst
   app_map2asc.rst
   app_map2col.rst
   app_mapattr.rst
   app_pcrcalc.rst
   app_resample.rst
   app_table.rst
