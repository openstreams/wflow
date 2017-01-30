Known bugs
----------
The following PCRaster Python operations are known to be deficient:

nominal, ordinal
   Some nominal and ordinal numbers larger than 198000000 are not correct represented. To avoid the problem use lower class values. In some cases large nominal or ordinal numbers are used in models where scalar numbers are more natural.

argorder and related functions
   A variable number of arguments is currently not supported.
