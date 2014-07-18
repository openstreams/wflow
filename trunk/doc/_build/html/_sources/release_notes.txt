Release notes
=============

Version 1.0 RC5
---------------
unsupported interim release

+ netcdf reading and writing added (filename should be configured in ini file, framework section: netcdfoutput, netcdfwritebuffer, netcdfinput)
+ summary sections (summary, summary_max, symmary_avg, ect) added to ini file to save maps at end of run
+ added option to save flow per subcatchment by setting pits at the end of each subcatchment in the ldd
+ added new tbl file for wflow_sbm (et_reftopot.tbl). Used to covert reference ET to potential ET. Set to 1 by default
+ better representation of open water ET in wflow_sbm
+ wflow_adapt can now convert log files to XML for Delft-FEWS

Version 1.0 RC4
---------------

unsupported interim release

+ tss (and csv) output refactored. The ini file can now hold multiple outputtss sections each with a diffrent maps for extracting/averaging