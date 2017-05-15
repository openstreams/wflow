Release notes
=============


2017.01
-------

+ wflow\_sbm renamed to wflow_sbm\_old
+ wflow\_sbm2 renamed to wflow\_sbm
+ Changes to hbv (quickflow)


2016.04
-------
.. note::

    Several none-backwards compatible changes will be part of this release. Use 2016.03 for older models

+ update soil names in sbm to the sbm2 names:
  + FirstZoneKsatVer -> KsatVer
  + FirstZoneMinCapacity -> SoilMinThickness
  + FirstZoneCapacity (FirstZoneThickness) -> SoilThickness
  + FirstZoneCapacity -> SoilWaterCapacity
  + FirstZoneFlux -> SatWaterFlux
  + FirstZoneDepth -> SatWaterDepth
+ [model]reinit moved to [run]reinit (same for fewsrun) in all models
+ added [rollingmean] section in framework
+ updates to topoflex
+ updates to BMI framework
+ updates to netcdf reader. [model] modeltype= can be specified
+ updated wflow_hbv to better resemble hbv-96 for the upper zone when kquickflow is determined (and not specified directly). This may break older calibrations


2016.03
-------
+ last tag before moving to new names in SBM
+ irrigation added to SBM/SBM2

2016.02
-------
+ added better BMI support
+ added bmi2runner and wflow_bmi_combined
+ updated date/time framework
+ added wflow_topoflex model
+ added reservoir support to wflow_routing
+ added support for other functions apart from averaging in output of time series
+ wflow_delwaq support netCDF files
+ wflow_sbm2 updated


2015.02
-------
+ Updated netcdf support -> now needs the pyproj module
+ Updated wflow_routing with two-layer (floodplain)
+ Redone multiplication of parameter framework. Now possible via ini and BMI
+ Added .mult postfix for tbl files to apply multiplication (untested)
+ Added bmi support
+ wflow_sbm s-curve sharpness of height distribution now based on upper and lower
  half of the distribution (average)
+ Added separate soil evap and RootZoneSoilMoisture to wflow_sbm
+ Added support in wflow_SBM to get vegetation parameters for Gash from LAI
+ Moved non-mature scripts to SandBox dir
+ Added unit tests for wflow_sbm, hbv and bmi

+ Started with wflow_sbm2 -> New version with new (more logical) variable names. This version will
  also have a better unsaturated zone (for G4INDO project). New Names:

  + FirstZoneKsatVer -> KsatVer
  + FirstZoneMinCapacity -> SoilMinThickness
  + FirstZoneCapacity (FirstZoneThickness) -> SoilThickness
  + FirstZoneCapacity -> SoilWaterCapacity
  + FirstZoneFlux -> SatWaterFlux
  + FirstZoneDepth -> SatWaterDepth

2015.01
-------
+ support for scalar input to wflow\_sbm/hbv has been removed.
+ added wf_updateparameters to framework. This allows the user to set parameters to update
  in the ini file but also list them in the parameters function in the model itself. This
  functionality should replace all manual reading of forcing data and static parameters


Version 1.0 RC7
---------------
unsupported interim release

+ added  HBV type lower zone to wflow\_sbm. Use MaxPercolation > 0 to use this zone. MaxLeakege > 0 will send
  water outside of the model
+ Test version of the wflow_W3RA model
+ Made two lateral flow options for sbm
+ Stopped support for pcraster version 3 and python 2.6
+ removed all the try/except from importing wflow. Now you
  NEED to install wflow as a package
+ Added seperate wflow\_routing module that includes the kinematic wave routing. This part will be removed from the
  wflow\_sbm and wflow\_hbv models
+ Added check in gash interception not to have more interception than available potential evap
+ Fixed capillary rise calculation to include timestep. This means that sub-daily models may need to be recalibrated

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