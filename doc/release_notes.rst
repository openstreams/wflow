Release notes
=============

2019.1
------

* Wflow\_sbm, redesign of lateral subsurface flow because of a bug that resulted in an overestimation of lateral subsurface flux:

    * The kinematic wave approach is used to route subsurface flow laterally.
    * Kinematic wave solution for surface (overland and river) and subsurface flow added to wflow\_funcs, and used by wflow\_sbm. 
      Numba (open-source JIT compiler) is used to accelerate these functions. This replaces the PCRaster kinematic wave solution.
    * A seperate kinematic wave reservoir for river (self.RiverRunoff) and overland flow (self.LandRunoff)
    * Lateral subsurface flow and overland flow feed into drainage network (channel) based on ratio slopes
      (slope upstream cell / (slope upstream cell + slope cell))
    * Option to provide land slope and river slope (staticmaps)
    * Removed sub-grid runoff generation
    * Removed re-infiltration of surface water

* Wflow\_hbv uses kinematic wave solution for surface flow available in wflow\_funcs (replacement of PCRaster kinematic wave solution)
* Option to enable iteration (smaller time steps) for kinematic wave solution for surface runoff (wflow\_hbv and wflow\_sbm)
* Added hydrological model wflow\_w3 to wflow framework, based on Australian Water Resources Assessment Landscape model (AWRA-L),
  an improved version compared to wflow\_W3RA.
* Added hydrological model wflow\_stream to wflow framework (STREAM (Spatial Tools for River Basins and Environment and Analysis of Management Options)).
* Added experimental version of wflow\_sediment to simulate soil erosion and delivery to the river system.
* Added the wflow\_emwaq module: provides a set of functions to create and link a Delft3D-WAQ or D-Emission model to a wflow model.
* Update of lake modelling for wflow\_hbv and wflow\_sbm (including also Modified Puls Approach) in wflow\_funcs
* Update of glacier model (used by wflow\_hbv and wflow\_sbm)
* Parameter Manning N of river can be linked to streamorder for wflow\_hbv
* setuptools_scm used for version wflow package



2018.1
-------

* -T and -S options now use date/time strings to better support netcdf
* date/time handling improved
* fews adapter and support simplified:

    * fewsrun removed
    * adapter need less arguments

* netcdf handling improved:

    * check for date/time in input and determines offset if there is a problem.
    * netcdf states now work (checks for date/time)

* wflow\_sbm:

    * Feddes reduction of root water uptake
    * Improved/fixed handling of open water and soil evaporation
    
* Added hydrological model PCR-GLOBWB ((PCRaster Global Water Balance) version v2.1.0_beta_1) (wflow\_pcrglobwb) to wflow framework
* Added hydrological model SPHY (Spatial Processes in HYdrology) version 2.1 to the wflow framework
* Simple reservoirs: possible to include reservoir areas as static map
* Parameter Manning N of river can be linked to streamorder for wflow\_sbm
* HBV lake functionality also available in wflow\_sbm


2017.01
-------

* wflow\_sbm renamed to wflow_sbm\_old
* wflow\_sbm2 renamed to wflow\_sbm:

    * This breaks older models. You will need to regenerate initial conditions. In the new model the unsaturated part of
      the soil can be split into several layers. With one layer the results should be the same but the name of the
      state file is now UStoreLayerDepth_0.map for the first layer (or for the whole layer if only one layer is
      defined). It used to be UStoreDepth.map

* Changes to wflow\_hbv:

    * Quickflow calculation
    * Addition of lake modelling (based on HBV96, complex reservoirs)
    * Simple reservoirs functionality also available in wflow\_hbv
    
* Added rice crop growth model LINTUL (wflow\_lintul) to wflow framework and coupled to wflow\_sbm (BMI)
* Added irrigation of rice fields to wflow\_sbm


2016.04
-------
.. note::

    Several none-backwards compatible changes will be part of this release. Use 2016.03 for older models

* update soil names in sbm to the sbm2 names:

    * FirstZoneKsatVer -> KsatVer
    * FirstZoneMinCapacity -> SoilMinThickness
    * FirstZoneCapacity (FirstZoneThickness) -> SoilThickness
    * FirstZoneCapacity -> SoilWaterCapacity
    * FirstZoneFlux -> SatWaterFlux
    * FirstZoneDepth -> SatWaterDepth
  
* [model]reinit moved to [run]reinit (same for fewsrun) in all models
* added [rollingmean] section in framework
* updates to topoflex
* updates to BMI framework
* updates to netcdf reader. [model] modeltype= can be specified
* updated wflow_hbv to better resemble hbv-96 for the upper zone when kquickflow is determined
  (and not specified directly). This may break older calibrations


2016.03
-------
* last tag before moving to new names in SBM
* irrigation added to SBM/SBM2


2016.02
-------
* added better BMI support
* added bmi2runner and wflow_bmi_combined
* updated date/time framework
* added wflow_topoflex model
* added reservoir support to wflow_routing
* added support for other functions apart from averaging in output of time series
* wflow_delwaq support netCDF files
* wflow_sbm2 updated


2015.02
-------
* Updated netcdf support -> now needs the pyproj module
* Updated wflow_routing with two-layer (floodplain)
* Redone multiplication of parameter framework. Now possible via ini and BMI
* Added .mult postfix for tbl files to apply multiplication (untested)
* Added bmi support
* wflow_sbm s-curve sharpness of height distribution now based on upper and lower
  half of the distribution (average)
* Added separate soil evap and RootZoneSoilMoisture to wflow_sbm
* Added support in wflow_SBM to get vegetation parameters for Gash from LAI
* Moved non-mature scripts to SandBox dir
* Added unit tests for wflow_sbm, hbv and bmi

* Started with wflow_sbm2 -> New version with new (more logical) variable names. This version will
  also have a better unsaturated zone (for G4INDO project). New Names:

    * FirstZoneKsatVer -> KsatVer
    * FirstZoneMinCapacity -> SoilMinThickness
    * FirstZoneCapacity (FirstZoneThickness) -> SoilThickness
    * FirstZoneCapacity -> SoilWaterCapacity
    * FirstZoneFlux -> SatWaterFlux
    * FirstZoneDepth -> SatWaterDepth

2015.01
-------
* support for scalar input to wflow\_sbm/hbv has been removed.
* added wf_updateparameters to framework. This allows the user to set parameters to update
  in the ini file but also list them in the parameters function in the model itself. This
  functionality should replace all manual reading of forcing data and static parameters


Version 1.0 RC7
---------------
unsupported interim release

* added  HBV type lower zone to wflow\_sbm. Use MaxPercolation > 0 to use this zone. MaxLeakege > 0 will send
  water outside of the model
* Test version of the wflow_W3RA model
* Made two lateral flow options for sbm
* Stopped support for pcraster version 3 and python 2.6
* removed all the try/except from importing wflow. Now you
  NEED to install wflow as a package
* Added seperate wflow\_routing module that includes the kinematic wave routing. This part will be removed from the
  wflow\_sbm and wflow\_hbv models
* Added check in gash interception not to have more interception than available potential evap
* Fixed capillary rise calculation to include timestep. This means that sub-daily models may need to be recalibrated

Version 1.0 RC5
---------------
unsupported interim release

* netcdf reading and writing added (filename should be configured in ini file, framework section: netcdfoutput, netcdfwritebuffer, netcdfinput)
* summary sections (summary, summary_max, symmary_avg, ect) added to ini file to save maps at end of run
* added option to save flow per subcatchment by setting pits at the end of each subcatchment in the ldd
* added new tbl file for wflow_sbm (et_reftopot.tbl). Used to covert reference ET to potential ET. Set to 1 by default
* better representation of open water ET in wflow_sbm
* wflow_adapt can now convert log files to XML for Delft-FEWS

Version 1.0 RC4
---------------

unsupported interim release

* tss (and csv) output refactored. The ini file can now hold multiple outputtss sections each with a diffrent maps for extracting/averaging