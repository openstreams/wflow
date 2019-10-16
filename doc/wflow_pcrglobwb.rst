The wflow_pcrglobwb Model
=========================

Introduction
------------

In 2018 the following PCR-GLOBWB ((PCRaster Global Water Balance) version was added to the wflow framework:

https://github.com/UU-Hydro/PCR-GLOBWB_model/tree/v2.1.0_beta_1


Changes made to the PCR-GLOBWB code
-----------------------------------

+ The original code was converted from Python 2 to Python 3.

+ For the different modules available in PCR-GLOBWB:

 + groundwater
 + landCover
 + landSurface
 + ncConverter
 + parameterSoilAndTopo
 + routing
 + virtualOS
 + waterBodies
 
 only the initialization function (including the states) was changed, to initialize these classes from wflow_pcrglobwb properly,
 so PCR-GLOBWB complies with the directory structure of wflow. Furthermore the checks in virtualOS whether the clone map
 and input maps (e.g. forcing) have the same attributes (cell size and domain), are switched off, including possible cropping and 
 resampling of input maps when clone map and input maps don't have the same attributes.

+ The wflow framework takes care of saving output, and the Reporting module of PCR-GLOBWB is not included. Also, the spinUp module of
  PCR-GLOBWB is not included in wflow. Initial conditions (cold state) can be set in the wflow_pcrglobwb.ini file as follows:
  
  ::
    
    [forestOptions]
    # initial conditions:
    interceptStorIni  = 0.0
    snowCoverSWEIni   = 0.0
    snowFreeWaterIni  = 0.0
    topWaterLayerIni  = 0.0
    storUpp000005Ini  = 0.0
    storUpp005030Ini  = 0.0
    storLow030150Ini  = 0.0
    interflowIni      = 0.0
    
  or default initial conditions (set in the code) are used when these are not set in the ini file. Warm states can be set in the ini file
  as follows:
  
  ::
  
    [forestOptions]
    # initial conditions:
    interceptStorIni    = landSurface.interceptStor_forest.map
    snowCoverSWEIni     = landSurface.snowCoverSWE_forest.map
    snowFreeWaterIni    = landSurface.snowFreeWater_forest.map
    topWaterLayerIni    = landSurface.topWaterLayer_forest.map
    storUppIni          = landSurface.storUpp005030_forest.map
    storLowIni          = landSurface.storLow030150_forest.map
    interflowIni        = landSurface.interflow_forest.map
    
  and should be available in the instate directory of the Case directory.
    
    
Ini file settings
-----------------

Below an example of a wflow_pcrglobwb.ini file:

::

    [framework]
    netcdfoutput = outmaps.nc
    netcdfinput = inmaps/forcing.nc
    netcdfwritebuffer=20
    EPSG = EPSG:4326

    [run]
    # either a runinfo file or a start and end-time are required
    starttime= 2002-01-01 00:00:00
    endtime= 2002-01-30 00:00:00
    reinit = 0
    timestepsecs = 86400
    runlengthdetermination=steps

    [model]
    modeltype = wflow_pcrglobwb

    [layout]
    # if set to zero the cell-size is given in lat/long (the default)
    sizeinmetres = 0

    [outputmaps]
    self.routing.subDischarge = Qro
    self.routing.waterBodyStorage = wbs
    self.landSurface.storUpp = su1
    self.landSurface.storLow = slo
    self.landSurface.actualET = aet
    self.landsurface.swAbstractionFractionData = swAbsF
    self.landSurface.totalPotET = pet  
    self.landSurface.gwRecharge = gwr
    self.landSurface.snowCoverSWE = swe
    self.groundwater.nonFossilGroundwaterAbs = nFAbs
    self.landSurface.fossilGroundwaterAbstr = FAbs
    self.landSurface.irrGrossDemand = IrrGD
    self.landSurface.nonIrrGrossDemand = nIrrGD

    [globalOptions]

    # Map of clone (must be provided in PCRaster maps)
    # - Spatial resolution and coverage are based on this map:
    cloneMap = wflow_clone.map

    # The area/landmask of interest 
    landmask = mask.map
    # If None, area/landmask is limited for cells with ldd value. 

    [landSurfaceOptions]
    debugWaterBalance = True

    numberOfUpperSoilLayers = 2

    topographyNC     = topoProperties5ArcMin.nc
    soilPropertiesNC = soilProperties5ArcMin.nc

    includeIrrigation = True

    # a pcraster map/value defining irrigation efficiency (dimensionless) - optional
    irrigationEfficiency = efficiency.map

    # netcdf time series for historical expansion of irrigation areas (unit: hectares). 
    # Note: The resolution of this map must be consisten with the resolution of cellArea. 
    historicalIrrigationArea = irrigationArea05ArcMin.nc

    includeDomesticWaterDemand = True
    includeIndustryWaterDemand = True
    includeLivestockWaterDemand = True

    # domestic and industrial water demand data (unit must be in m.day-1)
    domesticWaterDemandFile = domestic_water_demand_version_april_2015.nc
    industryWaterDemandFile = industry_water_demand_version_april_2015.nc
    livestockWaterDemandFile = livestock_water_demand_version_april_2015.nc

    # desalination water supply (maximum/potential/capacity)
    #desalinationWater = desalination_water_version_april_2015.nc # should be included
    # zone IDs (scale) at which allocations of groundwater and surface water (as well as desalinated water) are performed  
    allocationSegmentsForGroundSurfaceWater = uniqueIds60min.nom_5min.map

    # predefined surface water - groundwater partitioning for irrigation demand (based on Siebert, 2010/2013: Global Map of Irrigation Areas version 5):
    irrigationSurfaceWaterAbstractionFractionData        = AEI_SWFRAC_5min.map
    irrigationSurfaceWaterAbstractionFractionDataQuality = AEI_QUAL_5min.map


    # predefined surface water - groundwater partitioning for irrigation demand (based on McDonald, 2014):
    maximumNonIrrigationSurfaceWaterAbstractionFractionData = max_city_sw_fraction_5min.map

    [forestOptions]
    name = forest
    debugWaterBalance = True

    # snow module properties
    snowModuleType      =  Simple
    freezingT           = 0.0
    degreeDayFactor     =  0.0025
    snowWaterHoldingCap =  0.1
    refreezingCoeff     =  0.05

    # other paramater values
    minTopWaterLayer = 0.0
    minCropKC        = 0.2
    minInterceptCap  = 0.0002

    landCoverMapsNC = forestProperties5ArcMin.nc

    # Parameters for the Arno's scheme:
    arnoBeta = None
    # If arnoBeta is defined, the soil water capacity distribution is based on this.
    # If arnoBeta is NOT defined, maxSoilDepthFrac must be defined such that arnoBeta will be calculated based on maxSoilDepthFrac and minSoilDepthFrac.

    cropCoefficientNC = cropKC_forest_daily366.nc
    interceptCapNC    = interceptCap_forest_daily366.nc
    coverFractionNC   = coverFraction_forest_daily366.nc

    # initial conditions:
    interceptStorIni  = landSurface.interceptStor_forest.map
    snowCoverSWEIni   = landSurface.snowCoverSWE_forest.map
    snowFreeWaterIni  = landSurface.snowFreeWater_forest.map
    topWaterLayerIni  = landSurface.topWaterLayer_forest.map
    storUppIni  = landSurface.storUpp005030_forest.map
    storLowIni  = landSurface.storLow030150_forest.map
    interflowIni      = landSurface.interflow_forest.map

    [grasslandOptions]
    name = grassland
    debugWaterBalance = True

    # snow module properties
    snowModuleType      =  Simple
    freezingT           = 0.0
    degreeDayFactor     =  0.0025
    snowWaterHoldingCap =  0.1
    refreezingCoeff     =  0.05

    # other paramater values
    minTopWaterLayer = 0.0
    minCropKC        = 0.2
    minInterceptCap  = 0.0002

    landCoverMapsNC = grasslandProperties5ArcMin.nc
    #
    # Parameters for the Arno's scheme:
    arnoBeta = None
    # If arnoBeta is defined, the soil water capacity distribution is based on this.
    # If arnoBeta is NOT defined, maxSoilDepthFrac must be defined such that arnoBeta will be calculated based on maxSoilDepthFrac and minSoilDepthFrac.

    cropCoefficientNC = cropKC_grassland_daily366.nc
    interceptCapNC    = interceptCap_grassland_daily366.nc
    coverFractionNC   = coverFraction_grassland_daily366.nc

    # initial conditions:
    interceptStorIni  = landSurface.interceptStor_grassland.map
    snowCoverSWEIni   = landSurface.snowCoverSWE_grassland.map
    snowFreeWaterIni  = landSurface.snowFreeWater_grassland.map
    topWaterLayerIni  = landSurface.topWaterLayer_grassland.map
    #storUpp000005Ini  = landSurface.storUpp000005_grassland.map
    storUppIni  = landSurface.storUpp005030_grassland.map
    storLowIni  = landSurface.storLow030150_grassland.map
    interflowIni      = landSurface.interflow_grassland.map

    [irrPaddyOptions]
    name = irrPaddy
    debugWaterBalance = True

    # snow module properties
    snowModuleType      =  Simple
    freezingT           = -0.0
    degreeDayFactor     =  0.0025
    snowWaterHoldingCap =  0.1
    refreezingCoeff     =  0.05
    #
    landCoverMapsNC = irrPaddyProperties30min.nc
    maxRootDepth     = 0.5
    #
    # Parameters for the Arno's scheme:
    arnoBeta = None
    # If arnoBeta is defined, the soil water capacity distribution is based on this.
    # If arnoBeta is NOT defined, maxSoilDepthFrac must be defined such that arnoBeta will be calculated based on maxSoilDepthFrac and minSoilDepthFrac.
    #
    # other paramater values
    minTopWaterLayer = 0.05
    minCropKC        = 0.2
    minInterceptCap  = 0.0002
    cropDeplFactor   = 0.2

    cropCoefficientNC = cropKC_irrPaddy_daily366.nc

    # initial conditions:
    interceptStorIni  = landSurface.interceptStor_irrPaddy.map
    snowCoverSWEIni   = landSurface.snowCoverSWE_irrPaddy.map
    snowFreeWaterIni  = landSurface.snowFreeWater_irrPaddy.map
    topWaterLayerIni  = landSurface.topWaterLayer_irrPaddy.map
    #storUpp000005Ini  = landSurface.storUpp000005_irrPaddy.map
    storUppIni  = landSurface.storUpp005030_irrPaddy.map
    storLowIni  = landSurface.storLow030150_irrPaddy.map
    interflowIni      = landSurface.interflow_irrPaddy.map

    [irrNonPaddyOptions]
    name = irrNonPaddy
    debugWaterBalance = True

    # snow module properties
    snowModuleType      =  Simple
    freezingT           = -0.0
    degreeDayFactor     =  0.0025
    snowWaterHoldingCap =  0.1
    refreezingCoeff     =  0.05
    #
    landCoverMapsNC  = irrNonPaddyProperties30min.nc
    maxRootDepth     = 1.0
    #
    # Parameters for the Arno's scheme:
    arnoBeta = None
    # If arnoBeta is defined, the soil water capacity distribution is based on this.
    # If arnoBeta is NOT defined, maxSoilDepthFrac must be defined such that arnoBeta will be calculated based on maxSoilDepthFrac and minSoilDepthFrac.
    #
    # other paramater values
    minTopWaterLayer = 0.0
    minCropKC        = 0.2
    minInterceptCap  = 0.0002
    cropDeplFactor   = 0.5

    cropCoefficientNC = cropKC_irrNonPaddy_daily366.nc

    # initial conditions:
    interceptStorIni  = landSurface.interceptStor_irrNonPaddy.map
    snowCoverSWEIni   = landSurface.snowCoverSWE_irrNonPaddy.map
    snowFreeWaterIni  = landSurface.snowFreeWater_irrNonPaddy.map
    topWaterLayerIni  = landSurface.topWaterLayer_irrNonPaddy.map
    #storUpp000005Ini  = landSurface.storUpp000005_irrNonPaddy.map
    storUppIni  = landSurface.storUpp005030_irrNonPaddy.map
    storLowIni  = landSurface.storLow030150_irrNonPaddy.map
    interflowIni      = landSurface.interflow_irrNonPaddy.map

    [groundwaterOptions]

    debugWaterBalance = True
    groundwaterPropertiesNC = groundwaterProperties5ArcMin_5min.nc

    # minimum value for groundwater recession coefficient (day-1)
    minRecessionCoeff = 1.0e-4

    limitFossilGroundWaterAbstraction = True
    minimumTotalGroundwaterThickness       = 0.000
    estimateOfTotalGroundwaterThickness    = thickness_05min_5min.map
    estimateOfRenewableGroundwaterCapacity = 0.0

    # annual pumping capacity for each region (unit: billion cubic meter per year), should be given in a netcdf file
    pumpingCapacityNC = regional_abstraction_limit_5min.nc


    # initial conditions:
    storGroundwaterIni = groundwater.storGroundwater.map
    storGroundwaterFossilIni = groundwater.storGroundwaterFossil.map
    #
    avgNonFossilGroundwaterAllocationLongIni  = groundwater.avgNonFossilAllocation.map
    avgNonFossilGroundwaterAllocationShortIni = groundwater.avgNonFossilAllocationShort.map
    avgTotalGroundwaterAbstractionIni         = groundwater.avgAbstraction.map        
    avgTotalGroundwaterAllocationLongIni      = groundwater.avgAllocation.map  
    avgTotalGroundwaterAllocationShortIni     = groundwater.avgAllocationShort.map   

    allocationSegmentsForGroundwater = uniqueIds30min.nom_5min.map
    #~ allocationSegmentsForGroundwater = None

    [routingOptions]
    debugWaterBalance = True

    lddMap      = lddsound_05min.map
    cellAreaMap = cellarea05min.map
    gradient    = ChannelGradient_05min.map

    # manning coefficient
    manningsN   = 0.04

    routingMethod = accuTravelTime
    # TODO: including kinematicWave
    #~ # Maximum length of a sub time step in seconds (optional and only used if either kinematicWave or simplifiedKinematicWave is used)
    #~ # - Note that too long sub time step may create water balance errors.
    #~ # - Default values: 3600 seconds for 30 arcmin ; 720 seconds for 5 arcmin
    #~ maxiumLengthOfSubTimeStep = 3600.
    #~ maxiumLengthOfSubTimeStep = 720.

    # dynamic flood plain options
    dynamicFloodPlain = False

    # lake and reservoir parameters
    waterBodyInputNC = waterBodies5ArcMin.nc
    onlyNaturalWaterBodies = False

    # composite crop factors for WaterBodies: 
    cropCoefficientWaterNC = cropCoefficientForOpenWater.nc
    minCropWaterKC         = 0.20

    # number of days (timesteps) that have been performed for spinning up initial conditions in the routing module (i.e. channelStorageIni, avgDischargeLongIni, avgDischargeShortIni, etc.)
    timestepsToAvgDischargeIni     = routing.timestepsToAvgDischarge.map
    # Note that: 
    # - maximum number of days (timesteps) to calculate long term average flow values (default: 5 years = 5 * 365 days = 1825)
    # - maximum number of days (timesteps) to calculate short term average values (default: 1 month = 1 * 30 days = 30)

    # initial conditions:
    waterBodyStorageIni	       = routing.waterBodyStorage.map
    avgLakeReservoirInflowShortIni = routing.avgInflow.map
    avgLakeReservoirOutflowLongIni = routing.avgOutflow.map
    channelStorageIni              = routing.channelStorage.map
    readAvlChannelStorageIni       = routing.readAvlChannelStorage.map
    avgDischargeLongIni            = routing.avgDischarge.map
    m2tDischargeLongIni            = routing.m2tDischarge.map
    avgBaseflowLongIni             = routing.avgBaseflow.map
    riverbedExchangeIni            = routing.riverbedExchange.map
    avgDischargeShortIni           = routing.avgDischargeShort.map
    subDischargeIni                = routing.subDischarge.map


An example model is available in \\wflow\\examples\\wflow_RhineMeuse_pcrglobwb\\. 

wflow_pcrglobwb module documentation
-------------------------------

.. automodule:: wflow_pcrglobwb
    :members:
    :undoc-members:
    :show-inheritance:
