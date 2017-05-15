# variation on the runoff demo (runoff2.mod)

# that will create some data to show Aguila features

# model for simulation of runoff
# 24 timesteps of 6 hours => modelling time one week

binding
 RainStations = rainstat.map;          # map with location of rainstations
 RainTimeSeries = rain.tss;            # timeseries with rain at rainstations
 RainZones = rainzone.map;             # reported stack of maps with rain
 SurfaceWater = rainfall;              # reported maps with rain (mm/6hours)
 SoilInfiltrationTable = infilcap.tbl; # table with infiltr. cap. of soil types
 SoilType = soil.map;                  # soil map
 InfiltrationCapacity = infilcap.map;  # reported map with infiltr. cap.
 Dem = dem.map;                        # digital elevation map
 Ldd = ldd.map;                        # reported local drain direction map
 UpstreamArea = upstreamArea.map;      # report the upstream area
 ConvConst = 216000;                   # conversion mm/6hours => m3/s
 RunOff = runoff;                      # reported stack of maps with

 SamplePlaces = samples.map;           # map with runoff sampling locations
 RunoffTimeSeries = runoff.tss;        # reported timeseries with runoff 
                                       # at sampling locations

timer
 1 28 1;

initial
 # coverage of meteorological stations for the whole area
 RainZones = spreadzone(RainStations, 0, 1);
 # create an infiltration capacity map (mm/6 hours), based on the soil map
 InfiltrationCapacity = lookupscalar(SoilInfiltrationTable, SoilType);
 # generate the local drain direction map on basis of the elevation map
 report Ldd = lddcreate(Dem, 1e31, 1e31, 1e31, 1e31);

 report UpstreamArea = accuflux(Ldd, cellarea());

dynamic
 # calculate and report maps with rainfall at each timestep (mm/6 hours)
 SurfaceWater = timeinputscalar(RainTimeSeries, RainZones);
 # compute both runoff and actual infiltration
 RunoffPerTimestep,Infiltration = accuthresholdflux, accuthresholdstate(
         Ldd,SurfaceWater,InfiltrationCapacity);
 # output runoff, converted to m3/s, at each timestep
 report RunOff = RunoffPerTimestep / ConvConst;
 # output runoff (converted to m3/s) at each timestep for selected locations
 report RunoffTimeSeries = timeoutput(SamplePlaces, RunOff);
