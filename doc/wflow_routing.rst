THe wflow_routing Model
=======================


Introduction
------------
The wflow routing module uses the pcraster kinematic wave to route water over a DEM. By adding a bankflull level
and a floodplainwidth to the configuration the model can also include estimated flow over a floodplain.

Method
======

       """
        Updates the kinematic wave reservoir. Should be run after updates to Q

        WL = A * pow(Q,B)/Bw
        WL/A * Bw = pow(Q,B)
        Q = pow(WL/A * Bw,1/B)
        """
        # TODO: This is not correct, need to subtract Channel runoff
        self.Qbankfull = pow(self.bankFull/self.AlphaCh * self.Bw,1.0/self.Beta)
        self.Qchannel = min(self.SurfaceRunoff,self.Qbankfull)
        self.floodcells  = boolean(ifthenelse(self.WaterLevelCH > self.bankFull, boolean(1), boolean(0)))
        self.Qfloodplain = max(0.0,self.SurfaceRunoff - self.Qbankfull)

        self.WaterLevelCH = self.AlphaCh * pow(self.Qchannel, self.Beta) / (self.Bw)
        self.WaterLevelFP = self.AlphaFP * pow(self.Qfloodplain, self.Beta) / (self.Bw + self.Pfp)
        #self.WaterLevelFP = max(0.0,self.WaterLevelCH - self.bankFull)
        #self.Alpha * pow(self.SurfaceRunoff-self.Qchannel, self.Beta) / (self.Bw + self.Pch) + self.bankFull)
        self.Qtot = pow(self.WaterLevelCH/self.AlphaCh * self.Bw,1.0/self.Beta) + pow(self.WaterLevelFP/self.AlphaFP * self.Pfp,1.0/self.Beta)

        # wetted perimeter (m)
        self.Pch = self.wetPerimiterCH(self.WaterLevelCH,self.Bw)
        self.Pfp = self.wetPerimiterFP(self.WaterLevelFP,self.floodPlainWidth)
        # Alpha
        self.AlphaFP = self.AlpTerm * pow(self.Pfp, self.AlpPow)
        self.AlphaCh = self.AlpTerm * pow(self.Pch, self.AlpPow)
        self.Alpha = self.AlpTerm * pow(self.Pch + self.Pfp, self.AlpPow)
        self.OldKinWaveVolume = self.KinWaveVolume
        self.KinWaveVolume = self.WaterLevelCH * self.Bw * self.DCL


Dependencies
------------
T


Configuration
-------------


It needs a number of settings in the ini file. The default name for the file
is wflow\_routing.ini. it is also possible to insert this section in the
wflow\_sbm or wflow\_hbv ini file and point to that file.

See below for an example: 

::

    [inputmapstacks]
    # Name of the mapstack with specific discharge (mm/timestep output from the hydrological model)
    IW= inmaps/IW
    # of if you next this within say the wflow\_hv model:
    # IW = outmaps/IW





A description of the implementation of the kinematic wave is given on the pcraster website
 
In addition to the settings in the ini file you need to give the model additional maps
or lookuptables in the staticmaps or intbl directories:

Lookup tables
~~~~~~~~~~~~~

:N.tbl:
    Manning's N fro all no-river cells. Defaults to 0.072

:N_River.tbl:
    Manning's N for the river cells. Defaults to 0.036

:N_FloodPlain.tbl:
    Manning's N for the floodplain. A floodplain is always linked to a river cell. Defaults to 2* N of the river

As with all models the lookup tables can be replaced by a map with the same name (but with the .map extension) in the staticmaps directory.

staticmaps
~~~~~~~~~~

:wflow_subcatch.map:
    Map of the subcatchment in the area. Usually shared with the hydrological model

:wflow_dem.map:
    The digital elevation model. Usually shared with the hydrological model

:wflow_ldd.map:
    The D8 local drainage network.

:wflow_river.map:
    Definition of the river cells.

:wflow_riverlength.map:
    Optional map that defines the actual legth of the river in each cell.

:wflow_riverlength_fact.map:
    Optional map that defines a multiplication factor fro the river length in each cell.

:wflow_gauges.map:
    Map of river gauges that can be used in outputting timeseries

:wflow_inflow.map:
    Optional map of inflow points into the surface water. Limited testing.

:wflow_riverwidth.map:
    Optional map of the width of the river for each river cell.

:wflow_floodplainwidth.map:
    Optional map of the width of the floodplain for each river cell.

:wflow_bankfulldepth.map:
    Optional map of the level at which the river starts to flood and water will also be conducted over the floodplain.

:wflow_floodplaindist.map:
    Optional map that defines the relation between the water level in the floodplain

:wflow_landuse.map:
    Required map of landuse/land cover. This map is used in the lookup tables to relate parameters to landuse/landcover.
    Usually shared with the hydrological model

:wflow_soil.map:
    Required map of soil type. Usually shared with the hydrological model




wflow_routing module documentation
----------------------------------

.. automodule:: wflow_routing
    :members:
    :undoc-members:
    :show-inheritance:
