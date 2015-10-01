THe wflow_routing Model
=======================


.. warning::

    The documentation of this model is incomplete

Introduction
------------
The wflow routing module uses the pcraster kinematic wave to route water over a DEM. By adding a bankflull level
and a floodplainwidth to the configuration the model can also inclue estimate flow over a floodplain.

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





A description of the implementation of the kinematic wave is given on the pcraster website
 
In addition to the settings in the ini file you need to give the model additional maps
or lookuptables in the staticmaps or intbl directories:






wflow_routing module documentation
----------------------------------

.. automodule:: wflow_routing
    :members:
    :undoc-members:
    :show-inheritance:
