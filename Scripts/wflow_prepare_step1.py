#!python
"""
wflow_prepare_step1
===================

wflow data preparation script. Data preparation can be done by hand or using 
the two scripts. This script does the first step. The second script does 
the resampling. This scripts need the pcraster and gdal executables to be
available in you search path.


Usage::

    wflow_prepare_step1 [-W workdir][-f][-h] -I inifile 
    
    -f force recreation of ldd if it already exists
    -h show this information
    -W set the working directory, default is current dir
    -I name of the ini file with settings
            
$Id: $  


.. todo:
    # Add a legend to the maps
    # possibility to give a file with gauge names (legend)
    # add shape files and qgis config for easy viewing
    
.. todo:
    # calibration script
    
.. todo:
    # plotting functions in python
    # reporting in python
"""

import wflow.wflow_lib as tr


import os
import os.path
import getopt
import configparser
import sys
import gc
import numpy as np
import pcraster as pcr


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


def configget(config, section, var, default):
    """
    gets parameter from config file and returns a default value
    if the parameter is not found
    """
    try:
        ret = config.get(section, var)
    except:
        print("returning default (" + default + ") for " + section + ":" + var)
        ret = default

    return ret


def OpenConf(fn):
    config = configparser.ConfigParser()
    config.optionxform = str

    if os.path.exists(fn):
        config.read(fn)
    else:
        print("Cannot open config file: " + fn)
        sys.exit(1)

    return config


def main():
    """
        
    :ivar masterdem: digital elevation model
    :ivar dem: digital elevation model
    :ivar river: optional river map
    """

    # Default values
    strRiver = 8
    masterdem = "dem.map"
    step1dir = "step1"
    step2dir = "step2"
    workdir = "."
    inifile = "wflow_prepare.ini"
    recreate = False
    snapgaugestoriver = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "W:hI:f",['version'])
    except getopt.error as msg:
        usage(msg)

    for o, a in opts:
        if o == "-W":
            workdir = a
        if o == "-I":
            inifile = a
        if o == "-h":
            usage()
        if o == "-f":
            recreate = True
        if o == "--version":
            import wflow
            print("wflow version: ", wflow.__version__)
            sys.exit(0)

    pcr.setglobaloption("unitcell")
    os.chdir(workdir)

    config = OpenConf(workdir + "/" + inifile)

    masterdem = configget(config, "files", "masterdem", "dem.map")
    pcr.setclone(masterdem)

    strRiver = int(configget(config, "settings", "riverorder", "4"))

    try:
        gauges_x = config.get("settings", "gauges_x")
        gauges_y = config.get("settings", "gauges_y")
    except:
        print("gauges_x and  gauges_y are required entries in the ini file")
        sys.exit(1)

    step1dir = configget(config, "directories", "step1dir", "step1")
    step2dir = configget(config, "directories", "step2dir", "step2")
    # upscalefactor = float(config.get("settings","upscalefactor"))

    corevolume = float(configget(config, "settings", "corevolume", "1E35"))
    catchmentprecipitation = float(
        configget(config, "settings", "catchmentprecipitation", "1E35")
    )
    corearea = float(configget(config, "settings", "corearea", "1E35"))
    outflowdepth = float(configget(config, "settings", "lddoutflowdepth", "1E35"))

    initialscale = int(configget(config, "settings", "initialscale", "1"))
    csize = float(configget(config, "settings", "cellsize", "1"))

    snapgaugestoriver = bool(
        int(configget(config, "settings", "snapgaugestoriver", "1"))
    )
    lddglobaloption = configget(config, "settings", "lddglobaloption", "lddout")
    pcr.setglobaloption(lddglobaloption)
    lu_water = configget(config, "files", "lu_water", "")
    lu_paved = configget(config, "files", "lu_paved", "")

    # X/Y coordinates of the gauges the system    
    X = np.fromstring(gauges_x, sep=',')
    Y = np.fromstring(gauges_y, sep=',')

    tr.Verbose = 1

    # make the directories to save results in
    if not os.path.isdir(step1dir + "/"):
        os.makedirs(step1dir)
    if not os.path.isdir(step2dir):
        os.makedirs(step2dir)

    if initialscale > 1:
        print("Initial scaling of DEM...")
        os.system(
            "resample -r "
            + str(initialscale)
            + " "
            + masterdem
            + " "
            + step1dir
            + "/dem_scaled.map"
        )
        print("Reading dem...")
        dem = pcr.readmap(step1dir + "/dem_scaled.map")
        ldddem = dem
    else:
        print ("Reading dem...")
        dem = pcr.readmap(masterdem)
        ldddem = dem

    try:
        catchmask = config.get("files", "catchment_mask")
    except:
        print("No catchment mask...")
    else:
        print("clipping DEM with mask.....")
        mask = pcr.readmap(catchmask)
        ldddem = pcr.ifthen(pcr.boolean(mask), ldddem)
        dem = pcr.ifthen(pcr.boolean(mask), dem)

    # See if there is a shape file of the river to burn in
    try:
        rivshp = config.get("files", "river")
    except:
        print("no river file specified")
        outletpointX = float(configget(config, "settings", "outflowpointX", "0.0"))
        outletpointY = float(configget(config, "settings", "outflowpointY", "0.0"))
    else:
        print("river file specified.....")
        try:
            outletpointX = float(configget(config, "settings", "outflowpointX", "0.0"))
            outletpointY = float(configget(config, "settings", "outflowpointY", "0.0"))
        except:
            print(
                "Need to specify the river outletpoint (a point at the end of the river within the current map)"
            )
            exit(1)

        outletpointmap = tr.points_to_map(dem, outletpointX, outletpointY, 0.5)
        pcr.report(outletpointmap, step1dir + "/outletpoint.map")
        # rivshpattr = config.get("files","riverattr")
        pcr.report(dem * 0.0, step1dir + "/nilmap.map")
        thestr = (
            "gdal_translate -of GTiff "
            + step1dir
            + "/nilmap.map "
            + step1dir
            + "/riverburn.tif"
        )
        os.system(thestr)
        rivshpattr = os.path.splitext(os.path.basename(rivshp))[0]
        os.system(
            "gdal_rasterize -burn 1 -l "
            + rivshpattr
            + " "
            + rivshp
            + " "
            + step1dir
            + "/riverburn.tif"
        )
        thestr = (
            "gdal_translate -of PCRaster "
            + step1dir
            + "/riverburn.tif "
            + step1dir
            + "/riverburn.map"
        )
        os.system(thestr)
        riverburn = pcr.readmap(step1dir + "/riverburn.map")
        # Determine regional slope assuming that is the way the river should run
        # Determine regional slope assuming that is the way the river should run
        # pcr.setglobaloption("unitcell")
        # demregional=pcr.windowaverage(dem,100)
        ldddem = pcr.ifthenelse(riverburn >= 1.0, dem - 1000, dem)

    pcr.setglobaloption("unittrue")
    upscalefactor = int(csize / pcr.celllength())

    print("Creating ldd...")
    ldd = tr.lddcreate_save(
        step1dir + "/ldd.map",
        ldddem,
        recreate,
        outflowdepth=outflowdepth,
        corevolume=corevolume,
        catchmentprecipitation=catchmentprecipitation,
        corearea=corearea,
    )

    print("Determining streamorder...")
    stro = pcr.streamorder(ldd)
    pcr.report(stro, step1dir + "/streamorder.map")
    strdir = pcr.ifthen(stro >= strRiver, stro)
    pcr.report(strdir, step1dir + "/streamorderrive.map")
    pcr.report(pcr.boolean(pcr.ifthen(stro >= strRiver, stro)), step1dir + "/rivers.map")

    pcr.setglobaloption("unittrue")
    # outlet (and other gauges if given)
    # TODO: check is x/y set if not skip this
    print("Outlet...")

    outlmap = tr.points_to_map(dem, X, Y, 0.5)

    if snapgaugestoriver:
        print("Snapping gauges to nearest river cells...")
        pcr.report(outlmap, step1dir + "/orggauges.map")
        outlmap = tr.snaptomap(outlmap, strdir)

    # noutletmap = tr.points_to_map(dem,XX,YY,0.5)
    # pcr.report(noutletmap,'noutlet.map')

    pcr.report(outlmap, step1dir + "/gauges.map")

    # check if there is a pre-define catchment map
    try:
        catchmask = config.get("files", "catchment_mask")
    except:
        print("No catchment mask, finding outlet")
        # Find catchment (overall)
        outlet = tr.find_outlet(ldd)
        sub = tr.subcatch(ldd, outlet)
        pcr.report(sub, step1dir + "/catchment_overall.map")
    else:
        print("reading and converting catchment mask.....")
        os.system(
            "resample -r "
            + str(initialscale)
            + " "
            + catchmask
            + " "
            + step1dir
            + "/catchment_overall.map"
        )
        sub = pcr.readmap(step1dir + "/catchment_overall.map")

    print("Scatch...")
    sd = tr.subcatch(ldd, pcr.ifthen(outlmap > 0, outlmap))
    pcr.report(sd, step1dir + "/scatch.map")

    pcr.setglobaloption("unitcell")
    print("Upscalefactor: " + str(upscalefactor))

    if upscalefactor > 1:
        gc.collect()
        print("upscale river length1 (checkerboard map)...")
        ck = tr.checkerboard(dem, upscalefactor)
        pcr.report(ck, step1dir + "/ck.map")
        pcr.report(dem, step1dir + "/demck.map")
        print("upscale river length2...")
        fact = tr.area_riverlength_factor(ldd, ck, upscalefactor)
        pcr.report(fact, step1dir + "/riverlength_fact.map")

        # print("make dem statistics...")
        dem_ = pcr.areaaverage(dem, ck)
        pcr.report(dem_, step1dir + "/demavg.map")

        print("Create DEM statistics...")
        dem_ = pcr.areaminimum(dem, ck)
        pcr.report(dem_, step1dir + "/demmin.map")
        dem_ = pcr.areamaximum(dem, ck)
        pcr.report(dem_, step1dir + "/demmax.map")
        # calculate percentiles
        order = pcr.areaorder(dem, ck)
        n = pcr.areatotal(pcr.spatial(pcr.scalar(1.0)), ck)
        #: calculate 25 percentile
        perc = tr.area_percentile(dem, ck, n, order, 25.0)
        pcr.report(perc, step1dir + "/dem25.map")
        perc = tr.area_percentile(dem, ck, n, order, 10.0)
        pcr.report(perc, step1dir + "/dem10.map")
        perc = tr.area_percentile(dem, ck, n, order, 50.0)
        pcr.report(perc, step1dir + "/dem50.map")
        perc = tr.area_percentile(dem, ck, n, order, 33.0)
        pcr.report(perc, step1dir + "/dem33.map")
        perc = tr.area_percentile(dem, ck, n, order, 66.0)
        pcr.report(perc, step1dir + "/dem66.map")
        perc = tr.area_percentile(dem, ck, n, order, 75.0)
        pcr.report(perc, step1dir + "/dem75.map")
        perc = tr.area_percentile(dem, ck, n, order, 90.0)
        pcr.report(perc, step1dir + "/dem90.map")
    else:
        print("No fancy scaling done. Going strait to step2....")
        pcr.report(dem, step1dir + "/demavg.map")
        Xul = float(config.get("settings", "Xul"))
        Yul = float(config.get("settings", "Yul"))
        Xlr = float(config.get("settings", "Xlr"))
        Ylr = float(config.get("settings", "Ylr"))
        gdalstr = (
            "gdal_translate  -projwin "
            + str(Xul)
            + " "
            + str(Yul)
            + " "
            + str(Xlr)
            + " "
            + str(Ylr)
            + " -of PCRaster  "
        )
        # gdalstr = "gdal_translate  -a_ullr " + str(Xul) + " " + str(Yul) + " " +str(Xlr) + " " +str(Ylr) + " -of PCRaster  "
        print(gdalstr)
        pcr.report(pcr.cover(1.0), step1dir + "/wflow_riverlength_fact.map")
        # Now us gdat tp convert the maps
        os.system(
            gdalstr
            + step1dir
            + "/wflow_riverlength_fact.map"
            + " "
            + step2dir
            + "/wflow_riverlength_fact.map"
        )
        os.system(
            gdalstr + step1dir + "/demavg.map" + " " + step2dir + "/wflow_dem.map"
        )
        os.system(
            gdalstr + step1dir + "/demavg.map" + " " + step2dir + "/wflow_demmin.map"
        )
        os.system(
            gdalstr + step1dir + "/demavg.map" + " " + step2dir + "/wflow_demmax.map"
        )
        os.system(
            gdalstr + step1dir + "/gauges.map" + " " + step2dir + "/wflow_gauges.map"
        )
        os.system(
            gdalstr + step1dir + "/rivers.map" + " " + step2dir + "/wflow_river.map"
        )
        os.system(
            gdalstr
            + step1dir
            + "/streamorder.map"
            + " "
            + step2dir
            + "/wflow_streamorder.map"
        )
        os.system(
            gdalstr + step1dir + "/gauges.map" + " " + step2dir + "/wflow_outlet.map"
        )
        os.system(
            gdalstr + step1dir + "/scatch.map" + " " + step2dir + "/wflow_catchment.map"
        )
        os.system(gdalstr + step1dir + "/ldd.map" + " " + step2dir + "/wflow_ldd.map")
        os.system(
            gdalstr + step1dir + "/scatch.map" + " " + step2dir + "/wflow_subcatch.map"
        )

        if lu_water:
            os.system(gdalstr + lu_water + " " + step2dir + "/WaterFrac.map")

        if lu_paved:
            os.system(gdalstr + lu_paved + " " + step2dir + "/PathFrac.map")

        try:
            lumap = config.get("files", "landuse")
        except:
            print("no landuse map...creating uniform map")
            # clone=pcr.readmap(step2dir + "/wflow_dem.map")
            pcr.setclone(step2dir + "/wflow_dem.map")
            pcr.report(pcr.nominal(1), step2dir + "/wflow_landuse.map")
        else:
            os.system(
                "resample --clone "
                + step2dir
                + "/wflow_dem.map "
                + lumap
                + " "
                + step2dir
                + "/wflow_landuse.map"
            )

        try:
            soilmap = config.get("files", "soil")
        except:
            print("no soil map..., creating uniform map")
            pcr.setclone(step2dir + "/wflow_dem.map")
            pcr.report(pcr.nominal(1), step2dir + "/wflow_soil.map")
        else:
            os.system(
                "resample --clone "
                + step2dir
                + "/wflow_dem.map "
                + soilmap
                + " "
                + step2dir
                + "/wflow_soil.map"
            )


if __name__ == "__main__":
    main()
