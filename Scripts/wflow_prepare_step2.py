#!python
"""
wflow_prepare_step2
===================

wflow data preparation script. Data preparation can be done by hand or using 
the two scripts. This script does the resampling. This scripts need the pcraster and gdal executables to be
available in you search path.


Usage::

    wflow_prepare_step2 [-W workdir][-f][-h] -I inifile 

::    
    
    -f force recreation of ldd if it already exists
    -h show this information
    -W set the working directory, default is current dir
    -I name of the ini file with settings


$Id: $
"""


import wflow.wflow_lib as tr

import os
import os.path
import getopt
import configparser
import sys
import numpy as np
import pcraster as pcr

tr.Verbose = 1


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


def configget(config, section, var, default):
    """
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


def resamplemaps(step1dir, step2dir):
    """
    Resample the maps from step1 and rename them in the process
    """
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/dem10.map "
        + step2dir
        + "/wflow_dem10.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/dem25.map "
        + step2dir
        + "/wflow_dem25.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/dem33.map "
        + step2dir
        + "/wflow_dem33.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/dem50.map "
        + step2dir
        + "/wflow_dem50.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/dem66.map "
        + step2dir
        + "/wflow_dem66.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/dem75.map "
        + step2dir
        + "/wflow_dem75.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/dem90.map "
        + step2dir
        + "/wflow_dem90.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/demavg.map "
        + step2dir
        + "/wflow_dem.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/demmin.map "
        + step2dir
        + "/wflow_demmin.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/demmax.map "
        + step2dir
        + "/wflow_demmax.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/riverlength_fact.map "
        + step2dir
        + "/wflow_riverlength_fact.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/catchment_overall.map "
        + step2dir
        + "/catchment_cut.map"
    )
    os.system(
        "resample --clone "
        + step2dir
        + "/cutout.map "
        + step1dir
        + "/rivers.map "
        + step2dir
        + "/wflow_riverburnin.map"
    )


def main():
    """

    """
    workdir = "."
    inifile = "wflow_prepare.ini"

    try:
        opts, args = getopt.getopt(sys.argv[1:], "W:hI:f", ['version'])
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

    os.chdir(workdir)

    config = OpenConf(workdir + "/" + inifile)

    step1dir = configget(config, "directories", "step1dir", "step1")
    step2dir = configget(config, "directories", "step2dir", "step2")
    snapgaugestoriver = bool(
        int(configget(config, "settings", "snapgaugestoriver", "1"))
    )

    # make the directories to save results in
    if not os.path.isdir(step1dir + "/"):
        os.makedirs(step1dir)
    if not os.path.isdir(step2dir):
        os.makedirs(step2dir)

    ##first make the clone map
    try:
        Xul = float(config.get("settings", "Xul"))
        Yul = float(config.get("settings", "Yul"))
        Xlr = float(config.get("settings", "Xlr"))
        Ylr = float(config.get("settings", "Ylr"))
    except:
        print("Xul, Xul, Xlr and  Ylr are required entries in the ini file")
        sys.exit(1)

    csize = float(configget(config, "settings", "cellsize", "1"))
    try:
        gauges_x = config.get("settings", "gauges_x")
        gauges_y = config.get("settings", "gauges_y")
    except:
        print("gauges_x and  gauges_y are required entries in the ini file")
        sys.exit(1)

    strRiver = int(configget(config, "settings", "riverorder_step2", "4"))

    corevolume = float(configget(config, "settings", "corevolume", "1E35"))
    catchmentprecipitation = float(
        configget(config, "settings", "catchmentprecipitation", "1E35")
    )
    corearea = float(configget(config, "settings", "corearea", "1E35"))
    outflowdepth = float(configget(config, "settings", "lddoutflowdepth", "1E35"))
    lddmethod = configget(config, "settings", "lddmethod", "dem")
    lddglobaloption = configget(config, "settings", "lddglobaloption", "lddout")
    pcr.setglobaloption(lddglobaloption)

    nrrow = round(abs(Yul - Ylr) / csize)
    nrcol = round(abs(Xlr - Xul) / csize)
    mapstr = (
        "mapattr -s -S -R "
        + str(nrrow)
        + " -C "
        + str(nrcol)
        + " -l "
        + str(csize)
        + " -x "
        + str(Xul)
        + " -y "
        + str(Yul)
        + " -P yb2t "
        + step2dir
        + "/cutout.map"
    )

    os.system(mapstr)
    pcr.setclone(step2dir + "/cutout.map")

    lu_water = configget(config, "files", "lu_water", "")
    lu_paved = configget(config, "files", "lu_paved", "")

    if lu_water:
        os.system(
            "resample --clone "
            + step2dir
            + "/cutout.map "
            + lu_water
            + " "
            + step2dir
            + "/wflow_waterfrac.map"
        )

    if lu_paved:
        os.system(
            "resample --clone "
            + step2dir
            + "/cutout.map "
            + lu_paved
            + " "
            + step2dir
            + "/PathFrac.map"
        )

    #
    try:
        lumap = config.get("files", "landuse")
    except:
        print("no landuse map...creating uniform map")
        clone = pcr.readmap(step2dir + "/cutout.map")
        pcr.report(pcr.nominal(clone), step2dir + "/wflow_landuse.map")
    else:
        os.system(
            "resample --clone "
            + step2dir
            + "/cutout.map "
            + lumap
            + " "
            + step2dir
            + "/wflow_landuse.map"
        )

    try:
        soilmap = config.get("files", "soil")
    except:
        print("no soil map..., creating uniform map")
        clone = pcr.readmap(step2dir + "/cutout.map")
        pcr.report(pcr.nominal(clone), step2dir + "/wflow_soil.map")
    else:
        os.system(
            "resample --clone "
            + step2dir
            + "/cutout.map "
            + soilmap
            + " "
            + step2dir
            + "/wflow_soil.map"
        )

    resamplemaps(step1dir, step2dir)

    dem = pcr.readmap(step2dir + "/wflow_dem.map")
    demmin = pcr.readmap(step2dir + "/wflow_demmin.map")
    demmax = pcr.readmap(step2dir + "/wflow_demmax.map")
    # catchcut = pcr.readmap(step2dir + "/catchment_cut.map")
    catchcut = pcr.readmap(step2dir + "/cutout.map")
    # now apply the area of interest (catchcut) to the DEM
    # dem=pcr.ifthen(catchcut >=1 , dem)
    #

    # See if there is a shape file of the river to burn in
    try:
        rivshp = config.get("files", "river")
    except:
        print("no river file specified")
        riverburn = pcr.readmap(step2dir + "/wflow_riverburnin.map")
    else:
        print("river file speficied.....")
        # rivshpattr = config.get("files","riverattr")
        pcr.report(dem * 0.0, step2dir + "/nilmap.map")
        thestr = (
            "gdal_translate -of GTiff "
            + step2dir
            + "/nilmap.map "
            + step2dir
            + "/wflow_riverburnin.tif"
        )
        os.system(thestr)
        rivshpattr = os.path.splitext(os.path.basename(rivshp))[0]
        os.system(
            "gdal_rasterize -burn 1 -l "
            + rivshpattr
            + " "
            + rivshp
            + " "
            + step2dir
            + "/wflow_riverburnin.tif"
        )
        thestr = (
            "gdal_translate -of PCRaster "
            + step2dir
            + "/wflow_riverburnin.tif "
            + step2dir
            + "/wflow_riverburnin.map"
        )
        os.system(thestr)
        riverburn = pcr.readmap(step2dir + "/wflow_riverburnin.map")
        # ldddem = pcr.ifthenelse(riverburn >= 1.0, dem -1000 , dem)

    # Only burn within the original catchment
    riverburn = pcr.ifthen(pcr.scalar(catchcut) >= 1, riverburn)
    # Now setup a very high wall around the catchment that is scale
    # based on the distance to the catchment so that it slopes away from the
    # catchment
    if lddmethod != "river":
        print("Burning in highres-river ...")
        disttocatch = pcr.spread(pcr.nominal(catchcut), 0.0, 1.0)
        demmax = pcr.ifthenelse(
            pcr.scalar(catchcut) >= 1.0,
            demmax,
            demmax + (pcr.celllength() * 100.0) / disttocatch,
        )
        pcr.setglobaloption("unitcell")
        # demregional=pcr.windowaverage(demmin,100)
        demburn = pcr.cover(pcr.ifthen(pcr.boolean(riverburn), demmin - 100.0), demmax)
    else:
        print("using average dem..")
        demburn = dem

    ldd = tr.lddcreate_save(
        step2dir + "/wflow_ldd.map",
        demburn,
        True,
        outflowdepth=outflowdepth,
        corevolume=corevolume,
        catchmentprecipitation=catchmentprecipitation,
        corearea=corearea,
    )

    # Find catchment (overall)
    outlet = tr.find_outlet(ldd)
    sub = tr.subcatch(ldd, outlet)
    pcr.report(sub, step2dir + "/wflow_catchment.map")
    pcr.report(outlet, step2dir + "/wflow_outlet.map")

    # make river map
    strorder = pcr.streamorder(ldd)
    pcr.report(strorder, step2dir + "/wflow_streamorder.map")

    river = pcr.ifthen(pcr.boolean(strorder >= strRiver), strorder)
    pcr.report(river, step2dir + "/wflow_river.map")

    # make subcatchments
    # os.system("col2map --clone " + step2dir + "/cutout.map gauges.col " + step2dir + "/wflow_gauges.map")
    X = np.fromstring(gauges_x, sep=',')
    Y = np.fromstring(gauges_y, sep=',')

    pcr.setglobaloption("unittrue")

    outlmap = tr.points_to_map(dem, X, Y, 0.5)
    pcr.report(outlmap, step2dir + "/wflow_gauges_.map")

    if snapgaugestoriver:
        print("Snapping gauges to river")
        pcr.report(outlmap, step2dir + "/wflow_orggauges.map")
        outlmap = tr.snaptomap(outlmap, river)

    outlmap = pcr.ifthen(outlmap > 0, outlmap)
    pcr.report(outlmap, step2dir + "/wflow_gauges.map")

    scatch = tr.subcatch(ldd, outlmap)
    pcr.report(scatch, step2dir + "/wflow_subcatch.map")


if __name__ == "__main__":
    main()
