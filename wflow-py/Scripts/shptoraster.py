
import os
import os.path
import getopt
import configparser
import sys

"""
This scripts converst the hbv shape file to raster maps for the
distributed version of HBV.

rem rasterize needs to be done on tif file as .map files cannot be
rem handled directby by gdal
"""

shpfile = "subbasins_rhein_wgs1984.shp"
nname = os.path.splitext(os.path.basename(shpfile))[0]

# In this dictionary the fiedl in the dbf is linked to the filename in the .map
pars = {
    "BETA": "BetaSeepage",
    "CFMAX": "Cfmax_org",
    "ALFA": "AlphaNL",
    "FOCFMAX": "FOCFMAX",
    "TTI": "TTI",
    "WHC": "WHC",
    "TT": "TT",
    "PERC": "PERC",
    "K4": "K4",
    "FC": "FC",
    "KHQ": "KHQ",
    "LP": "LP",
    "HQ": "HQ",
    "CFR": "CFR",
    "CEVPFO": "CEVPFO",
    "EPF": "EPF",
    "RFCF": "RFCF",
    "CFLUX": "Cflux",
    "SFCF": "SFCF",
    "ICFI": "ICFI",
    "PCORR_": "Pcorr",
    "ECORR": "ECORR",
    "ICFO": "ICFO",
}

os.system('pcrcalc "nilmap.map=scalar(if(cutout.map>99,10))"')
os.system("gdal_translate -of GTiff nilmap.map " + " subcatch.tif")
os.system("gdal_rasterize -a ID -l " + nname + " " + shpfile + " subcatch.tif")
os.system("gdal_translate -of PCRaster  subcatch.tif subcatch.map")

for zz in pars:
    print(pars[zz])
    os.system("gdal_translate -of GTiff nilmap.map " + pars[zz] + ".tif")
    os.system(
        "gdal_rasterize -a "
        + zz
        + " -l "
        + nname
        + " "
        + shpfile
        + " "
        + pars[zz]
        + ".tif"
    )
    os.system("gdal_translate -of PCRaster " + pars[zz] + ".tif " + pars[zz] + ".map")

# now some specifics
forest = 3

os.system('pcrcalc "CEVPF.map=if(wflow_landuse.map == 3, CEVPFO.map, 1.0)"')
os.system('pcrcalc "ICF.map=if(wflow_landuse.map == 3, ICFO.map, ICFI.map)"')
os.system(
    'pcrcalc "Cfmax.map=if(wflow_landuse.map == 3, Cfmax_org.map * FOCFMAX.map, Cfmax_org.map)"'
)
os.system('pcrcalc "ECORR.map=ECORR.map * 10.0"')
