
import os
import os.path
import getopt
import ConfigParser
import sys
"""
rem rasterize needs to be done on tif file as .map files cannot be
rem handled directby by gdal
"""

shpfile = "subbasins_rhein_wgs1984.shp"
nname = os.path.splitext(os.path.basename(shpfile))[0]

# In this dictionary the fiedl in the dbf is linked to the filename in the .map
pars = {"BETA": "BetaSeepage",
        "CFMAX": "Cfmax",
        "ALPHA": "AlphaNL",
        "TTI" : "TTI",
        "TT" : "TT",
        "PERC" : "PERC",
        "K4" :"K4",
        "FC" : "FC",
        "KHQ" : "KHQ",
        "LP": "LP",
        "HQ" : "HQ",
        "CFR" : "CFR",
        "CEVPF" : "CEVPF"
}


os.system('pcrcalc "nilmap.map=scalar(if(scalar(cutout.map) >= 10.0,1.0))"')
for zz in pars:
    print pars[zz]
    os.system("gdal_translate -of GTiff nilmap.map " + pars[zz] + ".tif")
    os.system("gdal_rasterize -a " + zz + " -l " +  nname +  " "  + shpfile + " "  + pars[zz] + ".tif")
    os.system("gdal_translate -of PCRaster " + pars[zz] + ".tif " + pars[zz] + ".map")
    