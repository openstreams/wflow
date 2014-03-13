rem rasterize needs to be done on tif file as .map files cannot be
rem handled directby by gdal

#WHC.tbl
#TTI.tbl
#TT.tbl
#PERC.tbl
#LP.tbl
#KHQ.tbl
#K4.tbl
#ICF.tbl (was icfo and icfi need to fix...)
#HQ.tbl
#FoCfmax.tbl (only linked to LLU ().6 for forest), kan be table)
#FC.tbl
#CFR.tbl
#Cfmax.tbl
#CEVPF.tbl (is dat CEVPFO???)
#BetaSeepage.tbl
#AlphaNL.tbl

gdal_rasterize -a BETA -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp BetaSeepage.tif
gdal_rasterize -a CFMAX -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp Cfmax.tif
gdal_rasterize -a ALPHA -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp AlphaNL.tif
gdal_rasterize -a WHC -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp WHC.tif
gdal_rasterize -a TTI -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp TTI.tif
gdal_rasterize -a TT -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp TT.tif
gdal_rasterize -a PERC -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp PERC.tif
gdal_rasterize -a K4 -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp K4.tif
gdal_rasterize -a FC -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp FC.tif
gdal_rasterize -a KHQ -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp KHQ.tif
gdal_rasterize -a LP -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp LP.tif
gdal_rasterize -a HQ -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp HQ.tif
gdal_rasterize -a CFR -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp CFR.tif
gdal_rasterize -a CEVPFO -l subbasins_rhein_wgs1984 subbasins_rhein_wgs1984.shp CEVPF.tif