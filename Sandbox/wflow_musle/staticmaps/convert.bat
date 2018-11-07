gdal_translate -of PCRaster -ot Float32 silt.asc wflow_silt.map
gdal_translate -of PCRaster -ot Float32 clay.asc wflow_clay.map
gdal_translate -of PCRaster -ot Float32 lai.asc wflow_lai.map
gdal_translate -of PCRaster -ot Float32 canopy_height.asc wflow_canopy_height.map
gdal_translate -of PCRaster -ot Float32 imp.asc wflow_impervious.map
gdal_translate -of PCRaster -ot Float32 idplt.asc wflow_idplt.map
gdal_translate -of PCRaster -ot Float32 sol_cov.asc wflow_sol_cov.map
gdal_translate -of PCRaster -ot Float32 dem.asc wflow_dem.map

pause