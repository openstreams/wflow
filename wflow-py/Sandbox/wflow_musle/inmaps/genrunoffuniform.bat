gdal_translate -of PCRaster -ot Float32 clone.asc clone.map
rem pcrcalc RUN00000.001 = clone.map  + 10
rem pcrcalc RUN00000.002 = clone.map  + 10
rem pcrcalc RUN00000.003 = clone.map  + 40
rem pcrcalc RUN00000.004 = clone.map  + 75
rem pcrcalc RUN00000.005 = clone.map  + 25
rem pcrcalc RUN00000.006 = clone.map  + 200
rem pcrcalc RUN00000.007 = clone.map  + 125
rem pcrcalc RUN00000.008 = clone.map  + 75
rem pcrcalc RUN00000.009 = clone.map  + 25
rem pcrcalc RUN00000.010 = clone.map  + 0

pcrcalc RUN00000.001 = clone.map  + 1
pcrcalc RUN00000.002 = clone.map  + 1
pcrcalc RUN00000.003 = clone.map  + 1
pcrcalc RUN00000.004 = clone.map  + 1
pcrcalc RUN00000.005 = clone.map  + 1
pcrcalc RUN00000.006 = clone.map  + 1
pcrcalc RUN00000.007 = clone.map  + 1
pcrcalc RUN00000.008 = clone.map  + 1
pcrcalc RUN00000.009 = clone.map  + 1
pcrcalc RUN00000.010 = clone.map  + 1

pause