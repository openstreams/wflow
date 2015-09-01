hillslope_orien4.map = plateauArab.map + plateauPast.map + hillslope.map (same as HPPPA4 Hillslope, Plateau with Pasture, Plateau with Arable land)
perc_plateau_orien.map = grid_maps\perc_plateauArab_orien.map + grid_maps\perc_plateauPast_orien.map

wflow_velocity made by
pcrcalc wflow_velocity.map = if(wflow_dem.map ge -100, scalar(3600))