.. _history:

*******
History
*******

1.3
===

1.2
===
- X-axis of the probability graph view now scales when the cutoff values of the attribute values are adjusted. (`#2907978 <http://sourceforge.net/tracker/?func=detail&aid=2907978&group_id=205121&atid=992360>`_)
- Time graph view now correctly shows probability values when appropriate. (`#2865535 <http://sourceforge.net/tracker/?func=detail&aid=2865535&group_id=205121&atid=992360>`_ and `#2865536 <http://sourceforge.net/tracker/?func=detail&aid=2865536&group_id=205121&atid=992360>`_)
- Added support for starting Aguila with the probability graph view. Also added support for configuring the probability graph view in the Xml startup configuration file. (`#2796370 <http://sourceforge.net/tracker/?func=detail&aid=2796370&group_id=205121&atid=992360>`_)
- Added support for visualizing :ref:`vector data <vectorDataFormats>` (attributes with a direction and a magnitude).
- Added support for directories when naming attributes that are located somewhere else than the current directory::

  $ aguila --timesteps [1,100] MyProject/MyData/Dem

- Improved panning when a 2D map is zoomed.
- Vector and Ldd attributes are not drawn when there are only a few pixels available per cell (large number of cells and/or zoomed out a lot). In that case a transparant cell is drawn to show the area of non-missing value cells. The vector and ldd directions are shown again when more pixels become available (when the map is zoomed into). This is done to keep Aguila responsive when large data sets are loaded. In such cases vector and ldd directions where not visible anyway.
- Re-added 'zoom by rectangle', by popular request. See section :ref:`mapControls`.
- For all data read using GDAL, the color interpretation is taken into account when determining the value scale of a raster attribute. For example, satellite imagery values stored in a byte geoTiff raster are treated as scalar data when the color interpretation stored in the raster is reported by GDal as being Gray. (`#2873963 <http://sourceforge.net/tracker/index.php?func=detail&aid=2873963&group_id=205121&atid=992360>`_)
- For all data read using GDAL, the PCRASTER_VALUESCALE meta data item is taken into account. If it is set, it overrides the color interpretation and type id that are otherwise used to determine the value scale, see section :ref:`rasterDataFormats`. (`#2873989 <http://sourceforge.net/tracker/?func=detail&aid=2873989&group_id=205121&atid=992360>`_)
- ...

1.1
===
- Added support for visualizing feature data.

1.0
===
First release in the 1-series. Since we are about to add some new features we first made an official release of the current version.
