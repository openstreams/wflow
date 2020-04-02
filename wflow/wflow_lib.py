# Copyright (c) J. Schellekens 2005-2011
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
wflow_lib - terrain analysis and hydrological library
-----------------------------------------------------

The goal of this module is to make a series functions to upscale maps (DEM)
and to  maintain as much of the information in a detailled dem when upscaling
to a coarser DEM. These include:

    - river length (per cell)
    - river network location
    - elevation distribution
    - other terrain analysis

the wflow_prepare scripts use this library extensively.

$Author: schelle $
$Id: wflow_lib.py 808 2013-10-04 19:42:43Z schelle $
$Rev: 808 $
"""

import gzip
import os.path
import sys
import zipfile

import numpy as np
import osgeo.gdal as gdal
import pcraster as pcr
import pcraster.framework


def pt_flow_in_river(ldd, river):
    """
    Returns all points (True) that flow into the mak river (boolean map with river set to True)

    :param ldd: Drainage network
    :param river: Map of river (True River, False non-river)
    :return ifmap: map with infrlo points into the river (True)
    :return ctach: catchment of each of the inflow points
    """

    dspts = pcr.downstream(ldd, pcr.cover(river, 0))
    dspts = pcr.ifthenelse(pcr.cover(river, 0) == 1, 0, dspts)

    catch = pcr.subcatchment(ldd, pcr.nominal(pcr.uniqueid(dspts)))

    return dspts, catch


def sum_list_cover(list_of_maps, covermap):
    """
    Sums a list of pcrastermap using cover to fill in missing values

    :param list_of_maps: list of maps to sum
    :param covermap: maps/ value to use fro cover

    :return: sum of list of maps (single map)
    """
    sum_ = pcr.cover(0.0)
    for map in list_of_maps:
        sum_ = sum_ + pcr.cover(map, covermap)

    return sum_


def idtoid(sourceidmap, targetidmap, valuemap):
    """
    tranfer the values from valuemap at the point id's in sourceidmap to the areas in targetidmap.

    :param pointmap:
    :param areamap:
    :param valuemap:
    :return:
    """

    _area = pcr.pcr2numpy(targetidmap, 0.0).copy().astype(float)
    _pt = pcr.pcr2numpy(sourceidmap, 0.0).copy()
    _val = pcr.pcr2numpy(valuemap, 0.0).copy()

    for val in np.unique(_pt):
        if val > 0:  #
            _area[_area == val] = np.mean(_val[_pt == val])

    retmap = pcr.numpy2pcr(pcr.Scalar, _area, 0.0)

    return retmap


Verbose = 0


def lddcreate_save(
    lddname,
    dem,
    force,
    corevolume=1e35,
    catchmentprecipitation=1e35,
    corearea=1e35,
    outflowdepth=1e35,
):
    """
    Creates an ldd if a file does not exists or if the force flag is used

    input:
        - lddname (name of the ldd to create)
        - dem (actual dem)
        - force (boolean to force recreation of the ldd)
        - outflowdepth (set to 10.0E35 normally but smaller if needed)

    Output:
        - the LDD

    """
    if os.path.exists(lddname) and not force:
        if Verbose:
            print(("Returning existing ldd", lddname))
            return pcr.readmap(lddname)
    else:
        if Verbose:
            print(("Creating ldd", lddname))
            LDD = pcr.lddcreate(dem, 10.0e35, outflowdepth, 10.0e35, 10.0e35)
            pcr.report(LDD, lddname)
            return LDD


def configget(config, section, var, default):
    """

    Gets a string from a config file (.ini) and returns a default value if
    the key is not found. If the key is not found it also sets the value
    with the default in the config-file

    Input:
        - config - python ConfigParser object
        - section - section in the file
        - var - variable (key) to get
        - default - default string

    Returns:
        - string - either the value from the config file or the default value


    """
    Def = False
    try:
        ret = config.get(section, var)
    except:
        Def = True
        ret = default
        configset(config, section, var, default, overwrite=False)

    default = Def
    return ret


def configset(config, section, var, value, overwrite=False):
    """
    Sets a string in the in memory representation of the config object
    Deos NOT overwrite existing values if overwrite is set to False (default)

    Input:
        - config - python ConfigParser object
        - section - section in the file
        - var - variable (key) to set
        - value - the value to set
        - overwrite (optional, default is False)

    Returns:
        - nothing

    """

    if not config.has_section(section):
        config.add_section(section)
        config.set(section, var, value)
    else:
        if not config.has_option(section, var):
            config.set(section, var, value)
        else:
            if overwrite:
                config.set(section, var, value)


def configsection(config, section):
    """
    gets the list of keys in a section

    Input:
        - config
        - section

    Output:
        - list of keys in the section
    """
    try:
        ret = config.options(section)
    except:
        ret = []

    return ret


def getrows():
    """
    returns the number of rows in the current map

    Input:
        - -

    Output:
        - nr of rows in the current clonemap as a scalar
    """
    a = pcr.pcr2numpy(pcr.celllength(), np.nan).shape[0]

    return a


def getcols():
    """
    returns the number of columns in the current map

    Input:
        - -

    Output:
        - nr of columns in the current clonemap as a scalar
    """
    a = pcr.pcr2numpy(pcr.celllength(), np.nan).shape[1]

    return a


def getgridparams():
    """ return grid parameters in a python friendly way

    Output:
        [ Xul, Yul, xsize, ysize, rows, cols]

        - xul - x upper left centre
        - yul - y upper left centre
        - xsize - size of a cell in x direction
        - ysize - size of a cell in y direction
        - cols - number of columns
        - rows - number of rows
        - xlr -  x lower right centre
        - ylr -  y lower right centre
    """
    # This is the default, but add for safety...
    pcr.setglobaloption("coorcentre")
    # x and Y are the same for now
    xy = pcr.pcr2numpy(pcr.celllength(), np.nan)[0, 0]
    xu = pcr.pcr2numpy(pcr.xcoordinate(1), np.nan)[0, 0]
    yu = pcr.pcr2numpy(pcr.ycoordinate(1), np.nan)[0, 0]
    ylr = pcr.pcr2numpy(pcr.ycoordinate(1), np.nan)[getrows() - 1, getcols() - 1]
    xlr = pcr.pcr2numpy(pcr.xcoordinate(1), np.nan)[getrows() - 1, getcols() - 1]

    return [xu, yu, xy, xy, getrows(), getcols(), xlr, ylr]


def snaptomap(points, mmap):
    """
    Snap the points in _points_ to nearest non missing
    values in _mmap_. Can be used to move gauge locations
    to the nearest rivers.

    Input:
        - points - map with points to move
        - mmap - map with points to move to

    Return:
        - map with shifted points
    """
    points = pcr.cover(points, 0)
    # Create unique id map of mmap cells
    unq = pcr.nominal(pcr.cover(pcr.uniqueid(pcr.defined(mmap)), pcr.scalar(0.0)))
    # Now fill holes in mmap map with lues indicating the closes mmap cell.
    dist_cellid = pcr.scalar(pcr.spreadzone(unq, 0, 1))
    # Get map with values at location in points with closes mmap cell
    dist_cellid = pcr.ifthenelse(points > 0, dist_cellid, 0)
    # Spread this out
    dist_fill = pcr.spreadzone(pcr.nominal(dist_cellid), 0, 1)
    # Find the new (moved) locations
    npt = pcr.uniqueid(pcr.boolean(pcr.ifthen(dist_fill == unq, unq)))
    # Now recreate the original value in the points maps
    ptcover = pcr.spreadzone(pcr.cover(points, 0), 0, 1)
    # Now get the org point value in the pt map
    nptorg = pcr.ifthen(npt > 0, ptcover)

    return nptorg


def riverlength(ldd, order):
    """
    Determines the length of a river using the ldd.
    only determined for order and higher.

    Input:
        - ldd, order (streamorder)

    Returns:
        - totallength,lengthpercell, streamorder
    """
    strorder = pcr.streamorder(ldd)
    strorder = pcr.ifthen(strorder >= pcr.ordinal(order), strorder)
    dist = pcr.max(
        pcr.celllength(), pcr.ifthen(pcr.boolean(strorder), pcr.downstreamdist(ldd))
    )

    return pcr.catchmenttotal(pcr.cover(dist, 0), ldd), dist, strorder


def upscale_riverlength(ldd, order, factor):
    """
    Upscales the riverlength using 'factor'
    The resulting maps can be resampled (e.g. using resample.exe) by factor and should
    include the accurate length as determined with the original higher
    resolution maps.  This function is **depricated**,
    use are_riverlength instead as this version
    is very slow for large maps

    Input:
        - ldd
        - minimum streamorder to include

    Output:
        - distance per factor cells
    """

    strorder = pcr.streamorder(ldd)
    strorder = pcr.ifthen(strorder >= order, strorder)
    dist = pcr.cover(
        pcr.max(
            pcr.celllength(), pcr.ifthen(pcr.boolean(strorder), pcr.downstreamdist(ldd))
        ),
        0,
    )
    totdist = pcr.max(
        pcr.ifthen(
            pcr.boolean(strorder),
            pcr.windowtotal(
                pcr.ifthen(pcr.boolean(strorder), dist), pcr.celllength() * factor
            ),
        ),
        dist,
    )

    return totdist


def area_riverlength_factor(ldd, Area, Clength):
    """
    ceates correction factors for riverlength for
    the largest streamorder in each area

    Input:
        - ldd
        - Area
        - Clength (1d length of a cell (pcr.sqrt(Area))

    Output:
        - distance per area

    """
    strorder = pcr.streamorder(ldd)
    strordermax = pcr.areamaximum(strorder, Area)
    dist = pcr.downstreamdist(ldd)
    # count nr of strorder cells in each area
    nr = pcr.areatotal(pcr.ifthen(strorder == strordermax, dist), Area)
    # N = pcr.sqrt(pcr.areatotal(pcr.scalar(pcr.boolean(Area)),Area))
    N = Clength
    factor = nr / N

    return factor


def area_river_burnin(ldd, dem, order, Area):
    """
  Calculates the lowest values in as DEM for each erea in an area map for
  river of order *order*

  Input:
      - ldd
      - dem
      - order
      - Area map

  Output:
      - dem
  """
    strorder = pcr.streamorder(ldd)
    strordermax = pcr.areamaximum(strorder, Area)
    maxordcell = pcr.ifthen(strordermax > order, strordermax)
    riverdem = pcr.areaminimum(dem, Area)

    return pcr.ifthen(pcr.boolean(maxordcell), riverdem)


def area_percentile(inmap, area, n, order, percentile):
    """
  calculates percentile of inmap per area
  n is the number of points in each area,
  order, the sorter order of inmap per area (output of
  areaorder(inmap,area))
  n is the output of pcr.areatotal(pcr.spatial(pcr.scalar(1.0)),area)

  Input:
      - inmap
      - area map
      - n
      - order (riverorder)
      - percentile

  Output:
      - percentile map

  """
    i = pcr.rounddown((n * percentile) / 100.0 + 0.5)  # index in order map
    perc = pcr.ifthen(i == order, inmap)

    return pcr.areaaverage(perc, area)


def find_outlet(ldd):
    """
    Tries to find the outlet of the largest catchment in the Ldd

    Input:
        - Ldd

    Output:
        - outlet map (single point in the map)
    """
    largest = pcr.mapmaximum(pcr.catchmenttotal(pcr.spatial(pcr.scalar(1.0)), ldd))
    outlet = pcr.ifthen(
        pcr.catchmenttotal(1.0, ldd) == largest, pcr.spatial(pcr.scalar(1.0))
    )

    return outlet


def subcatch(ldd, outlet):
    """
    Determines a subcatchment map using LDD and outlet(s). In the resulting
    subcatchment map the i's of the catchment are determiend by the id's of
    the outlets.

    Input:
        - ldd
        - Outlet - maps with points for each outlet.

    Output:
        - map of subcatchments
    """
    subcatch = pcr.subcatchment(ldd, pcr.ordinal(outlet))

    return subcatch


def areastat(Var, Area):
    """
    Calculate several statistics of *Var* for each unique id in *Area*

    Input:
        - Var
        - Area

    Output:
        - Standard_Deviation,Average,Max,Min

    """
    Avg = pcr.areaaverage(Var, Area)
    Sq = (Var - Avg) ** 2
    N = pcr.areatotal(pcr.spatial(pcr.cellarea()), Area) / pcr.cellarea()
    Sd = (pcr.areatotal(Sq, Area) / N) ** 0.5
    Max = pcr.areamaximum(Var, Area)
    Min = pcr.areaminimum(Var, Area)

    return Sd, Avg, Max, Min


def checkerboard(mapin, fcc):
    """
    checkerboard create a checkerboard map with unique id's in a
    fcc*fcc cells area. The resulting map can be used
    to derive statistics for (later) upscaling of maps (using the fcc factor)

    .. warning: use with unitcell to get most reliable results!

    Input:
        - map (used to determine coordinates)
        - fcc (size of the areas in cells)

    Output:
        - checkerboard type map
    """
    msker = pcr.defined(mapin)
    ymin = pcr.mapminimum(pcr.ycoordinate(msker))
    yc = (pcr.ycoordinate((msker)) - ymin) / pcr.celllength()
    yc = pcr.rounddown(yc / fcc)
    # yc = yc/fcc
    xmin = pcr.mapminimum(pcr.xcoordinate((msker)))
    xc = (pcr.xcoordinate((msker)) - xmin) / pcr.celllength()
    xc = pcr.rounddown(xc / fcc)
    # xc = xc/fcc

    yc = yc * (pcr.mapmaximum(xc) + 1.0)

    xy = pcr.ordinal(xc + yc)

    return xy


def subcatch_stream(
    ldd,
    threshold,
    min_strahler=-999,
    max_strahler=999,
    assign_edge=False,
    assign_existing=False,
    up_area=None,
):
    """
    (From Deltares Hydrotools)

    Derive catchments based upon strahler threshold
    Input:
        ldd -- pcraster object direction, local drain directions
        threshold -- integer, strahler threshold, subcatchments ge threshold
            are derived
        min_strahler -- integer, minimum strahler threshold of river catchments
            to return
        max_strahler -- integer, maximum strahler threshold of river catchments
            to return
        assign_unique=False -- if set to True, unassigned connected areas at
            the edges of the domain are assigned a unique id as well. If set
            to False, edges are not assigned
        assign_existing=False == if set to True, unassigned edges are assigned
            to existing basins with an upstream weighting. If set to False,
            edges are assigned to unique IDs, or not assigned
    output:
        stream_ge -- pcraster object, streams of strahler order ge threshold
        subcatch -- pcraster object, subcatchments of strahler order ge threshold

    """
    # derive stream order

    stream = pcr.streamorder(ldd)
    stream_ge = pcr.ifthen(stream >= threshold, stream)
    stream_up_sum = pcr.ordinal(pcr.upstream(ldd, pcr.cover(pcr.scalar(stream_ge), 0)))
    # detect any transfer of strahler order, to a higher strahler order.
    transition_strahler = pcr.ifthenelse(
        pcr.downstream(ldd, stream_ge) != stream_ge,
        pcr.boolean(1),
        pcr.ifthenelse(
            pcr.nominal(ldd) == 5,
            pcr.boolean(1),
            pcr.ifthenelse(
                pcr.downstream(ldd, pcr.scalar(stream_up_sum)) > pcr.scalar(stream_ge),
                pcr.boolean(1),
                pcr.boolean(0),
            ),
        ),
    )
    # make unique ids (write to file)
    transition_unique = pcr.ordinal(pcr.uniqueid(transition_strahler))

    # derive upstream catchment areas (write to file)
    subcatch = pcr.nominal(pcr.subcatchment(ldd, transition_unique))

    if assign_edge:
        # fill unclassified areas (in pcraster equal to zero) with a unique id, above the maximum id assigned so far
        unique_edge = pcr.clump(pcr.ifthen(subcatch == 0, pcr.ordinal(0)))
        subcatch = pcr.ifthenelse(
            subcatch == 0,
            pcr.nominal(pcr.mapmaximum(pcr.scalar(subcatch)) + pcr.scalar(unique_edge)),
            pcr.nominal(subcatch),
        )
    elif assign_existing:
        # unaccounted areas are added to largest nearest draining basin
        if up_area is None:
            up_area = pcr.ifthen(
                pcr.boolean(pcr.cover(stream_ge, 0)), pcr.accuflux(ldd, 1)
            )
        riverid = pcr.ifthen(pcr.boolean(pcr.cover(stream_ge, 0)), subcatch)

        friction = 1.0 / pcr.scalar(
            pcr.spreadzone(pcr.cover(pcr.ordinal(up_area), 0), 0, 0)
        )  # *(pcr.scalar(ldd)*0+1)
        delta = pcr.ifthen(
            pcr.scalar(ldd) >= 0,
            pcr.ifthen(
                pcr.cover(subcatch, 0) == 0,
                pcr.spreadzone(pcr.cover(riverid, 0), 0, friction),
            ),
        )
        subcatch = pcr.ifthenelse(pcr.boolean(pcr.cover(subcatch, 0)), subcatch, delta)

    # finally, only keep basins with minimum and maximum river order flowing through them
    strahler_subcatch = pcr.areamaximum(stream, subcatch)
    subcatch = pcr.ifthen(
        pcr.ordinal(strahler_subcatch) >= min_strahler,
        pcr.ifthen(pcr.ordinal(strahler_subcatch) <= max_strahler, subcatch),
    )

    return stream_ge, pcr.ordinal(subcatch)


def subcatch_order_a(ldd, oorder):
    """
    Determines subcatchments using the catchment order

    This version uses the last cell BELOW order to derive the
    catchments. In general you want the _b version

    Input:
        - ldd
        - order - order to use

    Output:
        - map with catchment for the given streamorder
    """
    outl = find_outlet(ldd)
    large = pcr.subcatchment(ldd, pcr.boolean(outl))
    stt = pcr.streamorder(ldd)
    sttd = pcr.downstream(ldd, stt)
    pts = pcr.ifthen((pcr.scalar(sttd) - pcr.scalar(stt)) > 0.0, sttd)
    dif = pcr.upstream(
        ldd,
        pcr.cover(
            pcr.ifthen(
                large,
                pcr.uniqueid(pcr.boolean(pcr.ifthen(stt == pcr.ordinal(oorder), pts))),
            ),
            0,
        ),
    )
    dif = pcr.cover(pcr.scalar(outl), dif)  # Add catchment outlet
    dif = pcr.ordinal(pcr.uniqueid(pcr.boolean(dif)))
    sc = pcr.subcatchment(ldd, dif)

    return sc, dif, stt


def subcatch_order_b(
    ldd, oorder, sizelimit=0, fill=False, fillcomplete=False, stoporder=0
):
    """
    Determines subcatchments using the catchment order

    This version tries to keep the number op upstream/downstream catchment the
    small by first dederivingatchment connected to the major river(the order) given, and fill
    up from there.

    Input:
        - ldd
        - oorder - order to use
        - sizelimit - smallest catchments to include, default is all (sizelimit=0) in number of cells
        - if fill is set to True the higer order catchment are filled also
        - if fillcomplete is set to True the whole ldd is filled with catchments.


    :returns sc, dif, nldd; Subcatchment, Points, subcatchldd
    """
    # outl = find_outlet(ldd)
    # large = pcr.subcatchment(ldd,pcr.boolean(outl))

    if stoporder == 0:
        stoporder = oorder

    stt = pcr.streamorder(ldd)
    sttd = pcr.downstream(ldd, stt)
    pts = pcr.ifthen((pcr.scalar(sttd) - pcr.scalar(stt)) > 0.0, sttd)
    maxorder = pcraster.framework.getCellValue(pcr.mapmaximum(stt), 1, 1)
    dif = pcr.uniqueid(pcr.boolean(pcr.ifthen(stt == pcr.ordinal(oorder), pts)))

    if fill:
        for order in range(oorder, maxorder):
            m_pts = pcr.ifthen((pcr.scalar(sttd) - pcr.scalar(order)) > 0.0, sttd)
            m_dif = pcr.uniqueid(
                pcr.boolean(pcr.ifthen(stt == pcr.ordinal(order), m_pts))
            )
            dif = pcr.uniqueid(pcr.boolean(pcr.cover(m_dif, dif)))

        for myorder in range(oorder - 1, stoporder, -1):
            sc = pcr.subcatchment(ldd, pcr.nominal(dif))
            m_pts = pcr.ifthen((pcr.scalar(sttd) - pcr.scalar(stt)) > 0.0, sttd)
            m_dif = pcr.uniqueid(
                pcr.boolean(pcr.ifthen(stt == pcr.ordinal(myorder - 1), m_pts))
            )
            dif = pcr.uniqueid(
                pcr.boolean(pcr.cover(pcr.ifthen(pcr.scalar(sc) == 0, m_dif), dif))
            )

        if fillcomplete:
            sc = pcr.subcatchment(ldd, pcr.nominal(dif))
            cs, m_dif, stt = subcatch_order_a(ldd, stoporder)
            dif = pcr.uniqueid(
                pcr.boolean(
                    pcr.cover(
                        pcr.ifthen(pcr.scalar(sc) == 0, pcr.ordinal(m_dif)),
                        pcr.ordinal(dif),
                    )
                )
            )

    scsize = pcr.catchmenttotal(1, ldd)
    dif = pcr.ordinal(pcr.uniqueid(pcr.boolean(pcr.ifthen(scsize >= sizelimit, dif))))
    sc = pcr.subcatchment(ldd, dif)

    # Make pit ldd
    nldd = pcr.lddrepair(pcr.ifthenelse(pcr.cover(dif, 0) > 0, 5, ldd))

    return sc, dif, nldd


def getRowColPoint(in_map, xcor, ycor):
    """
    returns the row and col in a map at the point given.
    Works but is rather slow.

    Input:
        - in_map - map to determine coordinates from
        - xcor - x coordinate
        - ycor - y coordinate

    Output:
        - row, column
    """
    x = pcr.pcr2numpy(pcr.xcoordinate(pcr.boolean(pcr.scalar(in_map) + 1.0)), np.nan)
    y = pcr.pcr2numpy(pcr.ycoordinate(pcr.boolean(pcr.scalar(in_map) + 1.0)), np.nan)
    XX = pcr.pcr2numpy(pcr.celllength(), 0.0)
    tolerance = 0.5  # takes a single point

    diffx = x - xcor
    diffy = y - ycor
    col_ = np.absolute(diffx) <= (XX[0, 0] * tolerance)  # cellsize
    row_ = np.absolute(diffy) <= (XX[0, 0] * tolerance)  # cellsize
    point = col_ * row_

    return point.argmax(0).max(), point.argmax(1).max()


def getValAtPoint(in_map, xcor, ycor):
    """
    returns the value in a map at the point given.
    works but is rather slow.

    Input:
        - in_map - map to determine coordinates from
        - xcor - x coordinate
        - ycor - y coordinate

    Output:
        - value
    """
    x = pcr.pcr2numpy(pcr.xcoordinate(pcr.defined(in_map)), np.nan)
    y = pcr.pcr2numpy(pcr.ycoordinate(pcr.defined(in_map)), np.nan)
    XX = pcr.pcr2numpy(pcr.celllength(), 0.0)
    themap = pcr.pcr2numpy(in_map, np.nan)
    tolerance = 0.5  # takes a single point

    diffx = x - xcor
    diffy = y - ycor
    col_ = np.absolute(diffx) <= (XX[0, 0] * tolerance)  # cellsize
    row_ = np.absolute(diffy) <= (XX[0, 0] * tolerance)  # cellsize
    point = col_ * row_
    pt = point.argmax()

    return themap.ravel()[pt]


def points_to_map(in_map, xcor, ycor, tolerance):
    """
    Returns a map with non zero values at the points defined
    in X, Y pairs. It's goal is to replace the pcraster col2map program.

    tolerance should be 0.5 to select single points
    Performance is not very good and scales linear with the number of points


    Input:
        - in_map - map to determine coordinates from
        - xcor - x coordinate (array or single value)
        - ycor - y coordinate (array or single value)
        - tolerance - tolerance in cell units. 0.5 selects a single cell\
        10 would select a 10x10 block of cells

    Output:
        - Map with values burned in. 1 for first point, 2 for second and so on
    """
    point = in_map * 0.0

    x = pcr.pcr2numpy(pcr.xcoordinate(pcr.defined(in_map)), np.nan)
    y = pcr.pcr2numpy(pcr.ycoordinate(pcr.defined(in_map)), np.nan)
    cell_length = float(pcr.celllength())

    # simple check to use both floats and numpy arrays
    try:
        c = xcor.ndim
    except:
        xcor = np.array([xcor])
        ycor = np.array([ycor])

    # Loop over points and "burn in" map
    for n in range(0, xcor.size):
        if Verbose:
            print(n)
        diffx = x - xcor[n]
        diffy = y - ycor[n]
        col_ = np.absolute(diffx) <= (cell_length * tolerance)  # cellsize
        row_ = np.absolute(diffy) <= (cell_length * tolerance)  # cellsize
        point = point + pcr.numpy2pcr(pcr.Scalar, ((col_ * row_) * (n + 1)), np.nan)

    return pcr.ordinal(point)


def detdrainlength(ldd, xl, yl):
    """
    Determines the drainaige length (DCL) for a non square grid

    Input:
        - ldd - drainage network
        - xl - length of cells in x direction
        - yl - length of cells in y direction

    Output:
        - DCL
    """
    # take into account non-square cells
    # if ldd is 8 or 2 use Ylength
    # if ldd is 4 or 6 use Xlength
    draindir = pcr.scalar(ldd)
    slantlength = pcr.sqrt(xl ** 2 + yl ** 2)
    drainlength = pcr.ifthenelse(
        draindir == 2,
        yl,
        pcr.ifthenelse(
            draindir == 8,
            yl,
            pcr.ifthenelse(
                draindir == 4, xl, pcr.ifthenelse(draindir == 6, xl, slantlength)
            ),
        ),
    )

    return drainlength


def detdrainwidth(ldd, xl, yl):
    """
    Determines width of drainage over DEM for a non square grid

    Input:
        - ldd - drainage network
        - xl - length of cells in x direction
        - yl - length of cells in y direction

    Output:
        - DCL
    """
    # take into account non-square cells
    # if ldd is 8 or 2 use Xlength
    # if ldd is 4 or 6 use Ylength
    draindir = pcr.scalar(ldd)
    slantwidth = (xl + yl) * 0.5
    drainwidth = pcr.ifthenelse(
        draindir == 2,
        xl,
        pcr.ifthenelse(
            draindir == 8,
            xl,
            pcr.ifthenelse(
                draindir == 4, yl, pcr.ifthenelse(draindir == 6, yl, slantwidth)
            ),
        ),
    )
    return drainwidth


def classify(
    inmap, lower=[0, 10, 20, 30], upper=[10, 20, 30, 40], classes=[2, 2, 3, 4]
):
    """
    classify a scaler maps accroding to the boundaries given in classes.

    """

    result = pcr.ordinal(pcr.cover(-1))
    for l, u, c in zip(lower, upper, classes):
        result = pcr.cover(
            pcr.ifthen(inmap >= l, pcr.ifthen(inmap < u, pcr.ordinal(c))), result
        )

    return pcr.ifthen(result >= 0, result)


def derive_HAND(dem, ldd, accuThreshold, rivers=None, basin=None):
    """
    Function derives Height-Above-Nearest-Drain.
    See http://www.sciencedirect.com/science/article/pii/S003442570800120X
    Input:
        dem -- pcraster object float32, elevation data
        ldd -- pcraster object direction, local drain directions
        accuThreshold -- upstream amount of cells as threshold for river
            delineation
        rivers=None -- you can provide a rivers layer here. Pixels that are
                        identified as river should have a value > 0, other
                        pixels a value of zero.
        basin=None -- set a boolean pcraster map where areas with True are estimated using the nearest drain in ldd distance
                        and areas with False by means of the nearest friction distance. Friction distance estimated using the
                        upstream area as weight (i.e. drains with a bigger upstream area have a lower friction)
                        the spreadzone operator is used in this case.
    Output:
        hand -- pcraster bject float32, height, normalised to nearest stream
        dist -- distance to nearest stream measured in cell lengths
            according to D8 directions
    """
    if rivers is None:
        stream = pcr.ifthenelse(
            pcr.accuflux(ldd, 1) >= accuThreshold, pcr.boolean(1), pcr.boolean(0)
        )
    else:
        stream = pcr.boolean(pcr.cover(rivers, 0))

    height_river = pcr.ifthenelse(stream, pcr.ordinal(dem * 100), 0)
    if basin is None:
        up_elevation = pcr.scalar(pcr.subcatchment(ldd, height_river))
    else:
        drainage_surf = pcr.ifthen(rivers, pcr.accuflux(ldd, 1))
        weight = 1.0 / pcr.scalar(
            pcr.spreadzone(pcr.cover(pcr.ordinal(drainage_surf), 0), 0, 0)
        )
        up_elevation = pcr.ifthenelse(
            basin,
            pcr.scalar(pcr.subcatchment(ldd, height_river)),
            pcr.scalar(pcr.spreadzone(height_river, 0, weight)),
        )
        # replace areas outside of basin by a spread zone calculation.
    hand = pcr.max(pcr.scalar(pcr.ordinal(dem * 100)) - up_elevation, 0) / 100
    dist = pcr.ldddist(ldd, stream, 1)

    return hand, dist


def sCurve(X, a=0.0, b=1.0, c=1.0):
    """
    sCurve function:

    Input:
        - X input map
        - C determines the steepness or "stepwiseness" of the curve.
          The higher C the sharper the function. A negative C reverses the function.
        - b determines the amplitude of the curve
        - a determines the centre level (default = 0)

    Output:
        - result
    """
    s = 1.0 / (pcr.exp(-c * (X - a)) + b)
    return s


def sCurveSlope(X, a=0.0, b=1.0, c=1.0):
    """
    First derivative of the sCurve defined by a,b,c at point X

    Input:
        - X - value to calculate for
        - a
        - b
        - c

    Output:
        - first derivative (slope) of the curve at point X
    """
    sc = sCurve(X, a=a, b=b, c=c)
    slope = sc * (1 - sc)
    return slope


def Gzip(fileName, storePath=False, chunkSize=1024 * 1024):
    """
        Usage: Gzip(fileName, storePath=False, chunksize=1024*1024)
        Gzip the given file to the given storePath and then remove the file.
        A chunk size may be selected. Default is 1 megabyte
        Input:
            fileName:   file to be GZipped
            storePath:  destination folder. Default is False, meaning the file will be zipped to its own folder
            chunkSize:  size of chunks to write. If set too large, GZip will fail with memory problems
    """
    if not storePath:
        pathName = os.path.split(fileName)[0]
        fileName = os.path.split(fileName)[1]
        curdir = os.path.curdir
        os.chdir(pathName)
    # open files for reading / writing
    r_file = open(fileName, "rb")
    w_file = gzip.GzipFile(fileName + ".gz", "wb", 9)
    dataChunk = r_file.read(chunkSize)
    while dataChunk:
        w_file.write(dataChunk)
        dataChunk = r_file.read(chunkSize)
    w_file.flush()
    w_file.close()
    r_file.close()
    os.unlink(fileName)  # We don't need the file now
    if not storePath:
        os.chdir(curdir)


# These come from GLOFRIS_Utils


def zipFiles(fileList, fileTarget):
    """
    Usage: zipFiles(fileList, fileTarget)
    zip the given list of files to the given target file
    Input:
        fileList:   list of files to be zipped
        fileTarget: target zip-file
    """
    zout = zipfile.ZipFile(fileTarget, "w", compression=zipfile.ZIP_DEFLATED)
    for fname in fileList:
        zout.write(fname, arcname=os.path.split(fname)[1])
    zout.close()


def readMap(fileName, fileFormat):
    """
    Read geographical file into memory

    :param fileName:
    :param fileFormat:
    :return x, y, data, FillVal:
    """
    # Open file for binary-reading
    mapFormat = gdal.GetDriverByName(fileFormat)
    mapFormat.Register()
    ds = gdal.Open(fileName)
    if ds is None:
        print("Could not open " + fileName + ". Something went wrong!! Shutting down")
        sys.exit(1)
        # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX = geotrans[1]
    resY = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX + resX / 2, originX + resX / 2 + resX * (cols - 1), cols)
    y = np.linspace(originY + resY / 2, originY + resY / 2 + resY * (rows - 1), rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1)  # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0, 0, cols, rows)
    FillVal = RasterBand.GetNoDataValue()
    RasterBand = None
    ds = None
    return x, y, data, FillVal


def cutMapById(data, subcatchmap, id, x, y, FillVal):
    """

    :param data: 2d numpy array to cut
    :param subcatchmap: 2d numpy array with subcatch
    :param id: id (value in the array) to cut by
    :param x: array with x values
    :param y:  array with y values
    :return: x,y, data
    """

    if len(data.flatten()) == len(subcatchmap.flatten()):
        scid = subcatchmap == id
        data[np.logical_not(scid)] = FillVal
        xid, = np.where(scid.max(axis=0))
        xmin = xid.min()
        xmax = xid.max()
        if xmin >= 1:
            xmin = xmin - 1
        if xmax < len(x) - 1:
            xmax = xmax + 1

        yid, = np.where(scid.max(axis=1))
        ymin = yid.min()
        ymax = yid.max()
        if ymin >= 1:
            ymin = ymin - 1
        if ymax < len(y) - 1:
            ymax = ymax + 1

        return (
            x[xmin:xmax].copy(),
            y[ymin:ymax].copy(),
            data[ymin:ymax, xmin:xmax].copy(),
        )
    else:
        return None, None, None


def writeMap(fileName, fileFormat, x, y, data, FillVal):
    """ Write geographical data into file"""

    verbose = False
    gdal.AllRegister()
    driver1 = gdal.GetDriverByName("GTiff")
    driver2 = gdal.GetDriverByName(fileFormat)

    # Processing
    if verbose:
        print("Writing to temporary file " + fileName + ".tif")
    # Create Output filename from (FEWS) product name and data and open for writing

    if data.dtype == np.int32:
        TempDataset = driver1.Create(
            fileName + ".tif", data.shape[1], data.shape[0], 1, gdal.GDT_Int32
        )
    else:
        TempDataset = driver1.Create(
            fileName + ".tif", data.shape[1], data.shape[0], 1, gdal.GDT_Float32
        )
    # Give georeferences
    xul = x[0] - (x[1] - x[0]) / 2
    yul = y[0] + (y[0] - y[1]) / 2
    TempDataset.SetGeoTransform([xul, x[1] - x[0], 0, yul, 0, y[1] - y[0]])
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data, 0, 0)
    TempBand.FlushCache()
    TempBand.SetNoDataValue(FillVal)
    # Create data to write to correct format (supported by 'CreateCopy')
    if verbose:
        print("Writing to " + fileName + ".map")
    outDataset = driver2.CreateCopy(fileName, TempDataset, 0)
    TempDataset = None
    outDataset = None
    if verbose:
        print("Removing temporary file " + fileName + ".tif")
    os.remove(fileName + ".tif")

    if verbose:
        print("Writing to " + fileName + " is done!")
