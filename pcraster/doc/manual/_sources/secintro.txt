

.. _secintro:

**************************************************************
Introduction to the PCRaster Package: Concepts, Package Layout
**************************************************************


.. _secintrointro:

Introduction
============

PCRaster is a Geographical Information System which consists of a set of
computer tools for storing, manipulating, analyzing and retrieving
geographic information. It is a raster-based system that uses a strict data
type checking mechanism. This means that data type information is added
to all spatial data, based upon the kind of attribute that the data represent.
The use of data types controls the way the data are stored in the database
and the possibilities for manipulation and analysis of the data. This strict
data type checking mechanism has the advantage that it helps the user to
organize the data and prevents the execution of operations that will make
a nonsense result. The spatial data types implemented discriminate between
various types of continuous fields and classified objects.


.. _PCRastermod:

.. index::
   single: PCRaster; modules of

PCRaster has a relatively open database. The architecture of the system permits the integration of environmental modelling functions with classical GIS functions such as database maintenance, screen display and hard copy output.  The modules for Cartographic and Dynamic Modelling are integrated with the GIS at a high level, which means that the GIS functions and modelling functions are incorporated in a single GIS and modelling language for performing both GIS and modelling operations. 


The module for geostatistical modelling (gstat module) are
integrated at a medium level with the GIS part of the package: both 
use a separate set of functions for manipulating the data; the map format
of the central database is used and files can be automatically exchanged
between the modules and the GIS part of the package. Exchange of ascii
files with any other modelling or GIS package can be done using PCRaster
conversion operations.


.. _PCRasterConcepts:



.. _CellConcept:

.. index::
   single: PCRaster; concept of

.. index::
   single: cell; concept of

The central concept of PCRaster is a discretization of the landscape in space, resulting in cells of information.  Each cell can be regarded as a set of attributes defining its properties, but one which can receive and transmit information to and from neighbouring cells. This representation of the landscape is often referred to as 2.5 D: the lateral directions in a landscape are represented by a set of neighbouring cells making up a map; relations in vertical directions, for instance between soil layers, are implemented using several attributes stored in each cell. GIS operations or operations used in modelling can be regarded as functions that induce a change in the properties of the cells on the basis of the relations within cells (between attributes on one cell location) or between cells, see figure below. In PCRaster, each functional module of the package represents a group of operations that change the properties of the cells in a specific way. 

.. _fig2.2:

.. figure:: ../figures/block.png

   The cell: relations between attributes  within the cell and relations in lateral directions with attributes stored in neighbouring cells.



.. _secintrogis:

GIS and Cartographic Modelling
==============================
.. _CartModModIntro:



.. _GISModIntro:



.. _digitizing:



.. _scanning:

.. index::
   single: Cartographic Modelling; module for, introduction

.. index::
   single: GIS; module for, introduction

.. index::
   single: digitizing

.. index::
   single: scanning

The central module of the PCRaster system is the group of PCRaster operations where the operations for Cartographic Modelling are integrated at a high level with the GIS functions of the package. The main GIS functions supported are user interfaces (i.a. screen display, hard copy output), conversion of data with other GIS packages and database management. No digitizing  and scanning functionalities are implemented, but data transfer to and from other GIS packages that support these functionalities is simple. Spatial data are stored in the database as :emphasis:`PCRaster maps`, this is a binary format used for representation of raster maps in PCRaster. 

.. _MapAlgebraIntro:

.. index::
   single: Map Algebra; description of

The Cartographic Modelling part consist of operators for the static analysis of maps. This set of operators follows the concept of Map Algebra  and Cartographic Modelling. There are several versions of Map Algebra, all with different names, but the concept used in PCRaster is strongly related to the concept of the MAP package designed by Tomlin (:ref:`tomlin80 <bibliography>`, :ref:`tomlin90 <bibliography>`) and the algebra used by :ref:`berry87b <bibliography>`. The Cartographic Modelling part consists of a set of primitive operators that induce a change in the properties of the cells, where the change in properties is calculated on the basis of some kind of dependency within cells (point operations) or between cells (neighbourhood operations, area operations, map operations). An extensive set of operators is available in the PCRaster system: several point operators (analytical and arithmetical functions, Boolean operators, operators for relations, comparison, rounding, field generation etc.), neighbourhood operators for calculations in moving windows (highpass filtering, edge filtering, moving averages, etc.), area operators for calculations within specified areas (for instance soil groups), operators for the calculation of cost paths. In the PCRaster package a rich suite of geomorphological and hydrological functions is available that goes behind the range of operations generally considered as Map Algebra. These include functions for visibility analysis, catchment analysis and routing of transport (drainage) of material in a catchment using interactively generated local drain direction maps and transport (routing) operations. 


This set of operators is a computer language designed especially for
spatial and temporal analysis. It is an algebraic language, which means
that the PCRaster operations can be applied and combined in the same way
as algebraic calculations. In general an operation is done by typing:

  | pcrcalc Result = PCRasterOperator(PCRasterExpression)

where pcrcalc activates the PCRaster operation shell and PCRasterOperator is one of the PCRaster operations resulting in the :emphasis:`Result`. The :emphasis:`Result` may be a map or a non spatial value. The operation is done on the map :emphasis:`PCRasterExpression`. This map is called an :emphasis:`Expression` because it may be a map, but it may also be a PCRaster operation or a set of operations that result in a map or a non spatial value. This means that several PCRaster operations may be nested in one command. For instance the maximum slope on a slope map generated on basis of the elevation map Elevation can be calculated in one command line using the mapmaximum and slope operators by typing: 

.. parsed-literal::

   pcrcalc Result = mapmaximum(slope(Elevation)) 

As explained above, the separate commands are applied from the command line. Additionally, by linking these commands into PCRaster scripts or programs, it is possible to perform a theoretically unlimited number of commands consecutively. This is usually referred to as (static) Cartographic Modelling. Cartographic Modelling does not have a concept of time: several operations are performed consecutively, but they do not necessarily represent a process over time: the operations performed represent one, static change in the property of cells.  


Elements of GIS and Catrographic modelling, like the database of the PCRaster package 
(:ref:`secdatbase`), GIS functions (:ref:`secimport`) and Cartographic Modelling (:ref:`secstat`) are described in the next chapters. 



.. _secintrodyn:

Dynamic Modelling
=================

The central idea of Cartographic Modelling is the derivation of new cell
attributes from these attributes already present, or from attributes
of neighbouring cells. In Dynamic Modelling the principle of spatial
modelling is elaborated further by adding the concept of time: new
attributes are computed as a function of attribute changes over time.



The Dynamic Modelling module is integrated at a high level with the part
of the package for GIS functions and Cartographic Modelling. It provides
a meta-language within which the user can build a dynamic model with the
operators that are also used for Cartographic Modelling. Extra operators
are added for creation of iterations through time and the reading of
time series. The dynamic modelling language can be used for building a
wide range of models, from very simple (point) models up to conceptually
complicated or physically based models for environmental modelling (for
instance erosion models). A dynamic model developed in the PCRaster
Dynamic Modelling module can carry out all steps performed in modelling
with ordinary low level GIS integrated models such as MODFLOW, MICROFEM
(i.e. validating and calibrating). It has the advantage that the model
is integrated at a high level with the GIS: data exchange problems do
not exist, because the database of the GIS :emphasis:`is` the database of the model and vice versa. So it is very easy to run models using distributed data sets imported, analyzed and manipulated with the GIS and Cartographic Modelling part of the PCRaster package or created with the other modules of the package. Additionally the results of model runs can be visualised and analyzed without further data exchange.  :ref:`secdyn` covers the Dynamic Modelling module. 



.. _secintrogstat:

gstat module: Geostatistical Modelling
======================================

Gstat is the module of PCRaster for Geostatistical Modelling.
Gstat is a seperate application, for details look at www.gstat.org.
Most Likely your PCRaster installation came with a working version of gstat
and seperate manual.
Here we give a short description
of geostatistical modelling and the functionality of gstat.



In geostatistical modelling it is assumed that the property of the cells is a
realization of a spatial random function. It includes modelling the spatial
dependence on basis of the known cell values and spatial prediction where
values are predicted at cells with an unknown value using cells with a
known value. In gstat, modelling the spatial dependence is done by
estimating the variogram, the (pseudo) cross variogram, covariogram or
cross covariogram and fitting (nested) variogram models (with interactive
graphical display). Tools for spatial prediction are simple, ordinary or
universal, univariable or multivariable, point or block kriging or conditional
simulation. Simple inverse distance functions are also available.  



The gstat module is integrated at a medium level with the GIS part of
the package. It is a separate module, but conversion of data with the central
database is simple. For interpolating point data, it uses a point data column
format also used in Geo-EAS (:ref:`secdatbasepointform`). This point data format can easily be converted to PCRraster map format. The output from gstat of spatial data is in PCRaster map format: when performing interpolations in gstat PCRaster formatted maps can be used as a mask to specify the area over which interpolations are done and the resolution of the interpolated map. The resulting interpolated maps are in PCRaster map format and can be visualized and analyzed using the PCRaster operators. 

