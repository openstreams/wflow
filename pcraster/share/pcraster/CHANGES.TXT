Changes
=======

PCRaster 4.1.0
--------------
* On Windows, shared libraries are now installed in a seperate directory called ``lib``, just like on the other platforms. The path to this directory does not have to be listed in the ``$PATH`` environment variable. This has the advantage that the shared libraries shipped with PCRaster will never conflict with shared libraries shipped with other software, and visa versa.

Aguila
^^^^^^
It is now possible again to select a specific column from a timeseries file, like it used to (`Ticket #635 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/635/>`_)

PCRaster Python package
^^^^^^^^^^^^^^^^^^^^^^^
* PCRaster data types can now be serialized by the pickle module (`Ticket #593 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/593/>`_). See also the ``pickle`` Python module. More about this in a future release.

GDAL
^^^^
* The submitted patches for the GDAL PCRaster driver support the dynamic creation of PCRaster maps (`Ticket #664 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/664/>`_) and contain some updates (`Ticket #679 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/679/>`_). We expect the functionality to be integrated in one of the next official GDAL releases.

Documentation
^^^^^^^^^^^^^
* Minor updates to documentation of ``spatial``, ``transient`` and ``horizontan`` operations (`Ticket #620 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/620/>`_, `Ticket #677 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/677/>`_, `Ticket #678 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/678/>`_)


PCRaster 4.0.2
--------------

This release improves the reliability of PCRaster applications relating to large data sets, includes functional enhancements for the Modflow extension, and fixes several bugs.


PCRaster model engine
^^^^^^^^^^^^^^^^^^^^^
* An erroneous check was replaced to allow for processing maps with more than 2^31 - 1 cells on 64bit systems (`Ticket #648 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/648/>`_)

col2map
^^^^^^^
* We fixed the incorrect allocation of cells when using maps with more than 2^31 - 1 cells (`Ticket #661 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/661/>`_)
* We fixed the incorrect printout of cells with more than one record (`Ticket #660 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/660/>`_)


Modflow extension
^^^^^^^^^^^^^^^^^

*  Added functions to obtain cell-by-cell values for storage, constant heads, and front/right/lower flow face (`Ticket #643 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/643/>`_)
*  Fixed incorrect reporting of the Python getRiverLeakage (`Ticket #663 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/663/>`_)
*  Fixed activation of well package for time steps > 1 (`Ticket #658 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/658/>`_)
*  Fixed setting of rewetting thresholds for top layer type laycon 1 (`Ticket #657 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/657/>`_)
*  Fixed uѕage of reserved unit numbers for Modflow input files (`Ticket #662 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/662/>`_)
*  Maps now can be reported at time steps where Modflow fails to converge (`Ticket #669 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/669/>`_)
*  Minor performance improvements reducing the memory and I/O load

Aguila
^^^^^^
* We fixed the incorrect colour assignment of 2D directional data types (`Ticket #670 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/670/>`_)
* We fixed the incorrect colour assignment of 3D directional and scalar data types (`Ticket #641 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/641/>`_)

Documentation
^^^^^^^^^^^^^
* Minor updates (`Ticket #628 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/628/>`_, `Ticket #659 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/659/>`_)


PCRaster 4.0.1
--------------
This is a bug fix release for 4.0.

Global options ``chezy`` and ``manning`` for dynwavestate, dynwaveflux, dynamicwave (pcrcalc, PCRaster Python package)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We discovered a documentation error for the operations dynwavestate, dynwaveflux and dynamicwave.
The manual stated that the Chezy algorithm was the default algorithm to calculate the dynamic flow equation.
In fact, it was calculated by the Manning algorithm by default.

If you did not use any global option, your results were calculated by the Manning equation. From now on, without specifying global options, results will be calculated by the Manning equation as well.

If you used either ``chezy`` or ``manning`` as global option, the corresponding algorithms were used. This behaviour remains unchanged.

To obtain values calculated with the Chezy algorithm, you now need to specify explicitly either
``--chezy`` on the command line, ``#! --chezy`` in PCRcalc scripts, or ``setglobaloption("chezy")`` in Python scripts.

dynamicwave (pcrcalc, PCRaster Python package)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We discovered and fixed a bug in the dynamicwave operation while using the Manning algorithm (`Ticket #609 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/609/>`_).
As the Manning algorithm was used as default (see the remarks above) it is expected that your model results will change.

Aguila
^^^^^^
* Fixed '#624 cannot open .tss files in PCRASTER 4 aguila version' (`Ticket #624 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/624/>`_)

PCRaster Python package
^^^^^^^^^^^^^^^^^^^^^^^
* Fixed a wrong number of arguments in the base class for dynamic models (`Ticket #603 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/603/>`_)

resample
^^^^^^^^
* Fixed a regression that caused the generation of MV in all cells while using the crop option (`Ticket #485 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/485/>`_)

pcrcalc
^^^^^^^
* Fixed a redundant timestep output (`Ticket #625 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/625/>`_)

Documentation
^^^^^^^^^^^^^
* The manual pages include updates for the mapattr application and the lookupstate and lookuppotential operations (`Ticket #613 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/613/>`_, `Ticket #601 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/601/>`_)

Developer information
^^^^^^^^^^^^^^^^^^^^^
* Ported machine_status.py to newer apt_pkg, updated list of required applications for compiling PCRaster (`Ticket #610 <https://sourceforge.net/p/pcraster/bugs-and-feature-requests/610/>`_)



PCRaster 4.0.0
--------------
General
^^^^^^^
* Changed the license of all PCRaster source code to the `GPLv3 <http://www.gnu.org/licenses/gpl-3.0.html>`_ open source license. Moved all sourcecode to the `PCRaster Open Source Tools site <https://sourceforge.net/projects/pcraster/>`_ at SourceForge.
* The installation process of PCRaster has been simplified. On all platforms we distribute a zip file which can be unzipped at a preferred location. After setting two environment variables, PCRaster is ready to be used. The goal is to make it possible to install multiple versions of PCRaster at the same time. This has the advantage that older models can still be run with older installed versions of PCRaster. And it allows us to keep improving PCRaster, even if we break backwards compatibility (we prefer not to, but sometimes there is a good reason).
* Removed support for reading HDF4 formatted rasters. Maintaining support for this format proved to be too much of a hassle.

pcrcalc
^^^^^^^
* Removed support for encrypting models.
* Removed support for license specific functionality (like missing value compression). All features that used to require a commercial license are available for everybody now.

resample
^^^^^^^^
* Fixed the spurious creation of adjacent raster cells while using resample as cookie cutter (`Ticket #463 <http://sourceforge.net/p/pcraster/bugs-and-feature-requests/463/>`_)

PCRaster Python package
^^^^^^^^^^^^^^^^^^^^^^^
* Updated the code to allow the garbage collector to reclaim memory used by some of the framework class instanceѕ, after the last reference goes out of scope.
* Updated the code to prevent that the memory used by the PCRaster Python extension increases during a model run.
* PCRaster Python package now depends on Python 2.7.
* PCRaster Python package uses lower case names for package names. Update all PCRaster related imports and change them to lower case. See also the `Style Guide for Python Code <http://www.python.org/dev/peps/pep-0008/>`_.
* Removed ``pcraster.numpy`` sub-package. Numpy functionality is merged in the ``pcraster`` main package and available without an explicit import of the ``numpy`` sub-package. Remove any import of ``pcraster.numpy`` and rename any calls of ``pcraster.numpy.pcr2numpy`` and ``pcraster.numpy.numpy2pcr`` to ``pcraster.pcr2numpy`` and ``pcraster.numpy2pcr``.
* Removed ``pcr2numarray`` and ``numarray2pcr`` which were already deprecated. Use ``pcr2numpy`` and ``numpy2pcr``.
* Reimplemented ``numpy2pcr``. It is faster now.
* Added a `setclone` overload taking `nrRows`, `nrCols`, `cellSize`, `west`, `north`. No need to pass the name of an existing raster anymore.

MODFLOW extension
^^^^^^^^^^^^^^^^^
* Fixed a crash.
* Renamed extension from ``PCRasterModflow`` to ``pcraster_modflow``.
