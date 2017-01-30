Introduction
------------
PCRaster is a powerful package of software for environmental dynamic modelling. It is a dynamic modelling system for distributed simulation models. The main applications of PCRaster are found in environmental modelling: geography, hydrology, ecology to name a few. Examples include rainfall runoff models, vegetation competition models and slope stability models.

PCRaster offers a large number of operations for spatio-temporal analysis. These high-level operations allow the development of models that are substantial shorter than programs offering the same functionality, but developed with system programming languages like C or FORTRAN. This makes PCRaster model scripts easy to construct, to modify and to understand. More information about PCRaster, including software, manuals and course material can be found at the PCRaster website: http://www.pcraster.eu.

Python is an interpreted, object-oriented programming language. Python scripts are platform independent, they can be executed on different operating systems like UNIX, Linux, Windows, Macintosh and others. Python has modules, classes, exceptions, very high level dynamic data types, dynamic typing and a clean, understandable syntax. More information about Python, including software and manuals can be found at the Python website: http://www.python.org.

The PCRaster Python extension enables you to use the PCRaster modelling engine from Python. This means that you can combine and call more than 150 PCRaster operations from Python. Combined with the flexible nature of the Python programming language this gives a powerful environment for creating environmental models.

Requirements
^^^^^^^^^^^^
We assume that the user is familiar with creating models in PCRaster and has basic Python programming knowledge.

The minimum requirement to use PCRaster Python is a functioning install of `Python version 2.X <http://www.python.org>`_. If you want to use the conversion functions to/from NumPy, in addition `NumPy <http://www.numpy.org>`_ must be installed. The exact version of PCRaster Python (e.g. 2.7 or up) and NumPy matters. These versions are mentioned at the `PCRaster package download page <http://pcraster.geo.uu.nl/downloads>`_ or in the release notes.

The PCRaster Python extension is part of the standard PCRaster installation. See the separate installation note for further information.
