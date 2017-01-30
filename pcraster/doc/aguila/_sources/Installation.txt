.. _installation:

************
Installation
************

.. note::

  The information in this section is related to the stand-alone binary Aguila distribution and is not valid when Aguila is distributed as part of a PCRaster distribution. The source code and information about building Aguila can be obtained from the `PCRaster Open Source Tools project website`__ on SourceForge.


Aguila installation packages can be found at the `files section`__ of the SourceForge project page.

__ http://pcraster.sourceforge.net
__ http://sourceforge.net/projects/pcraster/files

.. _microsoftWindows:

Microsoft Windows
=================
The software is distributed as a Windows installer. Install Aguila by executing the installer. After installation there is an option to uninstall Aguila again.

The following directories are created within the installation directory::

  bin    # Contains the executable and the required shared libraries.
  demos  # Contains demo data.
  doc    # Contains this manual, in pdf format.
  share  # Contains support files.

When you install Aguila as a user with Administrator rights, the search path for executables is automatically updated to include the path to the Aguila binary. Otherwise, you need to update the search path yourself:

#. Right-click on ``My Computer`` and select ``Properties``.
#. Go to ``Advanced`` tab and select ``Environment Variables``.
#. Edit the ``PATH`` variable from the ``User variables for user`` list, or add it if it not already exists.
#. Insert the path to the Aguila binary into the edit field (separate paths using a semicolon).
#. Open a new command prompt which will have the new environment setting.

.. _gnuLinux:

GNU Linux
=========

.. warning::

  The binary distribution of Aguila for Linux is created on the PCRaster development machines. This version of Aguila will run on systems with similar properties (currently Debian__ testing__ on i686__ and x86_64__). On other machines this version may not run [#]_. Here is a list of Linux platforms for which Aguila is reported to run:

  * i686 / Debian testing
  * i686 / Ubuntu 8.04
  * i686 / Ubuntu 9.04
  * x86_64 / Debian testing
  * x86_64 / Ubuntu 8.04

__ http://www.debian.org
__ http://www.debian.org/releases/testing
__ http://en.wikipedia.org/wiki/Intel_P6
__ http://en.wikipedia.org/wiki/X86-64
.. [#] We are in the process of making Aguila LSB__-compliant.
__ http://www.linuxfoundation.org/en/LSB

The software is distributed as a compressed tar file which can be unzipped as follows::

  tar zxf Aguila-<architecture>-<version>-Linux.tar.gz

This will result in a directory with the following subdirectories::

  bin    # Contains the executable.
  demos  # Contains demo data.
  doc    # Contains this manual, in pdf format.
  lib    # Contains the required shared libraries.
  share  # Contains support files.

To be able to use Aguila you need to add the path to the ``bin`` directory to the search path for executables (``PATH``).

.. note::

  It is not necessary anymore to add the ``lib`` directory to the search path for shared libraries (``LD_LIBRARY_PATH``).

.. tip::

  To debug issues when Aguila fails to start you can look at the output of the following command::

    ldd path_to_aguila/aguila | less

  This command will report whether certain required shared libraries cannot be found.

Mac OS X
========
Aguila compiles on `Mac OS X`__ too, but since we currently do not have a development machine available with this operating system installed, we cannot provide binary packages, unfortunately.

__ http://www.apple.com/macosx

