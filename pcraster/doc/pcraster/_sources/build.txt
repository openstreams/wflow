Build PCRaster
==============

.. note::

   Building PCRaster yourself is something you want to prevent. It may take a lot of time and you may end up being frustrated. Are you sure you want to proceed? Is it really necessary?

The PCRaster software depends on many other software. For example, Qt is used in the implementation of GUI components, and GDAL is used to read various raster formats. Before being able to build PCRaster, these 3rd party software sources must be downloaded and built. After that, the PCRaster sources can be downloaded and built.


Check out all PCRaster project sources
--------------------------------------
This step is not necessary if you already have the source code. Download and execute the `clone_pcraster_sources.sh`_ script.

.. _clone_pcraster_sources.sh: https://sourceforge.net/p/pcraster/pcraster/ci/master/tree/environment/scripts/clone_pcraster_sources.sh

.. code-block:: bash

   cd $HOME/tmp
   wget https://sourceforge.net/p/pcraster/pcraster/ci/master/tree/environment/scripts/clone_pcraster_sources.sh?format=raw
   mv clone_pcraster_sources.sh\?format\=raw clone_pcraster_sources.sh
   cd $HOME
   mkdir -p development/{projects,objects}
   bash ~/tmp/clone_pcraster_sources.sh development/projects


Build and install 3rd party software
------------------------------------
Determine where to put the 3rd party stuff. Create this directory. This is the root directory of the directory that will be created when building the 3rd party stuff. Examples are: ``/home/pcrteam_extern``, ``$HOME/pcrteam_extern``, ``/mnt/pcrteam_extern``, ``C:\PcrTeamExtern``. The actual path name is not relevant.

Before executing the ``build_pcrteam_extern.sh`` script mentioned below, you may want to set ``CC``/``CXX`` to point to the correct compilers.

The ``build_pcrteam_extern.sh`` script depends on certain tools to be installed. On Debian-based systems you can run ``.../devenv/scripts/machine_status.py`` to tell you what you still need to install, if anything.

.. code-block:: bash

   # Create a temp location with lots of disk space.
   mkdir $HOME/tmp/3rd

   # Create a location where to install 3rd party software.
   mkdir $HOME/pcrteam_extern

   # Download/build/install 3rd party software.
   $HOME/development/projects/devenv/scripts/build_pcrteam_extern.sh $HOME/tmp/3rd $HOME/pcrteam_extern/master-`date +"%Y%m%d"`

   # Create a symbolic link that can be updated should a new version of
   # pcrteam_extern be installed.
   cd $HOME/pcrteam_extern
   ln -s master-`date +"%Y%m%d"` master


Build PCRaster
--------------
Build PCRaster projects
~~~~~~~~~~~~~~~~~~~~~~~
Before building the PCRaster projects, you need to set some environment variables.

.. code-block:: bash

   unset OBJECTS

   export DEVELOPMENT_ROOT=$HOME/development
   export PROJECTS=$DEVELOPMENT_ROOT/projects
   export PCRTEAM_EXTERN_ROOT=$HOME/pcrteam_extern

   source "$PROJECTS/devenv/configuration/profiles/Utils.sh"

Now, your environment is setup to build any of the projects that comprise PCRaster. Next you need to configure your environment for the specific project you want to build. There are bash scripts that can be sourced that do that for you. You can configure aliases to make sourcing these scripts easy:

.. code-block:: bash

   alias PCRaster-master="source $PROJECTS/pcraster/environment/configuration/PCRaster-master"

Now you can just type the name of the project (CamelCased) appended by the branch name, and optionally followed by the build type (release or debug, default is debug): ``<project_name>-<branch_name> [build_type]``.

To build everything, choose the ``PCRaster`` project. If you need to set ``CC``/``CXX`` explicitly, then do that *before* configuring the environment for the project.

.. code-block:: bash

   PCRaster-master
   rebuild_projects.py


Create PCRaster package
~~~~~~~~~~~~~~~~~~~~~~~
To create a PCRaster package for distribution you can use the ``make_pcraster_package.sh`` shell script from ``$PCRASTER/environment/scripts``. It builds the software in the current directory and creates a PCRaster package that can be copied to other locations. On Linux:

.. code-block:: bash

   PCRaster-master release
   cd $HOME/tmp
   make_pcraster_package.sh
