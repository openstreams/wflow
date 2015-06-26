# -*- coding: utf-8 -*-
"""
******************
pcraster.framework
******************
The ``pcraster.framework`` package contains classes for various model types, and framework types.

Model types are meant to be overriden by the user. By implementing these types, the user is certain that his model fulfills the requirements of the framework type used to execute the model. The goal is to enable the model to be as small as possible, containing hardly more than the statements making up the actual model.

Supported model types are:

* Static model
* Dynamic model
* Monte Carlo model
* Particle filter model
* Kalman filter model

Framework types are provided by the package to execute the model created by the user. They contain the boilerplate logic to call certain model methods in the right order. Apart from that, they contain some additional goodies, like a progress indicator.

Supported framework types are:

* Static model framework
* Dynamic model framework
* Monte Carlo framework
* Particle filter framework
* Kalman filter framework

Please find more detailed information in these next sections:

* :ref:`Static models <staticModels>`
* :ref:`Dynamic models <dynamicModels>`
* :ref:`Monte Carlo method <monteCarloMethod>`
* :ref:`Particle filter method <particleFilterMethod>`
* :ref:`Kalman filter method <kalmanFilterMethod>`

.. todo::

  Can we replace the term `static model` with something better? Same for dynamic.

.. todo::

  Get rid of the dependency on the PCRaster Python extension. At runtime, we can test whether the extension is installed. If not, than we can turn off support for certain features.

.. todo::

  Can we get rid of the support for "either 'run' or 'initial' and 'dynamic'" in the Dynamic Model Concept?

.. _staticModels:

Static models
=============

.. _staticModelConcept:

Static model concept
--------------------

.. todo::

  Document.

To use the static framework the user must implement:

* either a run or an initial method

Classes
-------
* :py:class:`pcraster.framework.staticBase.StaticBase`
* :py:class:`pcraster.framework.staticPCRasterBase.StaticModel`
* :py:class:`pcraster.framework.staticFramework.StaticFramework`

.. automodule:: pcraster.framework.staticBase
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.staticPCRasterBase
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.staticFramework
   :members:
   :undoc-members:
   :show-inheritance:

.. _dynamicModels:

Dynamic models
==============

.. _dynamicModelConcept:

Dynamic model concept
---------------------

.. todo::

  Document.

To use the dynamic framework the user must implement the following methods in this class:

* either "run" or "initial" and "dynamic"

Classes
-------
.. automodule:: pcraster.framework.dynamicBase
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.dynamicPCRasterBase
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.dynamicFramework
   :members:
   :undoc-members:
   :show-inheritance:

.. _monteCarloMethod:

Monte Carlo method
==================

.. _monteCarloModelConcept:

Monte Carlo model concept
-------------------------

.. todo::

  Document.

* premcloop
* postmcloop

Classes
-------
.. automodule:: pcraster.framework.mcBase
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.mcPCRasterBase
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.mcFramework
   :members:
   :undoc-members:
   :show-inheritance:

.. _particleFilterMethod:

Particle filter method
======================

.. _particleFilterModelConcept:

Particle filter model concept
-----------------------------

.. todo::

  Document.

Classes
-------
.. automodule:: pcraster.framework.pfBase
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.pfPCRasterBase
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.particleFilterFramework
   :members:
   :undoc-members:
   :show-inheritance:

.. _kalmanFilterMethod:

Kalman filter method
====================
Classes
-------
.. automodule:: pcraster.framework.kfBase
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.kfPCRasterBase
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.kalmanFilterFramework
   :members:
   :undoc-members:
   :show-inheritance:

Stuff
=====
.. automodule:: pcraster.framework.frameworkBase
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.generalfunctions
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.aggregationfunctions
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pcraster.framework.Timeoutput
   :members:
   :undoc-members:
   :show-inheritance:
"""

from staticPCRasterBase import *
from dynamicPCRasterBase import *
from mcPCRasterBase import *
from pfPCRasterBase import *
from kfPCRasterBase import *


from staticFramework import *
from dynamicFramework import *
from mcFramework import *
from particleFilterFramework import *
from kalmanFilterFramework import *
from generalfunctions import *
from aggregationfunctions import *
from Timeoutput import *
# for backwards compatibility: generateNameXY
from frameworkBase import *
