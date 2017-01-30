************
Introduction
************
The framework supports static models without a time component and temporals models wich simulate a dynamic processes. Both types can be used for Monte Carlo (:ref:`Doucet2001 <Doucet2001>`), particle filter (:ref:`Weerts2006 <Weerts2006>`) or Ensemble Kalman filter (:ref:`Evensen1994 <Evensen1994>`) simulations. With the PCRaster Python extension (:ref:`pcraster2008 <pcraster2008>`, :ref:`Karssenberg2007 <Karssenberg2007>`) models can be extended by a spatial component.

Two steps in the model development cycle are the conversion of the conceptual model structure into computer code and the assimilation or calibration of the model with observational data (:ref:`Karssenberg2006 <Karssenberg2006>`. This framework combines the tasks of model construction and optimisation in a single framework.

The framework does not require to use the PCRaster environmental modelling language even so the example scripts in the following sections solely make use of it.

Requirements
============
To use the framework you need to install `Python 2.7`_ and `NumPy`_. To generate the output graphs for the particle filter you need to install `Graphviz`_.

.. _Python 2.7: http://www.python.org
.. _NumPy: http://numpy.scipy.org
.. _Graphviz: http://www.graphviz.org

Furthermore the framework directory must be added to the `PYTHONPATH` environment variable.

The examples in the demo directory require the PCRaster Python extension.
