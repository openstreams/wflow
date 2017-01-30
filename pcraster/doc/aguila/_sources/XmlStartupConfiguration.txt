.. _xmlStartupConfiguration:

*************************
Xml Startup Configuration
*************************

While configuration of Aguila is possible with the use of command line options and a configuration file, a more sophisticated scheme is possible by configuring Aguila with an XML start up file. The XML interface is unique in the following ways:

- Specify draw properties per data item.
- Specify real dates (e.g. 2005-02-10T18:15:00) and real time intervals (e.g. 24 hours).
- Map data with different time dimensions to a single time line.

The XML interface is intended as an intermediate interface for other applications. For example, from a Python__ script it is relatively easy to generate the XML file. To use the XML interface a basic understanding of XML is required.

__ http://www.python.org

The XML Schema (``Aguila.xsd``) and samples of XML script files are distributed with a data set in the ``demos/xml`` distribution directory. References to ``.xml`` files in this chapter can be found in that directory. Note that the Aguila application does not require the file ``Aguila.xsd`` to be present. The correct version of ``Aguila.xsd`` is compiled into the Aguila application. For ease of reference, the schema is copied verbatim in section :ref:`aguilaXsd`

Basic structure
===============
In this section we only discuss the overall structure of the XML file. The detailed XML structure is documented in element specification within ``Aguila.xsd``. Each Aguila XML file should have the following structure:

.. highlight:: xml

::

  <?xml version="1.0" encoding="UTF-8"?>
  <aguila
     xmlns="http://www.pcraster.nl/pcrxml"
     xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
     xsi:schemaLocation="http://www.pcraster.nl/pcrxml Aguila.xsd">
    <visualisationGroup>
     <!-- contents consisting of elements:
          - cursorValueMonitorFile  0 or 1
          - searchSpace             0 or 1
          - data                    0 or more
          - view                    0 or more
       -->
    </visualisationGroup>
  </aguila>

The elements of the basic structure will be discussed in the reverse order of appearance.

``view``
  A number of views can be specified. All view contents is specified by the dataset name as an ``item`` element. That same name may optionally appear as a data element. ``Example1.xml`` illustrates a number of points on the view specification:

  - The first appearance yields the position in the cursor and values matrix. Use a ``valueOnly`` element to obtain the preferred order (``dem.map``).
  - Time series file and cursor value graphs can be combined. See ``runoff`` and ``runoff.tss`` combined in a view in ``example1.xml``.

``data``
  A ``data`` element is used to tag a dataset name used in a ``view`` element with additional specifications: its draw properties and ``dataSpace``. For example, in ``example1.xml`` the colour assignment of ``upstreamArea.map`` is set to ``shifted logarithmic``.

``searchSpace``
  With the ``searchSpace`` element one can describe the dimensions in which all data should be found if the data itself has no ``dataSpace`` sub element. In ``example1.xml`` the 28 model steps are set in the ``searchSpace`` and the timesteps are mapped to a real calendar date with a 6 hour time increment.

``cursorValueMonitorFile``
  See the section called :ref:`programOptions`.

``fileToGetCursorValue``
  This element defines what file is read to set the applications cursor at each time ``Get`` is pressed in the ``Cursor Value Window``. Currently only the x and y (space dimensions) are read. The file should be an XML file with the ``aguilaCursor`` element as root. See ``get.xml`` set in ``example1.xml``.

.. _aguilaXsd:

Aguila.xsd
==========

.. literalinclude:: Aguila.xsd
   :language: xml
