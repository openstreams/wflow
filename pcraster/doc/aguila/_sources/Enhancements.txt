.. _enhancements:

************
Enhancements
************

Some ideas for improving Aguila. Improvements can be new features or redesigns of existing features.

Attribute browser
=================
Attributes are stored in datasets. In a dataset one or more attributes can be stored. A dataset can be stored in many different ways: in a single file, in multiple files, in a database, by a WCS service, etc, etc.

Currently, users needs to know the names of the attributes, which are not always obvious. For example, attributes in a Shapefile are hidden in layers, which are hidden in Shapefiles. The user needs to use tools like gdalinfo to find out the names of the attributes.

Aguila needs an ``Open Attribute`` 'dialog' which is smart enough to find attributes in a certain location, whatever the storage mechanism. For this to work Dal needs to be updated. Each format driver needs some ``browse`` functionality that, given a location, returns information about the available attributes.

Dal's enhancement doc contains info on Dal's side of the implementation. Here, we only list some considerations for the user interface.

Requirements
------------
- User must be able to select different locations:

  - Directories of the local filesystem.
  - Url's (WFS, WCS).

- It must be obvious what the characteristics of the found attributes are (spatial, temporal, uncertain, scenarios).
- Found attributes may have the same name. It must be obvious which dataset each attribute is part of.
- Since browsing for attributes may be slow:

  - There must be some progress indicator.
  - The user must be able to cancel (or stop) the browse.
  - Browse results must be cached and there must be an option to force a rebrowse.

- It must be possible to drag found attributes onto running visualisations.
- The rest of the application must remain usable while the dialog scans for attributes.
- It must be possible to open more than one attribute browser at the same time.
- It must be possible to start a new visualisation from the context menu of the attribute. Only valid visualisations should be listed in the menu.
- Found attributes must be added to the window incrementally, not after the browse is finished.
- The attribute browser must be usable by other clients as well (Nutshell, Script editor, ...). It may even be useful on its own, as a smarter replacement of gdalinfo. It can be developed as an independent component and added to Aguila once it becomes usable. It probably should only depend on Dal, not on Aguila. Aguila depends on the attribute browser.
- The most recently visited locations must be available for selection within and between sessions.
- The attribute view mode must be selectable (icons, details, ...).

Design
------
From the requirements folows that the folowing graphical elements are part of the dialog:

- The location being browsed.
- A means to change the browse location, eg:

  - Push button which starts a dialog for selecting a directory in the local filesystem.
  - Push button which starts a dialog for entering urls of WCS and WFS servers.

- A push button to start the browse. In case chached data is used, this button may double as a 'Refresh' button.
- A push button to stop the browse.
- A table with information about the attributes found:

  - Name
  - Dataset
  - Dimensional properties
  - Format

- A progress indicator. Probably a progressbar. Amount of work is equal to number of drivers used. Current progress equals the number of drivers finished.

Implemenation
-------------
- The dialog must run in a seperate thread.
- The attribute browser is a Widget which can be put in a dialog or in any other widget type.
- When rebrowsing the a location another time, use the driver information of the previous browse to see whether that data can still be opened with it. Most probably it can. Drivers are tried in a sequence, and can be ordered based on previous experiences. Try to get a small set of files ASAP by using existing knowledge.

Signals for high-level browse class::

  attributeClicked(size_t id);
  attributeDoubleClicked(size_t id);
  customContextMenuRequested(size_t id);

2.5D View
=========
Current implementation makes use of OpenGL, which is hard to maintain if you don't use its API regularly. Switch to some higher level library (Ogre, ...).


