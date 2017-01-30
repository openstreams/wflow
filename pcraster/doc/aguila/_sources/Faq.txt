.. _faq:

***
FAQ
***

Why can't my dataset be opened?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Did you set PCRASTER_DAL_FORMATS (see section :ref:`environmentVariables`) and is the dataset formatted in a different format than the ones listed by this variable? Solution is to either add the name of this new format to the variable or to unset the variable.

Another possibility is that Aguila currently does not have support for reading the dataset. Aguila contains support for reading many data formats, but not for all data formats. See section :ref:`support` for information about how to request a new feature.

