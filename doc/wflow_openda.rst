Linking wflow to OpenDA
=======================

.. warning:: 

	This is experimental and incomplete


Prerequisites
-------------
+ download the Wflow openda kernel binary release.
+ install openda


Configuration
-------------

::


  <?xml version="1.0" encoding="UTF-8"?>
  <bmiModelFactoryConfig xmlns="http://www.openda.org" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.openda.org bmiModelFactoryConfig.xsd">
      <pythonModel>
          <pythonPath>../../../wflow_bin</pythonPath>
          <moduleName>wflow.wflow_bmi</moduleName>
          <className>wflowbmi_csdms</className>
          <!-- You must give an absolute path to the py.exe of the binary distribution -->
          <pythonExecutable>d:\2015.02\GLOFFIS_SA\Modules\OpenStreams\wflow_bin\py.exe</pythonExecutable>
      </pythonModel>
      <modelTemplateDirectory>../../combined</modelTemplateDirectory>
      <modelConfigFile>wflow_wr3a.ini</modelConfigFile>
     <bmiModelForcingsConfig>
          <dataObject>
              <className>org.openda.exchange.dataobjects.NetcdfDataObject</className>
              <file>Precipitationmapstack.nc</file>
              <arg>true</arg>
              <arg>false</arg>
          </dataObject>
      </bmiModelForcingsConfig>
  </bmiModelFactoryConfig>