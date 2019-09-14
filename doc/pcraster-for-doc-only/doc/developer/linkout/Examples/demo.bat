REM add path to find PCRasterModelEngine.dll and depends
PATH=..\..\..\..\bin;%PATH%

copy  data\*.*
mkdir allStatisticsResults

REM run dynamicModel.xml with memory Only transfer
python test.py

REM run other samples
python LinkOutTestDriver.py allStatistics.xml
python LinkOutTestDriver.py bilStatistics.xml
python LinkOutTestDriver.py lookup.xml
python LinkOutTestDriver.py staticScript.xml
python LinkOutTestDriver.py statisticsMask.xml
python LinkOutTestDriver.py compression.xml
python LinkOutTestDriver.py areaMap.xml

REM the ones with errors
python LinkOutTestDriver.py notValid.xml
python LinkOutTestDriver.py compressionNoAreaMap.xml
python LinkOutTestDriver.py dynamicNoTimer.xml
python LinkOutTestDriver.py memoryOnlyIO_1.xml
python LinkOutTestDriver.py memoryOnlyIO_2.xml
python LinkOutTestDriver.py memoryOnlyIO_3.xml
python LinkOutTestDriver.py memoryOnlyIO_4.xml
python LinkOutTestDriver.py memoryOnlyIO_5.xml
python LinkOutTestDriver.py memoryOnlyIO_6.xml

REM done in csharp unit test: memoryOnlyIO_7.xml
