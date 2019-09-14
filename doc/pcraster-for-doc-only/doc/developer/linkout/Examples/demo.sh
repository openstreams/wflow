echo incomplete linux sample

set -e
cp      data/*.* .
mkdir -p allStatisticsResults

python LinkOutTestDriver.py allStatistics.xml
python LinkOutTestDriver.py bilStatistics.xml

# the Memory Transfer only example
python test.py
