import csv # import the python csv library 
from Libraries.StandardFunctions import *
from Libraries.ChartFunctions import *

testfilepath = 'c:\\repos\\wflow-git\\examples\\wflow_rhine_sbm\\run_default\\prec.csv'

column = 4


#Read a CSV file
with open(testfilepath) as csvfile: # open the file test2.csv
    lines = csv.reader(csvfile, delimiter=',') # read lines as collection of arrays

    olines=[]
    # print all values
    for line in lines:
        oline = []
        thecol = 0
        for field in line:
        	thecol = thecol + 1
        	if thecol == 1:
        		oline.append(float(field))
        	if thecol == column:
        		oline.append(float(field))
            	
        olines.append(oline)

print olines


lineSeries = CreateLineSeries(olines)
# Configure the line series
lineSeries.Color = Color.Red
lineSeries.Width = 1
lineSeries.PointerVisible = True
lineSeries.PointerSize = 5
lineSeries.PointerColor = Color.Red
lineSeries.PointerLineVisible = True
lineSeries.PointerLineColor = Color.DarkRed
lineSeries.Transparency = 50

chart = CreateChart([lineSeries])
# Show the chart
OpenView(chart)