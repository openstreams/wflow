import csv # import the python csv library 
from Libraries.StandardFunctions import *
from Libraries.ChartFunctions import *
import WflowDeltashell.wflib as wfl
import os

# Needed if this .net thing is not loaded yet
import clr
clr.AddReference("System.Windows.Forms")
from System.Windows.Forms import OpenFileDialog, DialogResult


def getcsvname():
	themodel = wfl.GetModelByPartialName('wflow')
		
	dialog = OpenFileDialog()
	
	if themodel:
		dialog.InitialDirectory = os.path.join(themodel.DirectoryPath,themodel.DefaultOutputDirectoryName)
	else:
		dialog.InitialDirectory ="C:\\"
	 
	dialog.Filter = "csv files (*.csv) | *.csv"
	dialog.FilterIndex = 1
	dialog.RestoreDirectory = False
	dialog.Title = "Select a WFlow result csv file: "
	
	if (dialog.ShowDialog() == DialogResult.OK):
		thefile = dialog.FileName
		
	return thefile


def plotit(filename):
	colors = [Color.Red,Color.Black, Color.Green, Color.Black, Color.Orange, Color.YellowGreen, Color.Azure,\
			Color.Cyan, Color.Brown, Color.Aquamarine, Color.CornflowerBlue, Color.DarkBlue, Color.BlanchedAlmond,\
			Color.LightYellow, Color.SlateGray]
			
	runid = os.path.basename(os.path.dirname(filename))
	
	lineSeries = []
	
	olines = wfl.readwfcsv(filename)

	
	for yval in range(1,len(olines[0])): 
		x = [[col[0],col[yval]]  for col in olines]
	
		lineSeries.append(CreateLineSeries(x))
	
	
		lineSeries[-1].Color = colors[yval -1]
		lineSeries[-1].Width = 1
		lineSeries[-1].PointerVisible = False
		lineSeries[-1].PointerSize = 0
		lineSeries[-1].PointerColor = Color.Red
		lineSeries[-1].PointerLineVisible = True
		lineSeries[-1].PointerLineColor = colors[yval -1]
		lineSeries[-1].Transparency = 0
		lineSeries[-1].ShowInLegend= True
		lineSeries[-1].Title="Id: " + str(yval)
	
	chart = CreateChart(lineSeries)
	
	# Configure the chart
	chart.TitleVisible = True
	chart.Title = "run: " + runid + ": " + os.path.basename(filename)
	chart.BackGroundColor = Color.White
	chart.Legend.Visible = True
	chart.Name = chart.Title
	
	
	# Configure the bottom axis
	chart.BottomAxis.Automatic = True
	chart.BottomAxis.Title = "Timestep"
	
	# Configure the left axis
	chart.LeftAxis.Title = os.path.basename(filename)
	
	
	# Show the chart
	view = OpenView(chart)
	


def complot(dirnames,file,ids):
	colors = [Color.Red,Color.Black, Color.Green, Color.Black, Color.Orange, Color.YellowGreen, Color.Azure,\
			Color.Cyan, Color.Brown, Color.Aquamarine, Color.CornflowerBlue, Color.DarkBlue, Color.BlanchedAlmond,\
			Color.LightYellow, Color.SlateGray]
			

	cnt = 0
	lineSeries = []
	for dirname in dirnames:
		runid = os.path.dirname(dirname)
		
		ftoread = os.path.join(dirname,file)
		
		if os.path.exists(ftoread):
			olines = wfl.readwfcsv(os.path.join(dirname,file))
			
			for yval in ids:
				x = [[col[0],col[yval]]  for col in olines]
			
				lineSeries.append(CreateLineSeries(x))
			
			
				lineSeries[-1].Color = colors[cnt]
				lineSeries[-1].Width = 1
				lineSeries[-1].PointerVisible = False
				lineSeries[-1].PointerSize = 0
				lineSeries[-1].PointerColor = Color.Red
				lineSeries[-1].PointerLineVisible = True
				lineSeries[-1].PointerLineColor = colors[cnt]
				lineSeries[-1].Transparency = 0
				lineSeries[-1].ShowInLegend= True
				lineSeries[-1].Title="Col: " + str(yval) + " Run: " + str(os.path.basename(runid))
				cnt = cnt + 1
	
	casename = os.path.basename(os.path.dirname(os.path.dirname(dirname)))
	chart = CreateChart(lineSeries)
	
	# Configure the chart
	chart.TitleVisible = True
	chart.Title = casename.upper() + ":" + os.path.basename(file)
	chart.BackGroundColor = Color.White
	chart.Legend.Visible = True
	chart.Name = chart.Title
	
	
	# Configure the bottom axis
	chart.BottomAxis.Automatic = True
	chart.BottomAxis.Title = "Timestep"
	
	# Configure the left axis
	chart.LeftAxis.Title = os.path.basename(file)
	
	
	# Show the chart
	view = OpenView(chart)
	
