from WflowDeltashell.plotcsv import *
from WflowDeltashell.wflib import *
# Needed if this .net thing is not loaded yet
import clr
clr.AddReference("System.Windows.Forms")
from System.Windows.Forms import OpenFileDialog, DialogResult

one = 'c:\\repos\wflow-git\\examples\\wflow_rhine_sbm\\pmult\\store.csv'
two = 'c:\\repos\wflow-git\\examples\\wflow_rhine_sbm\\run_default\\store.csv'



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
	
casename = os.path.dirname(os.path.dirname(thefile))
csvfile = os.path.basename(thefile)

print casename

runs = getrunids(casename)

print runs

complot(runs,csvfile,[2])