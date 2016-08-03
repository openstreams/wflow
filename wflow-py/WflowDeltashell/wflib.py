# Library/utility functions for wflow/Deltashell
#Read a CSV file
import csv # import the python csv library 
import os
import glob


def GetItemByPartialName(list, name):
    """ Returns the first item in the list
        that has the provided name"""
    for item in list :
        if  name.upper() in item.Name.upper(): 
            return item
            
def GetModelByPartialName(modelName):
    """Searches for a model with the provided name"""
    return GetItemByPartialName(Application.ModelService.GetAllModels(RootFolder), modelName)
    

def readwfcsv(fname):
    """
    read csv file in a list of lists
    
    """
    with open(fname) as csvfile: # open the file test2.csv
        lines = csv.reader(csvfile, delimiter=',') # read lines as collection of arrays
    
        olines=[]
        # print all values
        for line in lines:
            oline = []
            thecol = 0
            if '#' not in line[0]:
                for field in line:              	
	                oline.append(float(field))
            	olines.append(oline)
            
    return olines


def getrunids(casedir):
	"""
	Ugly method to get the run ids. This is absolutely not failsave
	"""
	from glob import glob
	dirs =  glob(casedir + "/*/")
	
	ret = []
	
	for dir in dirs:
		dn = os.path.basename(os.path.dirname(dir))
		if dn not in "intbl staticmaps inmaps instate intss outstate":
			ret.append(dir)
			print dn
			
	return ret			
	