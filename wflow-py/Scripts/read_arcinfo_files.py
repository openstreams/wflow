#!/usr/bin/python
#
# Created on August 17, 2000
#  by Keith Cherkauer
#
# This python script contains library functions whic read and write
# standard Arc/Info ASCII grid files, returning a data dictionary
# with information from the file.
#
# Modified 1/16/07 so that the file is read in line by line, rather than all at once.  This
#  reduces memory usage, especially important when running in Windows.

# import string functions
from string import split
from string import atoi
from string import atof


def read_ARCINFO_ASCII_grid(
    filename,
    FLAG="UNFILTERED",
    INTflag=0,
    Extent={"North": 0, "South": 0, "East": 0, "West": 0},
):
    """This routine reads a standard Arc/Info ASCII grid file and 
    returns the coordinates and values in a data dictionary.  An 
    optional flag passed to the routine determines if values for 
    the entire grid are returned (UNFILTERED: default) or if only 
    the cells containing data are returned (FILTERED).  Filtering 
    is based of the NODATA_value defined in the file.  Extent defines
    the boundaries of a smaller sub-grid that will be extracted."""

    try:
	f = open(filename,"r")
    except IOError as E:
	print("ERROR: Unable to open or read the ArcInfo grid file %s" % filename)
        print(E)
	fileTable = { "Ncells" : 0 }
	return ( fileTable )

    fileTable = {
	"ncols"        : atoi(split(f.readline())[1]),
	"nrows"        : atoi(split(f.readline())[1]),
	"xllcorner"    : atof(split(f.readline())[1]),
	"yllcorner"    : atof(split(f.readline())[1]),
	"cellsize"     : atof(split(f.readline())[1]),
	"NODATA_value" : atoi(split(f.readline())[1]),
	}
    print(fileTable)

    if Extent["North"] and Extent["South"] and Extent["East"] and Extent["West"]:
        # redefine grid size for desired extent
        nrows = int((Extent["North"] - Extent["South"]) / fileTable["cellsize"])
        ncols = int((Extent["East"] - Extent["West"]) / fileTable["cellsize"])
    else:
        # define full extent
        Extent["North"] = (
            fileTable["yllcorner"] + fileTable["cellsize"] * fileTable["nrows"]
        )
        Extent["South"] = fileTable["yllcorner"]
        Extent["East"] = (
            fileTable["xllcorner"] + fileTable["cellsize"] * fileTable["ncols"]
        )
        Extent["West"] = fileTable["xllcorner"]
        ncols = fileTable["ncols"]
        nrows = fileTable["nrows"]
    # compute number of cells required
    fileTable["Ncells"] = nrows * ncols
    # allocate memory for all cells
    fileTable["cell"] = [ fileTable["NODATA_value"] ]*fileTable["Ncells"]

    print("North %f, South %f, East %f, West %f" % ( fileTable["yllcorner"] + fileTable["cellsize"]*fileTable["nrows"], fileTable["yllcorner"], fileTable["xllcorner"] + fileTable["cellsize"]*fileTable["ncols"], fileTable["xllcorner"] ))
    print(Extent)
    print(fileTable["nrows"], nrows)
    print(fileTable["ncols"], ncols)

    cellidx = 0
    tmprcnt = 0
    tmpccnt = 0
    for i in range(fileTable["nrows"]):
        # get current line from file
        line = f.readline()
        # compute current latitude
        lat = (
            fileTable["yllcorner"]
            + (fileTable["nrows"] - i) * fileTable["cellsize"]
            - fileTable["cellsize"] / 2.
        )
        if lat >= Extent["South"] and lat <= Extent["North"]:
            tmprcnt = tmprcnt + 1
            # if latitude falls within defined extent, split lines to get column values
            if ( INTflag ): tmpvals = list(map(atoi,split(line)))
            else: tmpvals = list(map(atof,split(line)))
            if ( len(tmpvals) != fileTable["ncols"] ):
                print("ERROR: Number of items in row %i (%i) does not match the number of defined columns (%i)." % (i, len(tmpvals), fileTable["ncols"]))
   	    for j in range(fileTable["ncols"]):
                # compute longitude for current cell
                lng = (
                    fileTable["xllcorner"]
                    + (j) * fileTable["cellsize"]
                    + fileTable["cellsize"] / 2.
                )
                if lng >= Extent["West"] and lng <= Extent["East"]:
                    if tmprcnt == 1:
                        tmpccnt = tmpccnt + 1
                    # if longitude within extent boundaries, store current location
                    try:
                        fileTable["cell"][cellidx] = {
                            "lat": lat,
                            "lng": lng,
                            "value": tmpvals[j],
                        }
                        cellidx = cellidx + 1
                    except IndexError as errstr:
                        # did not allocate enough memory, add additional cells
                        fileTable["cell"] = fileTable["cell"] + [
                            {"lat": lat, "lng": lng, "value": tmpvals[j]}
                        ]
                        cellidx = cellidx + 1
        del(line)

    print("Number of rows filled: %i of %i" % ( tmprcnt, nrows ))
    print("Number of cols filled: %i of %i" % ( tmpccnt, ncols ))
    print("Number of cells filled: %i of %i" % ( cellidx, fileTable["Ncells"] ))
    if tmprcnt != nrows: nrows = tmprcnt
    if tmpccnt != ncols: ncols = tmpccnt
    if cellidx < fileTable["Ncells"]:
        fileTable["cell"] = fileTable["cell"][:cellidx]
    if cellidx != fileTable["Ncells"]:
        fileTable["Ncells"] = cellidx

    if FLAG == "FILTERED":
	if fileTable["NODATA_value"] == 0:
	    fileTable["cell"] = [x for x in fileTable["cell"] if x["value"] != 0]
	else:
	    fileTable["cell"] = [x for x in fileTable["cell"] if x["value"] != -9999]
	fileTable["Ncells"] = len(fileTable["cell"])

        # reset grid boundaries to agree with defined extent
    fileTable["ncols"] = ncols
    fileTable["nrows"] = nrows
    fileTable["xllcorner"] = Extent["West"]
    fileTable["yllcorner"] = Extent["South"]

    f.close()
    return fileTable


def write_ARCINFO_ASCII_grid(filename, gridTable, INTflag=0):
    """This routine writes a standard Arc/Info ASCII grid file values
    of which are stored in a data dictionary."""

    try:
	f = open(filename,"w")
    except IOError as E:
	print("ERROR: Unable to open or write the ArcInfo grid file %s" % filename)
	return ( 0 )

    f.write("ncols\t%i\n" % gridTable["ncols"])
    f.write("nrows\t%i\n" % gridTable["nrows"])
    f.write("xllcorner\t%f\n" % gridTable["xllcorner"])
    f.write("yllcorner\t%f\n" % gridTable["yllcorner"])
    f.write("cellsize\t%f\n" % gridTable["cellsize"])
    f.write("NODATA_value\t%i\n" % gridTable["NODATA_value"])

    idx = 0
    CellsWritten = 0
    for i in range(gridTable["nrows"]):
        lat = (
            gridTable["yllcorner"]
            + (gridTable["nrows"] - i) * gridTable["cellsize"]
            - gridTable["cellsize"] / 2.
        )
        tmpstr = ""
        for j in range(gridTable["ncols"]):
            lng = (
                gridTable["xllcorner"]
                + (j) * gridTable["cellsize"]
                + gridTable["cellsize"] / 2.
            )
            if idx < gridTable["Ncells"] and (
                abs(lat - gridTable["cell"][idx]["lat"]) <= gridTable["cellsize"] / 2.
                and abs(lng - gridTable["cell"][idx]["lng"])
                <= gridTable["cellsize"] / 2.
            ):
                if (
                    INTflag
                    or gridTable["cell"][idx]["value"] == gridTable["NODATA_value"]
                ):
                    tmpstr = tmpstr + "%i " % gridTable["cell"][idx]["value"]
                else:
                    tmpstr = tmpstr + "%f " % gridTable["cell"][idx]["value"]
                if gridTable["cell"][idx]["value"] != gridTable["NODATA_value"]:
                    CellsWritten = CellsWritten + 1
                idx = idx + 1
            else:
                tmpstr = tmpstr + "%i " % gridTable["NODATA_value"]
        f.write(tmpstr[:-1] + "\n")

    f.close()

    print("%s: %i cells of %i contain data." % ( filename, CellsWritten,
                                                 gridTable["ncols"]*gridTable["nrows"] ))

    return 1
