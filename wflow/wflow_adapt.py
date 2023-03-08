# Test version of wflow Delft-FEWS adapter
#
# Wflow is Free software, see below:
#
# Copyright (c) J. Schellekens 2005-2011
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
wflow_adapt.py: Simple wflow Delft-FEWS adapter in python. This file can be run
as a script from the command-line or be used as a module that provides (limited)
functionality for converting PI-XML files to .tss and back.

*Usage pre adapter:*

**wflow_adapt** -M Pre -t InputTimeseriesXml -I inifile

*Usage postadapter:*
    
**wflow_adapt**-M Post -t InputTimeseriesXml -s inputStateFile -I inifile 
              -o outputStateFile -r runinfofile -w workdir -C case [-R runId]

Issues:
    
- Delft-Fews exports data from 0 to timestep. PCraster starts to count at 1.
  Renaming the files is not desireable. The solution is the add a delay of 1 
  timestep in the GA run that exports the mapstacks to wflow.
- Not tested very well.
- There is a considerable amount of duplication (e.g. info in the runinfo.xml and
  the .ini file that you need to specify again :-())

 .. todo::

     rewrite and simplify

$Author: schelle $
$Id: wflow_adapt.py 915 2014-02-10 07:33:56Z schelle $
$Rev: 915 $
     
"""

import configparser
import getopt
import logging
import logging.handlers
import os
import shutil
import sys
from datetime import *
from xml.etree.ElementTree import *

import numpy
import wflow.pcrut as pcrut
import wflow.wflow_lib as wflow_lib

outMaps = ["run.xml", "lev.xml"]
iniFile = "wflow_sbm.ini"
case = "not_set"
runId = "run_default"

logfile = "wflow_adapt.log"


def make_uniek(seq, idfun=None):
    # Order preserving
    return list(_f10(seq, idfun))


def _f10(seq, idfun=None):
    seen = set()
    if idfun is None:
        for x in seq:
            if x in seen:
                continue
            seen.add(x)
            yield x
    else:
        for x in seq:
            x = idfun(x)
            if x in seen:
                continue
            seen.add(x)
            yield x


fewsNamespace = "http://www.wldelft.nl/fews/PI"


def setlogger(logfilename, loggername, thelevel=logging.INFO):
    """
    Set-up the logging system and return a logger object. Exit if this fails
    """

    try:
        # create logger
        logger = logging.getLogger(loggername)
        if not isinstance(thelevel, int):
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(thelevel)
        ch = logging.FileHandler(logfilename, mode="w")
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        # create formatter
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(module)s - %(levelname)s - %(message)s"
        )
        # add formatter to ch
        ch.setFormatter(formatter)
        console.setFormatter(formatter)
        # add ch to logger
        logger.addHandler(ch)
        logger.addHandler(console)
        logger.debug("File logging to " + logfilename)
        return logger
    except IOError:
        print("ERROR: Failed to initialize logger with logfile: " + logfilename)
        sys.exit(2)


def log2xml(logfile, xmldiag):
    """
    Converts a wflow log file to a Delft-Fews XML diag file

    """

    trans = {"WARNING": "2", "ERROR": "1", "INFO": "3", "DEBUG": "4"}
    if os.path.exists(logfile):
        with open(logfile, "r") as fi:
            lines = fi.readlines()

        with open(xmldiag, "w") as fo:
            fo.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            fo.write('<Diag xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" \n')
            fo.write(
                'xmlns="http://www.wldelft.nl/fews/PI" xsi:schemaLocation="http://www.wldelft.nl/fews/PI \n'
            )
            fo.write(
                'http://fews.wldelft.nl/schemas/version1.0/pi-schemas/pi_diag.xsd" version="1.2">\n'
            )
            for aline in lines:
                translator = aline.maketrans("", "", "><&\"'")
                alineesc = aline.translate(translator)
                fo.write(alineesc)
            fo.write("</Diag>\n")


def pixml_state_updateTime(inxml, outxml, DT):
    """
    Reads the pi-state xml file inxml and updates the data/time of
    the state using datetime. Writes updated file to outxml

    - Can be use in scripts to set the date.time of the
      output state.xml that Delft-FEWS writes.

    .. warning::

            This function does not fully parse the xml file and will only work properly
        if the xml files date the dateTime element written on one line.

    """

    if os.path.exists(inxml):
        datestr = DT.strftime("%Y-%m-%d")
        timestr = DT.strftime("%H:%M:%S")

        with open(inxml, "r") as fi:
            lines = fi.readlines()

        with open(outxml, "w") as fo:
            for aline in lines:
                pos = aline.find("dateTime")
                if pos >= 0:
                    fo.write(
                        '<dateTime date="' + datestr + '" time="' + timestr + '"/>\n'
                    )
                else:
                    fo.write(aline)

    else:
        print((inxml + " does not exists."))


def pixml_totss_dates(nname, outputdir):
    """
    Gets Date/time info from XML file and creates .tss files with:

        - Day of year
        - Hour of day
        - Others may follow

    """

    if os.path.exists(nname):
        with open(nname, "r") as f:
            tree = parse(f)

        PItimeSeries = tree.getroot()
        series = PItimeSeries.findall(".//{" + fewsNamespace + "}series")

        events = series[0].findall(".//{" + fewsNamespace + "}event")
        with open(outputdir + "/YearDay.tss", "w") as f:
            with open(outputdir + "/Hour.tss", "w") as ff:
                # write the header
                f.write("Parameter YearDay taken from " + nname + "\n")
                ff.write("Parameter Hour taken from " + nname + "\n")
                f.write("2\n")
                ff.write("2\n")
                for i in range(1, 3):
                    f.write("Data column " + str(i) + "\n")
                    ff.write("Data column " + str(i) + "\n")
                i = 1
                for ev in events:
                    dt = datetime.strptime(
                        ev.attrib["date"] + ev.attrib["time"], "%Y-%m-%d%H:%M:%S"
                    )
                    f.write(str(i) + "\t" + dt.strftime("%j\n"))
                    ff.write(str(i) + "\t" + dt.strftime("%H\n"))
                    i += 1
    else:
        print((nname + " does not exists."))


def pixml_totss(nname, outputdir):
    """
    Converts and PI xml timeseries file to a number of tss files.

    The tss files are created using the following rules:

        - tss filename determined by the content of the parameter element with a ".tss" postfix
        - files are created in "outputdir"
        - multiple locations will be multiple columns in the tss file written in order
          of appearance in the XML file

    """

    if os.path.exists(nname):
        with open(nname, "r") as f:
            tree = parse(f)

        PItimeSeries = tree.getroot()
        seriesStationList = PItimeSeries.findall(
            ".//{" + fewsNamespace + "}stationName"
        )
        LocList = []
        for station in seriesStationList:
            LocList.append(station.text)

        Parameters = PItimeSeries.findall(".//{" + fewsNamespace + "}parameterId")
        ParList = []
        for par in Parameters:
            ParList.append(par.text)

        uniqueParList = make_uniek(ParList)

        colsinfile = len(ParList)

        series = PItimeSeries.findall(".//{" + fewsNamespace + "}series")

        # put whole lot in a dictionary
        val = {}
        parlocs = {}
        i = 0
        for par in uniqueParList:
            parlocs[par] = 1

        for thisS in series:
            par = thisS.find(".//{" + fewsNamespace + "}parameterId").text
            events = thisS.findall(".//{" + fewsNamespace + "}event")
            locs = thisS.findall(".//{" + fewsNamespace + "}locationId")

            i = 0
            for ev in events:
                parlocs[par] = 1
                if (i, par) in val:
                    theval = val[i, par] + "\t" + ev.attrib["value"]
                    val[i, par] = theval
                    parlocs[par] = parlocs[par] + 1
                else:
                    val[i, par] = ev.attrib["value"]
                i += 1
        nrevents = i

        for par in uniqueParList:
            with open(outputdir + "/" + par + ".tss", "w") as f:
                # write the header
                f.write("Parameter " + par + " taken from " + nname + "\n")
                f.write(str(colsinfile + 1) + "\n")
                f.write("Timestep\n")
                for i in range(0, colsinfile):
                    f.write("Data column " + str(i) + "\n")
                for i in range(0, nrevents):
                    f.write(str(i + 1) + "\t" + val[i, par] + "\n")

    else:
        print((nname + " does not exists."))


def tss_topixml(tssfile, xmlfile, locationname, parametername, Sdate, timestep):
    """
    Converts a .tss file to a PI-xml file

    """
    missval = "-999.0"

    # try:
    tss, header = pcrut.readtss(tssfile)

    # except:
    #    logger.error("Tss file not found or corrupt: ", tssfile)
    #    return

    # Add dummpy first timesteps
    if len(tss.shape) > 1:
        dumm = tss[0, :].copy()
        dumm[:] = -999.0
        tss = numpy.vstack((dumm, tss))
    else:
        dumm = tss.copy()
        dumm[:] = -999.0
        tss = numpy.vstack((dumm, tss))

    # replace np.nan with missing values
    tss[numpy.isnan(tss)] = missval

    trange = timedelta(seconds=timestep * (tss.shape[0]))

    extraday = timedelta(seconds=timestep)
    # extraday = timedelta(seconds=0)
    # Sdate = Sdate + extraday
    # Edate = Sdate + trange - extraday
    Sdate = Sdate + extraday
    Edate = Sdate + trange - extraday - extraday

    Sdatestr = Sdate.strftime("%Y-%m-%d")
    Stimestr = Sdate.strftime("%H:%M:%S")

    Edatestr = Edate.strftime("%Y-%m-%d")
    Etimestr = Edate.strftime("%H:%M:%S")
    with open(xmlfile, "w") as fo:
        fo.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fo.write(
            '<TimeSeries xmlns="http://www.wldelft.nl/fews/PI" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.wldelft.nl/fews/PI http://fews.wldelft.nl/schemas/version1.0/pi-schemas/pi_timeseries.xsd" version="1.2">\n'
        )
        fo.write("<timeZone>0.0</timeZone>\n")
        count = 0

        for col in tss.transpose():
            count = count + 1
            fo.write("<series>\n")
            fo.write("<header>\n")
            fo.write("<type>instantaneous</type>\n")
            fo.write("<locationId>" + header[count - 1] + "</locationId>\n")
            fo.write("<parameterId>" + parametername + "</parameterId>\n")
            fo.write('<timeStep unit="second" multiplier="' + str(timestep) + '"/>\n')
            fo.write('<startDate date="' + Sdatestr + '" time="' + Stimestr + '"/>\n')
            fo.write('<endDate date="' + Edatestr + '" time="' + Etimestr + '"/>\n')
            fo.write("<missVal>" + str(missval) + "</missVal>\n")
            fo.write("<stationName>" + header[count - 1] + "</stationName>\n")
            fo.write("</header>\n")
            # add data here
            xdate = Sdate
            xcount = 1
            for pt in col:
                if xcount > 1:
                    Ndatestr = xdate.strftime("%Y-%m-%d")
                    Ntimestr = xdate.strftime("%H:%M:%S")
                    fo.write(
                        '<event date="'
                        + Ndatestr
                        + '" time="'
                        + Ntimestr
                        + '" value="'
                        + str(pt)
                        + '" />\n'
                    )
                    xdate = xdate + timedelta(seconds=timestep)
                xcount = xcount + 1
            fo.write("</series>\n")

        fo.write("</TimeSeries>\n")

    return tss


def mapstackxml(
    mapstackxml, mapstackname, locationname, parametername, Sdate, Edate, timestepsecs
):
    """
    writes a mapstack xml file
    """
    Sdatestr = Sdate.strftime("%Y-%m-%d")
    Stimestr = Sdate.strftime("%H:%M:%S")
    Edatestr = Edate.strftime("%Y-%m-%d")
    Etimestr = Edate.strftime("%H:%M:%S")
    with open(mapstackxml, "w") as fo:
        fo.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fo.write(
            '<MapStacks xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns="http://www.wldelft.nl/fews/PI" xsi:schemaLocation="http://www.wldelft.nl/fews/PI http://fews.wldelft.nl/schemas/version1.0/pi-schemas/pi_mapstacks.xsd" version="1.2">\n'
        )
        fo.write("<geoDatum>WGS 1984</geoDatum>\n")
        fo.write("<timeZone>0.0</timeZone>\n")
        fo.write("<mapStack>\n")
        fo.write("<locationId>" + locationname + "</locationId>\n")
        fo.write("<parameterId>" + parametername + "</parameterId>\n")
        fo.write('<timeStep unit="second" multiplier="' + str(timestepsecs) + '"/>\n')
        fo.write('<startDate date="' + Sdatestr + '" time="' + Stimestr + '"/>\n')
        fo.write('<endDate date="' + Edatestr + '" time="' + Etimestr + '"/>\n')
        fo.write("<file>\n")
        fo.write('    <pcrgrid file="' + mapstackname + '"/>\n')
        fo.write("</file>\n")
        fo.write("</mapStack>\n")
        fo.write("</MapStacks>\n")


def getTimeStepsfromRuninfo(xmlfile, timestepsecs):
    """
    Gets the number of  timesteps from the FEWS runinfo file.
    """
    if os.path.exists(xmlfile):
        with open(xmlfile, "r") as f:
            tree = parse(f)

        runinf = tree.getroot()
        sdate = runinf.find(".//{" + fewsNamespace + "}startDateTime")
        ttime = sdate.attrib["time"]
        if len(ttime) == 12:  # Hack for milliseconds in testrunner runifo.xml...
            ttime = ttime.split(".")[0]

        edate = runinf.find(".//{" + fewsNamespace + "}endDateTime")
        sd = datetime.strptime(sdate.attrib["date"] + ttime, "%Y-%m-%d%H:%M:%S")
        ed = datetime.strptime(
            edate.attrib["date"] + edate.attrib["time"], "%Y-%m-%d%H:%M:%S"
        )
        diff = ed - sd

        if timestepsecs < 86400:  # assume hours
            return (diff.seconds + diff.days * 86400) / timestepsecs + 1
        else:
            return diff.days + 1  # Should actually be + 1 but fews starts at 0!
    else:
        print((xmlfile + " does not exists."))


def getEndTimefromRuninfo(xmlfile):
    """
    Gets the endtime of the run from the FEWS runinfo file
    """
    if os.path.exists(xmlfile):
        with open(xmlfile, "r") as f:
            tree = parse(f)
        runinf = tree.getroot()
        edate = runinf.find(".//{" + fewsNamespace + "}endDateTime")
        ed = datetime.strptime(
            edate.attrib["date"] + edate.attrib["time"], "%Y-%m-%d%H:%M:%S"
        )
    else:
        print((xmlfile + " does not exists."))
        ed = None

    return ed


def getStartTimefromRuninfo(xmlfile):
    """
    Gets the starttime from the FEWS runinfo file
    """
    if os.path.exists(xmlfile):
        with open(xmlfile, "r") as f:
            tree = parse(f)
        runinf = tree.getroot()
        edate = runinf.find(".//{" + fewsNamespace + "}startDateTime")
        ttime = edate.attrib["time"]
        if len(ttime) == 12:  # Hack for millisecons in testrunner runinfo.xml...
            ttime = ttime.split(".")[0]
        ed = datetime.strptime(edate.attrib["date"] + ttime, "%Y-%m-%d%H:%M:%S")
        # ed = pa
    else:
        return None

    return ed


def getMapStacksFromRuninfo(xmlfile):
    """
    Gets the list of mapstacks fews expect from the runinfo file and create those
    """

    if os.path.exists(xmlfile):
        with open(xmlfile, "r") as f:
            tree = parse(f)
        runinf = tree.getroot()
        edate = runinf.find(".//{" + fewsNamespace + "}startDateTime")
        ed = datetime.strptime(
            edate.attrib["date"] + edate.attrib["time"], "%Y-%m-%d%H:%M:%S"
        )
    else:
        print((xmlfile + " does not exists."))

    return ed


def pre_adapter(INxmlTimeSeries, logger):
    list_xmlTimeSeries = INxmlTimeSeries.split()
    for xmlTimeSeries in list_xmlTimeSeries:
        logger.info("Converting " + xmlTimeSeries + " ..... ")
        pixml_totss(xmlTimeSeries, case + "/intss/")
        pixml_totss_dates(xmlTimeSeries, case + "/intss/")
        # writeNrTimesteps()


def usage():
    print("wflow_adapter -M Pre -t InputTimeseriesXml -I inifile")
    print("wflow_adapter -M Post -t InputTimeseriesXml -s inputStateFile -I inifile")
    print("              -o outputStateFile -r runinfofile -w workdir -C case")


def main():
    """
    Main entry for using the module as a command line program (e.g. from the Delft-FEWS GA)
    """
    global case
    global runId
    timestepsecs = 86400
    xmldiagfname = "wflow_diag.xml"
    adaptxmldiagfname = "wflow_adapt_diag.xml"
    logfname = "wflow.log"
    netcdfoutput = False

    try:
        opts, _ = getopt.getopt(sys.argv[1:], "-M:-t:-s:-o:-r:-w:-C:-I:R:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))
        usage()
        sys.exit(2)

    if not opts:
        usage()
        sys.exit(2)

    xmlTimeSeries = ""
    stateFile = ""

    mode = "Pre"
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-t"):
            xmlTimeSeries = a
        elif o in ("-R"):
            runId = a
        elif o in ("-o"):
            stateFile = a
        elif o in ("-s"):
            inputStateFile = a
        elif o in ("-r"):
            runinfofile = a
        elif o in ("-w"):
            workdir = a
        elif o in ("-C"):
            case = a
        elif o in ("-I"):
            iniFile = a
        elif o in ("-M"):
            mode = a
        else:
            assert False, "unhandled option"

    # Try and read config file and set default options
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(workdir + "/" + case + "/" + iniFile)

    # get timestep from wflow ini
    timestepsecs = int(
        wflow_lib.configget(config, "run", "timestepsecs", str(timestepsecs))
    )
    netcdf = wflow_lib.configget(config, "framework", "netcdfoutput", "None")
    if netcdf != "None":
        netcdfoutput = True

    logger = setlogger(logfile, "wflow_adapt")

    if mode == "Pre":
        logger.info("Starting preadapter")
        pre_adapter(xmlTimeSeries, logger)
        logger.info("Ending preadapter")
        sys.exit(0)
    elif mode == "Post":
        logger.info("Starting postadapter")

        # Step1: update the state xml files
        pixml_state_updateTime(
            inputStateFile, stateFile, getEndTimefromRuninfo(runinfofile)
        )

        # Step 2: make XML files to go with the output mapstacks if the output is not in netcdf
        # Get outputmapstacks from wflow ini
        mstacks = config.options("outputmaps")
        # Create XML files for all mapstacks if not netcdf
        if not netcdfoutput:
            for a in mstacks:
                var = config.get("outputmaps", a)
                logger.debug(
                    "Creating mapstack xml: "
                    + workdir
                    + "/"
                    + case
                    + "/"
                    + runId
                    + "/"
                    + var
                    + ".xml"
                )
                mapstackxml(
                    workdir + "/" + case + "/" + runId + "/outmaps/" + var + ".xml",
                    var + "?????.???",
                    var,
                    var,
                    getStartTimefromRuninfo(runinfofile),
                    getEndTimefromRuninfo(runinfofile),
                    timestepsecs,
                )

            # Back hack to work around the 0 based FEWS problem and create a double timestep zo that we have connection between subsequent runs in FEWS
            try:
                shutil.copy(
                    workdir + "/" + case + "/instate/SurfaceRunoff.map",
                    workdir + "/" + case + "/" + runId + "/outmaps/run00000.000",
                )
                shutil.copy(
                    workdir + "/" + case + "/instate/WaterLevel.map",
                    workdir + "/" + case + "/" + runId + "/outmaps/lev00000.000",
                )
                shutil.copy(
                    workdir + "/" + case + "/instate/WSO.map",
                    workdir + "/" + case + "/" + runId + "/outmaps/WSO00000.000",
                )
                shutil.copy(
                    workdir + "/" + case + "/instate/LAI.map",
                    workdir + "/" + case + "/" + runId + "/outmaps/LAI00000.000",
                )
            except:
                logger.warning("Cannot copy Surfacerunoff and/or level")

        # Step 3:
        # now check for tss files in the ini file and convert to XML
        stop = 0
        secnr = 0
        while stop == 0:
            if stop == 1:
                break
            try:
                thissection = "outputtss_" + str(secnr)
                tssfiles = config.options(thissection)
                secnr = secnr + 1
                sDate = getStartTimefromRuninfo(runinfofile)

                for aa in tssfiles:
                    if aa not in "samplemap":
                        tssFile = (
                            workdir
                            + "/"
                            + case
                            + "/"
                            + runId
                            + "/"
                            + config.get(thissection, aa)
                        )
                        logger.debug(
                            "Creating xml from tss: "
                            + tssFile
                            + "==> "
                            + tssFile
                            + ".xml"
                        )
                        tss_topixml(
                            tssFile,
                            tssFile + ".xml",
                            "wflow",
                            config.get(thissection, aa),
                            sDate,
                            timestepsecs,
                        )
            except:
                stop = 1

        # Convert log file of model code
        log2xml(case + "/" + runId + "/" + logfname, xmldiagfname)
        logger.info("Ending postadapter")
        # convert logfile of adapter
        log2xml(logfile, adaptxmldiagfname)
    else:
        sys.exit(2)

    # ...


if __name__ == "__main__":
    main()
