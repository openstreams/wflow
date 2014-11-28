#!/usr/bin/python

# Wflow is Free software, see below:
#
# Copyright (c) J. Schellekens/Deltares 2005-2014
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



# $Rev:: 904           $:  Revision of last commit
# $Author:: schelle    $:  Author of last commit
# $Date:: 2014-01-13 1#$:  Date of last commit
"""

    Usage:
    wflow_prepare_rad -D DEM [-S start day][-E end day][-M][-x lon][-y lat][-h][-l loglevel][-T minutes]

    -D DEM Filename of the digital elevation model
    -S Startday - Start day of the simulation (1 Jan is 1)
    -E EndDay - End day of the simulation
    -T minuts - timeresolution in minutss (60 default is 1 hour)
    -M The DEM xy units are in metres (instead of lat/lon)
    -x longitute of the map left (if map xy in metres)
    -y lattitude of the map bottom (if map xy in metres)
    -l loglevel Set loglevel (DEBUG, INFO, WARNING,ERROR)
    -s start hour (per day) of the calculations (default =1)
    -e end hour (per day) of the calculations (default = 23)
    -h This information

"""


from pcraster import *
import sys
import logging
import wflow.pcrut as pcr
import getopt
import numpy as np






def correctrad(Day,Hour,Lat,Lon,Slope,Aspect,Altitude,Altitude_UnitLatLon):
    """ 
    Determines radiation over a DEM assuming clear sky for a specified hour of
    a day
    
    :var Day: Day of the year (1-366)
    :var Hour: Hour of the day (0-23)
    :var Lat: map with latitudes for each grid cell
    :var Lon: map with longitudes for each grid cell
    :var Slope: Slope in degrees
    :var Aspect: Aspect in degrees relative to north for each cell
    :var Altitude: Elevation in metres
    :var Altitude_Degree: Elevation in degrees. If the actual pcraster maps
                          are in lat lon this maps should hold the Altitude converted
                          to degrees. If the maps are in metres this maps should also
                          be in metres

    :return Stot: Total radiation on the dem, shadows not taken into account
    :return StotCor: Total radiation on the dem taking shadows into acount
    :return StotFlat: Total radiation on the dem assuming a flat surface
    :return Shade: Map with shade (0) or no shade (1) pixels
    """

    Sc  = 1367.0          # Solar constant (Gates, 1980) [W/m2]
    Trans   = 0.6             # Transmissivity tau (Gates, 1980)    
    pi = 3.1416
    a = pow(100,5.256)
    #report(a,"zz.map")

    AtmPcor = pow(((288.0-0.0065*Altitude)/288.0),5.256)
    #Lat = Lat * pi/180
    ##########################################################################
    # Calculate Solar Angle and correct radiation ############################
    ##########################################################################
    # Solar geometry
    # ----------------------------
    # SolDec  :declination sun per day  between +23 & -23 [deg]
    # HourAng :hour angle [-] of sun during day
    # SolAlt  :solar altitude [deg], height of sun above horizon
    # SolDec  = -23.4*cos(360*(Day+10)/365);
    # Now added a new function that should work on all latitudes! 
    #theta    =(Day-1)*2 * pi/365  # day expressed in radians
    theta    =(Day-1)*360.0/365.0  # day expressed in degrees
     
    SolDec =180/pi *  (0.006918-0.399912 * cos(theta)+0.070257 * sin(theta) -  0.006758 * cos(2*theta)+0.000907 * sin(2*theta) -  0.002697 *  cos(3*theta)+0.001480 * sin(3*theta))
    
    #HourAng = 180/pi * 15*(Hour-12.01)
    HourAng = 15.0*(Hour-12.01) 
    SolAlt  = scalar(asin(scalar(sin(Lat)*sin(SolDec)+cos(Lat)*cos(SolDec)*cos(HourAng))))
    
    # Solar azimuth                    
    # ----------------------------
    # SolAzi  :angle solar beams to N-S axes earth [deg]
    SolAzi = scalar(acos((sin(SolDec)*cos(Lat)-cos(SolDec)* sin(Lat)*cos(HourAng))/cos(SolAlt)))
    SolAzi = ifthenelse(Hour <= 12, SolAzi, 360 - SolAzi)
    

    # Surface azimuth
    # ----------------------------
    # cosIncident :cosine of angle of incident; angle solar beams to angle surface
    cosIncident = sin(SolAlt)*cos(Slope)+cos(SolAlt)*sin(Slope)*cos(SolAzi-Aspect)
    # Fro flat surface..  
    FlatLine = spatial(scalar(0.00001))
    FlatSpect = spatial(scalar(0.0000))
    cosIncidentFlat = sin(SolAlt)*cos(FlatLine)+cos(SolAlt)*sin(FlatLine)*cos(SolAzi-FlatSpect)
    # Fro flat surface..    
    #cosIncident = sin(SolAlt) + cos(SolAzi-Aspect)


    # Critical angle sun
    # ----------------------------
    # HoriAng  :tan maximum angle over DEM in direction sun, 0 if neg 
    # CritSun  :tan of maximum angle in direction solar beams
    # Shade    :cell in sun 1, in shade 0
    # NOTE: for a changing DEM in time use following 3 statements and put a #
    #       for the 4th CritSun statement
    HoriAng   = cover(horizontan(Altitude_UnitLatLon,directional(SolAzi)),0)
    #HoriAng   = horizontan(Altitude,directional(SolAzi))
    HoriAng   = ifthenelse(HoriAng < 0, scalar(0), HoriAng)
    CritSun   = ifthenelse(SolAlt > 90, scalar(0), scalar(atan(HoriAng)))
    Shade   = SolAlt > CritSun
    #Shade = spatial(boolean(1))
    # Radiation outer atmosphere
    # ----------------------------
    #report(HoriAng,"hor.map")

    OpCorr = Trans**((sqrt(1229+(614*sin(SolAlt))**2) -614*sin(SolAlt))*AtmPcor)    # correction for air masses [-] 
    Sout   = Sc*(1+0.034*cos(360*Day/365.0)) # radiation outer atmosphere [W/m2]
    Snor   = Sout*OpCorr                   # rad on surface normal to the beam [W/m2]

    # Radiation at DEM
    # ----------------------------
    # Sdir   :direct sunlight on dem surface [W/m2] if no shade
    # Sdiff  :diffuse light [W/m2] for shade and no shade
    # Stot   :total incomming light Sdir+Sdiff [W/m2] at Hour
    # Radiation :avg of Stot(Hour) and Stot(Hour-HourStep)
    # NOTE: PradM only valid for HourStep & DayStep = 1

    
    SdirCor   = ifthenelse(Snor*cosIncident*scalar(Shade)<0,0.0,Snor*cosIncident*scalar(Shade))
    Sdir   = ifthenelse(Snor*cosIncident<0,0.0,Snor*cosIncident)
    SdirFlat   = ifthenelse(Snor*cosIncidentFlat<0,0.0,Snor*cosIncidentFlat)
    Sdiff  = ifthenelse(Sout*(0.271-0.294*OpCorr)*sin(SolAlt)<0, 0.0, Sout*(0.271-0.294*OpCorr)*sin(SolAlt))
    #AtmosDiffFrac = ifthenelse(Sdir > 0, Sdiff/Sdir, 1)
    Shade = ifthenelse(Sdir <=0, 0,Shade)

    # Stot   = cover(Sdir+Sdiff,windowaverage(Sdir+Sdiff,3));     # Rad [W/m2]
    Stot   = Sdir + Sdiff                                             # Rad [W/m2]
    StotCor   = SdirCor + Sdiff                                   # Rad [W/m2]
    StotFlat = SdirFlat + Sdiff
    
    
     
    return Stot, StotCor, StotFlat, Shade, 


def GenRadMaps(SaveDir,Lat,Lon,Slope,Aspect,Altitude,DegreeDem,logje,start=1,end=2,interval=60,shour=1,ehour=23):
    """ 
    Generates daily radiation maps for a whole year.
    It does so by running correctrad for a whole year with hourly
    steps and averaging this per day.
    """

    Intperday = 1440./interval
    Starthour = shour
    EndHour = ehour
    Calcsteps = Intperday/24 * 24
    calchours = np.arange(Starthour,EndHour,24/Intperday)

    for Day in range(start,end+1):
        avgrad = 0.0 * Altitude
        _avgrad = 0.0 * Altitude
        _flat = 0.0 * Altitude
        avshade = 0.0 * Altitude
        id = 1
        logje.info("Calulations for day: " + str(Day))
        for Hour in calchours:
            logje.info("Hour: " + str(Hour))
            cradnodem, crad,  flat, shade = correctrad(Day,float(Hour),Lat,Lon,Slope,Aspect,Altitude,DegreeDem)
            avgrad=avgrad + crad
            _flat = _flat + flat
            _avgrad=_avgrad + cradnodem
            avshade=avshade + scalar(shade)
            nrr = "%03d" % id
            #report(crad,"tt000000." + nrr)
            #report(shade,"sh000000." + nrr)
            #report(cradnodem,"ttr00000." + nrr)
            id = id + 1
        
        nr = "%0.3d" % Day
        report(avgrad/Calcsteps,SaveDir + "/RAD00000." + nr)
        report(_avgrad/Calcsteps,SaveDir + "/_RAD0000." + nr)
        report(avshade,SaveDir + "/SHADE000." + nr)
        report(_flat/Calcsteps,SaveDir + "/FLAT0000." + nr)
        report(ifthen((Altitude + 300) > 0.0, cover(avgrad/_flat,1.0)),SaveDir + "/RATI0000." + nr)

def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)


def main(argv=None):
    """
    Perform command line execution of the model.
    """

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return


    try:
        opts, args = getopt.getopt(argv, 'hD:Mx:y:l:O:S:E:T:s:e:')
    except getopt.error, msg:
        usage(msg)


    thedem = "mydem.map"
    xymetres = False
    lat = 52
    lon = 10
    loglevel = logging.DEBUG
    outputdir="output_rad"
    startday = 1
    endday = 2
    calc_interval = 60
    shour=1
    ehour=23

    for o, a in opts:
        if o == '-h': usage()
        if o == '-O': outputdir = a
        if o == '-D': thedem = a
        if o == '-M': xymetres = true
        if o == '-x': lat = int(a)
        if o == '-y': lon = int(a)
        if o == '-S': startday = int(a)
        if o == '-E': endday = int(a)
        if o == '-T': calc_interval = int(a)
        if o == '-l': exec "thelevel = logging." + a
        if o == '-s': shour = int(a)
        if o == '-e': ehour = int(a)


    logger = pcr.setlogger("wflow_prepare_rad.log","wflow_prepare_rad",thelevel=loglevel)
    if not os.path.exists(thedem):
        logger.error("Cannot find dem: " + thedem + " exiting.")
        sys.exit(1)
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    logger.debug("Reading dem: " " thedem")
    setclone(thedem)
    dem = readmap(thedem)

    logger.debug("Calculating slope and aspect...")
    if xymetres:
        LAT = spatial(scalar(lat))
        LON= spatial(scalar(lon))
        Slope = max(0.00001,slope(dem))
        DEMxyUnits = dem
    else:
        LAT= ycoordinate(boolean(dem))
        LON = xcoordinate(boolean(dem))
        Slope = slope(dem)
        xl, yl, reallength = pcr.detRealCellLength(dem * 0.0, 0)
        Slope = max(0.00001, Slope * celllength() / reallength)
        DEMxyUnits = dem * celllength() / reallength

    # Get slope in degrees
    Slope = scalar(atan(Slope))
    Aspect = cover(scalar(aspect(dem)),0.0)

    GenRadMaps(outputdir,LAT,LON,Slope,Aspect,dem,DEMxyUnits,logger,start=startday,end=endday,interval=calc_interval,shour=shour,ehour=ehour)




if __name__ == "__main__":
    main()
