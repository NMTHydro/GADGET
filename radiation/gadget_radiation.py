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
Determine topographically corrected downscaled solar radiation from the NLDAS-2 forcing

usage:

    gadget_radiation.py -I inifile [-S start][-E end][-l loglevel]

    -I inifile - ini file with settings which data to get
    -S start - start timestep (default is 1)
    -E end - end timestep default is last one defined in ini file (from date)
    -l loglevel - DEBUG, INFO, WARN, ERROR
"""
#from pcraster import *
import getopt, sys, os
#import pcraster as pcr
import datetime
import numpy as np
import numexpr as ne
from e2o_utils import *
#import gc
#import psutil
import time
#from gadget_lib import *
from memory_profiler import profile


nthreads=3
start_time = time.time()
ne.set_num_threads(nthreads)

def usage(*args):
    """
    Print usage information

    -  *args: command line arguments given
    """
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)


def save_as_mapsstack(lat,lon,data,times,directory,prefix="GADGET",oformat="Gtiff"):

    cnt = 0
    if not os.path.exists(directory):
        os.mkdir(directory)
    for a in times:
            mapname = getmapname(cnt,prefix,date)
            #print "saving map: " + os.path.join(directory,mapname)
            #writeMap(os.path.join(directory,mapname),oformat,lon,lat[::-1],flipud(data[cnt,:,:]),-999.0)
            writeMap(os.path.join(directory,mapname),oformat,lon,lat,data[cnt,:,:],-999.0)
            cnt = cnt + 1

def save_as_mapsstack_per_day(lat,lon,data,ncnt,date,directory,prefix="GADGET",oformat="Gtiff",FillVal=1E31,gzip=True):
    import platform

    if not os.path.exists(directory):
        os.mkdir(directory)
    mapname = getmapname(ncnt,prefix,date)
    #print "saving map: " + os.path.join(directory,mapname)
    #writeMap(os.path.join(directory,mapname),oformat,lon,lat[::-1],flipud(data[:,:]),-999.0)
    writeMap(os.path.join(directory,mapname),oformat,lon,lat,data[:,:],FillVal)
    if gzip:
        if 'Linux' in platform.system():
            os.system('gzip ' + os.path.join(directory,mapname))


#@profile
def correctRsin(Rsin,currentdate,radiationCorDir,highResDEM,resLowResDEM,lat,longitude, logger):
    """
    Corrects incoming radiation using outputed solar incidence angles maps (r.sun) at each hour for each day of the year
    (here defined from r.sun as angle between the solar beam and plane of surface at each DEM pixel)

    :param Rsin:
    :param currentdate:
    :param radiationCorDir:
    :param logger:
    :return:  corrected incoming radiation, Ra, Rdif, Rbeam, Rtilt, AI, SolAlt
    """

    #get day of year
    logger.info("Correcting incoming radiation with DEM...")
    tt  = currentdate.timetuple()
    JULDAY = tt.tm_yday

    #Current for GMT to MST
    mst = datetime.timedelta(hours=7) #7 hour shift between NLDAS noon and local noon in data, 7 hours during daylight savings
    day = currentdate - mst
    #day = currentdate
    hrs = day.hour
    correct = 0  #Switch for calculating direct beam topographic correction 1 = calculate, 0 = no calculation
    correct2 = 0
    
    FillVal = 0
    Rsin[Rsin < 0.0] = FillVal
    Rsin[Rsin == 0.0] = 0.00001
    Rsinmax = np.max(Rsin)
    
    if Rsinmax > 1:

        if hrs >= 5 and hrs < 23:

            #Read in incidence angle if SW rad data is > 1 and during daylight hours calculated
            path = 'I:\NMIncAngLatLong'
            date = 'NMIncang' + '_' + day.strftime('%j') + '_' + day.strftime('%H') + '.tif'
            resX, resY, cols, rows, LinkeLong, LinkeLat, incang, FillVal = readMap((os.path.join(path,date)),'Gtiff',logger)
            logger.info("Reading Incident Angle Map:" + date)
            incang[incang < 0] = 0
            correct = 1

        #If past last incidence angle map set current incidence angle to all 0s
        if hrs >= 23:
            incang = np.zeros_like(Rsin)
            correct = 1

        #Set time to previous timestep, read in incidence angle
        minus1 = datetime.timedelta(hours=1)
        day = day - minus1
        hrs = day.hour

        if hrs >= 5 and hrs < 23:
                path = 'I:\NMIncAngLatLong'
                date = 'NMIncang' + '_' + day.strftime('%j') + '_' + day.strftime('%H') + '.tif'
                resX, resY, cols, rows, LinkeLong, LinkeLat, incang1, FillVal = readMap((os.path.join(path,date)),'Gtiff',logger)
                logger.info("Reading Incident Angle Map:" + date)
                incang1[incang1 < 0] = 0
                correct2 = 1

        #If previous time step before incidence angle map set previous incidence angle to all 0s
        if hrs <= 5:
                incang1 = np.zeros_like(Rsin)
                correct2 = 1
                        
        if correct == 1:

            #Get extraterrestrial radiation for hourly timestep
            Ra, Opcorr, SolAlt, SolAlt1 = Ra_rad(currentdate, lat, longitude, highResDEM, resLowResDEM, correct2)

            #Apply elevation adjustment based on optical path length with altitude
            Rsin = ne.evaluate("Rsin*Opcorr")
            
            #Calculate clearness index using Ra
            kt = np.maximum(ne.evaluate("Rsin/Ra"),0.0001)
                      
            #Partition radiation following Ruiz-Ariaz, 2010
            kd = np.maximum(0.952 - (1.041*np.exp(-1*np.exp(2.300-4.702*kt))),0.0001) 


            #Read Sky View Factor map for diffuse radiation calculation
            path = 'I:\\StatewideWater\\GADGET\\DEM'
            SVFmap = 'SAGA_NMBuff_50k_SKVF_sm.tif'
            resX, resY, cols, rows, SVFLong, SVFLat, SVF, FillVal = readMap((os.path.join(path,SVFmap)),'Gtiff',logger)

            Incrad = np.radians(incang)

            #Apply direct beam correction based on average of previous hour and current hour incidence angle
            if correct2 == 1:
                Incrad1 = np.radians(incang1)
                Rtilt = np.maximum(ne.evaluate("sin(Incrad)/sin(SolAlt)"),0.0001)
                Rtilt1 = np.maximum(ne.evaluate("sin(Incrad1)/sin(SolAlt1)"),0.0001)
                
                Rbeam = ne.evaluate("Rtilt*((1-kd)*Rsin)")
                Rbeam1 = ne.evaluate("Rtilt1*((1-kd)*Rsin)")
                Rbeam = (Rbeam + Rbeam1) / 2
                Rtilt = (Rtilt + Rtilt1) / 2

            if correct2 != 1:
                #Calculate magnitude of topographic correction for beam radiation
                # Incidence angle (solar beam to plane of tilted surface) divided by solar altitude angle (solar beam to plane of horizontal surface)
                Rtilt = np.maximum(ne.evaluate("sin(Incrad)/sin(SolAlt)"),0.0001)

                #Direct beam component based on beam proportion (1-kd) of decomposed global radiation
                Rbeam = ne.evaluate("Rtilt*((1-kd)*Rsin)")

            #Apply Hay & Davies' model for diffuse radiation on tilted surface
            #Anisotropy index: Horizontal Irradiance / Horizontal Extraterrestrial Irradiance
            AI =  np.minimum(np.maximum(ne.evaluate("((1-kd)*(Rsin / Ra))"),0.0001),1)   
            #Rdif = ne.evaluate('SVF*kd*Rsin')  #Alternative diffuse radiation calc using SVF and assuming isotropic radiation
            Rdif = ne.evaluate('kd*Rsin*((AI*Rtilt)+(SVF*(1-AI)))')
            Rdif[Rdif < 0] = 0
            Rdif[Rdif > 1000] = 0

            #Sum beam and diffuse rad components to reaggregate incoming total SW rad
            Rtot = Rbeam + Rdif

            #Convert SolAlt to degrees to compare with incidence angle r.sun values
            SolAlt *= 180
            SolAlt /= 3.1415

        else:
            Rtot = np.zeros_like(Rsin)
            SolAlt = np.zeros_like(Rsin)

    else:
        Rtot = np.zeros_like(Rsin)
        SolAlt = np.zeros_like(Rsin)


    
    

##    kt = np.zeros_like(Rsin)  
##    kd = np.zeros_like(Rsin)
##    Ra = np.zeros_like(Rsin)
##    Rbeam = np.zeros_like(Rsin)
##    Rtilt = np.zeros_like(Rsin)
##    Rdif = np.zeros_like(Rsin)
##    AI = np.zeros_like(Rsin)
##    Opcorr = np.zeros_like(Rsin)    


    FillVal = 0 
    Rtot[Rtot > 2000] = FillVal
    Rtot[Rtot < 0] = FillVal

    return Rtot, SolAlt
#, Ra, Rdif, Rbeam, Rtilt, AI, SolAlt, Opcorr

def Ra_rad(currentdate, lat, longitude, Altitude, AltAltitude, correct2):
    """
    Calculates extraterrestrial radiation using ASCE method based on latitude(s) for hourly time periods

    :param date:
    :param Lat:
    :param timestep:
    :return:
    Ra for period
    Optical adjustment to global incoming shortwave radiation as function of elevation
    """
    Gsc = 4.92  #Solar constant hourly timestep [MJ m-2 hr-1]
    tstep = 1
    Trans = 0.6 #Thau transmissivity for Air Mass Optical correction
    AtmPcor = np.power(((288.0-0.0065*Altitude)/288.0),5.256)
    AtmPcorAlt = np.power(((288.0-0.0065*AltAltitude)/288.0),5.256)

    #Determine fractional day of year
    tt  = currentdate.timetuple()
    JULDAY = tt.tm_yday

    #Current for GMT to MST (UTC-7) ignoring daylight savings
    mst = datetime.timedelta(hours=7)
    currentdate = currentdate - mst
    t1 = float(currentdate.hour - 1)
    t2 = float(currentdate.hour + 1)
    
    #t = (t1 + t2 )/ 2 #mid-point of period 
    #if correct2 == 0:
    t = float(currentdate.hour)
    
    #CALCULATE EXTRATERRESTRIAL RADIATION
    #Latitude radians
    LatRad = ne.evaluate("lat*pi/180.0")
    # declination (rad)
    #declin = ne.evaluate("0.4093*(sin(((2.0*pi*JULDAY)/365.0)-1.39))")
    #R.sun version
    jday = (2.0*pi*JULDAY)/365.25
    declin = ne.evaluate("arcsin(0.3978*sin(jday-1.4+0.0355*sin(jday-0.0489)))")

    # sunset hour angle
    arccosInput = ne.evaluate("(-(tan(LatRad))*(tan(declin)))")
    arccosInput = np.minimum(1,arccosInput)
    arccosInput = np.maximum(-1,arccosInput)
    sunangle = ne.evaluate("arccos(arccosInput)")
    negsun = -1 * sunangle
    
    #Seasonal correction for solar time
    b = 2 * pi * (JULDAY - 81 ) /  364.0
    Sc = 0.1645 * math.sin(2*b) - 0.1255 * math.cos(b) - 0.025 * math.sin(b)
    
    #solar time angle midpoint of period
    midsun = ne.evaluate("3.14 /12 * ((t + 0.06667*(105.0-abs(longitude))+Sc)-12.0)")  #105 is longitude of center of local time zone, here MST
    midsun1 = ne.evaluate("3.14 /12 * ((t1 + 0.06667*(105.0-abs(longitude))+Sc)-12.0)")  #105 is longitude of center of local time zone, here MST

    #solar time angle previous time step
    #sunminus = ne.evaluate("3.14 /12 * ((t1 + 0.06667*(105.0-abs(longitude))+Sc)-12.0)")
    #solar time angle next time step
    #sunplus = ne.evaluate("3.14 /12 * ((t2 + 0.06667*(105.0-abs(longitude))+Sc)-12.0)")

    #solar time angles end and beginning of periods
    #modify integration limits around sunrise / sunset
    omega1 = np.minimum(sunangle[0,:],np.maximum(ne.evaluate("midsun - (pi*tstep / 24.0)"), negsun[0,:]))
    omega2 = np.minimum(sunangle[0,:],np.maximum(ne.evaluate("midsun + (pi*tstep / 24.0)"), negsun[0,:]))
    omega1 = np.minimum(omega1, omega2)


    #Optical correction for elevation using relative air mass calculation
    #Solar altitude or solar elevation angle (Rads)
    SolAlt = ne.evaluate("arcsin((sin(LatRad)*sin(declin)+cos(LatRad)*cos(declin)*cos(midsun)))")
    SolAlt1 = ne.evaluate("arcsin((sin(LatRad)*sin(declin)+cos(LatRad)*cos(declin)*cos(midsun1)))")
    

#    #Take average solar altitutde with angle from previous timestep if value is low (around sunset) 
#    SolAltminus = ne.evaluate("arcsin((sin(LatRad)*sin(declin)+cos(LatRad)*cos(declin)*cos(sunminus)))")
#    #Take average solar altitutde with angle from next timestep if value is low (around sunrise) 
#    SolAltplus = ne.evaluate("arcsin((sin(LatRad)*sin(declin)+cos(LatRad)*cos(declin)*cos(sunplus)))")
#    SolAltmin = np.min(SolAlt)
#    
#    if (t > 16) and (SolAltmin < 0.05):
#        SolAlt += SolAltminus
#        SolAlt /= 2
#    if (t < 9 ) and (SolAltmin < 0.05):
#        SolAlt += SolAltplus 
#        SolAlt /= 2
    
    OpCorr = ne.evaluate("Trans**((sqrt(1229.0+(614.0*sin(SolAlt))**2)-614.0*sin(SolAlt))*AtmPcor)")
    AltOpCorr = ne.evaluate("Trans**((sqrt(1229.0+(614.0*sin(SolAlt))**2)-614.0*sin(SolAlt))*AtmPcorAlt)")

    #Optical correction between low-res DEM to high-res DEM
    OpDEM = ne.evaluate("OpCorr/AltOpCorr")
    OpDEM[np.isnan(OpDEM)]=0
    OpDEM[OpDEM < 0] = 1
    
    # distance of earth to sun
    distsun = ne.evaluate("1+0.033*(cos((2*pi*JULDAY)/365.0))")
    # Ra = water equivalent extra terrestiral radiation in MJ day-1
    Ra = ne.evaluate("((12*Gsc) / 3.14) * distsun * ((omega2-omega1)*(sin(LatRad))*(sin(declin))+(cos(LatRad))*(cos(declin))*(sin(omega2)-sin(omega1)))")
    Ra /=  0.003600   #Convert MJ hr to Watts m-2 (divided by 3600 * 1E-6)
    Ra[Ra < 0] = 0.0001
    
    return Ra, OpDEM, SolAlt, SolAlt1



#### MAIN ####

#@profile
def main(argv=None):
    # Set all sorts of defaults.....
    serverroot = 'F:\\StatewideWater'
    wrrsetroot = '\\NLDAS'

    #available variables with corresponding file names and standard_names as in NC files
    variables = ['SurfaceIncidentShortwaveRadiation']
    filenames = ["NLDAS_FORA0125_H_"]
    standard_names = ['SW radiation flux downwards (surface)']


    #tempdir

    #defaults in absence of ini file
    startyear = 1980
    endyear= 1980
    startmonth = 1
    endmonth = 12
    latmin = 51.25
    latmax = 51.75
    lonmin = 5.25
    lonmax = 5.75
    startday = 1
    endday = 1
    getDataForVar = True
    calculateEvap = True
    evapMethod = None
    downscaling = None
    resampling = None
    StartStep = 1
    EndStep = 0

    nrcalls = 0
    loglevel=logging.INFO

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            exit()
    try:
        opts, args = getopt.getopt(argv, 'I:l:S:E:')
    except getopt.error, msg:
        usage(msg)

    for o, a in opts:
        if o == '-I': inifile = a
        if o == '-E': EndStep = int(a)
        if o == '-S': StartStep = int(a)
        if o == '-l': exec "loglevel = logging." + a

    logger = setlogger("gadget_radiation.log","gadget_radiation",level=loglevel)
    #logger, ch = setlogger("e2o_getvar.log","e2o_getvar",level=loglevel)
    logger.info("Reading settings from ini: " + inifile)
    theconf = iniFileSetUp(inifile)


    # Read period from file
    startyear = int(configget(logger,theconf,"selection","startyear",str(startyear)))
    endyear = int(configget(logger,theconf,"selection","endyear",str(endyear)))
    endmonth = int(configget(logger,theconf,"selection","endmonth",str(endmonth)))
    startmonth = int(configget(logger,theconf,"selection","startmonth",str(startmonth)))
    endday = int(configget(logger,theconf,"selection","endday",str(endday)))
    startday = int(configget(logger,theconf,"selection","startday",str(startday)))
    start = datetime.datetime(startyear,startmonth,startday)
    end = datetime.datetime(endyear,endmonth,endday)

    #read remaining settings from in file
    lonmax = float(configget(logger,theconf,"selection","lonmax",str(lonmax)))
    lonmin = float(configget(logger,theconf,"selection","lonmin",str(lonmin)))
    latmax = float(configget(logger,theconf,"selection","latmax",str(latmax)))
    latmin = float(configget(logger,theconf,"selection","latmin",str(latmin)))
    BB = dict(
           lon=[ lonmin, lonmax],
           lat= [ latmin, latmax]
           )
    serverroot = configget(logger,theconf,"url","serverroot",serverroot)
    wrrsetroot = configget(logger,theconf,"url","wrrsetroot",wrrsetroot)
    oformat = configget(logger,theconf,"output","format","PCRaster")
    odir = configget(logger,theconf,"output","directory","output/")
    oprefix = configget(logger,theconf,"output","prefix","E2O")
    radcordir = configget(logger,theconf,"downscaling","radiationcordir","output_rad")
    FNhighResDEM = configget(logger,theconf,"downscaling","highResDEM","downscaledem.map")
    FNlowResDEM = configget(logger,theconf,"downscaling","lowResDEM","origdem.map")
    saveAllData = int(configget(logger,theconf,"output","saveall","0"))
    
    #Set clone for DEM
    #pcr.setclone(FNhighResDEM)

    # Check whether downscaling should be applied
    resamplingtype   = configget(logger,theconf,"selection","resamplingtype","linear")
    downscaling   = configget(logger,theconf,"selection","downscaling",downscaling)
    resampling   = configget(logger,theconf,"selection","resampling",resampling)

    if downscaling == 'True' or resampling =="True":
        # get grid info
        resX, resY, cols, rows, highResLon, highResLat, highResDEM, FillVal = readMap(FNhighResDEM,'GTiff',logger)
        LresX, LresY, Lcols, Lrows, lowResLon, lowResLat, lowResDEM, FillVal = readMap(FNlowResDEM,'GTiff',logger)
        #writeMap("DM.MAP","PCRaster",highResLon,highResLat,highResDEM,FillVal)
        #elevationCorrection, highResDEM, resLowResDEM = resampleDEM(FNhighResDEM,FNlowResDEM,logger)
        demmask=highResDEM != FillVal
        mismask=highResDEM == FillVal
        Ldemmask=lowResDEM != FillVal
        Lmismask=lowResDEM == FillVal
        # Fille gaps in high res DEM with Zeros for ineterpolation purposes
        lowResDEM[Lmismask] = 0.0
        resLowResDEM = resample_grid(lowResDEM,lowResLon, lowResLat,highResLon, highResLat,method=resamplingtype,FillVal=0.0)

        lowResDEM[Lmismask] = FillVal
        elevationCorrection = highResDEM - resLowResDEM



    #Check whether evaporation should be calculated
    relevantVars = ['Temperature','SurfaceIncidentShortwaveRadiation']
        
    currentdate = start
    ncnt = 0
    if EndStep == 0:
        EndStep = (end - start).days + 1

    logger.info("Will start at step: " + str(StartStep) + " date/time: " + str(start + datetime.timedelta(hours=StartStep)))
    logger.info("Will stop at step: " + str(EndStep) + " date/time: " + str(start + datetime.timedelta(hours=EndStep)))


    while currentdate <= end:
        ncnt += 1
        start_steptime = time.time()
        
        #Correct for GMT to MST
        mst = datetime.timedelta(hours=7)
        day = currentdate - mst
        #day = currentdate
        hrs = day.hour

        #Only correct for times with incidence angle map and during daylight
        while (hrs < 5) or (hrs >= 23):
            ncnt +=1
            currentdate += datetime.timedelta(hours=1)
            day = currentdate - mst
            #print(currentdate)
            #print(hrs)
            hrs = day.hour

        if ncnt > 0 and ncnt >= StartStep and ncnt <= EndStep:
            # Get all daily datafields needed and aad to list
            relevantDataFields = []
            # Get all data fro this timestep
            mapname = os.path.join(odir,getmapname(ncnt,oprefix,currentdate))
            if os.path.exists(mapname) or os.path.exists(mapname + ".gz") or os.path.exists(mapname + ".zip"):
                logger.info("Skipping map: " + mapname)
            else:
                for i in range (0,len(variables)):
                    if variables[i] in relevantVars:

                        filename = filenames[i]
                        logger.info("Getting data field: " + filename)
                        standard_name = standard_names[i]
                        # logger.info("Get file list..")
                        #  tlist, timelist = get_times_daily(currentdate,currentdate,serverroot, wrrsetroot, filename,logger)
                        # logger.info("Get dates..")
                        #  ncstepobj = getstepdaily(tlist,BB,standard_name,logger)

                        timestepSeconds = 3600
                        logger.info("Get file list..")
                        tlist, timelist = get_times(currentdate,currentdate,serverroot, wrrsetroot, filename,timestepSeconds,logger )
                        logger.info("Get dates..")
                        ncstepobj = getstep(tlist,BB,standard_name,timestepSeconds,logger)
                        
                        logger.info("Get data...: " + str(timelist))
                        mstack = ncstepobj.getdates(timelist)
#                        print('Mstack shape from get dates')
#                        print mstack.shape
                        mean_as_map = flipud(mstack.mean(axis=2))  #Time dimension is number 3(2) instead of 1st (0)
#                        print('Meanasmap shape from flipud')
#                        print mstack.shape
                        
                        #tmp = resample_grid(mean_as_map,ncstepobj.lon, ncstepobj.lat,highResLon, highResLat,method='nearest',FillVal=FillVal)

                        #save_as_mapsstack_per_day(highResLat,highResLon,tmp,int(ncnt),odir,prefix=str(i),oformat=oformat,FillVal=FillVal)
                                                

                        logger.info("Get data body...")
                        if downscaling == 'True' or resampling =="True":
                            logger.info("Downscaling..." + variables[i])
                            #save_as_mapsstack_per_day(ncstepobj.lat,ncstepobj.lon,mean_as_map,int(ncnt),'temp',prefixes[i],oformat='GTiff')
                            #mean_as_map = resample(FNhighResDEM,prefixes[i],int(ncnt),logger)

                            if downscaling == "True":
                                if variables[i]     == 'SurfaceIncidentShortwaveRadiation':
                                    mean_as_map = resample_grid(mean_as_map,ncstepobj.lon, ncstepobj.lat,highResLon,highResLat,method=resamplingtype,FillVal=FillVal)
                                    
                                    #only needed once
                                    if nrcalls == 0:
                                        nrcalls = nrcalls + 1
                                        latitude = ncstepobj.lat[:]
                                        #assuming a resolution of 0.5 degrees
                                        LATITUDE = np.ones(((8*(latmax-latmin)),(8*(lonmax-lonmin)))) #Instead of 2 (0.5 degrees), use 8 (1/8 degree)
                                        for i in range (0,int((8*(lonmax-lonmin)))):
                                            LATITUDE[:,i]=LATITUDE[:,i]*latitude
                                        #save_as_mapsstack_per_day(ncstepobj.lat,ncstepobj.lon,LATITUDE,int(ncnt),'temp','lat',oformat=oformat)
                                        #LATITUDE = resample(FNhighResDEM,'lat',int(ncnt),logger)
                                        LATITUDE = np.zeros_like(highResDEM)
                                        for i in range(0,LATITUDE.shape[1]):
                                            LATITUDE[:,i] = highResLat
                                    
                                    mean_as_map[mismask] = FillVal
                                    mean_as_map, SolAlt   = correctRsin(mean_as_map,currentdate,radcordir,highResDEM,resLowResDEM,LATITUDE,highResLon,logger)
                                    #mean_as_map, Ra, Rd, Rb, Rtilt, AI, SolAlt, Opcorr   = correctRsin(mean_as_map,currentdate,radcordir,highResDEM,resLowResDEM,LATITUDE,highResLon,logger)
                                    mean_as_map[mismask] = FillVal
                                    #Rd[mismask]          = FillVal

                        relevantDataFields.append(mean_as_map)
                        
                                            
                        

                #assign longitudes and lattitudes grids
                if downscaling == 'True' or resampling == "True":
                    lons = highResLon
                    lats = highResLat
                else:
                    lons = ncstepobj.lon
                    lats = ncstepobj.lat
                    
                    
                logger.info("Saving solar radiation data for: " +str(currentdate))
                #logger.info("Diffuse partioning coefficient: " +str(currentdate))
                save_as_mapsstack_per_day(lats,lons,relevantDataFields[0],int(ncnt),currentdate,odir,prefix='RSIN',oformat=oformat,FillVal=FillVal)
                #save_as_mapsstack_per_day(lats,lons,Opcorr,int(ncnt),currentdate,odir,prefix='Opcorr',oformat=oformat,FillVal=FillVal)
                #save_as_mapsstack_per_day(lats,lons,Rd,int(ncnt),currentdate,odir,prefix='Rdif',oformat=oformat,FillVal=FillVal)
                #save_as_mapsstack_per_day(lats,lons,Rb,int(ncnt),currentdate,odir,prefix='Rb',oformat=oformat,FillVal=FillVal)
                #save_as_mapsstack_per_day(lats,lons,Rtilt,int(ncnt),currentdate,odir,prefix='Rtilt',oformat=oformat,FillVal=FillVal)
                #save_as_mapsstack_per_day(lats,lons,AI,int(ncnt),currentdate,odir,prefix='AI',oformat=oformat,FillVal=FillVal)
                #save_as_mapsstack_per_day(lats,lons,SolAlt,int(ncnt),currentdate,odir,prefix='SunAlt',oformat=oformat,FillVal=FillVal)
                #save_as_mapsstack_per_day(lats,lons,Ra,int(ncnt),currentdate,odir,prefix='Ra',oformat=oformat,FillVal=FillVal)

                
                Ra = []
                Rd = []
                Rb = []
                Rtilt = []
                AI = []
                SolAlt = []
                
        else:
            pass
        
        relevantDataFields=[]
        gc.collect()
        compsteptime = (time.time() - start_steptime)
        a = open("gadget_rad_comptime.txt", "a")
        a.write(str(currentdate) + ' Computation time: '+ str(compsteptime)+ ' seconds' + '\n')
        a.close()        
        
        currentdate += datetime.timedelta(hours=1)

    logger.info("Done.")
    comptime = (time.time() - start_time)
    logger.info("Computation time : " +str(comptime) +' seconds')

if __name__ == "__main__":
    main()
