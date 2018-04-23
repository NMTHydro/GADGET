
"""
Created on Sun Dec 27 16:22:20 2015

@author: PETERS
"""

import numpy as np
import gdal
import os
import datetime
from dateutil import rrule
import matplotlib.pyplot as plt
import time
import numexpr as ne
import osr
from numba import jit


def readMap(fileName, fileFormat):
    """
    read geographical file into memory

    :param fileName: Name of the file to read
    :param fileFormat: Gdal format string
    :param logger: logger object
    :return: resX, resY, cols, rows, x, y, data, FillVal
    """

    #Open file for binary-reading

    mapFormat = gdal.GetDriverByName(fileFormat)
    mapFormat.Register()
    ds = gdal.Open(fileName)
    prj = ds.GetProjection()
#    if ds is None:
#        logger.error('Could not open ' + fileName + '. Something went wrong!! Shutting down')
#        sys.exit(1)
    # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX    = geotrans[1]
    resY    = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX+resX/2,originX+resX/2+resX*(cols-1),cols)
    y = np.linspace(originY+resY/2,originY+resY/2+resY*(rows-1),rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1) # there's only 1 band, starting from 1

    data = RasterBand.ReadAsArray(0,0,cols,rows)
    FillVal = RasterBand.GetNoDataValue()
    RasterBand = None
    del ds
    return resX, resY, cols, rows, x, y, data, prj, FillVal
    
def writeMap(fileName, fileFormat, x, y, data,  prj, FillVal):
    """
    Write geographical data into file. Also replave NaN bu FillVall

    :param fileName:
    :param fileFormat:
    :param x:
    :param y:
    :param data:
    :param FillVal:
    :return:
    """


    verbose = False
    gdal.AllRegister()
    driver1 = gdal.GetDriverByName('GTiff')
    driver2 = gdal.GetDriverByName(fileFormat)

    data[np.isnan(data)] = FillVal
    # Processing
    if verbose:
        print 'Writing to temporary file ' + fileName + '.tif'
        print "Output format: " + fileFormat
    # Create Output filename from (FEWS) product name and date and open for writing
    TempDataset = driver1.Create(fileName + '.tif',data.shape[1],data.shape[0],1,gdal.GDT_Float32)
    # Give georeferences
    xul = x[0]-(x[1]-x[0])/2
    yul = y[0]+(y[0]-y[1])/2

    TempDataset.SetGeoTransform( [ xul, x[1]-x[0], 0, yul, 0, y[1]-y[0] ] )
    TempDataset.SetProjection(prj)
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data,0,0)
    TempBand.FlushCache()
    TempBand.SetNoDataValue(FillVal)
    # Create data to write to correct format (supported by 'CreateCopy')
    if verbose:
        print 'Writing to ' + fileName 
    outDataset = driver2.CreateCopy(fileName, TempDataset, 0)
    TempDataset = None
    outDataset = None
    if verbose:
        print 'Removing temporary file ' + fileName + '.tif'
    os.remove(fileName + '.tif');

    if verbose:
        print 'Writing to ' + fileName + ' is done!'
        
        



#@jit
def main():
    
    start_time = time.time()

    #path = r'I:\StatewideWater\Linux\processing\UTM\PM2000'
    #os.chdir(path)
    
    #resX, resY, cols, rows, x, y, elev, prj, FillVal  = readMap("PM_NM_2000_001.tif",'Gtiff')
    #Rcorday = np.zeros_like(elev)
    #Rcorday = np.zeros((2525,2272)) #UTM Array Size
    #Rcorday = np.zeros((2331,2585)) #Lat LONG WGS84 Array Size 
    Rcorday = np.zeros((5,16)) #Lat LONG WGS84 Array Size for METDATA
    #start = datetime.datetime(2001,1,1,0)
    #end = datetime.datetime(2013,12,31,23)
    
    startday = datetime.datetime(2000,1,1,0)
    endday = datetime.datetime(2000,12,31,23)
    
    
    
    #yearit = 1
    #endyear = endday.year
    outpath = 'G:\\Walnut\\FlatRadWm2\\'
    
    #while endyear < 2008:
        #for years in rrule.rrule(rrule.YEARLY, dtstart=start, until=end):
        #Rcorday = np.zeros((2525,2272)) #UTM size
        
    
            #print(path2)
    for day in rrule.rrule(rrule.DAILY, dtstart=startday, until=endday):
        #path2 = 'B:\\GADGET\\PM_2007arid'  
        path2 = 'G:\\Walnut\\rsun_METDATA'
        
        dates = 'Rbflat_'+ day.strftime('%j') + '.tif'
        fullfile = os.path.join(path2,dates)
        print(fullfile)
        resX, resY, cols, rows, Lon, Lat, Rcor, prj, FillVal = readMap(os.path.join(path2,dates), 'Gtiff')
##        Rcor[np.isnan(Rcor)] = 0
##        Rcor[np.isinf(Rcor)] = 0
##        Rcor[Rcor < 0] = 0
##        Rcor[Rcor > 50] = 0
##        Rcor[Rcor == FillVal] = 0
        Rcorday = Rcor / 24.0 #Divide by 24 hrs to go from Wh.m-2.day-1 to W m-2
       # testInc = "IncAng" + count3d
       # report(nmfill,testInc)
       #aguila(Rcorday)
        
        Rcorday[np.isnan(Rcorday)] = FillVal
        Rcorday[np.isinf(Rcorday)] = FillVal
        Rcorday[Rcorday < 0] = FillVal
        Rcorday[Rcorday > 5000] = FillVal
    
        srs = osr.SpatialReference(prj)
        sr_wkt = srs.ExportToWkt()
        output_file_name = dates
        outname = outpath + '\\' + output_file_name
        print('Saving converted daily radiation for ' + str(output_file_name))
        
        if not os.path.exists(outpath):
            os.makedirs(outpath)    
        
        writeMap(outname,'Gtiff', Lon, Lat, Rcorday, sr_wkt, FillVal)
            
        #Add year until 2014
##        dtstart = startday
##        newstart = dtstart.year + yearit
##        startday = startday.replace(year=newstart)
##        until = endday
##        newend = until.year + yearit
##        endday = endday.replace(year=newend)
##        endyear = newend
        
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()
