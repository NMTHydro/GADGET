
"""
The ``clearsky`` module contains several methods
to calculate clear sky GHI, DNI, and DHI.
"""

from __future__ import division

import os
from collections import OrderedDict

import numpy as np
import pandas as pd

#import gdal
#from osgeo import osr, gdal
import osr, gdal
import datetime
from dateutil import rrule
from datetime import timedelta
import time
import numexpr as ne
from pip import logger

from scipy import interpolate
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt

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
    if ds is None:
        logger.error('Could not open ' + fileName + '. Something went wrong!! Shutting down')
        gdal.sys.exit(1) #sys.exit(1) - changed GELP
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
    #print(x)
    #print(y)

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
        print('Writing to temporary file ' + fileName + '.tif')
        print("Output format: " + fileFormat)
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
        print('Writing to ' + fileName)
    outDataset = driver2.CreateCopy(fileName, TempDataset, 0)
    TempDataset = None
    outDataset = None
    if verbose:
        print('Removing temporary file ' + fileName + '.tif')
    os.remove(fileName + '.tif');

    if verbose:
        print('Writing to ' + fileName + ' is done!')


def lookup_linke_turbidity(time, latitude, longitude, filepath=None,
                           interp_turbidity=True):
    """
    Look up the Linke Turibidity from the ``LinkeTurbidities.mat``
    data file supplied with pvlib.
    Parameters
    ----------
    time : pandas.DatetimeIndex
    latitude : float
    longitude : float
    filepath : string
        The path to the ``.mat`` file.
    interp_turbidity : bool
        If ``True``, interpolates the monthly Linke turbidity values
        found in ``LinkeTurbidities.mat`` to daily values.
    Returns
    -------
    turbidity : Series
    """

    # The .mat file 'LinkeTurbidities.mat' contains a single 2160 x 4320 x 12
    # matrix of type uint8 called 'LinkeTurbidity'. The rows represent global
    # latitudes from 90 to -90 degrees; the columns represent global longitudes
    # from -180 to 180; and the depth (third dimension) represents months of
    # the year from January (1) to December (12). To determine the Linke
    # turbidity for a position on the Earth's surface for a given month do the
    # following: LT = LinkeTurbidity(LatitudeIndex, LongitudeIndex, month).
    # Note that the numbers within the matrix are 20 * Linke Turbidity,
    # so divide the number from the file by 20 to get the
    # turbidity.

    try:
        import scipy.io
    except ImportError:
        raise ImportError('The Linke turbidity lookup table requires scipy. ' +
                          'You can still use clearsky.ineichen if you ' +
                          'supply your own turbidities.')

    if filepath is None:
        pvlib_path = os.path.dirname(os.path.abspath(__file__))
        filepath = os.path.join(pvlib_path, 'data', 'LinkeTurbidities.mat')
    print(filepath)

    #mat = scipy.io.loadmat(filepath)
    linke_turbidity_table = filepath

    latitude_index = np.arange(0,2160,1)
    longitude_index = np.arange(0,4320,1)

    #g = linke_turbidity_table[latitude_index][longitude_index]
    g = filepath
    g3 = g[1,:,:]
    print(g3.shape)
    x,y = np.meshgrid(longitude_index,latitude_index)
    f = interpolate.interp2d(x,y,g,kind='bilinear')
    

    if interp_turbidity:
        # Data covers 1 year.
        # Assume that data corresponds to the value at
        # the middle of each month.
        # This means that we need to add previous Dec and next Jan
        # to the array so that the interpolation will work for
        # Jan 1 - Jan 15 and Dec 16 - Dec 31.
        # Then we map the month value to the day of year value.
        # This is approximate and could be made more accurate.
        g2 = np.concatenate([[g[-1,:,:]], g[:,:,:], [g[0,:,:]]])
        days = np.linspace(-15, 380, num=14)
        doy = int(time.strftime('%j'))
        print('day of the year:',doy)
        print(g2.shape)

        # TODO - Figure out why this was commented...
   #     linke_turbidity = pd.Series(np.interp(doy, days, g2),
   #                                 index=time)
   #
   # else:
   #     linke_turbidity = pd.DataFrame(time.month, index=time)
   #     # apply monthly data
   #     linke_turbidity = linke_turbidity.apply(lambda x: g[x[0]-1], axis=1)
   #
   #
   #
   #  linke_turbidity /= 20.
   #
   #  return linke_turbidity



def time_interpolation(day, Lat, Lon):

    # Peter's old paths
    # base_dir = 'G:\\GADGET\\LinkeTurbidity\\crop_Walnut'
    # output_linke_write_NDimage = 'G:\\GADGET\\LinkeTurbidity\\linke_METDATA_by_day_Walnut'

    # New Paths
    base_dir = '/Volumes/SeagateBackupPlusDrive/Gadget_2014_2015/Linke_Turbidity/proj'
    output = '/Volumes/SeagateBackupPlusDrive/Gadget_2014_2015/Linke_Turbidity/linke_metdata_by_day'

    this_month = day.month
    this_monthobj = datetime.datetime(day.year,this_month,day.day,0)

    if day.day < 16:
        in_one_month = this_monthobj.replace(day=1)
        thisyear = day.year

        if day.month == 1:
            thisyear -= 1

        this_month = subtract_one_month(this_monthobj)
        this_month = this_month.month
        this_monthobj = datetime.datetime(thisyear, this_month, day.day, 0)
        # this_monthobj.replace(year=thisyear)

        # strmonth = int(this_month.strftime('%m'))

    elif day.day >= 16:
        thisyear = day.year
        this_monthobj = datetime.datetime(thisyear, this_month, day.day, 0)

        if day.month == 12:
            thisyear += 1

        in_one_month = add_one_month(this_monthobj)
        in_one_month.replace(year=thisyear)


    # day.replace(year=thisyear)
    print('day year', day.year)
    print('this year', thisyear)
    print('this month', this_month)
    print('in one month', in_one_month)

    month = '%02d' % this_month
    str_nextmonth = int(in_one_month.strftime('%m'))
    mid_thismonth = this_monthobj.replace(day=1) + datetime.timedelta(days=14)
    print('middle of this month',mid_thismonth)
    mid_nextmonth = in_one_month + datetime.timedelta(days=14)
    print('middle of next month', mid_nextmonth)
    next_month = '%02d' % str_nextmonth

    dates = '{}_longlat_wgs84.tif'.format(month)
    print(dates)
    next_date = '{}_longlat_wgs84.tif'.format(next_month)
    print(next_date)
    resX, resY, cols, rows, Lon, Lat, Linke, prj, FillVal = readMap(os.path.join(base_dir,dates), 'Gtiff')
    # print(Linke.shape)

    Linke = Linke / 20.

    resX, resY, cols, rows, Lon, Lat, LinkeNext, prj, FillVal = readMap(os.path.join(base_dir,next_date), 'Gtiff')

    LinkeNext = LinkeNext / 20.

    latitude_index = np.arange(5)
    longitude_index = np.arange(16)
    arr1 = Linke
    arr2 = LinkeNext

    # rejoin Linke, LinkeNext into a single array of shape (2, 2160, 4320)
    arr = np.r_['0,3', Linke, LinkeNext]
    # print('arr.shape',arr.shape)

    # define the grid coordinates where you want to interpolate
    Y, X = np.meshgrid(longitude_index,latitude_index)
    # print('X',X)
    # print(X.shape)
    # print('Y',Y)
    # print(Y.shape)

#        a = np.linspace(15, 365, num=12)

    #Setup time variables for interpolation
    # days_dif = mid_nextmonth - day
    days_dif = mid_nextmonth - day
    days_dif = days_dif.days
    max_days_diff = mid_nextmonth - mid_thismonth
    max_days_diff = max_days_diff.days
    print('day',day)
    print('days difference from mid next month',days_dif)
    print('out of max days difference from mid this month', max_days_diff)
    interp = 1 - (days_dif / max_days_diff) #1 = weight completely next month linke values, 0 = previous month
    print('interp ratio between monthly images',interp)


#        
#     def interp(m1, t1, m2, t2, t_interp):
#         (t2-t1) * (t_interp-t1)

    # 0.5 corresponds to half way between arr1 and arr2
    coordinates = np.ones((5,16))*interp, X, Y
    coordones = np.ones((5,16))*interp
    #print('coordinates',coordinates)
    #print(coordones.shape)

    # given arr interpolate at coordinates
    newarr = ndimage.map_coordinates(arr, coordinates, order=2)

#        fig, ax = plt.subplots(nrows=3)
#        cmap = plt.get_cmap('Greys')
#
#        vmin = np.min([arr1.min(), newarr.min(), arr2.min()])
#        vmax = np.max([arr1.max(), newarr.max(), arr2.max()])
#        ax[0].imshow(arr1, interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
#        ax[1].imshow(newarr, interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
#        ax[2].imshow(arr2, interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
#        ax[0].set_xlabel('arr1')
#        ax[1].set_xlabel('interpolated')
#        ax[2].set_xlabel('arr2')
#        plt.show()

    return newarr

def add_one_month(dt0):
    dt1 = dt0.replace(day=1)       
    dt2 = dt1 + timedelta(days=32)
    dt3 = dt2.replace(day=1)
    return dt3

def subtract_one_month(dt0):
    dt1 = dt0.replace(day=1)
    dt2 = dt1 - timedelta(days=2)  
    dt3 = dt2.replace(day=1)
    return dt3

##count = 0
##for this_month in rrule.rrule(rrule.MONTHLY, dtstart=startday, until=endday):
##    month = this_month.month
##    strmonth = int(this_month.strftime('%m'))
##    month = '%02d' % strmonth
##    next_month = add_one_month(this_month)
##    dates = 'linke{}_longlat_wgs84.tif'.format(month)
##    next_date = 'linke{}_longlat_wgs84.tif'.format(month)
##    resX, resY, cols, rows, Lon, Lat, Linke, prj, FillVal = readMap(os.path.join(base_dir,dates), 'Gtiff')
##    print(Linke.shape)
##
##    if count == 0:
##        resX, resY, cols, rows, Lon, Lat, LinkeNext, prj, FillVal = readMap(os.path.join(base_dir,next_date), 'Gtiff')
##        agg = np.array([Linke,LinkeNext])
##        agg = np.vstack([Linke[np.newaxis,...], LinkeNext[np.newaxis,...] ])
##
##    if count > 0:
##        #agg = np.concatenate((agg[...,np.newaxis],Linke[...,np.newaxis]),axis=2)
##        agg = np.vstack([agg, Linke[np.newaxis,...] ])
##        print(agg.shape)
##
##    count += 1

#FillVal = 0
#agg[np.isinf(agg)] = FillVal

#Write stacked Linke values to new Raster
##stacked = 'LinkeAgg.tif'
##outname = os.path.join(output_linke_write_NDimage, stacked)
##print(agg.shape)
##print(agg)
##
##writeMap(outname,'Gtiff', Lon, Lat, agg, sr_wkt, FillVal)

# base_dir='G:\\linketurbidity-master\\linketurbidity\\crop_METDATA'
# output_linke_write_NDimage='G:\\linketurbidity-master\\linketurbidity\\linke_METDATA_by_day'

# TODO - if you are a new user change this for your own computer.
# base_dir = 'G:\\GADGET\\LinkeTurbidity\\crop_Walnut'
base_dir = '/Volumes/SeagateBackupPlusDrive/Gadget_2014_2015/Linke_Turbidity/proj'
# output_linke_write_NDimage = 'G:\\GADGET\\LinkeTurbidity\\linke_METDATA_by_day_Walnut'
output = '/Volumes/SeagateBackupPlusDrive/Gadget_2014_2015/Linke_Turbidity/linke_metdata_by_day'

# Don't change these?
startday = datetime.datetime(2000,1,1,0)
endday = datetime.datetime(2000,12,31,0)

# the ref map is just an example map from the previous step output_linke_write_NDimage from GdalProjLinkeLatLong
ref_map = 'April.tif' #
resX, resY, cols, rows, Lon, Lat, Linke, prj, FillVal = readMap(os.path.join(base_dir,ref_map), 'Gtiff')

srs = osr.SpatialReference(prj)
sr_wkt = srs.ExportToWkt()

for day in rrule.rrule(rrule.DAILY, dtstart=startday, until=endday):
    nr = day.strftime('%j')
    dates = 'linkewgs84_{}.tif'.format(nr)
    fullfile = os.path.join(output,dates)
    print(fullfile)

    #Interpolate daily values of linke turbidity
    #linke_daily = lookup_linke_turbidity(day, Lat, Lon, agg, interp_turbidity=True)
    linke_daily = time_interpolation(day, Lat, Lon)

    #Write daily values to new daily rasters
    daily_doy = 'linkewgs84_{}.tif'.format(nr)
    outname = os.path.join(output, daily_doy)
    writeMap(outname,'Gtiff', Lon, Lat, linke_daily, sr_wkt, FillVal)
    

