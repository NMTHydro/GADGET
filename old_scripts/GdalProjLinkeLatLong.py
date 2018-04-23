import os, fnmatch
import datetime
from dateutil import rrule
import gdal
from subprocess import call

# home = 'I:\\StatewideWater\\NM_IncAng'
INPUT_FOLDER = '/Volumes/SeagateBackupPlusDrive/Gadget_2014_2015/Linke_Turbidity/tifs'
OUTPUT_FOLDER = '/Volumes/SeagateBackupPlusDrive/Gadget_2014_2015/Linke_Turbidity/proj'
#extent='F:\\StatewideWater\\Python\\e2o_dstools_org\\DEM\\Zfactor\\PenWGSDomain.shp'

os.chdir(INPUT_FOLDER)

startday = datetime.datetime(2000, 7, 1, 0)
endday = datetime.datetime(2000, 12, 31, 22)

starthour = datetime.datetime(2000, 1, 1, 5)
endhour = datetime.datetime(2000, 1, 1, 22)

# os.chdir(home)

# filter = '*.tif'

def findRasters (dir, filter):
    for root, dirs, files in os.walk(dir):
        print("root", root)
        print('dirs', dirs)
        print("files", files)
        for file in fnmatch.filter(files, filter):
            print("file", file)
            yield file
print('does this even matter?')
for raster in findRasters(INPUT_FOLDER, '*.tif'): # '**_longlat_wgs84.tif'
    print("raster", raster)
    inRaster = INPUT_FOLDER + '/' + str(raster)
    outRaster = OUTPUT_FOLDER + '/' + str(raster)

    gcp ='-gcp 0.0 0.0 -180.0 90.0 -gcp 0.0 2160.0 -180.0 -90.0 -gcp 4320.0 0.0 180.0 90.0'
##    gcp ='-gcp 0.0 0.0 -180.0 90.0 -gcp -2160.0 0.0 -180.0 -90.0 -gcp 0.0 4320.0 180.0 90.0'

    trans = 'gdal_translate -of GTiff -a_srs EPSG:4326 %s -r "near" %s %s' % (gcp, inRaster, outRaster)
    print(trans)
    call(trans, shell=True)

##    warp = 'gdalwarp -overwrite -multi -s_srs EPSG:4326 -t_srs EPSG:4326 -r "near" -te -180 -90 180 90 -dstnodata "-999" %s %s' % (inRaster, outRaster)
##    call(warp)


