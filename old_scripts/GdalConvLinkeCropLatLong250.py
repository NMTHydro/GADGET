import os, fnmatch
import gdal
from subprocess import call
import datetime
from dateutil import rrule

extent='I:\\StatewideWater\\GIS_shapefiles\\NM_ExtentLatLong.shp'

##os.chdir (INPUT_FOLDER)

def findRasters (path, filter):
    for root, dirs, files in os.walk(path):
        for file in fnmatch.filter(files, filter):
            yield file

INPUT_FOLDER = 'H:\\GADGET\\LinkeTurbidity\\proj'
OUTPUT_FOLDER = 'H:\\GADGET\\LinkeTurbidity\\crop_250m_Walnut_E'
print(INPUT_FOLDER)

for raster in findRasters(INPUT_FOLDER, '*.tif'):
    (FOLDER, rastname)= os.path.split (raster)

    inRaster = INPUT_FOLDER + '/' + rastname
    temp = OUTPUT_FOLDER + '/' + 'temp.tif'
    outRaster = OUTPUT_FOLDER + '/' + rastname
    print(inRaster)
    print(outRaster)

    warp = 'gdalwarp -overwrite -s_srs EPSG:4326\
    -t_srs EPSG:4326 -te -110.1900273268889379 31.6392870650104605 -109.5463334551972707 31.8112862372033760 -tr 0.00260605 0.00260605\
     -r "cubic" -multi -srcnodata "-3.40282346639e+038" -dstnodata -999 %s %s' % (inRaster, temp)
    print(warp)
    call(warp)

    warp2 = 'gdalwarp -overwrite -s_srs EPSG:4326\
     -t_srs EPSG:4326 -te -110.1900273268889379 31.6392870650104605 -109.5463334551972707 31.8112862372033760 -tr 0.00260605 0.00260605\
     -r "cubic" -multi -srcnodata -999 -dstnodata -999 %s %s' % (temp, outRaster)
    call(warp2)


