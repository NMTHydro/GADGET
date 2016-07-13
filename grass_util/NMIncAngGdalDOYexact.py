import grass.script as grass
import grass.script.setup as gsetup
import datetime
import os, fnmatch
from dateutil import rrule

search = 'F:\\Grass_newdem\\NewMexico\\NM_buffer\\cell'
path = 'F:\\StatewideWater\\Grass_scripts\\NM_wide'
path2 = 'I:\\StatewideWater\\NM_IncAng'
dest = path2 + '\\UTM'

os.chdir(path)


startday = datetime.datetime(2000,1,1,0)
endday = datetime.datetime(2000,12,31,22)

starthour = datetime.datetime(2000,12,31,5)
endhour = datetime.datetime(2000,12,31,22)


a = open("NMoutIncang.txt", "wb")
count = 0

for day in rrule.rrule(rrule.DAILY, dtstart=startday, until=endday):
    for hour in rrule.rrule(rrule.HOURLY, dtstart=starthour, until=endhour):
        strcount = str(count)
        strdoy = int(day.strftime('%j'))
        strhour = int(hour.strftime('%H'))
        day3 = '%03d' % strdoy
        hour4 = '%04d' % strhour
        
        grass_file = 'nmfull_incangle' + day3 + '_' + hour4 
        dates = 'NMIncang' + '_' + day.strftime('%j') + '_' + hour.strftime('%H') + '.tif'
        a.write(str(grass_file + '\n'))


        inRaster = grass_file
        #convert file to str to add destination, file extension
        out = dest + '/' + dates
        grass.run_command('r.out.gdal', flags="c", input=inRaster, output=out, nodata="0", format='GTiff', type='Float32')
a.close()


