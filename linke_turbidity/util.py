# ===============================================================================
# Copyright 2018 ross
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ===============================================================================
import datetime
import os
from subprocess import call

from dateutil import rrule
from dateutil.relativedelta import relativedelta
from numpy import meshgrid, r_, ones
from numpy.ma import arange
from scipy import ndimage

from config import cfg
from map_utils import read_raster, write_raster, find_rasters


def extract_dates(day, verbose=False):
    this_month = day.month
    this_monthobj = datetime.datetime(day.year, this_month, day.day, 0)

    if day.day < 16:
        in_one_month = this_monthobj.replace(day=1)
        thisyear = day.year

        if day.month == 1:
            thisyear -= 1

        this_month = this_monthobj - relativedelta(months=1)
        this_month = this_month.month
        this_monthobj = datetime.datetime(thisyear, this_month, day.day, 0)

    elif day.day >= 16:
        thisyear = day.year
        this_monthobj = datetime.datetime(thisyear, this_month, day.day, 0)

        if day.month == 12:
            thisyear += 1

        in_one_month = this_monthobj + relativedelta(months=1)
        in_one_month.replace(year=thisyear)
    if verbose:
        print('this month', this_month)
        print('in one month', in_one_month)

    fmt = '{:02d}'
    month = fmt.format(this_month)
    str_nextmonth = int(in_one_month.strftime('%m'))
    mid_thismonth = this_monthobj.replace(day=1) + datetime.timedelta(days=14)

    mid_nextmonth = in_one_month + datetime.timedelta(days=14)
    if verbose:
        print('middle of this month', mid_thismonth)
        print('middle of next month', mid_nextmonth)
    next_month = fmt.format(str_nextmonth)
    return month, next_month, mid_thismonth, mid_nextmonth


def get_paths(month, next_month):
    current_month_name = '{}_longlat_wgs84.tif'.format(month)
    next_month_name = '{}_longlat_wgs84.tif'.format(next_month)

    cp = os.path.join(cfg['linke_turbidity_dir'], current_month_name)
    np = os.path.join(cfg['linke_turbidity_dir'], next_month_name)
    return cp, np


def time_interpolation(day, verbose=False):
    month, next_month, mid_month, mid_next_month = extract_dates(day)
    current_path, next_path = get_paths(month, next_month)

    current = read_raster(current_path)
    next_map = read_raster(next_path)

    linke = current[6] / 20.
    nlinke = next_map[6] / 20.
    c, r = 247, 66

    latitude_index = arange(r)
    longitude_index = arange(c)

    # rejoin Linke, LinkeNext into a single array of shape (2, 2160, 4320)
    arr = r_['0,3', linke, nlinke]
    # print('arr.shape',arr.shape)

    # define the grid coordinates where you want to interpolate
    y, x = meshgrid(longitude_index, latitude_index)

    # Setup time variables for interpolation
    days_dif = mid_next_month - day
    days_dif = days_dif.days
    max_days_diff = mid_next_month - mid_month
    max_days_diff = max_days_diff.days

    interp = 1 - (days_dif / max_days_diff)  # 1 = weight completely next month linke values, 0 = previous month
    if verbose:
        print('day', day)
        print('days difference from mid next month', days_dif)
        print('out of max days difference from mid this month', max_days_diff)
        print('interp ratio between monthly images', interp)

    # 0.5 corresponds to half way between arr1 and arr2
    coordinates = ones((r, c)) * interp, x, y

    # given arr interpolate at coordinates
    newarr = ndimage.map_coordinates(arr, coordinates, order=2)

    return newarr


def do_linke_turbidity():
    print 'Doing linke turbidity'

    for day in rrule.rrule(rrule.DAILY,
                           dtstart=cfg.get('startday'),
                           until=cfg['endday']):
        bm = cfg['basemap_dict']

        lat, lon = bm['lat'], bm['lon']
        linke_daily = time_interpolation(day)

        nr = day.strftime('%j')
        daily_doy = 'linkewgs84_{}.tif'.format(nr)
        outname = os.path.join(cfg['linke_output_dir'], daily_doy)
        write_raster(outname, 'Gtiff', lon, lat, linke_daily, cfg['linke_sr_wkt'], bm['fill'], verbose=False)

    print 'Linke turbidity calculation finished.'


def project_lat_long():
    root = cfg['linke_output_dir']
    out_root = cfg['linke_output_proj_dir']

    for raster in find_rasters(root):
        out = os.path.join(out_root, raster)
        raster = os.path.join(root, raster)
        gcp = '-gcp 0.0 0.0 -180.0 90.0 -gcp 0.0 2160.0 -180.0 -90.0 -gcp 4320.0 0.0 180.0 90.0'
        cmd = 'gdal_translate -of GTiff -a_srs EPSG:4326 {} -r "near" {} {}' % (gcp, raster, out)

        call(cmd, shell=True)


def do_warp():
    root = cfg['linke_output_proj_dir']
    out_root = cfg['linke_output_warp_dir']

    ssrs = 'EPSG: 4326'
    tsrs = 'EPSG: 4326'
    te = '-110.2098918542965151 31.6254675894567470 -109.5432356466024828 31.8337976622254928'
    tr = '0.041666 0.041666'
    ts = '6 18'
    r = 'cubic'
    srcnodata = -3.40282346639e+038

    warp_cmd = 'gdalwarp -overwrite -s_srs {} -t_srs {} -te {} -r "{}"'.format(ssrs, tsrs, te, r)
    for raster in find_rasters(root):
        temp = os.path.join(root, 'temp.tif')
        in_raster = os.path.join(root, raster)
        out_raster = os.path.join(out_root, raster)

        cmd = '{} -ts {} -multi -srcnodata "{}" -dstnodata -999 {} {}'.format(warp_cmd, ts, srcnodata, in_raster, temp)
        call(cmd)

        cmd = '{} -tr {} -multi -srcnodata -999 -dstnodata -999 {} {}'.format(warp_cmd, tr, temp, out_raster)
        call(cmd)

    os.remove(temp)


def do_downsize():
    pass

# ============= EOF =============================================
