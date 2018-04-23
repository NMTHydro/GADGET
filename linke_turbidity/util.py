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

from dateutil import rrule
from dateutil.relativedelta import relativedelta
from numpy import meshgrid, r_, ones
from numpy.ma import arange
from scipy import ndimage

from config import cfg
from map_utils import read_map, write_map


def extract_dates(day):

    this_month = day.month
    this_monthobj = datetime.datetime(day.year, this_month, day.day, 0)

    if day.day < 16:
        in_one_month = this_monthobj.replace(day=1)
        thisyear = day.year

        if day.month == 1:
            thisyear -= 1

        this_month = this_monthobj-relativedelta(months=1)
        this_month = this_month.month
        this_monthobj = datetime.datetime(thisyear, this_month, day.day, 0)

    elif day.day >= 16:
        thisyear = day.year
        this_monthobj = datetime.datetime(thisyear, this_month, day.day, 0)

        if day.month == 12:
            thisyear += 1

        in_one_month = this_monthobj+relativedelta(months=1)
        in_one_month.replace(year=thisyear)

    print('this month', this_month)
    print('in one month', in_one_month)

    fmt = '{:02d}'
    month = fmt.format(this_month)
    str_nextmonth = int(in_one_month.strftime('%m'))
    mid_thismonth = this_monthobj.replace(day=1) + datetime.timedelta(days=14)
    print('middle of this month', mid_thismonth)
    mid_nextmonth = in_one_month + datetime.timedelta(days=14)
    print('middle of next month', mid_nextmonth)
    next_month = fmt.format(str_nextmonth)
    return month, next_month, mid_thismonth, mid_nextmonth


def get_paths(month, next_month):
    current_month_name = '{}_longlat_wgs84.tif'.format(month)
    next_month_name = '{}_longlat_wgs84.tif'.format(next_month)

    cp = os.path.join(cfg['linke_turbidity_dir'], current_month_name)
    np = os.path.join(cfg['linke_turbidity_dir'], next_month_name)
    return cp, np


def time_interpolation(day, lat, lon):
    month, next_month, mid_month, mid_next_month = extract_dates(day)
    current_path, next_path = get_paths(month, next_month)

    current = read_map(current_path)
    next_map = read_map(next_path)

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
    print('day', day)
    print('days difference from mid next month', days_dif)
    print('out of max days difference from mid this month', max_days_diff)
    interp = 1 - (days_dif / max_days_diff)  # 1 = weight completely next month linke values, 0 = previous month
    print('interp ratio between monthly images', interp)

    # 0.5 corresponds to half way between arr1 and arr2
    coordinates = ones((r, c)) * interp, x, y

    # given arr interpolate at coordinates
    newarr = ndimage.map_coordinates(arr, coordinates, order=2)

    return newarr


def do_linke_turbidity():
    for day in rrule.rrule(rrule.DAILY,
                           dtstart=cfg.get('startday'),
                           until=cfg['endday']):
        bm = cfg['basemap_dict']

        lat, lon = bm['lat'], bm['lon']
        linke_daily = time_interpolation(day, lat, lon)

        nr = day.strftime('%j')
        daily_doy = 'linkewgs84_{}.tif'.format(nr)
        outname = os.path.join(cfg['linke_output_dir'], daily_doy)
        write_map(outname, 'Gtiff', lon, lat, linke_daily, cfg['linke_sr_wkt'], bm['fill'])


def project_lat_long():
    pass
# ============= EOF =============================================
