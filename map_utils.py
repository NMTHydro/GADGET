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
import os
import sys
from math import isnan

import gdal
from numpy import linspace
gdal.AllRegister()


def read_map(path, file_format):
    """
    read geographical file into memory

    :param path: Name of the file to read
    :param file_format: Gdal format string
    :param logger: logger object
    :return: resX, resY, cols, rows, x, y, data, FillVal
    """

    # Open file for binary-reading

    if not os.path.isfile(path):
        print '{} is not a valid file. Exiting'.format(path)
        sys.exit(1)

    mapFormat = gdal.GetDriverByName(file_format)
    mapFormat.Register()
    ds = gdal.Open(path)
    prj = ds.GetProjection()

    # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX = geotrans[1]
    resY = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = linspace(originX + resX / 2, originX + resX / 2 + resX * (cols - 1), cols)
    y = linspace(originY + resY / 2, originY + resY / 2 + resY * (rows - 1), rows)

    # Retrieve raster
    raster_band = ds.GetRasterBand(1)  # there's only 1 band, starting from 1
    # print(x)
    # print(y)

    data = raster_band.ReadAsArray(0, 0, cols, rows)
    fill_value = raster_band.GetNoDataValue()
    del ds

    return resX, resY, cols, rows, x, y, data, prj, fill_value


def add_extension(path, ext='.tif'):
    if not path.endswith('.tif'):
        if not ext.startswith('.'):
            ext = '.{}'.format(ext)

        path = '{}{}'.format(path, ext)

    return path


def write_map(path, file_format, x, y, data, prj, fill, verbose=False):
    """
    Write geographical data into file. Also replave NaN bu FillVall

    :param path:
    :param file_format:
    :param x:
    :param y:
    :param data:
    :param fill:
    :return:
    """

    driver1 = gdal.GetDriverByName('GTiff')
    # driver2 = gdal.GetDriverByName(file_format)

    data[isnan(data)] = fill
    # Processing
    if verbose:
        print 'Writing to temporary file {}.tif'.format(path)
        print "Output format: {}".format(file_format)

    path = add_extension(path)
    # Create Output filename from (FEWS) product name and date and open for writing
    temp_dataset = driver1.Create(path, data.shape[1], data.shape[0], 1, gdal.GDT_Float32)

    # Give georeferences
    xul = x[0] - (x[1] - x[0]) / 2
    yul = y[0] + (y[0] - y[1]) / 2

    temp_dataset.SetGeoTransform([xul, x[1] - x[0], 0, yul, 0, y[1] - y[0]])
    temp_dataset.SetProjection(prj)
    # get rasterband entry
    temp_band = temp_dataset.GetRasterBand(1)
    # fill rasterband with array
    temp_band.WriteArray(data, 0, 0)
    temp_band.FlushCache()
    temp_band.SetNoDataValue(fill)
    # Create data to write to correct format (supported by 'CreateCopy')
    if verbose:
        print 'Writing to {}'.format(path)

    # outDataset = driver2.CreateCopy(path, temp_dataset, 0)
    # temp_dataset = None
    # outDataset = None
    #
    # if verbose:
    #     print 'Removing temporary file ' + path + '.tif'
    # os.remove(path + '.tif');
    #
    # if verbose:
    #     print 'Writing to ' + path + ' is done!'

# ============= EOF =============================================
