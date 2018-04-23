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
from datetime import datetime

import osr
import yaml

from map_utils import read_map

cfg = {
    'use_linke_turbidity': True,
    'startday': datetime(2000, 1, 1, 0),
    'endday': datetime(2000, 12, 31, 0),
    'basemap': 'data/01_longlat_wgs84.tif',
    'linke_turbidity_dir': 'data',
    'linke_output_dir': 'output'
}


def load_configuration(path=None):
    if path is None:
        path = os.path.join(os.path.expanduser('~'), '.gadget.yml')

    if os.path.isfile(path):
        with open(path, 'r') as rfile:
            yd = yaml.load(rfile)

            cfg.update(yd)


def load_basemap():
    keys = ('resX', 'resY', 'cols', 'rows', 'lon', 'lat', 'linke', 'prj', 'fill')
    vals = read_map(cfg['basemap'], 'Gtiff')
    cfg['basemap_dict'] = bm = dict(zip(keys, vals))

    srs = osr.SpatialReference(bm['prj'])
    sr_wkt = srs.ExportToWkt()
    cfg['linke_sr_wkt'] = sr_wkt
# ============= EOF =============================================
