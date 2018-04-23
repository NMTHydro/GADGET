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
from config import load_configuration, cfg


def welcome():

    welcome_str = '''===========================================================================================
 _______  _______  ______   _______  _______  _______ 
|       ||   _   ||      | |       ||       ||       |
|    ___||  |_|  ||  _    ||    ___||    ___||_     _|
|   | __ |       || | |   ||   | __ |   |___   |   |  
|   ||  ||       || |_|   ||   ||  ||    ___|  |   |  
|   |_| ||   _   ||       ||   |_| ||   |___   |   |  
|_______||__| |__||______| |_______||_______|  |___|  


Developed by Peter Revelle 2016
Updated by Jake Ross, Gabe Parrish 2018
New Mexico Tech/New Mexico Bureau of Geology
==========================================================================================
'''

    print welcome_str


def main():
    welcome()
    load_configuration()

    # step 1

    if cfg['use_linke_turbidity']:
        from linke_turbidity.util import do_linke_turbidity
        do_linke_turbidity()

    # step 2


if __name__ == '__main__':
    main()
# ============= EOF =============================================
