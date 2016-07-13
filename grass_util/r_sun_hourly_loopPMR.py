from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import imagery as i
from grass.script.core import run_command

dem="nmfullDemUtm"
aspect="nmfull_aspect"
slope="nmfull_slope"
horizon="NMhorRadians"
horizonsteps="10"


def format_day(day):
    return '%03d' % day

day = range(1,366)

for d in day:
    binarybeam="nmfull_shade" + format_day(d)
    incangle="nmfull_incangle" + format_day(d)

    run_command("r.sun.hourlyPMR",overwrite = "True",flags="b",
					 elevation=dem,
					 aspect=aspect,
					 slope=slope,
                          albedo_value="0.23",
					 start_time="5",
					 end_time="22",
                          civil_time="-7",
					 day=d,
					 year=2000,
					 beam_rad_basename=binarybeam,
                                         horizon_basename=horizon,
                                         horizon_step=horizonsteps,
                                         incidout_basename=incangle,
                                         nprocs=4)

