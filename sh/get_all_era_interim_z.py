#!/usr/bin/env python

#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "1979-01-01/to/2017-05-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "levelist": "500",
    "levtype": "pl",
    "param": "129.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "format" : "netcdf",
    "type": "an",
    "target": "/barnes-scratch/sbrey/era_interim_nc_6_hourly/z_all.nc",
})
