s111
======
Python API and Utilities for Working with IHO S-111 Data Formats

Overview
--------

S-111 is an IHO standard outlining formats for storing and sending surface
water current data and metadata.

-   S-111 Data Coding Formats (DCF):

        1.  Time-series at fixed station
        2.  Regularly-gridded arrays
        3.  Ungeorectified Grid (i.e. irregular grid)
        4.  Time series for moving platform

*NOTE: Only DCF 2 & 3 are currently supported*

Example Usage
-------------

**Create an S-111 Data Coding Format 2 File:**
```python
import numpy
import datetime
from s100py import s111
speed = numpy.array([[0.34, 0.35, 0.35, 0.36, 0.37],
                     [0.4 , 0.41, 0.41, 0.42, 0.43],
                     [0.48, 0.49, 0.49, 0.49, 0.5],
                     [0.53, 0.54, 0.54, 0.55, 0.56],
                     [0.59, 0.6 , 0.61, 0.62, 0.63],
                     [0.34, 0.34, 0.35, 0.35, 0.37]])

direction = numpy.array([[189.1, 188.6, 188.1, 188.8, 189.7],
                         [192.6, 192.9, 192.8, 192.9, 193.1],
                         [193.6, 193.8, 193.9, 193.8, 193.6],
                         [192.8, 192.7, 192.9, 193.1, 193.3],
                         [195.0, 195.4, 195.7, 195.8, 195.7],
                         [194.3, 194.1, 194.0, 194.0, 194.0]])

grid_properties = {
        'maxx': -76.064,
        'minx': -76.03,
        'miny': 37.034,
        'maxy': 37.062,
        'cellsize_x': 0.005695343,
        'cellsize_y': 0.0056991577,
        'nx': 6,
        'ny': 5
}

# Example metadata
metadata = {
    'horizontalDatumReference': 'EPSG',
    'horizontalDatumValue': 4326,
    'metadata': f'MD_test_s111.XML',
    'epoch': 'G1762',
    'geographicIdentifier': 'RegionName',
    'speedUncertainty': -1.0, # Default or Unknown values
    'directionUncertainty': -1.0, # Default or Unknown values
    'verticalUncertainty': -1.0, # Default or Unknown values
    'horizontalPositionUncertainty': -1.0, # Default or Unknown values
    'timeUncertainty': -1.0, # Default or Unknown values
    'surfaceCurrentDepth': 0, 
    'depthTypeIndex': 2, # 2:Sea surface
    'commonPointRule': 3, # 3:high
    'interpolationType': 10, # 10:discrete
    'typeOfCurrentData': 6, # 6:Hydrodynamic model forecast (F)
    'methodCurrentsProduct': 'ROMS_Hydrodynamic_Model_Forecasts',
    'datetimeOfFirstRecord': '2021-01-07T13:00:00'

}

datetime_value = datetime.datetime(2021, 1, 7, 12, 0, 0)

data_coding_format = 2

update_meta = {
        'dateTimeOfLastRecord': '2021-01-07T13:00:00',
        'numberOfGroups': 1,
        'numberOfTimes': 1,
        'timeRecordInterval': 0,
        'num_instances': 1
    }

data_file = s111.utils.create_s111("test_s111.h5")

s111.utils.add_metadata(metadata, data_file)
s111.utils.add_data_from_arrays(speed, direction, data_file, grid_properties, datetime_value, data_coding_format)
s111.utils.update_metadata(data_file, grid_properties, update_meta)

s111.utils.write_data_file(data_file)
```

For S-111 Developers
--------------------
- [S-111 Module Documentation](https://s100py.readthedocs.io/en/latest/s111.html#s111-module-docs)
- [S-100 Module Documentation](https://s100py.readthedocs.io/en/latest/s100.html)

Authors
-------

-   Erin Nagel (UCAR), <erin.nagel@noaa.gov>
-   Barry Gallagher (NOAA), <barry.gallagher@noaa.gov>
-   Jason Greenlaw 



