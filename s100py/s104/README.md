s104
======
Python API and Utilities for Working with IHO S-104 Data Formats

Overview
--------

S-104 is an IHO standard outlining formats for storing and sending
water level data and metadata for surface navigation.

-   S-104 Data Coding Formats (DCF):

        1.  Time-series at fixed station
        2.  Regularly-gridded arrays
        3.  Ungeorectified Grid
        7.  TIN
        8.  Time Series at fixed stations (stationwise)

*NOTE: Only DCF 2 is currently supported*

Example Usage
-------------

**Create an S-104 DCF2 File:**
```python
import numpy
import datetime
from s100py import s104
water_level_height = numpy.array([[ 0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,
                                    0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34, 0.34,  0.34],
                                  [ 0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34, 0.34,  0.34,  0.34,
                                    0.34,  0.34,  0.34,  0.34,  0.34,  0.34, 0.34,  0.34,  0.34,  0.34],
                                  [ 0.34,  0.34,  0.34, 0.34,  0.35,  0.35,  0.35,  0.35,  0.35,  0.35,
                                    0.35,  0.35, 0.35,  0.35,  0.35,  0.35,  0.35,  0.35,  0.35,  0.35]])

water_level_trend = numpy.array([[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

grid_properties = {
        'maxx': 134.175,
        'minx': 134.216,
        'miny': 6.982,
        'maxy': 7.254,
        'cellsize_x': 0.013595581,
        'cellsize_y': 0.013597488,
        'nx': 3,
        'ny': 20
}

# Example metadata
metadata = {
    'horizontalCRS': 4362, #EPSG code
    'metadata': f'MD_test_s104.XML',
    'geographicIdentifier': 'RegionName',
    'waterLevelHeightUncertainty': -1.0, # Default or Unknown values
    'verticalUncertainty': -1.0, # Default or Unknown values
    'horizontalPositionUncertainty': -1.0, # Default or Unknown values
    'timeUncertainty': -1.0, # Default or Unknown values
    'waterLevelTrendThreshold': 0.2,
    'verticalCS': 6499, # EPSG code
    'verticalCoordinateBase': 2, # 2:Vertical Datum
    'verticalDatumReference': 2, # 2:EPSG
    'verticalDatum': 1027, # EPSG code
    'commonPointRule': 4, # 4:all
    'interpolationType': 10, # 10:discrete
    'typeOfWaterLevelData': 5, # 5:Hydrodynamic model forecast (F)
    'methodWaterLevelProduct': 'ADCIRC_Hydrodynamic_Model_Forecasts',
    'datetimeOfFirstRecord': '2020-09-26T16:00:00'

}

datetime_value = datetime.datetime(2020, 9, 26, 15, 0, 0)

data_coding_format = 2

update_meta = {
        'dateTimeOfLastRecord': '2020-09-26T16:00:00',
        'numberOfGroups': 1,
        'numberOfTimes': 1,
        'timeRecordInterval': 0,
        'num_instances': 1
    }

data_file = s104.utils.create_s104("test_s104.h5")

s104.utils.add_metadata(metadata, data_file)
s104.utils.add_data_from_arrays(water_level_height, water_level_trend, data_file, grid_properties, datetime_value, data_coding_format)
s104.utils.update_metadata(data_file, grid_properties, update_meta)

s104.utils.write_data_file(data_file)
```

For S-104 Developers
--------------------
- [S-104 Module Documentation](https://s100py.readthedocs.io/en/latest/s104.html#s104-module-docs)
- [S-100 Module Documentation](https://s100py.readthedocs.io/en/latest/s100.html)

Authors
-------

-   Erin Nagel (UCAR), <erin.nagel@noaa.gov>
-   Barry Gallagher (NOAA), <barry.gallagher@noaa.gov>
-   Jason Greenlaw



