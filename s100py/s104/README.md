s104
======
Python Utilities and API for Working with IHO S-104 Data Formats

Overview
--------

S-104 is an IHO standard outlining formats for storing and sending
water level data and metadata for surface navigation.

Example Usage
-------------

**Create an S-104 DCF2 File:**
```python
import numpy
import datetime
from s100py import s104
water_level_height = numpy.array([[ 0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,
                                     0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,
                                     0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,
                                     0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,
                                     0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,
                                     0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,
                                     0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,
                                     0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,  0.34,
                                     0.34,  0.35,  0.35,  0.35,  0.35,  0.35,  0.35,  0.35,  0.35,
                                     0.35,  0.35,  0.35,  0.35,  0.35,  0.35,  0.35,  0.35,  0.35]])

water_level_trend = numpy.array([[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3,
                                    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]])

grid_properties = {
        'maxx': -75.30278,
        'minx': -75.59722,
        'miny': 37.202778,
        'maxy': 37.497223,
        'cellsize_x': 0.005554199,
        'cellsize_y': 0.005558014,
        'nx': 1,
        'ny': 220
}

# Example metadata
metadata = {
    'horizontalCRS': 3855,
    'metadata': f'MD_test_s104.XML',
    'geographicIdentifier': 'Region',
    'waterLevelHeightUncertainty': -1.0, # Default or Unknown values
    'verticalUncertainty': -1.0, # Default or Unknown values
    'horizontalPositionUncertainty': -1.0, # Default or Unknown values
    'timeUncertainty': -1.0, # Default or Unknown values
    'waterLevelTrendThreshold': 0.02,
    'verticalCS': 6499,
    'verticalCoordinateBase': 2,
    'verticalDatumReference': 2,
    'verticalDatum': 1027,
    'commonPointRule': 4, # 'all'
    'interpolationType': 10, #'discrete'
    'typeOfCurrentData': 5, # Hydrodynamic model forecast (F)
    'methodCurrentsProduct': 'ADCIRC_Hydrodynamic_Model_Forecasts',
    'datetimeOfFirstRecord': '2020-09-26T16:00:00'

}

datetime_value = datetime.datetime(2020, 9, 26, 15, 0, 0)

data_coding_format = 2

update_meta = {
        'dateTimeOfLastRecord': '2020-09-26T17:00:00',
        'numberOfGroups': 2,
        'numberOfTimes': 2,
        'timeRecordInterval': 3600,
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
-   Jason Greenlaw (ERT), <jason.greenlaw@noaa.gov>
-   Barry Gallagher (NOAA), <barry.gallagher@noaa.gov>



