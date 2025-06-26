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
        8.  Time Series at fixed stations (stationwise)
- 
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
        'maxx': -75.30278,
        'minx': -75.59722,
        'miny': 37.202778,
        'maxy': 37.497223,
        'cellsize_x': 0.005695343,
        'cellsize_y': 0.0056991577,
        'nx': 5,
        'ny': 6
}

datetime_forecast_issuance = datetime.datetime(2025, 6, 5, 12, 0, 0)


datetime_interval = datetime.timedelta(seconds=0)

# Example metadata
metadata = {
    'horizontalCRS': 4326,
    'geographicIdentifier': 'RegionName',
    'speedUncertainty': -1.0, # Default or Unknown values
    'directionUncertainty': -1.0, # Default or Unknown values
    'verticalUncertainty': -1.0, # Default or Unknown values
    'horizontalPositionUncertainty': -1.0, # Default or Unknown values
    'surfaceCurrentDepth': 0, 
    'depthTypeIndex': 1, # 1:heightOrDepth
    'commonPointRule': 3, # 3:high
    'interpolationType': 10, # 10:discrete
    'dataDynamicity': 5, # 5:Hydrodynamic model forecast (F)
    'methodCurrentsProduct': 'Hydrodynamic_Model_Forecasts',
    'dateTimeOfFirstRecord': '20250605T130000Z',
    'verticalCS': 6498,
    'verticalDatumReference': 1,
    'verticalCoordinateBase': 2, # Only allowed valued 2:verticalDatum
    'verticalDatum': 48,
    'datasetDeliveryInterval': 'PT6H',
    'issueDate': datetime_forecast_issuance,
    'issueTime': datetime_forecast_issuance
}

data_coding_format = 2

data_file = s111.utils.create_s111("test_s111.h5", 2)

s111.utils.add_metadata(metadata, data_file)
s111.utils.add_surface_current_instance(data_file)
data_series_time_001 = datetime_forecast_issuance + datetime.timedelta(hours=1)
s111.utils.add_data_from_arrays(speed, direction, data_file, grid_properties, data_series_time_001, data_coding_format)
s111.utils.update_metadata(data_file, grid_properties, metadata)

s111.utils.write_data_file(data_file)
```

For S-111 Developers
--------------------
- [S-111 Module Documentation](https://s100py.readthedocs.io/en/latest/s111.html#s111-module-docs)
- [S-100 Module Documentation](https://s100py.readthedocs.io/en/latest/s100.html)

Authors
-------

-   Erin Nagel (ERT), <erin.nagel@noaa.gov>
-   Barry Gallagher (NOAA), <barry.gallagher@noaa.gov>
-   Jason Greenlaw 



