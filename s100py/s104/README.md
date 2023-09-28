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

water_level_height_001 = numpy.array([
                                    [    0.2 ,     0.2 ,     0.2 ,     0.19,     0.19,     0.19,
                                         0.19,     0.19,     0.19,     0.19,     0.19,     0.18,  0.18],
                                    [    0.2 ,     0.2 ,     0.2 ,     0.2 ,     0.19,     0.19,
                                         0.19,     0.19,     0.19,     0.19,     0.19,     0.18,  0.18],
                                    [    0.22,     0.2 ,     0.2 ,     0.2 ,     0.19,     0.19,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [    0.23,     0.21,     0.2 ,     0.2 ,     0.2 ,     0.19,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [    0.24,     0.21,     0.21,     0.2 ,     0.2 ,     0.19,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [-9999.  ,     0.21,     0.21,     0.2 ,     0.2 ,     0.19,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [-9999.  ,     0.23,     0.21,     0.21,     0.2 ,     0.2 ,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [-9999.  ,     0.25,     0.21,     0.21,     0.2 ,     0.2 ,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [-9999.  , -9999.  ,     0.22,     0.21,     0.2 ,     0.2 ,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [-9999.  , -9999.  ,     0.23,     0.21,     0.2 ,     0.19,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [-9999.  ,     0.27,     0.22,     0.2 ,     0.2 ,     0.19,
                                         0.19,     0.18,     0.18,     0.18,     0.18,     0.18,  0.18],
                                    [    0.26,     0.23,     0.2 ,     0.2 ,     0.19,     0.19,
                                         0.19,     0.18,     0.18,     0.18,     0.18,     0.18,  0.18],
                                    [    0.23,     0.2 ,     0.2 ,     0.2 ,     0.19,     0.19,
                                         0.19,     0.18,     0.18,     0.18,     0.18,     0.18,  0.18],
                                    [    0.21,     0.2 ,     0.2 ,     0.19,     0.19,     0.19,
                                         0.18,     0.18,     0.18,     0.18,     0.18,     0.18,  0.18],
                                    [    0.2 ,     0.2 ,     0.19,     0.19,     0.19,     0.19,
                                         0.18,     0.18,     0.18,     0.18,     0.18,     0.18,  0.18],
                                    [    0.19,     0.19,     0.19,     0.19,     0.19,     0.18,
                                         0.18,     0.18,     0.18,     0.18,     0.18,     0.18, 0.18]])

water_level_trend_001 = numpy.array([
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]])

water_level_height_002 = numpy.array([
                                    [    0.18,     0.17,     0.17,     0.17,     0.17,     0.17,
                                         0.16,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [    0.18,     0.18,     0.17,     0.17,     0.17,     0.17,
                                         0.16,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [    0.2 ,     0.18,     0.18,     0.17,     0.17,     0.17,
                                         0.16,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [    0.21,     0.18,     0.18,     0.18,     0.17,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [    0.21,     0.19,     0.18,     0.18,     0.17,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [-9999.  ,     0.19,     0.18,     0.18,     0.18,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [-9999.  ,     0.21,     0.19,     0.18,     0.18,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [-9999.  ,     0.22,     0.19,     0.18,     0.18,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [-9999.  , -9999.  ,     0.2 ,     0.18,     0.18,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [-9999.  , -9999.  ,     0.2 ,     0.18,     0.18,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [-9999.  ,     0.25,     0.2 ,     0.18,     0.17,     0.17,
                                         0.16,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [    0.24,     0.21,     0.18,     0.18,     0.17,     0.17,
                                         0.16,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [    0.21,     0.18,     0.18,     0.17,     0.17,     0.17,
                                         0.16,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [    0.19,     0.18,     0.17,     0.17,     0.17,     0.16,
                                         0.16,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [    0.18,     0.17,     0.17,     0.17,     0.17,     0.16,
                                         0.16,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [    0.17,     0.17,     0.17,     0.17,     0.16,     0.16,
                                         0.16,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16]])

grid_properties = {
        'maxx': -169.642699,
        'minx': -169.812302,
        'miny': -19.165940,
        'maxy': -18.944464,
        'cellsize_x': 0.01384226425,
        'cellsize_y': 0.0130463781538463,
        'nx': 16,
        'ny': 13
}

datetime_forecast_issuance = datetime.datetime(2021, 9, 1, 0, 0, 0)

datetime_interval = datetime.timedelta(seconds=3600)

# Example metadata
metadata = {
    'horizontalCRS': 4362, #EPSG code
    'geographicIdentifier': 'RegionName',
    'waterLevelHeightUncertainty': -1.0, # Default or Unknown values
    'verticalUncertainty': -1.0, # Default or Unknown values
    'horizontalPositionUncertainty': -1.0, # Default or Unknown values
    'waterLevelTrendThreshold': 0.2,
    'verticalCS': 6499, # EPSG code
    'verticalDatumReference': 1, # 2:EPSG
    'verticalDatum': 12, # EPSG code
    'commonPointRule': 4, # 4:all
    'interpolationType': 10, # 10:discrete
    'dataDynamicity': 5, # 5:Hydrodynamic model forecast (F)
    'methodWaterLevelProduct': 'ADCIRC_Hydrodynamic_Model_Forecasts',
    'datetimeOfFirstRecord': '20210901T010000Z',
    'trendInterval': 60, # minutes
    'datasetDeliveryInterval': 'PT6H',
    'issueDateTime': datetime_forecast_issuance
}

data_coding_format = 2

update_meta = {
        'dateTimeOfLastRecord': '20210901T020000Z',
        'numberOfGroups': 2,
        'numberOfTimes': 2,
        'timeRecordInterval': 3600,
        'num_instances': 1
    }

data_file = s104.utils.create_s104("test_s104.h5")

s104.utils.add_metadata(metadata, data_file)
data_series_time_001 = datetime_forecast_issuance + datetime_interval
s104.utils.add_data_from_arrays(water_level_height_001, water_level_trend_001, data_file, grid_properties, data_series_time_001, data_coding_format)
data_series_time_002 = data_series_time_001 + datetime_interval

trend = numpy.round((water_level_height_002 - water_level_height_001), decimals=2)

water_level_trend_002 = numpy.where(( -1 * metadata['waterLevelTrendThreshold'] < trend) &
                                (trend < metadata['waterLevelTrendThreshold']),3,
                                numpy.where(trend >= metadata['waterLevelTrendThreshold'], 2,
                                    numpy.where( trend <= -1 * metadata['waterLevelTrendThreshold'], 1,numpy.any(trend))))

s104.utils.add_data_from_arrays(water_level_height_002, water_level_trend_002, data_file, grid_properties, data_series_time_002, data_coding_format)

s104.utils.update_metadata(data_file, grid_properties, update_meta)

s104.utils.write_data_file(data_file)
```

For S-104 Developers
--------------------
- [S-104 Module Documentation](https://s100py.readthedocs.io/en/latest/s104.html#s104-module-docs)
- [S-100 Module Documentation](https://s100py.readthedocs.io/en/latest/s100.html)

Authors
-------

-   Erin Nagel (ERT), <erin.nagel@noaa.gov>
-   Barry Gallagher (NOAA), <barry.gallagher@noaa.gov>
-   Jason Greenlaw



