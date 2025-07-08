s104
======
Python API and Utilities for Working with IHO S-104 Data Formats

Overview
--------

S-104 is an IHO standard outlining formats for storing and sending
water level data and metadata for surface navigation.

-   S-104 Data Coding Formats (DCF):

        2.  Regularly-gridded arrays

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
                                    [-9999.0,     0.21,     0.21,     0.2 ,     0.2 ,     0.19,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [-9999.0,     0.23,     0.21,     0.21,     0.2 ,     0.2 ,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [-9999.0,     0.25,     0.21,     0.21,     0.2 ,     0.2 ,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [-9999.0, -9999.0,     0.22,     0.21,     0.2 ,     0.2 ,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [-9999.0, -9999.0,     0.23,     0.21,     0.2 ,     0.19,
                                         0.19,     0.19,     0.19,     0.19,     0.18,     0.18,  0.18],
                                    [-9999.0,     0.27,     0.22,     0.2 ,     0.2 ,     0.19,
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
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

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
                                    [-9999.0,     0.19,     0.18,     0.18,     0.18,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [-9999.0,     0.21,     0.19,     0.18,     0.18,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [-9999.0,     0.22,     0.19,     0.18,     0.18,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [-9999.0, -9999.0,     0.2 ,     0.18,     0.18,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [-9999.0, -9999.0,     0.2 ,     0.18,     0.18,     0.17,
                                         0.17,     0.16,     0.16,     0.16,     0.16,     0.16, 0.16],
                                    [-9999.0 ,     0.25,     0.2 ,     0.18,     0.17,     0.17,
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
        'nx': 13,
        'ny': 16
}

datetime_forecast_issuance = datetime.datetime(2021, 9, 1, 0, 0, 0)

datetime_interval = datetime.timedelta(seconds=3600)

# Example metadata
# Optional attributes are denoted below and can be removed
metadata = {
    'horizontalCRS': 4326, # EPSG code
    'waterLevelHeightUncertainty': -1.0, # -1.0 (unknown) or positive value (m)
    'horizontalPositionUncertainty': -1.0, # -1.0 (unknown) or positive value (m)
    'verticalUncertainty': -1.0, # -1.0 (unknown) or positive value (m)
    'waterLevelTrendThreshold': 0.2,
    'verticalCS': 6499, # Height-Meters–Orientation Up
    'verticalDatumReference': 1, # 1:S100_VerticalAndSoundingDatum
    'verticalCoordinateBase': 2, # 2:verticalDatum (Only value allowed)
    'verticalDatum': 12, # 12:MLLW
    'commonPointRule': 4, # 4:all
    'interpolationType': 1, # 1:nearestneighbor (Only value allowed)
    'dataDynamicity': 5, # 5:Hydrodynamic model forecast (F)
    'issueDateTime': datetime_forecast_issuance, # All times are in UTC, DateTime format
    'dateTimeOfFirstRecord': '20210901T010000Z', # All times are in UTC, DateTime format
    'geographicIdentifier': 'Add Region Name', # Optional
    'methodWaterLevelProduct': 'ADCIRC_Hydrodynamic_Model_Forecasts', # Optional
    'trendInterval': 60, # Optional (Minutes)
    'datasetDeliveryInterval': 'PT6H', #Optional (ISO 8601 duration, format `YnMnDTnHnMnS`)
    'epoch': '2005.0', # Optional
    'verticalDatumEpoch': 'NOAA_NTDE_1983-2001', # Optional
    'timeUncertainty': -1.0 # Optional
}

data_coding_format = 2 # Only value allowed in S-104 Ed 2.0

data_file = s104.utils.create_s104("test_s104.h5", 2)

s104.utils.add_metadata(metadata, data_file)
s104.utils.add_water_level_instance(data_file)
data_series_time_001 = datetime_forecast_issuance + datetime_interval
s104.utils.add_data_from_arrays(water_level_height_001, water_level_trend_001, data_file, grid_properties, data_series_time_001, data_coding_format)
data_series_time_002 = data_series_time_001 + datetime_interval

water_level_difference = numpy.round((water_level_height_002 - water_level_height_001), decimals=2)
water_level_difference = numpy.ma.masked_where((water_level_height_002 == -9999.0) | (water_level_height_001 == -9999.0), water_level_difference)

#Trend is determined using S-104 Edition 2.0.0 Annex A-2 Feature Attributes Section 2: Water Level Trend (waterLevelTrend)
water_level_trend_002 = numpy.where((-metadata['waterLevelTrendThreshold'] < water_level_difference) & (water_level_difference < metadata['waterLevelTrendThreshold'] ), 3,
                 numpy.where(water_level_difference >= metadata['waterLevelTrendThreshold'] , 2,
                         numpy.where(water_level_difference <= -metadata['waterLevelTrendThreshold'] , 1, 0)))

water_level_trend_002 = numpy.ma.filled(water_level_trend_002, 0)

s104.utils.add_data_from_arrays(water_level_height_002, water_level_trend_002, data_file, grid_properties, data_series_time_002, data_coding_format)

s104.utils.update_metadata(data_file, grid_properties, metadata)

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



