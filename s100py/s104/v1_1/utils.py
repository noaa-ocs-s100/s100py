""" Functions to create S104 files and populate with data from other sources

"""
import logging
import sys
import datetime

import numpy

from ...s1xx import s1xx_sequence
from .api import S104File, FILLVALUE_HEIGHT, FILLVALUE_TREND, S104Exception


def _get_S104File(output_file):
    """
    Small helper function to convert the
    output_file parameter into a S104File
    """
    if isinstance(output_file, S104File):
        data_file = output_file
    else:
        try:
            data_file = S104File(output_file, "w")
        except TypeError as typeerr:
            msg = "Failed to create S104File using {}".format(str(output_file))
            logging.error(msg)
            raise type(typeerr)(msg).with_traceback(sys.exc_info()[2])

    return data_file


def create_s104(output_file, dcf) -> S104File:
    """ Creates or updates an S104File object.
    Default values are set for any data that doesn't have options or are mandatory to be filled in the S104 spec.

    Parameters
    ----------
    output_file
        S104File object
    dcf
       S100 Data Coding Format (Int)

    Returns
    -------
    data_file
        The S104File object created or updated by this function.


    """
    data_file = _get_S104File(output_file)
    root = data_file.root
    root.water_level = data_file.make_container_for_dcf(dcf)
    root.water_level.water_level_create()

    root.feature_information_create()
    group_f = root.feature_information
    group_f.feature_code_create()
    group_f.water_level_feature_dataset_create()

    water_level_feature_dataset = root.feature_information.water_level_feature_dataset

    water_level_height_info = water_level_feature_dataset.append_new_item()
    water_level_height_info.code = "waterLevelHeight"
    water_level_height_info.name = "Water Level Height"
    water_level_height_info.unit_of_measure = "metre"
    water_level_height_info.datatype = "H5T_FLOAT"
    water_level_height_info.fill_value = f"{FILLVALUE_HEIGHT:0.02f}"
    water_level_height_info.lower = "-99.99"
    water_level_height_info.upper = "99.99"
    water_level_height_info.closure = "closedInterval"

    water_level_trend_info = water_level_feature_dataset.append_new_item()
    water_level_trend_info.code = "waterLevelTrend"
    water_level_trend_info.name = "Water Level Trend"
    water_level_trend_info.unit_of_measure = ""
    water_level_trend_info.datatype = "H5T_ENUM"
    water_level_trend_info.fill_value = FILLVALUE_TREND
    water_level_trend_info.lower = ""
    water_level_trend_info.upper = ""
    water_level_trend_info.closure = ""

    return data_file


def add_metadata(metadata: dict, data_file) -> S104File:
    """  Updates an S104File object based on input metadata.

    Parameters
    ----------
    data_file
        S104File object
    metadata
        a dictionary of metadata describing the data passed in,
        metadata should have the following key/value pairs:
            - "productSpecification": The product specification used to create
            this dataset.
            - "horizontalCRS": Horizontal Datum EPSG code.
            - "metadata": File name for the associated discovery metadata (xml)
            - "geographicIdentifier": Location of the data, ex: "Tampa Bay".
                An empty string ("") is the default.
            - "waterLevelHeightUncertainty": In (meters) arises from the
            hydrodynamic model, and the spatial interpolation method.
            The default, denoting a missing value, is -1.0.
            - "verticalUncertainty": Accuracy of vertical datum
                The default, denoting a missing value, is -1.0.
            - "horizontalPositionUncertainty": Accuracy of geolocation
            techniques, model grid accuracy. The default, denoting a missing
            value, is -1.0.
            - "timeUncertainty": Sensor accuracy, data time tagging accuracy
                The default, denoting a missing value, is -1.0.
            - "waterLevelTrendThreshold": Critical value used to determine
            steady water level trend. Units are meters/hour (m/hr).
            - "verticalCS": Vertical datum EPSG Code.
            - "verticalDatumReference": For verticalCoordinateBase(2) only
                - 'S-100 vertical datum': 1
                - 'EPSG': 2
            - "verticalDatum":
                - 'meanLowWaterSprings': 1
                - 'meanLowerLowWaterSprings': 2
                - 'meanSeaLevel': 3
                - 'lowestLowWater': 4
                - 'meanLowWater': 5
                - 'lowestLowWaterSprings': 6
                - 'approximateMeanLowWaterSprings': 7
                - 'indianSpringLowWater': 8
                - 'lowWaterSprings': 9
                - 'approximateLowestAstronomicalTide': 10
                - 'nearlyLowestLowWater': 11
                - 'meanLowerLowWater': 12
                - 'lowWater': 13
                - 'approximateMeanLowWater': 14
                - 'approximateMeanLowerLowWater': 15
                - 'meanHighWater': 16
                - 'meanHighWaterSprings': 17
                - 'highWater': 18
                - 'approximateMeanSeaLevel': 19
                - 'highWaterSprings': 20
                - 'meanHigherHighWater': 21
                - 'equinoctialSpringLowWater': 22
                - 'lowestAstronomicalTide': 23
                - 'localDatum': 24
                - 'internationalGreatLakesDatum1985': 25
                - 'meanWaterLevel': 26
                - 'lowerLowWaterLargeTide': 27
                - 'higherHighWaterLargeTide': 28
                - 'nearlyHighestHighWater': 29
                - 'highestAstronomicalTide': 30
                - 'balticSeaChartDatum2000': 44
                - 'internationalGreatLakesDatum2020: 46
                - 'seaFloor': 47
                - 'seaSurface': 48
                - 'hydrographicZero: 49
            - "verticalDatumReference":
                - 'S-100 Vertical datum': 1
                - 'EPSG code': 2
           - "commonPointRule":
                - 'average': 1
                - 'low': 2
                - 'high': 3
                - 'all': 4
           - "interpolationType": Interpolation method recommended for
           evaluation of the S100_GridCoverage.
                    - 'nearestneighbor': 1
                    - 'linear': 2
                    - 'quadratic': 3
                    - 'cubic': 4
                    - 'bilinear': 5
                    - 'biquadratic': 6
                    - 'bicubic': 7
                    - 'lostarea': 8
                    - 'barycentric': 9
                    - 'discrete': 10
            - "typeOfWaterLevelData:
                    - 'observation': 1
                    - 'astronomicalPrediction': 2
                    - 'analysisOrHybrid': 3
                    - 'hydrodynamicHindcast': 4
                    - 'hydrodynamicForecast': 5
                    - 'observedMinusPredicted': 6
                    - 'observedMinusAnalysis': 7
                    - 'observedMinusHindcast': 8
                    - 'observedMinusForecast': 9
            - "methodWaterLevelProduct": Brief description of tide gauge type,
            forecast method or model, etc.
            - "dateTimeOfFirstRecord": Valid time of earliest value, 'YYYYMMDDTHHMMSSZ'
            - "datasetDeliveryInterval": The expected time interval between availability of successive
                datasets for time-varying data. Must be formatted as 'PnYnMnDTnHnMnS' (ISO 8601 duration)
            - "trendInterval": The interval over which trend at a particular time is calculated.
                Unit: minutes.

    Returns
    -------
    data_file
        An S104File object updated by this function.

    """
    root = data_file.root
    water_level_feature = root.water_level

    water_level_feature_instance_01 = water_level_feature.water_level.append_new_item()

    water_level_feature_instance_01.water_level_group_create()

    water_level_feature_instance_01.uncertainty_dataset_create()
    water_level_height_uncertainty = water_level_feature_instance_01.uncertainty_dataset.append_new_item()
    water_level_height_uncertainty.name = "waterLevelHeight"
    water_level_height_uncertainty.value = metadata["waterLevelHeightUncertainty"]

    water_level_feature_instance_01.time_record_interval = 0

    utc_now = datetime.datetime.now(datetime.timezone.utc)

    if "issueDate" in metadata:
        root.issue_date = metadata["issueDate"]
    else:
        root.issue_date = utc_now.date()
    if "issueTime" in metadata:
        root.issue_time = metadata["issueTime"]
    else:
        root.issue_time = utc_now.time()

    root.product_specification = S104File.PRODUCT_SPECIFICATION
    root.metadata = ""
    root.dataset_delivery_interval = metadata["datasetDeliveryInterval"]
    root.trend_interval = metadata["trendInterval"]
    root.horizontal_crs = metadata["horizontalCRS"]
    root.geographic_identifier = metadata["geographicIdentifier"]
    root.water_level_trend_threshold = metadata["waterLevelTrendThreshold"]
    root.vertical_cs = metadata["verticalCS"]
    root.vertical_datum_reference = metadata["verticalDatumReference"]
    root.vertical_datum = metadata["verticalDatum"]
    water_level_feature.common_point_rule = metadata["commonPointRule"]
    water_level_feature.interpolation_type = metadata["interpolationType"]
    water_level_feature.vertical_uncertainty = metadata["verticalUncertainty"]
    water_level_feature.horizontal_position_uncertainty = metadata["horizontalPositionUncertainty"]
    water_level_feature.method_water_level_product = metadata["methodWaterLevelProduct"]
    water_level_feature_instance_01.data_dynamicity = metadata["dataDynamicity"]
    water_level_feature_instance_01.date_time_of_first_record = metadata["dateTimeOfFirstRecord"]

    return data_file


def add_data_from_arrays(height: s1xx_sequence, trend, data_file, grid_properties: dict, datetime_value, data_coding_format) -> S104File:
    """  Updates an S104File object based on numpy array/h5py datasets.
        Calls :any:`create_s104` then fills in the HDF5 datasets with the
        supplied water level height and trend numpy.arrays.

        Raises an S104Exception if the shapes of the water level height and
        trend (if not None) grids are not equal.

        Parameters
        ----------
        height
            1d or 2d array containing water level heights
        trend
            1d or 2d array containing water level trends
        data_file
            S104File object
        datetime_value
            datetime object
        grid_properties
            a dictionary of metadata describing the grids passed in,
            metadata can have the following key/value pairs:
                - "minx": West bound longitude
                - "maxx": East bound longitude
                - "miny": South bound latitude
                - "maxy": North bound latitude
                - "cellsize_x": Only for DCF2, grid spacing longitude
                - "cellsize_y": Only for DCF2, grid spacing latitude
                - "nx": Only for DCF2, number of points longitudinal
                - "ny": Only for DCF2, number of points latitudinal
                - "latitude": Only for DCF3, latitude of nodes
                - "longitude": Only for DCF3, longitudes of nodes
                - "nodes": Only for DCF3, number of nodes
        data_coding_format
            - 'Time Series at fixed stations': 1
            - 'Regularly-Gridded arrays': 2
            - 'Ungeorectified Grid': 3
            - 'TIN': 7
            - 'Time Series at fixed stations (stationwise)': 8

        Returns
        -------
        data_file
            An S104File object updated by this function.

        """
    root = data_file.root
    water_level_feature = root.water_level
    water_level_feature_instance_01 = root.water_level.water_level[0]

    if data_coding_format == 2:
        water_level_feature.data_coding_format = data_coding_format
        water_level_feature_instance_01.start_sequence = "0,0"
        water_level_feature.sequencing_rule_scan_direction = "longitude,latitude"
        water_level_feature.sequencing_rule_type = 1
        water_level_feature_instance_01.grid_origin_longitude = grid_properties['minx']
        water_level_feature_instance_01.grid_origin_latitude = grid_properties['miny']
        water_level_feature_instance_01.grid_spacing_longitudinal = grid_properties['cellsize_x']
        water_level_feature_instance_01.grid_spacing_latitudinal = grid_properties['cellsize_y']

        water_level_feature_instance_01.num_points_latitudinal = grid_properties['ny']
        water_level_feature_instance_01.num_points_longitudinal = grid_properties['nx']

    elif data_coding_format == 3:
        water_level_feature.data_coding_format = data_coding_format
        water_level_feature_instance_01.number_of_nodes = grid_properties['nodes']

        water_level_feature_instance_01.positioning_create()
        positioning = water_level_feature_instance_01.positioning
        positioning.geometry_values_create()
        positioning.geometry_values.longitude = grid_properties['longitude']
        positioning.geometry_values.latitude = grid_properties['latitude']

    elif data_coding_format == 7:
        water_level_feature.data_coding_format = data_coding_format
        water_level_feature_instance_01.number_of_nodes = grid_properties['nodes']
        water_level_feature_instance_01.number_of_triangles = grid_properties['num_triangles']

        water_level_feature_instance_01.positioning_create()
        positioning = water_level_feature_instance_01.positioning
        positioning.geometry_values_create()
        positioning.geometry_values.longitude = grid_properties['longitude']
        positioning.geometry_values.latitude = grid_properties['latitude']
        positioning.adjacency = grid_properties['adjacency']
        positioning.triangles = grid_properties['triangles']

    water_level_feature_instance_01.east_bound_longitude = grid_properties['maxx']
    water_level_feature_instance_01.west_bound_longitude = grid_properties['minx']
    water_level_feature_instance_01.south_bound_latitude = grid_properties['miny']
    water_level_feature_instance_01.north_bound_latitude = grid_properties['maxy']

    water_level_feature.axis_names = numpy.array(["longitude", "latitude"])
    root.water_level.dimension = len(water_level_feature.axis_names)

    min_height = numpy.min(height[numpy.where(height != FILLVALUE_HEIGHT)])
    max_height = numpy.max(height[numpy.where(height != FILLVALUE_HEIGHT)])
    water_level_feature.min_dataset_height = numpy.round(min_height, decimals=2)
    water_level_feature.max_dataset_height =  numpy.round(max_height, decimals=2)

    if min_height < water_level_feature.min_dataset_height and min_height != FILLVALUE_HEIGHT:
        water_level_feature.min_dataset_height = min_height

    if max_height > water_level_feature.max_dataset_height and max_height != FILLVALUE_HEIGHT:
        water_level_feature.max_dataset_height = max_height

    if numpy.ma.is_masked(height):
        height = height.filled(FILLVALUE_HEIGHT)

    height = numpy.round(height, decimals=2)
    trend.astype(int)

    if height.shape != trend.shape:
        raise S104Exception("Water level height & trend grids have different shapes")

    water_level_group_object = water_level_feature_instance_01.water_level_group.append_new_item()
    water_level_group_object.time_point = datetime_value

    water_level_group_object.values_create()
    grid = water_level_group_object.values
    grid.water_level_height = height
    grid.water_level_trend = trend

    return data_file


def update_metadata(data_file, grid_properties: dict, update_meta: dict) -> S104File:
    """  Updates an S104File object based on dynamic metadata.

          Parameters
          ----------
          data_file
              S104File object
          grid_properties
              a dictionary of metadata describing the dynamic data passed in,
              metadata can have the following key/value pairs:
                  - "minx": West bound longitude
                  - "maxx": East bound longitude
                  - "miny": South bound latitude
                  - "maxy": North bound latitude
                  - "cellsize_x": Only for DCF2, grid spacing longitude
                  - "cellsize_y": Only for DCF2, grid spacing latitude
                  - "nx": Only for DCF2, number of points longitudinal
                  - "ny": Only for DCF2, number of points latitudinal
                  - "latitude": Only for DCF3, latitude of nodes
                  - "longitude": Only for DCF3, longitudes of nodes
                  - "nodes": Only for DCF3, number of nodes
          update_meta
              a dictionary of dynamic metadata, metadata can have the following
              key/value pairs:
                  - "dateTimeOfLastRecord": Valid ISO 8601 time of latest value
                  - "numberOfGroups": Number of forecasts
                  - "numberOfTimes": Number of valid times
                  - "timeRecordInterval": Time between forecasts in seconds
                  - "num_instances": Number of water level feature instances

          Returns
          -------
          data_file
              An S104File object updated by this function.

    """
    root = data_file.root
    water_level_feature = root.water_level
    water_level_feature.num_instances = update_meta["num_instances"]
    water_level_feature_instance_01 = root.water_level.water_level[0]

    water_level_feature_instance_01.date_time_of_last_record = update_meta['dateTimeOfLastRecord']
    water_level_feature_instance_01.num_grp = update_meta['numberOfGroups']
    water_level_feature_instance_01.number_of_times = update_meta['numberOfTimes']
    water_level_feature_instance_01.time_record_interval = update_meta['timeRecordInterval']

    root.east_bound_longitude = grid_properties["maxx"]
    root.west_bound_longitude = grid_properties["minx"]
    root.south_bound_latitude = grid_properties["miny"]
    root.north_bound_latitude = grid_properties["maxy"]

    return data_file


def write_data_file(data_file):
    """  Writes file structure, metadata, data and closes S104File object."""

    data_file.write()
    data_file.flush()
    data_file.close()

