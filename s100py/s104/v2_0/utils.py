""" Functions to create S104 files and populate with data from other sources

"""
import logging
import sys
import datetime

import numpy

from ...s1xx import s1xx_sequence
from .api import S104File, FILLVALUE_HEIGHT, FILLVALUE_TREND, FILLVALUE_UNCERTAINTY, S104Exception


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


def create_s104(output_file, dcf, uncertainty=False) -> S104File:
    """ Creates or updates an S104File object.
    Default values are set for any data that doesn't have options or are mandatory to be filled in the S104 spec.

    Parameters
    ----------
    output_file
        S104File object
    dcf
       S100 Data Coding Format (Int)
    uncertainty
        (Bool, optional, default is False) Feature attribute uncertainty,
        which represents the uncertainty at a particular grid point,
        may be omitted if the uncertainty is unknown or the same
        value at all grid points, if True the feature attribute
        with mandatory values are added to Group_F.

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
    data_file.set_feature_information_defaults(water_level_feature_dataset)

    if uncertainty:
        data_file.set_water_level_uncertainty_defaults(water_level_feature_dataset)

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
            - "verticalCoordinateBase": Only allowed value 2: verticalDatum
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
            - "dataDynamicity":
                    - 'observation': 1
                    - 'astronomicalPrediction': 2
                    - 'analysisOrHybrid': 3
                    - 'hydrodynamicForecast': 5
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

    utc_now = datetime.datetime.now(datetime.timezone.utc)

    try:
        # General metadata
        root.product_specification = S104File.PRODUCT_SPECIFICATION
        root.issue_date = (metadata["issueDate"] if "issueDate" in metadata else utc_now.date())
        root.horizontal_crs = int(metadata["horizontalCRS"])

        # Optional general metadata
        if "geographicIdentifier" in metadata:
            root.geographic_identifier = metadata["geographicIdentifier"]
        if "horizontalCRS" in metadata and metadata["horizontalCRS"] == -1:
            root.name_of_horizontal_crs = metadata["nameOfHorizontalCRS"]
            root.type_of_horizontal_crs = metadata["typeOfHorizontalCRS"]
            root.horizontal_cs = metadata["horizontalCS"]
            root.horizontal_datum = metadata["horizontalDatum"]
            if "horizontalDatum" in metadata and metadata["horizontalDatum"] == -1:
                root.name_of_horizontal_datum = metadata["nameOfHorizontalDatum"]
                root.prime_meridian = metadata["primeMedian"]
                root.spheroid = metadata["spheroid"]
                if "typeOfHorizontalCRS" in metadata and metadata["typeOfHorizontalCRS"] == 2:
                    root.projection_method = metadata["projectionMethod"]
                    if "projectionMethod" in metadata:
                        root.projection_parameter_1 = metadata["projectionParameter1"]
                        root.projection_parameter_2 = metadata["projectionParameter2"]
                        root.projection_parameter_3 = metadata["projectionParameter3"]
                        root.projection_parameter_4 = metadata["projectionParameter4"]
                        root.projection_parameter_5 = metadata["projectionParameter5"]
                        root.false_northing = metadata["falseNorthing"]
                        root.false_easting = metadata["falseEasting"]
        if "epoch" in metadata:
            root.epoch = metadata["epoch"]

        # Additional general metadata
        root.water_level_trend_threshold = metadata["waterLevelTrendThreshold"]
        # Optional additional general metadata
        if "datasetDeliveryInterval" in metadata:
            root.dataset_delivery_interval = metadata["datasetDeliveryInterval"]
        if "trendInterval" in metadata:
            root.trend_interval = metadata["trendInterval"]
        if "verticalDatumEpoch" in metadata:
            root.vertical_datum_epoch = metadata["verticalDatumEpoch"]

        # Additional restrictions on core general metadata for S-104
        root.issue_time = (metadata["issueTime"] if "issueTime" in metadata else utc_now.time())
        root.vertical_cs = metadata["verticalCS"]
        root.vertical_coordinate_base = metadata["verticalCoordinateBase"]
        root.vertical_datum_reference = metadata["verticalDatumReference"]
        root.vertical_datum = metadata["verticalDatum"]

        # WaterLevel Feature type metadata
        water_level_feature.common_point_rule = metadata["commonPointRule"]
        water_level_feature.vertical_uncertainty = metadata["verticalUncertainty"]
        water_level_feature.horizontal_position_uncertainty = metadata["horizontalPositionUncertainty"]
        # Optional feature type metadata
        if "timeUncertainty" in metadata:
            water_level_feature.time_uncertainty = metadata["timeUncertainty"]
        # Additional feature type metadata
        water_level_feature.method_water_level_product = metadata["methodWaterLevelProduct"]
        # Feature type metadata dataCodingFormat = 2 (regGrid) feature type metadata
        water_level_feature.interpolation_type = metadata["interpolationType"]
        # Optional, Allowed values 1: XMin, YMin (“Lower left”) or 5:Barycenter (centroid) of cell
        if "dataOffsetCode" in metadata:
            water_level_feature.data_offset_code = metadata["dataOffsetCode"]

    except KeyError as e:
        raise S104Exception(f"AttributeError: Mandatory S-104 attribute {e} not found in the metadata dictionary")

    return data_file


def add_water_level_instance(data_file):
    """ Adds the water_level object container to the S104File object,
        a feature instance is created each time the function is called
        (e.g.`WaterLevel.01`, `WaterLevel.02). A series of feature
        instances should be implemented in the same dataset for any
        data collection that varies by extent, location, time, or grid
        size.

        Parameters
        ----------
        data_file
            S104File object

        Returns
        -------
        data_file
            An S104File object updated by this function.
    """
    root = data_file.root
    water_level_feature = root.water_level

    water_level_feature_instance = water_level_feature.water_level.append_new_item()
    water_level_feature_instance.water_level_group_create()


    return data_file


def add_data_from_arrays(height: s1xx_sequence, trend, data_file, grid_properties: dict, datetime_value,
                         data_coding_format, uncertainty=None) -> S104File:
    """ Updates an S104File object with the supplied water level height and trend numpy.arrays.

        Raises an S104Exception if the shapes of the water level height and
        trend (if not None) grids are not equal.

        Parameters
        ----------
        height
            2d array containing water level heights
        trend
            2d array containing water level trends
        data_file
            S104File object
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
        datetime_value
            datetime object
        data_coding_format
            - 'Regularly-Gridded arrays': 2
        uncertainty
            (Optional, default None) The value can be a 2d array
            containing water level uncertainty, a positive fixed float
            single value, or None, a 2d array will be inserted into the
            water level values compound dataset, a single value will be
            inserted into the water level features uncertainty dataset
            and if None an uncertainty dataset will be created with the
            fill value

        Returns
        -------
        data_file
            An S104File object updated by this function.

        """
    root = data_file.root
    water_level_feature = root.water_level

    feature_instance_object_id = len(root.water_level.water_level) - 1
    water_level_feature_instance = None
    try:
        water_level_feature_instance = root.water_level.water_level[feature_instance_object_id]
    except IndexError as e:
        raise S104Exception(f"IndexError: Water Level Feature Instance {e}, feature instance does not exist, "
                            f"use 'add_water_level_instance()' to add a water level feature instance")

    if isinstance(uncertainty, float):
        water_level_feature_instance.uncertainty_dataset_create()
        water_level_height_uncertainty = water_level_feature_instance.uncertainty_dataset.append_new_item()
        water_level_height_uncertainty.name = "waterLevelHeight"
        water_level_height_uncertainty.value = uncertainty

    if uncertainty is None:
        water_level_feature_instance.uncertainty_dataset_create()
        water_level_height_uncertainty = water_level_feature_instance.uncertainty_dataset.append_new_item()
        water_level_height_uncertainty.name = "waterLevelHeight"
        water_level_height_uncertainty.value = FILLVALUE_UNCERTAINTY

    if data_coding_format == 2:
        water_level_feature.data_coding_format = data_coding_format
        water_level_feature_instance.start_sequence = "0,0"
        water_level_feature.sequencing_rule_scan_direction = "longitude,latitude"
        water_level_feature.sequencing_rule_type = 1
        water_level_feature_instance.grid_origin_longitude = grid_properties['minx']
        water_level_feature_instance.grid_origin_latitude = grid_properties['miny']
        water_level_feature_instance.grid_spacing_longitudinal = grid_properties['cellsize_x']
        water_level_feature_instance.grid_spacing_latitudinal = grid_properties['cellsize_y']

        water_level_feature_instance.num_points_latitudinal = grid_properties['ny']
        water_level_feature_instance.num_points_longitudinal = grid_properties['nx']

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

    height[numpy.isnan(height)] = FILLVALUE_HEIGHT
    trend[numpy.isnan(trend)] = FILLVALUE_TREND

    if not numpy.ma.is_masked(height):
        height = numpy.ma.masked_equal(height, FILLVALUE_HEIGHT)

    if not numpy.ma.is_masked(trend):
        trend =  numpy.ma.masked_array(trend, mask=height.mask)
        trend =  numpy.ma.filled(trend, FILLVALUE_TREND)

    height = numpy.round(height, decimals=2)
    trend.astype(int)

    if height.shape != trend.shape:
        raise S104Exception("Water level height & trend grids have different shapes")

    water_level_group_object = water_level_feature_instance.water_level_group.append_new_item()
    water_level_group_object.time_point = datetime_value

    water_level_group_object.values_create()
    grid = water_level_group_object.values
    grid.water_level_height = height
    grid.water_level_trend = trend

    feature_info = root.feature_information.water_level_feature_dataset

    if isinstance(uncertainty, numpy.ndarray):
        uncertainty_feature_info = False
        for i in range(len(feature_info)):
            if feature_info[i].name == 'Uncertainty':
                uncertainty_feature_info = True
                uncertainty = numpy.round(uncertainty, decimals=2)
                grid.water_level_uncertainty = uncertainty
                break
        if not uncertainty_feature_info:
            raise S104Exception("AttributeError: Water level uncertainty is not present in Group_F, values grid must"
                                " conform to the feature information group, see create_s104()")

    return data_file

def update_metadata(data_file, grid_properties: dict, metadata: dict) -> S104File:
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
          metadata
              a dictionary of metadata describing the data, metadata must have
              the following key/value pairs:
                  - "dateTimeOfFirstRecord": Valid time of earliest value, 'YYYYMMDDTHHMMSSZ'
                  - "dataDynamicity": Classification of data according to the relationship between the time of its
                  collection, generation, or calculation of generation parameters, in relation to the time of
                  publication of the dataset
                      - 'observation': 1
                      - 'astronomicalPrediction': 2
                      - 'analysisOrHybrid': 3
                      - 'hydrodynamicForecast': 5

          Returns
          -------
          data_file
              An S104File object updated by this function.

    """
    try:
        root = data_file.root
        water_level_feature = root.water_level
        num_feature_instances = len(root.water_level.water_level)
        feature_instance_object_id = num_feature_instances - 1
        water_level_feature_instance = root.water_level.water_level[feature_instance_object_id]

        num_groups = len(water_level_feature_instance.water_level_group)
        num_groups_id = num_groups - 1

        last_time_point = water_level_feature_instance.water_level_group[num_groups_id].time_point
        last_datetime_record = last_time_point.strftime("%Y%m%dT%H%M%SZ")
        time_record_interval = 0
        if num_groups >= 2:
            first_timestamp = water_level_feature_instance.water_level_group[0].time_point
            second_timestamp = water_level_feature_instance.water_level_group[1].time_point

            interval = second_timestamp - first_timestamp
            time_record_interval = interval.total_seconds()

        water_level_feature_instance.date_time_of_first_record = metadata["dateTimeOfFirstRecord"]
        # Additional feature instance metadata
        water_level_feature_instance.data_dynamicity = metadata["dataDynamicity"]

        water_level_feature.num_instances = num_feature_instances
        water_level_feature_instance.date_time_of_last_record = last_datetime_record
        water_level_feature_instance.num_grp = num_groups
        water_level_feature_instance.number_of_times = num_groups
        water_level_feature_instance.time_record_interval = time_record_interval

        root.east_bound_longitude = grid_properties["maxx"]
        root.west_bound_longitude = grid_properties["minx"]
        root.south_bound_latitude = grid_properties["miny"]
        root.north_bound_latitude = grid_properties["maxy"]

    except KeyError as e:
        raise S104Exception(f"KeyError: S-104 attribute {e} not found in the metadata dictionary")

    return data_file


def write_data_file(data_file):
    """  Writes file structure, metadata, data and closes S104File object."""

    data_file.write()
    data_file.flush()
    data_file.close()
