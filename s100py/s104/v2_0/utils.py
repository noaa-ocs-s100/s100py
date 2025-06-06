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
        (Optional, default is False) Features attribute uncertainty,
        which represents the uncertainty at a particular grid point,
        may be omitted if the uncertainty is unknown or the same
        value at all grid points ('Bool')

    Returns
    -------
    data_file
        The S104File object created or updated by this function.


    """
    data_file = _get_S104File(output_file)
    root = data_file.root
    # Set Water Level default structure and values as directed in
    # S-104 Ed 2.0 production specification mandatory requirements
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

    if uncertainty:
        water_level_uncertainty_info = water_level_feature_dataset.append_new_item()
        water_level_uncertainty_info.code = "uncertainty"
        water_level_uncertainty_info.name = "Uncertainty"
        water_level_uncertainty_info.unit_of_measure = "metre"
        water_level_uncertainty_info.datatype = "H5T_FLOAT"
        water_level_uncertainty_info.fill_value = f"{FILLVALUE_UNCERTAINTY:0.02f}"
        water_level_uncertainty_info.lower = "0.00"
        water_level_uncertainty_info.upper = "99.99"
        water_level_uncertainty_info.closure = "closedInterval"

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
            - "datetimeOfFirstRecord": Valid time of earliest value, 'YYYYMMDDTHHMMSSZ'
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

    water_level_feature_instance_01.time_record_interval = 0

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

        # Feature Instance Metadata
        water_level_feature_instance_01.date_time_of_first_record = metadata["datetimeOfFirstRecord"]
        # Additional feature instance metadata
        water_level_feature_instance_01.data_dynamicity = metadata["dataDynamicity"]

    except KeyError as e:
        raise S104Exception(f"Error: Mandatory S-104 attribute {e} not found in the metadata dictionary")

    return data_file


def add_data_from_arrays(height: s1xx_sequence, trend, data_file, grid_properties: dict, datetime_value,
                         data_coding_format, uncertainty=None) -> S104File:
    """  Updates an S104File object based on numpy array/h5py datasets.
        Calls :any:`create_s104` then fills in the HDF5 datasets with the
        supplied water level height and trend numpy.arrays.

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
        domain_extent_polygon
            (Optional) Spatial extent of the domain coverage only to be used
            if different vertical coverages exits, dataset contains coordinates
            of bounding polygon vertices as a closed ring, the first and last
            coordinates will contain the same values)

        Returns
        -------
        data_file
            An S104File object updated by this function.

        """
    root = data_file.root
    water_level_feature = root.water_level
    water_level_feature_instance_01 = root.water_level.water_level[0]

    if isinstance(uncertainty, float):
        water_level_feature_instance_01.uncertainty_dataset_create()
        water_level_height_uncertainty = water_level_feature_instance_01.uncertainty_dataset.append_new_item()
        water_level_height_uncertainty.name = "waterLevelHeight"
        water_level_height_uncertainty.value = uncertainty

    if uncertainty is None:
        water_level_feature_instance_01.uncertainty_dataset_create()
        water_level_height_uncertainty = water_level_feature_instance_01.uncertainty_dataset.append_new_item()
        water_level_height_uncertainty.name = "waterLevelHeight"
        water_level_height_uncertainty.value = FILLVALUE_UNCERTAINTY

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
        height = height.filled(f"{FILLVALUE_HEIGHT:0.02f}")
        trend = trend.filled(f"{FILLVALUE_TREND}")

    height = numpy.round(height, decimals=2)
    trend.astype(int)

    if isinstance(uncertainty, numpy.ndarray):
        uncertainty = numpy.round(uncertainty, decimals=2)

    if height.shape != trend.shape:
        raise S104Exception("Water level height & trend grids have different shapes")

    water_level_group_object = water_level_feature_instance_01.water_level_group.append_new_item()
    water_level_group_object.time_point = datetime_value
    # Optional values group metadata
    # water_level_group_object.water_level_trend_threshold = root.water_level_trend_threshold
    # water_level_group_object.trend_interval = root.trend_interval

    water_level_group_object.values_create()
    grid = water_level_group_object.values
    grid.water_level_height = height
    grid.water_level_trend = trend

    if isinstance(uncertainty, numpy.ndarray):
        grid.water_level_uncertainty = uncertainty

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
