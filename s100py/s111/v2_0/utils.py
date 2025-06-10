""" Functions to create S111 data from other sources

"""
import logging
import sys
import datetime
import numpy
import warnings
from typing import Union, Optional

from ...s1xx import s1xx_sequence
from .api import S111File, FILLVALUE, S111Exception, VERTICAL_DATUM, VERTICAL_DATUM_REFERENCE


def _get_S111File(output_file):
    """ Small helper function to convert the output_file parameter into a S111File"""
    if isinstance(output_file, S111File):
        data_file = output_file
    else:
        try:
            data_file = S111File(output_file, "w")
        except TypeError as typeerr:
            msg = "Failed to create S111File using {}".format(str(output_file))
            logging.error(msg)
            raise type(typeerr)(msg).with_traceback(sys.exc_info()[2])

    return data_file

def create_s111(output_file, dcf, speed_uncertainty=False, direction_uncertainty=False) -> S111File:
    """ Creates or updates an S111File object.
    Default values are set for any data that doesn't have options
    or are mandatory to be filled in the S111 spec.

    Parameters
    ----------
    output_file
        S111File object
    dcf
       S100 Data Coding Format (Int)
    speed_uncertainty
       (Bool, optional, default is False) Feature attribute characterising
       the accuracy of a speed value, or of the magnitude component of
       a velocity. The estimate is as defined within a particular confidence
       level and expressed as a positive value, if True the feature attribute
       and default values are added to Group_F.
    direction_uncertainty
       (Bool, optional, default is False) Feature attribute characterising
       the best estimate of the accuracy of a bearing, if True the feature
       attribute and default values are added to Group_F.

    Returns
    -------
    data_file
        The S111File object created or updated by this function.


    """
    data_file = _get_S111File(output_file)

    root = data_file.root
    root.surface_current = data_file.make_container_for_dcf(dcf)
    root.surface_current.surface_current_create()

    root.feature_information_create()
    group_f = root.feature_information
    group_f.feature_code_create()
    group_f.surface_current_feature_dataset_create()

    surface_current_feature_dataset = root.feature_information.surface_current_feature_dataset
    data_file.set_feature_information_defaults(surface_current_feature_dataset)

    if speed_uncertainty:
        data_file.set_speed_uncertainty_defaults(surface_current_feature_dataset)
    if direction_uncertainty:
        data_file.set_direction_uncertainty_defaults(surface_current_feature_dataset)

    return data_file


def add_metadata(metadata: dict, data_file) -> S111File:
    """  Updates an S111File object based on input metadata.

    Parameters
    ----------
    data_file
        S111File object
    metadata
        a dictionary of metadata describing the grids passed in,
        metadata should have the following key/value pairs:
            - "productSpecification": The product specification used to create this dataset.
            - "horizontalCRS":  Horizontal EPSG code or -1.
            - "metadata": File name for the associated discovery metadata (xml)
            - "geographicIdentifier": Location of the data, ex: "Chesapeake Bay".
                    An empty string ("") is the default.
                - "speedUncertainty": In (knots) arises from the current meter data on which the model is verified,
                    the hydrodynamic model, and the spatial interpolation method.
                    The default, denoting a missing value, is -1.0.
                - "directionUncertainty": In (degrees) arises from the current meter data on which the model is verified,
                    the hydrodynamic model, and the spatial interpolation method.
                    The default, denoting a missing value, is -1.0.
                - "verticalUncertainty": Accuracy of vertical datum
                    The default, denoting a missing value, is -1.0.
                - "horizontalPositionUncertainty": Accuracy of geolocation techniques, model grid accuracy
                    The default, denoting a missing value, is -1.0.
                - "timeUncertainty": Sensor accuracy, data time tagging accuracy
                    The default, denoting a missing value, is -1.0.
                - "surfaceCurrentDepth": Depth or height (depthTypeIndex=1) or layer thickness (depthTypeIndex=2) (m)
                - "depthTypeIndex":
                    - 'heightOrDepth': 1
                    - 'layerAverage': 2
               - "commonPointRule":
                    - 'average': 1
                    - 'low': 2
                    - 'high': 3
                    - 'all': 4
               - "interpolationType":
                    Interpolation method recommended for evaluation of the S100_GridCoverage.
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
                    - 'hydrodynamicHindcast': 4
                    - 'hydrodynamicForecast': 5
                - "methodCurrentsProduct": Brief description of current meter type, forecast method or model, etc.
                - "datetimeOfFirstRecord": Valid time of the earliest value, 'YYYYMMDDTHHMMSSZ'
                - "datasetDeliveryInterval": The expected time interval between availability of successive datasets
                    for time-varying data. Must be formatted as 'PnYnMnDTnHnMnS' (ISO 8601 duration)

        Returns
        -------
        data_file
            An S111File object updated by this function.

        """
    root = data_file.root
    surface_current_feature = root.surface_current

    surface_current_feature_instance_01 = surface_current_feature.surface_current.append_new_item()

    surface_current_feature_instance_01.surface_current_group_create()

    surface_current_feature_instance_01.uncertainty_dataset_create()
    speed_uncertainty = surface_current_feature_instance_01.uncertainty_dataset.append_new_item()
    speed_uncertainty.name = "surfaceCurrentSpeed"
    speed_uncertainty.value = metadata["speedUncertainty"]
    direction_uncertainty = surface_current_feature_instance_01.uncertainty_dataset.append_new_item()
    direction_uncertainty.name = "surfaceCurrentDirection"
    direction_uncertainty.value = metadata["directionUncertainty"]

    surface_current_feature.min_dataset_current_speed = 0
    surface_current_feature.max_dataset_current_speed = 0
    surface_current_feature_instance_01.time_record_interval = 0

    utc_now = datetime.datetime.now(datetime.timezone.utc)
    try:
        root.product_specification = S111File.PRODUCT_SPECIFICATION
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
        if "datasetDeliveryInterval" in metadata:
            root.dataset_delivery_interval = metadata["datasetDeliveryInterval"]

        root.depth_type_index = metadata["depthTypeIndex"]
        root.surface_current_depth = metadata["surfaceCurrentDepth"]

        # Additional restrictions on core general metadata for S-111
        root.issue_time = (metadata["issueTime"] if "issueTime" in metadata else utc_now.time())
        if metadata["verticalCoordinateBase"] != 2:
            warnings.warn("Warning the only allowed value for verticalCoordinateBase is verticalDatum (2)")
            metadata["verticalCoordinateBase"] = 2
        else:
            root.vertical_coordinate_base = metadata["verticalCoordinateBase"]

        if "depthTypeIndex" in metadata:
            if metadata["depthTypeIndex"] == 1:
                root.vertical_cs = metadata["verticalCS"]
                root.vertical_datum_reference = metadata.get('verticalDatumReference', VERTICAL_DATUM_REFERENCE.s100VerticalDatum)
                if metadata["verticalDatumReference"] == 1:
                    root.vertical_datum = metadata.get("verticalDatum", VERTICAL_DATUM.seaSurface)

        # SurfaceCurrent Feature Type Metadata
        surface_current_feature.common_point_rule = metadata["commonPointRule"]
        surface_current_feature.horizontal_position_uncertainty = metadata["horizontalPositionUncertainty"]
        surface_current_feature.vertical_uncertainty = metadata["verticalUncertainty"]

        # Optional SurfaceCurrent Feature type metadata
        if "timeUncertainty" in metadata:
            surface_current_feature.time_uncertainty = metadata["timeUncertainty"]
        if "dataOffsetCode" in metadata:
            surface_current_feature.data_offset_code = metadata["dataOffsetCode"]

        # SurfaceCurrent Feature Type Additional Metadata
        if "methodCurrentsProduct" in metadata:
            surface_current_feature.method_currents_product = metadata["methodCurrentsProduct"]

        # Optional SurfaceCurrent.NN Feature Instance Metadata
        if "datetimeOfFirstRecord" in metadata:
            surface_current_feature_instance_01.date_time_of_first_record = metadata["datetimeOfFirstRecord"]

        # SurfaceCurrent.NN Feature Instance Additional Metadata
        surface_current_feature_instance_01.data_dynamicity = metadata["dataDynamicity"]

    except KeyError as e:
        raise S111Exception(f"Error: Mandatory S-111 attribute {e} not found in the metadata dictionary")

    return data_file


def add_data_from_arrays(speed: s1xx_sequence, direction: s1xx_sequence, data_file,
                         grid_properties: dict, datetime_value, data_coding_format,
                         speed_uncertainty: Optional[s1xx_sequence]=None, direction_uncertainty: Optional[s1xx_sequence]=None ) -> S111File:
    """  Updates an S111File object based on numpy array/h5py datasets.
        Calls :any:`create_s111` then fills in the HDF5 datasets with the supplied speed and direction numpy.arrays.

        Raises an S11Exception if the shapes of the speed and direction (if not None) grids are not equal.

        Parameters
        ----------
        speed
            1d or 2d array containing surface current speeds.
        direction
            1d or 2d array containing surface current directions.
        data_file
            S111File object
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
            - 'Time series at fixed stations': 1,
            - 'Regularly-gridded arrays': 2,
            - 'Ungeorectified gridded arrays': 3,
            - 'Time series data for one moving platform': 4
            - 'Time Series at fixed stations (stationwise)': 8
        speed_uncertainty (Optional)
            1d or 2d array containing surface current speed uncertainty.
        direction_uncertainty (Optional)
            1d or 2d array containing surface current direction uncertainty.

        Returns
        -------
        data_file
            An S111File object updated by this function.

        """

    root = data_file.root
    surface_current_feature = root.surface_current
    surface_current_feature_instance_01 = root.surface_current.surface_current[0]

    if speed.shape != direction.shape:
        raise S111Exception("Speed and Direction grids have different shapes")

    if data_coding_format == 2:
        surface_current_feature.data_coding_format = data_coding_format
        surface_current_feature_instance_01.start_sequence = "0,0"
        surface_current_feature.sequencing_rule_scan_direction = "longitude,latitude"
        surface_current_feature.sequencing_rule_type = 1
        surface_current_feature_instance_01.grid_origin_longitude = grid_properties['minx']
        surface_current_feature_instance_01.grid_origin_latitude = grid_properties['miny']
        surface_current_feature_instance_01.grid_spacing_longitudinal = grid_properties['cellsize_x']
        surface_current_feature_instance_01.grid_spacing_latitudinal = grid_properties['cellsize_y']

        surface_current_feature_instance_01.num_points_latitudinal = grid_properties['ny']
        surface_current_feature_instance_01.num_points_longitudinal = grid_properties['nx']

    elif data_coding_format == 3:
        surface_current_feature.data_coding_format = data_coding_format
        surface_current_feature_instance_01.number_of_nodes = grid_properties['nodes']

        surface_current_feature_instance_01.positioning_create()
        positioning = surface_current_feature_instance_01.positioning
        positioning.geometry_values_create()
        geometry_values = positioning.geometry_values
        geometry_values.longitude = grid_properties['longitude']
        geometry_values.latitude = grid_properties['latitude']

    if data_coding_format == 2 or data_coding_format == 3:
        surface_current_feature.interpolation_type = 10

    surface_current_feature_instance_01.east_bound_longitude = grid_properties['maxx']
    surface_current_feature_instance_01.west_bound_longitude = grid_properties['minx']
    surface_current_feature_instance_01.south_bound_latitude = grid_properties['miny']
    surface_current_feature_instance_01.north_bound_latitude = grid_properties['maxy']

    surface_current_feature.axis_names = numpy.array(["longitude", "latitude"])
    root.surface_current.dimension = len(surface_current_feature.axis_names)

    min_speed = numpy.min(speed[numpy.where(speed != FILLVALUE)])
    max_speed = numpy.max(speed[numpy.where(speed != FILLVALUE)])
    surface_current_feature.min_dataset_current_speed = numpy.round(min_speed, decimals=2)
    surface_current_feature.max_dataset_current_speed =  numpy.round(max_speed, decimals=2)

    if min_speed < surface_current_feature.min_dataset_current_speed and min_speed != FILLVALUE:
        surface_current_feature.min_dataset_current_speed = min_speed

    if max_speed > surface_current_feature.max_dataset_current_speed and max_speed != FILLVALUE:
        surface_current_feature.max_dataset_current_speed = max_speed

    if numpy.ma.is_masked(speed):
        speed = speed.filled(FILLVALUE)
        direction = direction.filled(FILLVALUE)

    speed = numpy.round(speed, decimals=2)
    direction = numpy.round(direction, decimals=1)

    surface_current_group_object = surface_current_feature_instance_01.surface_current_group.append_new_item()
    surface_current_group_object.time_point = datetime_value

    surface_current_group_object.values_create()
    grid = surface_current_group_object.values
    grid.surface_current_speed = speed
    grid.surface_current_direction = direction

    feature_info = root.feature_information.surface_current_feature_dataset

    if isinstance(speed_uncertainty, numpy.ndarray):
        speed_uncertainty_feature_info = False
        for i in range(len(feature_info)):
            if feature_info[i].name == 'Speed Uncertainty':
                speed_uncertainty_feature_info = True
                grid.speed_uncertainty = speed_uncertainty
                break
        if not speed_uncertainty_feature_info:
            raise S111Exception("AttributeError: Speed uncertainty is not present in Group_F, values grid must"
                                " conform to the feature information group, see create_s111()")

    if isinstance(direction_uncertainty, numpy.ndarray):
        direction_uncertainty_feature_info = False
        for i in range(len(feature_info)):
            if feature_info[i].name == 'Direction Uncertainty':
                direction_uncertainty_feature_info = True
                grid.direction_uncertainty = direction_uncertainty
                break
        if not direction_uncertainty_feature_info:
            raise S111Exception("AttributeError: Direction uncertainty is not present in Group_F, values grid must"
                                " conform to the feature information group, see create_s111()")

    return data_file


def update_metadata(data_file, grid_properties: dict, update_meta: dict) -> S111File:
    """  Updates an S111File object based on dynamic metadata.

          Parameters
          ----------
          data_file
              S111File object
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
          update_meta
              a dictionary of dynamic metadata, metadata can have the following
              key/value pairs:
                  - "dateTimeOfLastRecord": Valid ISO 8601 time of latest value
                  - "numberOfGroups": Number of forecasts
                  - "numberOfTimes": Number of valid times
                  - "timeRecordInterval": Time between forecasts in seconds
                  - "num_instances": Number of surface current feature instances

          Returns
          -------
          data_file
              An S111File object updated by this function.

          """
    root = data_file.root
    surface_current_feature = root.surface_current
    surface_current_feature.num_instances = update_meta["num_instances"]
    surface_current_feature_instance_01 = root.surface_current.surface_current[0]

    surface_current_feature_instance_01.date_time_of_last_record = update_meta['dateTimeOfLastRecord']
    surface_current_feature_instance_01.num_grp = update_meta['numberOfGroups']
    surface_current_feature_instance_01.number_of_times = update_meta['numberOfTimes']
    surface_current_feature_instance_01.time_record_interval = update_meta['timeRecordInterval']

    root.east_bound_longitude = grid_properties["maxx"]
    root.west_bound_longitude = grid_properties["minx"]
    root.south_bound_latitude = grid_properties["miny"]
    root.north_bound_latitude = grid_properties["maxy"]

    return data_file


def write_data_file(data_file):
    """  Writes file structure, metadata, data and closes S111File object."""

    data_file.write()
    data_file.flush()
    data_file.close()

