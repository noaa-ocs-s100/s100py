""" Functions to create S111 data from other sources

"""
import logging
import sys
import datetime

import numpy


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


def create_s111(output_file, dcf) -> S111File:
    """ Creates or updates an S111File object.
    Default values are set for any data that doesn't have options or are mandatory to be filled in the S111 spec.

    Parameters
    ----------
    output_file
        S111File object
    dcf
       S100 Data Coding Format (Int)

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

    surface_current_speed_info = surface_current_feature_dataset.append_new_item()
    surface_current_speed_info.code = "surfaceCurrentSpeed"
    surface_current_speed_info.name = "Surface Current Speed"
    surface_current_speed_info.unit_of_measure = "knot"
    surface_current_speed_info.datatype = "H5T_FLOAT"
    surface_current_speed_info.fill_value = f"{FILLVALUE:0.02f}"
    surface_current_speed_info.lower = "0.00"
    surface_current_speed_info.upper = ""
    surface_current_speed_info.closure = "geSemiInterval"

    surface_current_direction_info = surface_current_feature_dataset.append_new_item()
    surface_current_direction_info.code = "surfaceCurrentDirection"
    surface_current_direction_info.name = "Surface Current Direction"
    surface_current_direction_info.unit_of_measure = "degree"
    surface_current_direction_info.datatype = "H5T_FLOAT"
    surface_current_direction_info.fill_value = f"{FILLVALUE:0.01f}"
    surface_current_direction_info.lower = "0.0"
    surface_current_direction_info.upper = "359.9"
    surface_current_direction_info.closure = "closedInterval"

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
            - "timeUncertainty": Sensor accuracy, ata time tagging accuracy
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
            - "dateTimeOfFirstRecord": Valid time of the earliest value, 'YYYYMMDDTHHMMSSZ'
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

    if metadata["issueDateTime"]:
        root.issue_date = metadata["issueDateTime"]
    else:
        root.issue_date = utc_now
    if metadata["issueDateTime"]:
        root.issue_time = metadata["issueDateTime"]
    else:
        root.issue_time = utc_now

    root.product_specification = S111File.PRODUCT_SPECIFICATION
    root.metadata = ""
    root.horizontal_crs = metadata["horizontalCRS"]
    root.dataset_delivery_interval = metadata["datasetDeliveryInterval"]
    root.geographic_identifier = metadata["geographicIdentifier"]
    root.surface_current_depth = metadata["surfaceCurrentDepth"]
    root.depth_type_index = metadata["depthTypeIndex"]
    surface_current_feature.common_point_rule = metadata["commonPointRule"]
    surface_current_feature.vertical_uncertainty = metadata["verticalUncertainty"]
    surface_current_feature.horizontal_position_uncertainty = metadata["horizontalPositionUncertainty"]
    surface_current_feature_instance_01.data_dynamicity = metadata["dataDynamicity"]
    surface_current_feature.method_currents_product = metadata["methodCurrentsProduct"]
    surface_current_feature_instance_01.date_time_of_first_record = metadata["dateTimeOfFirstRecord"]

    # Optional
    # surface_current_feature.time_uncertainty = metadata["timeUncertainty"]

    if "depthTypeIndex" in metadata:
        if metadata["depthTypeIndex"] == 1:
            root.vertical_cs = metadata["verticalCS"]
            root.vertical_datum_reference = metadata.get('verticalDatumReference', VERTICAL_DATUM_REFERENCE.s100VerticalDatum)
            if metadata["verticalDatumReference"] == 1:
                root.vertical_datum = metadata.get("verticalDatum", VERTICAL_DATUM.seaSurface)

    return data_file


def add_data_from_arrays(speed: s1xx_sequence, direction: s1xx_sequence, data_file, grid_properties: dict, datetime_value, data_coding_format) -> S111File:
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

    min_speed = numpy.round(numpy.nanmin(speed), decimals=2)
    max_speed = numpy.round(numpy.nanmax(speed), decimals=2)

    if min_speed < surface_current_feature.min_dataset_current_speed:
        surface_current_feature.min_dataset_current_speed = min_speed

    if max_speed > surface_current_feature.max_dataset_current_speed:
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
