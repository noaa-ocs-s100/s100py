import os
import numpy
import datetime
import fire
import re
# import logging

from thyme.model import model, roms, fvcom, pom, hycom
from s100py.s111 import S111File

# """
# Default logging level.
# """
# LOGGING_LEVEL = logging.DEBUG
#
# # Use global/module-level logging
# logging.basicConfig(level=LOGGING_LEVEL)
# logger = logging.getLogger(__name__)

# Default fill value for NetCDF variables
FILLVALUE = -9999.0

# Default depth in meters
DEFAULT_TARGET_DEPTH = 4.5
PRODUCER_CODE = "US"
DATA_CODING_FORMAT = 2
TYPE_OF_CURRENT_DATA = 6

MODELTYPE_FVCOM = 'fvcom'
MODELTYPE_HYCOM = 'hycom'
MODELTYPE_POM = 'pom'
MODELTYPE_ROMS = 'roms'

MODEL_FILE_CLASS = {
    MODELTYPE_FVCOM: fvcom.FVCOMFile,
    MODELTYPE_HYCOM: hycom.HYCOMFile,
    MODELTYPE_POM: pom.POMFile,
    MODELTYPE_ROMS: roms.ROMSFile
}

MODEL_INDEX_CLASS = {
    MODELTYPE_FVCOM: fvcom.FVCOMIndexFile,
    MODELTYPE_HYCOM: hycom.HYCOMIndexFile,
    MODELTYPE_POM: pom.POMIndexFile,
    MODELTYPE_ROMS: roms.ROMSIndexFile
}


PRODUCT_DESCRIPTION_FVCOM = 'FVCOM_Hydrodynamic_Model_Forecasts'
PRODUCT_DESCRIPTION_HYCOM = 'HYCOM_Hydrodynamic_Model_Forecasts'
PRODUCT_DESCRIPTION_POM = 'POM_Hydrodynamic_Model_Forecasts'
PRODUCT_DESCRIPTION_ROMS = 'ROMS_Hydrodynamic_Model_Forecasts'

MODELS = {
    'cbofs': {
        'region': 'Chesapeake_Bay',
        'product': PRODUCT_DESCRIPTION_ROMS,
        'model_type': MODELTYPE_ROMS,
    },
    'nyofs': {

        'region': 'Port_of_New_York_and_New_Jersey',
        'product': PRODUCT_DESCRIPTION_POM,
        'model_type': MODELTYPE_POM,
    }
}


class CLI:
    """
    Container for methods exposed through CLI via Python Fire.
    """
    def __init__(self):
        pass

    def convert(self, model_file_path, model_index_path, output_path):
        """Convert NetCDF hydrodynamic model to S111 format.

        Args:
            model_file_path: Path to native ofs model file.
            model_index_path: Path to model index file.
            output_path:  Path to output hdf5 file.

        Returns:
            List of paths to HDF5 files created.
        """
        s111_path_prefix = os.path.split(output_path)[0]
        model_filename = os.path.split(model_file_path)[-1]
        model_name = model_filename.split('.')[1]

        model_index = MODEL_INDEX_CLASS[MODELS[model_name]['model_type']](model_index_path)
        model_file = MODEL_FILE_CLASS[MODELS[model_name]['model_type']](model_file_path)

        # Path format/prefix for output S111 files. Forecast initialization (reference).
        if os.path.isdir(s111_path_prefix):
            if not s111_path_prefix.endswith('/'):
                s111_path_prefix += '/'

            file_date = re.findall(r'\d{8}', model_filename)[0]
            cycletime = re.findall(r'(?<=t)[^t:]+(?=:?z)', model_filename)[0]
            file_datetime = file_date + cycletime
            file_issuance = f"{file_date}T{cycletime}Z"

            dt = datetime.datetime.strptime(file_datetime, '%Y%m%d%H')

            transit = datetime.datetime.strptime('19870101', '%Y%m%d')
            g730 = datetime.datetime.strptime('19940629', '%Y%m%d')
            g873 = datetime.datetime.strptime('19970129', '%Y%m%d')
            g1150 = datetime.datetime.strptime('20020120', '%Y%m%d')
            g1674 = datetime.datetime.strptime('20120208', '%Y%m%d')
            g1762 = datetime.datetime.strptime('20131016', '%Y%m%d')

            if dt < g730:
                epoch = 'TRANSIT'
            elif dt < g873:
                epoch = 'G730'
            elif dt < g1150:
                epoch = 'G873'
            elif dt < g1674:
                epoch = 'G1150'
            elif dt < g1762:
                epoch = 'G1674'
            elif dt >= g1762:
                epoch = 'G1762'

            s111_path_prefix += f'S111{PRODUCER_CODE}_{file_issuance}_{model_name.upper()}_TYP{DATA_CODING_FORMAT}'

        try:
            model_index.open()
            model_file.open()

            with S111File(f'{s111_path_prefix}.h5', "w", driver=None) as s111_file:

                s111_file.create_empty_metadata()
                root = s111_file.root

                root.product_specification = "INT.IHO.S-111.1.0"
                root.metadata = 'MD_{}.XML'.format(os.path.split(s111_path_prefix)[-1])
                root.horizontal_datum_reference = "EPSG"
                root.horizontal_datum_value = 4326
                root.epoch = epoch
                root.geographic_identifier = MODELS[model_name]['region']
                now = datetime.datetime.now()
                root.issue_date = now.strftime('%Y%m%d')
                root.issue_time = now.strftime('%H%M%SZ')

                # Additional S-111 root metadata
                root.surface_current_depth = -1 * DEFAULT_TARGET_DEPTH
                root.depth_type_index = 2

                del root.vertical_datum
                del root.extent_type_code
                del root.meta_features

                # Add Group_F, featureCode, and SurfaceCurrent feature information dataset
                surf_current_feature = root.feature_information.surface_current_dataset

                surface_current_speed_info = surf_current_feature.append_new_item()
                surface_current_speed_info.initialize_properties(True)
                surface_current_speed_info.code = "surfaceCurrentSpeed"
                surface_current_speed_info.name = "Surface current speed"
                surface_current_speed_info.unit_of_measure = "knots"
                surface_current_speed_info.lower = "0.0"
                surface_current_speed_info.upper = ""
                surface_current_speed_info.closure = "geSemiInterval"

                surface_current_direction_info = surf_current_feature.append_new_item()
                surface_current_direction_info.initialize_properties(True)
                surface_current_direction_info.code = "surfaceCurrentDirection"
                surface_current_direction_info.name = "Surface current direction"
                surface_current_direction_info.unit_of_measure = "arc-degrees"
                surface_current_speed_info.lower = "0"
                surface_current_direction_info.upper = "360"
                surface_current_direction_info.closure = "geLtInterval"

                root.surface_current.axis_names = numpy.array(["longitude", "latitude"])
                root.surface_current.common_point_rule = 3
                root.surface_current.data_coding_format = 2
                root.surface_current.dimension = 2
                root.surface_current.sequencing_rule_scan_direction = "Longitude, Latitude"
                root.surface_current.interpolation_type = 10
                root.surface_current.num_instances = 1
                root.surface_current.sequencing_rule_type = 1
                root.surface_current.time_uncertainty = -1.0
                root.surface_current.vertical_uncertainty = -1.0
                root.surface_current.horizontal_position_uncertainty = -1.0

                # Additional SurfaceCurrentContainer Class Metadata
                root.surface_current.type_of_current_data = 6
                root.surface_current.method_currents_product = MODELS[model_name]['product']

                # Add Feature Instance
                surface_current_01 = root.surface_current.surface_current.append_new_item()
                surface_current_01.initialize_properties(True)

                # Feature Instance Uncertainty Dataset
                speed_uncertainty = surface_current_01.uncertainty_dataset.append_new_item()
                speed_uncertainty.name = "surfaceCurrentSpeed"
                speed_uncertainty.value = -1.0
                direction_uncertainty = surface_current_01.uncertainty_dataset.append_new_item()
                direction_uncertainty.name = "surfaceCurrentDirection"
                direction_uncertainty.value = -1.0

                del surface_current_01.grid_spacing_vertical
                del surface_current_01.grid_origin_vertical
                del surface_current_01.extent_type_code
                del surface_current_01.num_points_vertical
                del surface_current_01.vertical_extent_maximum_z
                del surface_current_01.vertical_extent_minimum_z

                # Additional Feature Instance Attributes
                surface_current_01.start_sequence = "0,0"
                time_str = model_file.datetime_values[0].strftime('%Y%m%dT%H%M%SZ')
                surface_current_01.datetime_first_record = numpy.string_(time_str)

                for time_index in range(len(model_file.datetime_values)):
                    # Get native-grid output with invalid/masked values removed
                    reg_grid_u, reg_grid_v = model_file.uv_to_regular_grid(model_index, time_index,
                                                                           DEFAULT_TARGET_DEPTH)

                    reg_grid_u = numpy.ma.masked_array(reg_grid_u, model_index.var_mask.mask)
                    reg_grid_v = numpy.ma.masked_array(reg_grid_v, model_index.var_mask.mask)

                    # Convert currents at regular grid points from u/v to speed/direction
                    speed, direction = model.regular_uv_to_speed_direction(reg_grid_u, reg_grid_v)

                    # Apply mask
                    direction = numpy.ma.masked_array(direction, model_index.var_mask.mask)
                    speed = numpy.ma.masked_array(speed, model_index.var_mask.mask)

                    # If any valid data points fall outside of the scipy griddata convex hull
                    # nan values will be used, if nan values are present
                    # add nan values to the original mask
                    if numpy.isnan(speed).any():
                        nan_mask_speed = numpy.ma.masked_invalid(speed)
                        nan_mask_direction = numpy.ma.masked_invalid(direction)
                        speed_mask = numpy.ma.mask_or(model_index.var_mask.mask, nan_mask_speed.mask)
                        direction_mask = numpy.ma.mask_or(model_index.var_mask.mask, nan_mask_direction.mask)

                        speed = numpy.ma.masked_array(speed, speed_mask)
                        direction = numpy.ma.masked_array(direction, direction_mask)

                    min_speed = numpy.round(numpy.nanmin(speed), decimals=2)
                    max_speed = numpy.round(numpy.nanmax(speed), decimals=2)

                    if numpy.ma.is_masked(speed):
                        speed = speed.filled(FILLVALUE)
                        direction = direction.filled(FILLVALUE)

                    # Format speed/direction
                    speed = numpy.round(speed, decimals=2)
                    direction = numpy.round(direction, decimals=1)

                    cellsize_x = model_index.var_x[1] - model_index.var_x[0]
                    cellsize_y = model_index.var_y[1] - model_index.var_y[0]

                    nx = model_index.dim_x.size
                    ny = model_index.dim_y.size

                    minx = numpy.nanmin(numpy.round(model_index.var_x, 7))
                    maxx = numpy.nanmax(numpy.round(model_index.var_x, 7))
                    miny = numpy.nanmin(numpy.round(model_index.var_y, 7))
                    maxy = numpy.nanmax(numpy.round(model_index.var_y, 7))

                    root.east_bound_longitude = minx
                    root.west_bound_longitude = maxx
                    root.south_bound_latitude = miny
                    root.north_bound_latitude = maxy

                    if min_speed < root.surface_current.min_dataset_current_speed:
                        root.surface_current.min_dataset_current_speed = min_speed

                    if max_speed > root.surface_current.max_dataset_current_speed:
                        root.surface_current.max_dataset_current_speed = max_speed

                    surface_current_01.datetime_last_record = model_file.datetime_values[time_index].strftime('%Y%m%dT%H%M%SZ')

                    surface_current_01.num_grp = len(model_file.datetime_values)
                    surface_current_01.number_of_times = len(model_file.datetime_values)

                    interval = model_file.datetime_values[1] - model_file.datetime_values[0]
                    surface_current_01.time_record_interval = interval.total_seconds()

                    surface_current_01.east_bound_longitude = minx
                    surface_current_01.west_bound_longitude = maxx
                    surface_current_01.south_bound_latitude = miny
                    surface_current_01.north_bound_latitude = maxy
                    surface_current_01.grid_origin_latitude = miny

                    surface_current_01.grid_origin_longitude = minx
                    surface_current_01.grid_origin_latitude = miny
                    surface_current_01.grid_spacing_longitudinal = cellsize_x
                    surface_current_01.grid_spacing_latitudinal = cellsize_y

                    surface_current_01.num_points_latitudinal = ny
                    surface_current_01.num_points_longitudinal = nx

                    # Add Group(s)
                    surface_current_group_object = surface_current_01.surface_current_group.append_new_item()
                    surface_current_group_object.values_create()
                    grid = surface_current_group_object.values
                    grid.surface_current_speed = speed
                    grid.surface_current_direction = direction
                    # Add Group(s) Attribute
                    surface_current_group_object.time_point = model_file.datetime_values[time_index].strftime('%Y%m%dT%H%M%SZ')

                    # TODO: Determine group values dataset chunk sizes(e.g grid.chunks)
                    # chunks = values_dset.chunks
                    # chunking_str = ','.join(str(x) for x in chunks)
                    # groupF_dset.attrs.create('chunking', chunking_str, dtype=h5py.special_dtype(vlen=str))
                    # feature_instance.attrs.create('instanceChunking', numpy.string_(chunking_str))

                s111_file.write()

        finally:
            model_index.close()
            model_file.close()


def main():
    fire.Fire(CLI())


if __name__ == "__main__":
    main()
