import os
import numpy
import datetime
import fire
# import logging

from thyme.model import model, roms
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

PRODUCT_DESCRIPTION_FVCOM = 'FVCOM_Hydrodynamic_Model_Forecasts'
PRODUCT_DESCRIPTION_HYCOM = 'HYCOM_Hydrodynamic_Model_Forecasts'
PRODUCT_DESCRIPTION_POM = 'POM_Hydrodynamic_Model_Forecasts'
PRODUCT_DESCRIPTION_ROMS = 'ROMS_Hydrodynamic_Model_Forecasts'


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
        model_file = roms.ROMSFile(model_file_path)
        model_index = roms.ROMSIndexFile(model_index_path)
        s111_path_prefix = os.path.split(output_path)[0]

        # Path format/prefix for output S111 files. Forecast initialization (reference).
        if os.path.isdir(s111_path_prefix):
            if not s111_path_prefix.endswith('/'):
                s111_path_prefix += '/'

            file_issuance = "20200129T12Z"

            s111_path_prefix += (
                f'S111{PRODUCER_CODE}_{file_issuance}_CBOFS_TYP{DATA_CODING_FORMAT}')

        try:
            model_index.open()
            model_file.open()
            for time_index in range(len(model_file.datetime_values)):
                # Get native-grid output with invalid/masked values removed
                reg_grid_u, reg_grid_v = model_file.uv_to_regular_grid(model_index, time_index, DEFAULT_TARGET_DEPTH)

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

                with S111File(f'{s111_path_prefix}.h5', "w", driver=None) as s111_file:

                    s111_file.create_empty_metadata()
                    root = s111_file.root

                    root.product_specification = "INT.IHO.S-111.1.0.1"
                    root.metadata = 'MD_{}.XML'.format(os.path.split(s111_path_prefix)[-1])
                    root.horizontal_datum_reference = "EPSG"
                    root.horizontal_datum_value = 4326
                    root.epoch = "G1762"
                    root.geographic_identifier = "Chesapeake Bay"
                    now = datetime.datetime.now()
                    root.issue_date = now.strftime('%Y%m%d')
                    root.issue_time = now.strftime('%H%M%SZ')
                    root.east_bound_longitude = minx
                    root.west_bound_longitude = maxx
                    root.south_bound_latitude = miny
                    root.north_bound_latitude = maxy

                    # Additional S-111 root metadata
                    root.surface_current_depth = -4.5
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

                    # Add Feature
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
                    root.surface_current.min_dataset_current_speed = min_speed
                    root.surface_current.max_dataset_current_speed = max_speed
                    root.surface_current.method_currents_product = PRODUCT_DESCRIPTION_ROMS

                    # Add Feature Instance
                    surface_current_01 = root.surface_current.surface_current.append_new_item()
                    surface_current_01.initialize_properties(True)

                    surface_current_01.east_bound_longitude = minx
                    surface_current_01.west_bound_longitude = maxx
                    surface_current_01.south_bound_latitude = miny
                    surface_current_01.north_bound_latitude = maxy
                    surface_current_01.grid_origin_latitude = miny

                    surface_current_01.grid_origin_longitude = minx
                    surface_current_01.grid_origin_latitude = miny
                    surface_current_01.grid_spacing_longitudinal = cellsize_x
                    surface_current_01.grid_spacing_latitudinal = cellsize_y
                    surface_current_01.num_grp = 1
                    surface_current_01.start_sequence = "0,0"
                    surface_current_01.num_points_latitudinal = ny
                    surface_current_01.num_points_longitudinal = nx

                    # TODO: Add uncertainty dataset
                    # surface_current_01.uncertainty_dataset.append_new_item()

                    del surface_current_01.grid_spacing_vertical
                    del surface_current_01.grid_origin_vertical
                    del surface_current_01.extent_type_code
                    del surface_current_01.num_points_vertical
                    del surface_current_01.vertical_extent_maximum_z
                    del surface_current_01.vertical_extent_minimum_z

                    # Additional Feature Instance Attributes
                    surface_current_01.number_of_times = 1
                    surface_current_01.time_record_interval = 0
                    surface_current_01.datetime_first_record = numpy.string_(file_issuance)
                    surface_current_01.datetime_last_record = numpy.string_(file_issuance)

                    # Add Group_001
                    surface_current_group_object = surface_current_01.surface_current_group.append_new_item()
                    surface_current_group_object.values_create()
                    grid = surface_current_group_object.values

                    # Add Group Attribute
                    surface_current_group_object.time_point = numpy.string_(file_issuance)

                    # Add data
                    grid.surface_current_speed = speed
                    grid.surface_current_direction = direction

                    # TODO: Determine group values dataset chunk sizes(e.g grid.chunks)
                    # chunks = values_dset.chunks
                    # chunking_str = ','.join(str(x) for x in chunks)
                    # groupF_dset.attrs.create('chunking', chunking_str, dtype=h5py.special_dtype(vlen=str))
                    # feature_instance.attrs.create('instanceChunking', numpy.string_(chunking_str))

                    # TODO: Add additional groups
                    # surface_current_group_object_2 = surface_current_01.surface_current_group.append_new_item()
                    # surface_current_group_object_2.values_create()
                    # grid_2 = surface_current_group_object_2.values
                    #
                    # grid_2.surface_current_speed = speed
                    # grid_2.surface_current_direction = direction

                    s111_file.write()

        finally:
            model_index.close()
            model_file.close()


def main():
    fire.Fire(CLI())


if __name__ == "__main__":
    main()
