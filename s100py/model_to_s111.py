import os
import numpy
import datetime
import fire
import re

from thyme.model import model, roms, fvcom, pom, hycom
from thyme.util import dateutil
from s100py.s111 import S111File

# Default fill value for NetCDF variables
FILLVALUE = -9999.0

# Default depth in meters
DEFAULT_TARGET_DEPTH = 4.5
PRODUCER_CODE = "US"
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
        'datetime_rounding': None
    },
    'nyofs': {

        'region': 'Port_of_New_York_and_New_Jersey',
        'product': PRODUCT_DESCRIPTION_POM,
        'model_type': MODELTYPE_POM,
        'datetime_rounding': dateutil.DatetimeRounding.NEAREST_HOUR
    }
}


class CLI:
    """
    Container for methods exposed through CLI via Python Fire.
    """
    def __init__(self):
        pass

    @staticmethod
    def determine_epoch(dt):
        """Convert NetCDF hydrodynamic model to S111 format.

        Args:
            dt: Date data was collected.

        Returns:
            Str: Epoch attribute.
        """
        epoch = None
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

        return epoch

    @staticmethod
    def convert_regular(model_file, model_index, time_index, idx):
        try:
            model_index.open()
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

            cellsize_x = model_index.var_x[1] - model_index.var_x[0]
            cellsize_y = model_index.var_y[1] - model_index.var_y[0]

            if model_index.dim_subgrid is not None and model_index.var_subgrid_id is not None:

                x_min = model_index.var_subgrid_x_min[idx]
                x_max = model_index.var_subgrid_x_max[idx]
                y_min = model_index.var_subgrid_y_min[idx]
                y_max = model_index.var_subgrid_y_max[idx]
                speed = speed[y_min:y_max + 1, x_min:x_max + 1]
                direction = direction[y_min:y_max + 1, x_min:x_max + 1]

                nx = 1 + model_index.var_subgrid_x_max[idx] - model_index.var_subgrid_x_min[idx]
                ny = 1 + model_index.var_subgrid_y_max[idx] - model_index.var_subgrid_y_min[idx]
                minx = model_index.var_x[model_index.var_subgrid_x_min[idx]]
                miny = model_index.var_y[model_index.var_subgrid_y_min[idx]]
                maxx = model_index.var_x[model_index.var_subgrid_x_max[idx]]
                maxy = model_index.var_y[model_index.var_subgrid_y_max[idx]]

            else:

                nx = model_index.dim_x.size
                ny = model_index.dim_y.size

                minx = numpy.nanmin(numpy.round(model_index.var_x, 7))
                maxx = numpy.nanmax(numpy.round(model_index.var_x, 7))
                miny = numpy.nanmin(numpy.round(model_index.var_y, 7))
                maxy = numpy.nanmax(numpy.round(model_index.var_y, 7))

        finally:
            model_index.close()

        return speed, direction, cellsize_x, cellsize_y, nx, ny, minx, maxx, miny, maxy

    @staticmethod
    def convert_irregular(model_file, time_index):
        # Get native-grid output with invalid/masked values removed
        u, v, latitude, longitude = model_file.output_native_grid(time_index, DEFAULT_TARGET_DEPTH)

        # Convert currents from u/v to speed/direction
        speed, direction = model.irregular_uv_to_speed_direction(u, v)

        minx = numpy.nanmin(numpy.round(longitude, 7))
        maxx = numpy.nanmax(numpy.round(longitude, 7))
        miny = numpy.nanmin(numpy.round(latitude, 7))
        maxy = numpy.nanmax(numpy.round(latitude, 7))

        return speed, direction, longitude, latitude, minx, maxx, miny, maxy

    def convert(self, model_file_path, output_path, data_coding_format, model_index_path=None):
        """Convert NetCDF hydrodynamic model to S111 format.

        Args:
            model_file_path: Path to native ofs model file.
            output_path:  Path to output hdf5 file.
            data_coding_format: 1: Time series at fixed stations
                                2: Regularly-gridded arrays
                                3: Ungeorectified gridded arrays
                                4: Moving platform
                                5. Irregular grid
                                6. Variable cell size
                                7. TIN
            model_index_path: Optional: Path to model index file.

        Returns:
            List of paths to HDF5 files created.
        """
        s111_path_prefix = os.path.split(output_path)[0]
        model_filename = os.path.split(model_file_path)[-1]
        model_name = model_filename.split('.')[1]

        s111_filenames = []
        subgrid_index = []

        if os.path.isdir(s111_path_prefix):
            if not s111_path_prefix.endswith('/'):
                s111_path_prefix += '/'

            file_date = re.findall(r'\d{8}', model_filename)[0]
            cycletime = re.findall(r'(?<=t)[^t:]+(?=:?z)', model_filename)[0]
            file_datetime = file_date + cycletime
            file_issuance = f"{file_date}T{cycletime}Z"

            dt = datetime.datetime.strptime(file_datetime, '%Y%m%d%H')
            epoch = self.determine_epoch(dt)

            if model_index_path is not None:
                try:
                    model_index = MODEL_INDEX_CLASS[MODELS[model_name]['model_type']](model_index_path)
                    model_index.open()

                    if model_index.dim_subgrid is not None and model_index.var_subgrid_id is not None:
                        for i in range(model_index.dim_subgrid.size):

                            if model_index.var_subgrid_name is not None:
                                x_min = model_index.var_subgrid_x_min[i]
                                x_max = model_index.var_subgrid_x_max[i]
                                y_min = model_index.var_subgrid_y_min[i]
                                y_max = model_index.var_subgrid_y_max[i]

                                subgrid_valid = model_index.var_mask[y_min:y_max + 1, x_min:x_max + 1]
                                valid_occurrences = numpy.count_nonzero(subgrid_valid == 1)

                                if valid_occurrences >= 20:
                                    subgrid_index.append(i)

                                if model_index.var_subgrid_name is not None:
                                    s111_subgrid_filename = f'S111{PRODUCER_CODE}_{file_issuance}_{model_name.upper()}_TYP{data_coding_format}_{model_index.var_subgrid_name[i]}'
                                    s111_filenames.append(s111_subgrid_filename)

                                else:
                                    s111_subgrid_filename = f'S111{PRODUCER_CODE}_{file_issuance}_{model_name.upper()}_TYP{data_coding_format}_FID_{model_index.var_subgrid_name[i]}'
                                    s111_filenames.append(s111_subgrid_filename)

                    else:
                        s111_filename = f'S111{PRODUCER_CODE}_{file_issuance}_{model_name.upper()}_TYP{data_coding_format}'
                        s111_filenames.append(s111_filename)

                finally:
                    model_index.close()

            else:
                s111_subgrid_filename = f'S111{PRODUCER_CODE}_{file_issuance}_{model_name.upper()}_TYP{data_coding_format}'
                s111_filenames.append(s111_subgrid_filename)

        model_file = MODEL_FILE_CLASS[MODELS[model_name]['model_type']](model_file_path, datetime_rounding=MODELS[model_name]['datetime_rounding'])

        if data_coding_format == 2:
            model_index = MODEL_INDEX_CLASS[MODELS[model_name]['model_type']](model_index_path)

        try:
            model_file.open()

            for idx, file in enumerate(s111_filenames):
                if idx in subgrid_index or not subgrid_index:

                    with S111File(f'{s111_path_prefix}{file}.h5', "w") as s111_file:

                        root = s111_file.root

                        root.surface_current_create()
                        surface_current_feature = root.surface_current

                        surface_current_feature.surface_current_create()
                        surface_current_feature_instance_01 = surface_current_feature.surface_current.append_new_item()

                        surface_current_feature_instance_01.surface_current_group_create()

                        surface_current_feature_instance_01.uncertainty_dataset_create()
                        speed_uncertainty = surface_current_feature_instance_01.uncertainty_dataset.append_new_item()
                        speed_uncertainty.name = "surfaceCurrentSpeed"
                        speed_uncertainty.value = -1.0
                        direction_uncertainty = surface_current_feature_instance_01.uncertainty_dataset.append_new_item()
                        direction_uncertainty.name = "surfaceCurrentDirection"
                        direction_uncertainty.value = -1.0

                        root.feature_information_create()
                        group_f = root.feature_information
                        group_f.feature_code_create()
                        group_f.surface_current_feature_dataset_create()

                        surface_current_feature_dataset = root.feature_information.surface_current_feature_dataset

                        surface_current_speed_info = surface_current_feature_dataset.append_new_item()
                        surface_current_speed_info.code = "surfaceCurrentSpeed"
                        surface_current_speed_info.name = "Surface current speed"
                        surface_current_speed_info.unit_of_measure = "knots"
                        surface_current_speed_info.datatype = "H5T_FLOAT"
                        surface_current_speed_info.fill_value = FILLVALUE
                        surface_current_speed_info.lower = "0.0"
                        surface_current_speed_info.upper = ""
                        surface_current_speed_info.closure = "geSemiInterval"

                        surface_current_direction_info = surface_current_feature_dataset.append_new_item()
                        surface_current_direction_info.code = "surfaceCurrentDirection"
                        surface_current_direction_info.name = "Surface current direction"
                        surface_current_direction_info.unit_of_measure = "arc-degrees"
                        surface_current_direction_info.datatype = "H5T_FLOAT"
                        surface_current_direction_info.fill_value = FILLVALUE
                        surface_current_direction_info.lower = "0"
                        surface_current_direction_info.upper = "360"
                        surface_current_direction_info.closure = "geLtInterval"

                        root.product_specification = "INT.IHO.S-111.1.0"
                        root.metadata = f'MD_{file}.XML'
                        root.horizontal_datum_reference = "EPSG"
                        root.horizontal_datum_value = 4326
                        root.epoch = epoch
                        root.geographic_identifier = MODELS[model_name]['region']
                        utc_now = datetime.datetime.utcnow()
                        root.issue_date = utc_now.strftime('%Y%m%d')
                        root.issue_time = utc_now.strftime('%H%M%SZ')
                        root.surface_current_depth = -1 * DEFAULT_TARGET_DEPTH
                        root.depth_type_index = 2

                        surface_current_feature.axis_names = numpy.array(["longitude", "latitude"])
                        surface_current_feature.common_point_rule = 3
                        surface_current_feature.data_coding_format = data_coding_format
                        surface_current_feature.interpolation_type = 10
                        surface_current_feature.num_instances = 1
                        surface_current_feature.time_uncertainty = -1.0
                        surface_current_feature.vertical_uncertainty = -1.0
                        surface_current_feature.horizontal_position_uncertainty = -1.0
                        surface_current_feature.type_of_current_data = 6
                        surface_current_feature.method_currents_product = MODELS[model_name]['product']

                        first_record = model_file.datetime_values[0].strftime('%Y%m%dT%H%M%SZ')
                        surface_current_feature_instance_01.date_time_of_first_record = numpy.string_(first_record)
                        surface_current_feature.min_dataset_current_speed = 0
                        surface_current_feature.max_dataset_current_speed = 0

                        for time_index in range(len(model_file.datetime_values)):
                            if data_coding_format == 2:
                                speed, direction, cellsize_x, cellsize_y, nx, ny, minx, maxx, miny, maxy = self.convert_regular(model_file, model_index, time_index, idx)

                                surface_current_feature_instance_01.start_sequence = "0,0"
                                surface_current_feature.sequencing_rule_scan_direction = "longitude, latitude"
                                surface_current_feature.sequencing_rule_type = 1
                                surface_current_feature_instance_01.grid_origin_longitude = minx
                                surface_current_feature_instance_01.grid_origin_latitude = miny
                                surface_current_feature_instance_01.grid_spacing_longitudinal = cellsize_x
                                surface_current_feature_instance_01.grid_spacing_latitudinal = cellsize_y

                                surface_current_feature_instance_01.num_points_latitudinal = ny
                                surface_current_feature_instance_01.num_points_longitudinal = nx

                            if data_coding_format == 3:
                                speed, direction, longitude, latitude, minx, maxx, miny, maxy = self.convert_irregular(model_file, time_index)

                                surface_current_feature_instance_01.number_of_nodes = longitude.size

                                surface_current_feature_instance_01.positioning_group_create()
                                positioning = surface_current_feature_instance_01.positioning_group
                                positioning.geometry_values_create()
                                geometry_values = positioning.geometry_values
                                geometry_values.longitude = longitude
                                geometry_values.latitude = latitude

                            root.east_bound_longitude = minx
                            root.west_bound_longitude = maxx
                            root.south_bound_latitude = miny
                            root.north_bound_latitude = maxy
                            root.surface_current.dimension = speed.ndim

                            min_speed = numpy.round(numpy.nanmin(speed), decimals=2)
                            max_speed = numpy.round(numpy.nanmax(speed), decimals=2)

                            if numpy.ma.is_masked(speed):
                                speed = speed.filled(FILLVALUE)
                                direction = direction.filled(FILLVALUE)

                            speed = numpy.round(speed, decimals=2)
                            direction = numpy.round(direction, decimals=1)

                            if min_speed < surface_current_feature.min_dataset_current_speed:
                                surface_current_feature.min_dataset_current_speed = min_speed

                            if max_speed > surface_current_feature.max_dataset_current_speed:
                                surface_current_feature.max_dataset_current_speed = max_speed

                            last_record = model_file.datetime_values[time_index].strftime('%Y%m%dT%H%M%SZ')
                            surface_current_feature_instance_01.date_time_of_last_record = numpy.string_(last_record)

                            surface_current_feature_instance_01.num_grp = len(model_file.datetime_values)
                            surface_current_feature_instance_01.number_of_times = len(model_file.datetime_values)

                            if len(model_file.datetime_values) == 1:
                                surface_current_feature_instance_01.time_record_interval = 0
                            else:
                                interval = model_file.datetime_values[1] - model_file.datetime_values[0]
                                surface_current_feature_instance_01.time_record_interval = int(interval.total_seconds())

                            surface_current_feature_instance_01.east_bound_longitude = minx
                            surface_current_feature_instance_01.west_bound_longitude = maxx
                            surface_current_feature_instance_01.south_bound_latitude = miny
                            surface_current_feature_instance_01.north_bound_latitude = maxy

                            surface_current_group_object = surface_current_feature_instance_01.surface_current_group.append_new_item()
                            surface_current_group_object.values_create()
                            grid = surface_current_group_object.values
                            grid.surface_current_speed = speed
                            grid.surface_current_direction = direction

                            surface_current_group_object.time_point = model_file.datetime_values[time_index].strftime('%Y%m%dT%H%M%SZ')

                        s111_file.write()

        finally:
            model_file.close()


def main():
    fire.Fire(CLI())


if __name__ == "__main__":
    main()
