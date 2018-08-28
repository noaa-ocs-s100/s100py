"""Utilities for creation and modification of S-111 compliant HDF-5 files.

S-111 is an IHO standard outlining formats for storing and sending surface
water current data and metadata.
"""
import contextlib
import datetime
import os
import numpy
import numpy.ma as ma

from .model import model

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


# Default fill value for NetCDF variables
FILLVALUE = -9999.0
# Default depth in meters
DEFAULT_TARGET_DEPTH = 4.5

# Dict lookup for HDF5 data type class names
# See https://github.com/h5py/h5py/blob/master/h5py/api_types_hdf5.pxd#L509
H5T_CLASS_T = {
    h5py.h5t.NO_CLASS: "H5T_NO_CLASS",
    h5py.h5t.INTEGER: "H5T_INTEGER",
    h5py.h5t.FLOAT: "H5T_FLOAT",
    h5py.h5t.TIME: "H5T_TIME",
    h5py.h5t.STRING: "H5T_STRING",
    h5py.h5t.BITFIELD: "H5T_BITFIELD",
    h5py.h5t.OPAQUE: "H5T_OPAQUE",
    h5py.h5t.COMPOUND: "H5T_COMPOUND",
    h5py.h5t.REFERENCE: "H5T_REFERENCE",
    h5py.h5t.ENUM: "H5T_ENUM",
    h5py.h5t.VLEN: "H5T_VLEN",
    h5py.h5t.ARRAY: "H5T_ARRAY"
}


class S111File:
    """Create and manage S-111 files.

    This class implements the context manager pattern and should thus be used
    similar to the following:

        with S111File("myfile.h5") as f:
            ...

    Attributes:
        path: Path (relative or absolute) to the file.
    """

    def __init__(self, path, clobber=False):
        """Initializes S111File object and opens h5 file at specified path.

        If `path` has an extension other than ".h5", it is replaced with
        ".h5"

        Args:
            path: Path of target hdf5 file. Must end in ".h5", otherwise its
                extension will be replaced with ".h5".
            clobber: (Optional, default False) If True, existing h5 file at
                specified path, if any, will be deleted and the new file will
                be opened in write mode.
        """
        filepath, file_extension = os.path.splitext(path)
        self.path = filepath + ".h5"

        if not os.path.exists(self.path) or clobber:
            # File doesn't exist, open in create (write) mode and add metadata
            self.h5_file = h5py.File(self.path, "w")
            # Create s111 structure, feature group, feature type container, and initial feature instance
            self.groupF = self.h5_file.create_group('Group_F')
            self.feature = self.h5_file.create_group('SurfaceCurrent')
            self.feature_instance = self.feature.create_group('SurfaceCurrent.01')
            self.groupF_dset = None
            # Add s111 structure and metadata
            self.add_structure_and_metadata()
        else:
            # File already exists, open in append mode
            self.h5_file = h5py.File(self.path, "r+")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5_file.close()

    def close(self):
        self.h5_file.flush()
        self.h5_file.close()

    def add_structure_and_metadata(self):
        """Add S100/S111 structure and empty global metadata attributes to HDF5 file."""

        # Add root group carrier metadata
        # String types
        self.h5_file.attrs.create('productSpecification', '', dtype=h5py.special_dtype(vlen=str))
        self.h5_file.attrs.create('issueTime', '', dtype=h5py.special_dtype(vlen=str))
        self.h5_file.attrs.create('issueDate', '', dtype=h5py.special_dtype(vlen=str))
        self.h5_file.attrs.create('epoch', '', dtype=h5py.special_dtype(vlen=str))
        self.h5_file.attrs.create('geographicIdentifier', '', dtype=h5py.special_dtype(vlen=str))
        self.h5_file.attrs.create('metadata', '', dtype=h5py.special_dtype(vlen=str))
        self.h5_file.attrs.create('horizontalDatumReference', '', dtype=h5py.special_dtype(vlen=str))
        # Integer types
        self.h5_file.attrs.create('horizontalDatumValue', 0, dtype=numpy.int32)
        # Enumeration types
        self.h5_file.attrs.create('depthTypeIndex', 0, dtype=numpy.int32)
        # Real types
        self.h5_file.attrs.create('surfaceCurrentDepth', 0, dtype=numpy.float32)
        self.h5_file.attrs.create('westBoundLongitude', 0, dtype=numpy.float32)
        self.h5_file.attrs.create('eastBoundLongitude', 0, dtype=numpy.float32)
        self.h5_file.attrs.create('southBoundLatitude', 0, dtype=numpy.float32)
        self.h5_file.attrs.create('northBoundLatitude', 0, dtype=numpy.float32)

        # Add feature container metadata
        # Integer types
        self.feature.attrs.create('dimension', 0, dtype=numpy.int32)
        self.feature.attrs.create('numInstances', 0, dtype=numpy.int32)
        # Real types
        self.feature.attrs.create('horizontalPositionUncertainty', -1.0, dtype=numpy.float32)
        self.feature.attrs.create('verticalUncertainty', -1.0, dtype=numpy.float32)
        self.feature.attrs.create('timeUncertainty', -1.0, dtype=numpy.float32)
        self.feature.attrs.create('speedUncertainty', -1.0, dtype=numpy.float32)
        self.feature.attrs.create('directionUncertainty', -1.0, dtype=numpy.float32)
        self.feature.attrs.create('minDatasetCurrentSpeed', 0, dtype=numpy.float32)
        self.feature.attrs.create('maxDatasetCurrentSpeed', 0, dtype=numpy.float32)
        # Enumeration types
        self.feature.attrs.create('dataCodingFormat', 0, dtype=numpy.int32)
        self.feature.attrs.create('commonPointRule', 0, dtype=numpy.int32)
        self.feature.attrs.create('interpolationType', 0, dtype=numpy.int32)
        self.feature.attrs.create('typeOfCurrentData', 0, dtype=numpy.int32)
        self.feature.attrs.create('sequenceRule.type', 0, dtype=numpy.int32)
        # String types
        self.feature.attrs.create('methodCurrentsProduct', '', dtype=h5py.special_dtype(vlen=str))
        self.feature.attrs.create('sequenceRule.scanDirection', '', dtype=h5py.special_dtype(vlen=str))

        # Add feature instance metadata
        # String types
        self.feature_instance.attrs.create('startSequence', '', dtype=h5py.special_dtype(vlen=str))
        self.feature_instance.attrs.create('dateTimeOfFirstRecord', '', dtype=h5py.special_dtype(vlen=str))
        self.feature_instance.attrs.create('dateTimeOfLastRecord', '', dtype=h5py.special_dtype(vlen=str))
        self.feature_instance.attrs.create('instanceChunking', '', dtype=h5py.special_dtype(vlen=str))
        # Integer types
        self.feature_instance.attrs.create('numPointsLongitudinal', 0, dtype=numpy.int32)
        self.feature_instance.attrs.create('numPointsLatitudinal', 0, dtype=numpy.int32)
        self.feature_instance.attrs.create('timeRecordInterval', 0, dtype=numpy.int32)
        self.feature_instance.attrs.create('numberOfTimes', 0, dtype=numpy.int32)
        self.feature_instance.attrs.create('numGRP', 0, dtype=numpy.int32)
        # Real types
        self.feature_instance.attrs.create('gridOriginLongitude', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('gridOriginLatitude', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('gridSpacingLongitudinal', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('gridSpacingLatitudinal', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('westBoundLongitude', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('eastBoundLongitude', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('southBoundLatitude', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('northBoundLatitude', 0, dtype=numpy.float32)

    def update_attributes(self, model_index, ofs_metadata, target_depth, subgrid_index=None):
        """Update HDF5 attributes based on grid properties.
        
        Args:
            model_index: `ModelIndexFile` instance representing model index file.
            subgrid_index: (Optional, default None) Index of subgrid, if any,
                that this S111File represents. Corresponds with index into
                subgrid dimension of model index file.
            ofs_metadata: `S111Metadata` instance describing metadata for geographic
                identifier and description of current meter type, forecast method,
                or model.
            target_depth: The water current at a specified target depth below
                the sea surface in meters, default target depth is 4.5 meters,
                target interpolation depth must be greater or equal to 0.
        """

        # Width between first two cells, grid spacing is uniform
        cellsize_x = model_index.var_x[1] - model_index.var_x[0]
        cellsize_y = model_index.var_y[1] - model_index.var_y[0]

        if subgrid_index is not None:
            if subgrid_index < 0 or subgrid_index >= model_index.dim_subgrid.size:
                raise Exception("Subgrid index [{}] out of model index subgrid dimension range [0-{}]".format(subgrid_index, model_index.dim_subgrid.size - 1))
            num_points_lon = 1 + model_index.var_subgrid_x_max[subgrid_index] - model_index.var_subgrid_x_min[
                subgrid_index]
            num_points_lat = 1 + model_index.var_subgrid_y_max[subgrid_index] - model_index.var_subgrid_y_min[
                subgrid_index]
            min_lon = model_index.var_x[model_index.var_subgrid_x_min[subgrid_index]]
            min_lat = model_index.var_y[model_index.var_subgrid_y_min[subgrid_index]]
            max_lon = model_index.var_x[model_index.var_subgrid_x_max[subgrid_index]]
            max_lat = model_index.var_y[model_index.var_subgrid_y_max[subgrid_index]]
        else:
            num_points_lon = model_index.dim_x.size
            num_points_lat = model_index.dim_y.size
            min_lon = numpy.nanmin(model_index.var_x)
            min_lat = numpy.nanmin(model_index.var_y)
            max_lon = numpy.nanmax(model_index.var_x)
            max_lat = numpy.nanmax(model_index.var_y)

        # X/Y coordinates are located at the center of each grid cell
        # Adjust X/Y coordinates to the lower left/upper right corner for compliance
        # with S111 product specification by adding/subtracting 1/2 pixel
        # size from min/max X/Y coordinates

        min_lon -= cellsize_x / 2
        min_lat -= cellsize_y / 2
        max_lon += cellsize_x / 2
        max_lat += cellsize_y / 2

        min_lon = numpy.round(min_lon, 7)
        min_lat = numpy.round(min_lat, 7)
        max_lon = numpy.round(max_lon, 7)
        max_lat = numpy.round(max_lat, 7)

        num_nodes = num_points_lon * num_points_lat
        current_depth = target_depth * -1

        # Update carrier metadata
        self.h5_file.attrs.modify('geographicIdentifier', ofs_metadata.region)
        self.h5_file.attrs.modify('productSpecification', S111Metadata.PRODUCT_SPECIFICATION)
        self.h5_file.attrs.modify('epoch', S111Metadata.EPOCH)
        self.h5_file.attrs.modify('horizontalDatumReference', S111Metadata.HORIZONTAL_DATUM_REFERENCE)
        self.h5_file.attrs.modify('horizontalDatumValue', S111Metadata.HORIZONTAL_DATUM_VALUE)
        self.h5_file.attrs.modify('depthTypeIndex', S111Metadata.DEPTH_TYPE_INDEX)
        self.h5_file.attrs.modify('surfaceCurrentDepth', current_depth)
        self.h5_file.attrs.modify('westBoundLongitude', min_lon)
        self.h5_file.attrs.modify('eastBoundLongitude', max_lon)
        self.h5_file.attrs.modify('southBoundLatitude', min_lat)
        self.h5_file.attrs.modify('northBoundLatitude', max_lat)
        self.h5_file.attrs.modify('metadata', S111Metadata.XML_REFERENCE)

        # Update feature container metadata
        self.feature.attrs.modify('dataCodingFormat', S111Metadata.DATA_CODING_FORMAT)
        self.feature.attrs.modify('interpolationType', S111Metadata.INTERPOLATION_TYPE)
        self.feature.attrs.modify('typeOfCurrentData', S111Metadata.TYPE_OF_CURRENT_DATA)
        self.feature.attrs.modify('commonPointRule', S111Metadata.COMMON_POINT_RULE)
        self.feature.attrs.modify('dimension', S111Metadata.DIMENSION)
        self.feature.attrs.modify('sequenceRule.type', S111Metadata.SEQUENCE_RULE_TYPE)
        self.feature.attrs.modify('methodCurrentsProduct', ofs_metadata.product)
        self.feature.attrs.modify('sequenceRule.scanDirection', S111Metadata.SEQUENCE_RULE_SCAN_DIRECTION)

        # Update feature instance metadata
        if self.feature.attrs['dataCodingFormat'] == S111Metadata.DATA_CODING_FORMAT:
            self.feature_instance.attrs.modify('gridOriginLongitude', min_lon)
            self.feature_instance.attrs.modify('gridOriginLatitude', min_lat)
            self.feature_instance.attrs.modify('gridSpacingLongitudinal', cellsize_x)
            self.feature_instance.attrs.modify('gridSpacingLatitudinal', cellsize_y)
            self.feature_instance.attrs.modify('numPointsLongitudinal', num_points_lon)
            self.feature_instance.attrs.modify('numPointsLatitudinal', num_points_lat)
            self.feature_instance.attrs.modify('startSequence', S111Metadata.START_SEQUENCE)
            self.feature_instance.attrs.modify('westBoundLongitude', min_lon)
            self.feature_instance.attrs.modify('eastBoundLongitude', max_lon)
            self.feature_instance.attrs.modify('southBoundLatitude', min_lat)
            self.feature_instance.attrs.modify('northBoundLatitude', max_lat)

    def add_data(self, time_value, reg_grid_speed, reg_grid_direction, cycletime):
        """Add data to the S111 file.
        
        As data is added, new Groups will be added and relevant attributes will
        be updated.

        Args:
            time_value: `datetime.datetime` instance representing valid time of
                the nowcast/forecast data being added.
            reg_grid_speed: `numpy.ma.masked_array` representing current speed
                data after interpolating to a regular grid and converting from
                u/v components.
            reg_grid_direction: `numpy.ma.masked_array` representing current
                direction data after interpolating to a regular grid and
                converting from u/v components.
            cycletime: `datetime.datetime` representing model cycle time.
        """

        # Add datasets to group_f and feature container once
        if len(self.groupF) == 0:
            # Add group_f compound dataset
            dtype = numpy.dtype([('code', h5py.special_dtype(vlen=str)),
                                 ('name', h5py.special_dtype(vlen=str)),
                                 ('uom.name', h5py.special_dtype(vlen=str)),
                                 ('fillValue', h5py.special_dtype(vlen=str)),
                                 ('dataType', h5py.special_dtype(vlen=str)),
                                 ('lower', h5py.special_dtype(vlen=str)),
                                 ('upper', h5py.special_dtype(vlen=str)),
                                 ('closure', h5py.special_dtype(vlen=str))])

            fdata = numpy.zeros((2,), dtype=dtype)
            fdata['code'][0] = 'surfaceCurrentSpeed'
            fdata['name'][0] = 'Surface current speed'
            fdata['uom.name'][0] = 'knots'
            fdata['fillValue'][0] = str(FILLVALUE)
            fdata['dataType'][0] = H5T_CLASS_T[h5py.h5t.FLOAT]
            fdata['lower'][0] = 0.0
            fdata['upper'][0] = ''
            fdata['closure'][0] = 'geSemiInterval'

            fdata['code'][1] = 'surfaceCurrentDirection'
            fdata['name'][1] = 'Surface current direction'
            fdata['uom.name'][1] = 'arc-degrees'
            fdata['fillValue'][1] = str(FILLVALUE)
            fdata['dataType'][1] = H5T_CLASS_T[h5py.h5t.FLOAT]
            fdata['lower'][1] = 0.0
            fdata['upper'][1] = 360
            fdata['closure'][1] = 'geLtInterval'

            self.groupF_dset = self.groupF.create_dataset('SurfaceCurrent', (2,), dtype=dtype)
            self.groupF_dset[...] = fdata

            # Add group_f feature code dataset
            fc_data = numpy.zeros((1,), dtype=h5py.special_dtype(vlen=str))
            fc_data[0] = 'SurfaceCurrent'
            feature_code = self.groupF.create_dataset('featureCode', (1,), dtype=h5py.special_dtype(vlen=str))
            feature_code[...] = fc_data

            # Add feature container dataset
            axis_names = numpy.zeros((2,), dtype=h5py.special_dtype(vlen=str))
            axis_names[0] = 'longitude'
            axis_names[1] = 'latitude'
            axis_dset = self.feature.create_dataset('axisNames', (2,), dtype=h5py.special_dtype(vlen=str))
            axis_dset[...] = axis_names

        # Create a list of all feature instance objects and groups
        feature_instance_objs = []
        self.feature_instance.visit(feature_instance_objs.append)
        feature_instance_groups = [obj for obj in feature_instance_objs if
                                   isinstance(self.feature_instance[obj], h5py.Group)]

        # Convert time value to string
        time_str = time_value.strftime('%Y%m%dT%H%M%SZ')

        if len(feature_instance_groups) == 0:
            new_group = self.feature_instance.create_group('Group_001')
            print("Creating", "Group_001", "dataset.")

            # Time attributes updated once
            issuance_time = cycletime.strftime('%H%M%SZ')
            issuance_date = cycletime.strftime('%Y%m%d')
            self.h5_file.attrs.modify('issueTime', numpy.string_(issuance_time))
            self.h5_file.attrs.modify('issueDate', numpy.string_(issuance_date))
            self.feature.attrs.modify('numInstances', len(self.feature_instance))
            self.feature_instance.attrs.modify('dateTimeOfFirstRecord', numpy.string_(time_str))
            self.feature_instance.attrs.modify('dateTimeOfLastRecord', numpy.string_(time_str))

        else:
            num_grps = len(feature_instance_groups)
            new_grp = num_grps + 1
            new_group = self.feature_instance.create_group('Group_{:03d}'.format(new_grp))
            print("Creating", "Group_{:03d}".format(new_grp), "dataset.")
            self.feature_instance.attrs.modify('dateTimeOfLastRecord', numpy.string_(time_str))

        # Update attributes from datasets added
        min_speed = numpy.nanmin(reg_grid_speed)
        max_speed = numpy.nanmax(reg_grid_speed)
        min_speed = numpy.round(min_speed, decimals=2)
        max_speed = numpy.round(max_speed, decimals=2)

        speed = reg_grid_speed.filled(FILLVALUE)
        direction = reg_grid_direction.filled(FILLVALUE)
        speed = numpy.round(speed, decimals=2)
        direction = numpy.round(direction, decimals=1)

        # Update attributes from datasets added
        new_group.attrs.create('timePoint', numpy.string_(time_str), None, h5py.special_dtype(vlen=str))
        self.feature_instance.attrs.modify('numberOfTimes', len(self.feature_instance))
        self.feature_instance.attrs.modify('numGRP', len(self.feature_instance))

        # Write data to empty feature instance group
        values_dtype = numpy.dtype([('surfaceCurrentSpeed', numpy.float32),
                                    ('surfaceCurrentDirection', numpy.float32)])

        values = numpy.zeros(speed.shape, dtype=values_dtype)
        values['surfaceCurrentSpeed'] = speed
        values['surfaceCurrentDirection'] = direction
        values_dset = new_group.create_dataset('values', speed.shape, dtype=values_dtype, chunks=True,
                                               compression='gzip', compression_opts=9)
        values_dset[...] = values

        # Update group_f attributes
        self.groupF_dset.attrs.create('chunking', str(values_dset.chunks), dtype=h5py.special_dtype(vlen=str))
        self.feature_instance.attrs.modify('instanceChunking', numpy.string_(str(values_dset.chunks)))

        if len(self.feature_instance) == 2:
            # Time record interval is the same through out the forecast.
            first_time = datetime.datetime.strptime((self.feature_instance.attrs['dateTimeOfFirstRecord']),
                                                    '%Y%m%dT%H%M%SZ')
            last_time = datetime.datetime.strptime((self.feature_instance.attrs['dateTimeOfLastRecord']),
                                                   '%Y%m%dT%H%M%SZ')
            interval = last_time - first_time
            time_interval = interval.total_seconds()
            self.feature_instance.attrs.modify('timeRecordInterval', time_interval)
            self.feature.attrs.modify('minDatasetCurrentSpeed', min_speed)
            self.feature.attrs.modify('maxDatasetCurrentSpeed', max_speed)

        else:
            # Update min and max speed attributes each time data is added 
            prior_min_speed = self.feature.attrs['minDatasetCurrentSpeed']
            prior_max_speed = self.feature.attrs['maxDatasetCurrentSpeed']
            if min_speed < prior_min_speed:
                self.feature.attrs.modify('minDatasetCurrentSpeed', min_speed)
            if max_speed > prior_max_speed:
                self.feature.attrs.modify('maxDatasetCurrentSpeed', max_speed)


class S111Metadata:
    """Contains s111 metadata to pass to s111File.

    PRODUCT_SPECIFICATION: The product specification used to create this dataset.
    EPOCH: Code denoting the epoch of the geodetic datum used by the CRS.
    HORIZONTAL_DATUM_REFERENCE: Reference to the register from which the horizontal datum value is taken.
    HORIZONTAL_DATUM_VALUE: Horizontal Datum of the entire dataset.
    DATA_CODING_FORMAT: 1:Time series at fixed stations, 2:Regularly gridded arrays, 3:Ungeorectified
        gridded arrays, 4:Time series for one moving platform.
    DEPTH_TYPE_INDEX: 1:Layer average, 2:Sea surface, 3:Vertical datum, 4:Sea bottom.
    TYPE_OF_CURRENT_DATA: 1:Historical, 2:Real-time, 3:Astronomical , 4:Hybrid, 5:Hydrodynamic model hindcast,
        6:Hydrodynamic Model forecast.
    METAFEATURES: Name of metafeatures file.
    METADATA: Name of XML metadata file.
    INTERPOLATION_TYPE: Interpolation method recommended for evaluation of the S100_GridCoverage.
    COMMON_POINT_RULE: The procedure used for evaluating geometric objects that overlap or lie fall on boundaries.
    DIMENSION: The dimension of the feature instance.
    SEQUENCE_RULE_TYPE: Method to assign values from the sequence of values to the grid coordinates (e.g. "linear").
    SEQUENCE_RULE_SCAN_DIRECTION: AxisNames, comma-separated (e.g. "latitude,longitude").
    START_SEQUENCE: Starting location of the scan.

    """
    PRODUCT_SPECIFICATION = numpy.string_('INT.IHO.S-111.1.0.0')
    EPOCH = numpy.string_('G1762')
    HORIZONTAL_DATUM_REFERENCE = numpy.string_('EPSG')
    HORIZONTAL_DATUM_VALUE = 4326
    DATA_CODING_FORMAT = 2
    DEPTH_TYPE_INDEX = 2
    TYPE_OF_CURRENT_DATA = 6
    METAFEATURES = None
    METADATA = None
    INTERPOLATION_TYPE = 10
    COMMON_POINT_RULE = 4
    DIMENSION = 2
    SEQUENCE_RULE_TYPE = 1
    SEQUENCE_RULE_SCAN_DIRECTION = numpy.string_('longitude,latitude')
    START_SEQUENCE = numpy.string_('0,0')
    XML_REFERENCE = numpy.string_('METADATA.XML')

    def __init__(self, region, product):
        self.region = numpy.string_(region)
        self.product = numpy.string_(product)


def convert_to_s111(model_index_file, model_files, s111_path_prefix, cycletime, ofs_model, ofs_metadata, target_depth=None):
    """Convert NetCDF model to regular grid in S111 format.

    If the supplied model index NetCDF contains information identifying
    subgrids, one S111 file will be generated for each subgrid. Otherwise, a
    single S111 file will be created for the entire domain.

    Note: Only a single time per model file is currently supported. If a model
    NetCDF includes more than one time/forecast, only the first will be
    extracted.

    Args:
        model_index_file: Path to model index NetCDF file containing
            pre-calculated grid and mask information.
        model_files: List of paths to one or more NetCDF model files.
            Files should be provided in ascending chronological order,
            as this order will be maintained when appending subsequent
            nowcasts/forecasts to each other in individual S111 files.
        s111_path_prefix: Path prefix for desired output location for generated
            S-111 files. If specified path is a directory, file(s) will be
            output to specified directory with autogenerated names. Otherwise,
            generated file(s) will be placed at specified file path, but with a
            filename suffix appended based on the properties of the target
            output grid and the model file.
        cycletime: `datetime.datetime` instance representing target cycle time
            of model forecast(s) being processed.
        ofs_model: Model identifier (e.g. "cbofs").
        ofs_metadata: `S111Metadata` instance describing metadata for geographic
            identifier and description of current meter type, forecast method,
            or model.
        target_depth: The water current at a specified target depth below
            the sea surface in meters, default target depth is 4.5 meters,
            target interpolation depth must be greater or equal to 0.

    Returns:
        List of paths to HDF5 files created.
    """

    # Path format/prefix for output S111 files. Forecast initialization (reference).
    if os.path.isdir(s111_path_prefix):
        if not s111_path_prefix.endswith("/"):
            s111_path_prefix += "/"
        file_issuance = cycletime.strftime('%Y%m%dT%HZ')
        s111_path_prefix += ("S111US_{}_{}_TYP2".format(file_issuance, str.upper(ofs_model)))

    if target_depth is None:
        target_depth = DEFAULT_TARGET_DEPTH

    s111_file_paths = []

    try:
        model_index_file.open()
        if model_index_file.dim_subgrid is not None and model_index_file.var_subgrid_id is not None:
            # Output to subgrids
            with contextlib.ExitStack() as stack:
                s111_files = []
                for i in range(model_index_file.dim_subgrid.size):
                    if model_index_file.var_subgrid_name is not None:
                        s111_file = S111File("{}_{}.h5".format(s111_path_prefix, model_index_file.var_subgrid_name[i]), clobber=True)
                    else:
                        s111_file = S111File("{}_FID_{}.h5".format(s111_path_prefix, model_index_file.var_subgrid_id[i]), clobber=True)

                    s111_file_paths.append(s111_file.path)
                    stack.enter_context(s111_file)
                    s111_file.update_attributes(model_index_file, ofs_metadata, target_depth, i)
                    s111_files.append(s111_file)

                for model_file in model_files:
                    try:
                        model_file.open()

                        # Call model method and convert and interpolate u/v to regular grid
                        reg_grid_u, reg_grid_v = model_file.uv_to_regular_grid(model_index_file, target_depth)

                        reg_grid_u = ma.masked_array(reg_grid_u, model_index_file.var_mask.mask)
                        reg_grid_v = ma.masked_array(reg_grid_v,  model_index_file.var_mask.mask)

                        # Convert currents at regular grid points from u/v to speed
                        # and direction
                        direction, speed = model.uv_to_speed_direction(reg_grid_u, reg_grid_v)

                        direction = ma.masked_array(direction, model_index_file.var_mask.mask)
                        speed = ma.masked_array(speed, model_index_file.var_mask.mask)

                        for subgrid_index, s111_file in enumerate(s111_files):
                            if os.path.isfile(s111_file.path):
                                x_min = model_index_file.var_subgrid_x_min[subgrid_index]
                                x_max = model_index_file.var_subgrid_x_max[subgrid_index]
                                y_min = model_index_file.var_subgrid_y_min[subgrid_index]
                                y_max = model_index_file.var_subgrid_y_max[subgrid_index]
                                subgrid_speed = speed[y_min:y_max + 1, x_min:x_max + 1]
                                subgrid_direction = direction[y_min:y_max + 1, x_min:x_max + 1]
                                if ma.count(subgrid_speed) >= 20:
                                    s111_file.add_data(model_file.time_val, subgrid_speed, subgrid_direction, cycletime)
                                else:
                                    s111_file.close()
                                    os.remove("{}".format(s111_file.path))
                                    s111_file_paths.remove(s111_file.path)
                    finally:
                        model_file.close()
        else:
            # Output to default grid (no subgrids)
            with S111File("{}.h5".format(s111_path_prefix), clobber=True) as s111_file:
                s111_file_paths.append(s111_file.path)
                s111_file.update_attributes(model_index_file, ofs_metadata, target_depth)
                for model_file in model_files:
                    try:
                        model_file.open()

                        # Call model method and convert and interpolate u/v to regular grid
                        reg_grid_u, reg_grid_v = model_file.uv_to_regular_grid(model_index_file, target_depth)

                        reg_grid_u = ma.masked_array(reg_grid_u, model_index_file.var_mask.mask)
                        reg_grid_v = ma.masked_array(reg_grid_v,  model_index_file.var_mask.mask)

                        # Convert currents at regular grid points from u/v to speed
                        # and direction
                        direction, speed = model.uv_to_speed_direction(reg_grid_u, reg_grid_v)

                        direction = ma.masked_array(direction, model_index_file.var_mask.mask)
                        speed = ma.masked_array(speed, model_index_file.var_mask.mask)

                        s111_file.add_data(model_file.time_val, speed, direction, cycletime)
                    finally:
                        model_file.close()
    finally:
        model_index_file.close()

    return s111_file_paths

