"""Utilities for creation and modification of S-111 compliant HDF-5 files.

S-111 is an IHO standard outlining formats for storing and sending surface
water current data and metadata.
"""
import contextlib
import datetime
import os
import numpy
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings('ignore', category=FutureWarning)
    import h5py

from thyme.model import model

# Default fill value for NetCDF variables
FILLVALUE = -9999.0

# Default depth in meters
DEFAULT_TARGET_DEPTH = 4.5

# Dict lookup for HDF5 data type class names
# See https://github.com/h5py/h5py/blob/master/h5py/api_types_hdf5.pxd#L509
H5T_CLASS_T = {
    h5py.h5t.NO_CLASS: 'H5T_NO_CLASS',
    h5py.h5t.INTEGER: 'H5T_INTEGER',
    h5py.h5t.FLOAT: 'H5T_FLOAT',
    h5py.h5t.TIME: 'H5T_TIME',
    h5py.h5t.STRING: 'H5T_STRING',
    h5py.h5t.BITFIELD: 'H5T_BITFIELD',
    h5py.h5t.OPAQUE: 'H5T_OPAQUE',
    h5py.h5t.COMPOUND: 'H5T_COMPOUND',
    h5py.h5t.REFERENCE: 'H5T_REFERENCE',
    h5py.h5t.ENUM: 'H5T_ENUM',
    h5py.h5t.VLEN: 'H5T_VLEN',
    h5py.h5t.ARRAY: 'H5T_ARRAY'
}


class S111File:
    """Create and manage S-111 files.

    This class implements the context manager pattern and should thus be used
    similar to the following:

        with S111File('myfile.h5') as f:
            ...

    If the file does not yet exist, it will be created. An existing file will
    be opened in append mode unless ``clobber=True``, in which case it will be
    overwritten.

    Attributes:
        path: Path (relative or absolute) to the .h5 file, including filename.
        filename: Name of the .h5 file.
        input_metadata: ``S111Metadata`` instance describing metadata for geographic
            identifier and description of current meter type, forecast method,
            or model.
        model_index: (Optional, default None) ``ModelIndexFile`` instance
            representing model index file.
        subgrid_index: (Optional, default None) Index of subgrid, if any,
            that this S111File represents. Corresponds with index into
            subgrid dimension of model index file.
        h5_file: Handle to the underlying ``h5py.File`` instance.
        feature: Handle to the underlying ``h5py.Group`` instance for the
            Feature.
        feature_instance: Handle to the underlying ``h5py.Group`` instance for
            the Feature Instance.
        groupF: Handle to the underlying ``h5py.Group`` instance for Group_F
            metadata.
        groupF_dset: Handle to the underlying ``h5py.Dataset`` for Group_F 
            metadata.
    """

    def __init__(self, path, input_metadata, data_coding_format, model_index=None, subgrid_index=None, clobber=False):
        """Initializes S111File object and opens h5 file at specified path.

        If ``path`` has an extension other than '.h5', it is replaced with
        '.h5'

        Args:
            path: Path of target hdf5 file. Must end in '.h5', otherwise its
                extension will be replaced with '.h5'.
            input_metadata: ``S111Metadata`` instance describing metadata for
                geographic identifier and description of current meter type,
                forecast method, or model.
            data_coding_format: 1:Time series at fixed stations, 2:Regularly gridded arrays,
                3:Ungeorectified gridded arrays, 4:Time series for one moving platform.
            model_index: (Optional, default None) ``ModelIndexFile`` instance
                representing model index file.
            subgrid_index: (Optional, default None) Index of subgrid, if any,
                that this S111File represents. Corresponds with index into
                subgrid dimension of model index file.
            clobber: (Optional, default False) If True, existing h5 file at
                specified path, if any, will be deleted and the new file will
                be opened in write mode.
        """
        prefix, extension = os.path.splitext(path)
        self.path = prefix + '.h5'
        self.filename = os.path.basename(self.path)
        self.model_index = model_index
        self.input_metadata = input_metadata
        self.data_coding_format = data_coding_format
        self.subgrid_index = subgrid_index

        if not os.path.exists(self.path) or clobber:
            # File doesn't exist, open in create (write) mode and add metadata
            self.h5_file = h5py.File(self.path, 'w')

            # Create s111 structure
            self.groupF = self.h5_file.create_group('Group_F')
            self.feature = self.h5_file.create_group('SurfaceCurrent')
            self.feature_instance = self.feature.create_group('SurfaceCurrent.01')
            self.feature_instance_groups = None
            self.groupF_dset = None

            # Add feature content
            self.add_feature_codes_content()
            self.add_feature_type_content()
            self.add_feature_instance_content()

            # Add metadata
            self.add_metadata()

        else:
            # File already exists, open in append mode
            self.h5_file = h5py.File(self.path, 'r+')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5_file.close()

    def close(self):
        self.h5_file.flush()
        self.h5_file.close()

    def add_feature_codes_content(self):
        """Add feature codes.

        This group specifies the S-100 feature to which the data applies.
        The group has no attributes and consists of two components:
        Feature name, a dataset with the name of the S-100 feature contained
        in the data product and a reference to the name.
        """
        # Add a feature name compound dataset
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

        # Add a feature code dataset
        fc_data = numpy.zeros((1,), dtype=h5py.special_dtype(vlen=str))
        fc_data[0] = 'SurfaceCurrent'
        feature_code = self.groupF.create_dataset('featureCode', (1,), dtype=h5py.special_dtype(vlen=str))
        feature_code[...] = fc_data

    def add_feature_type_content(self):
        """Add feature type content to the S111 file."""
        # Add horizontal and vertical axis names in feature type
        axis_names = numpy.zeros((2,), dtype=h5py.special_dtype(vlen=str))
        axis_names[0] = 'longitude'
        axis_names[1] = 'latitude'
        axis_dset = self.feature.create_dataset('axisNames', (2,), dtype=h5py.special_dtype(vlen=str))
        axis_dset[...] = axis_names

    def add_feature_instance_content(self):
        """Add feature instance content to the S111 file."""
        # Add feature instance uncertainty compound dataset
        u_dtype = numpy.dtype([('name', h5py.special_dtype(vlen=str)), ('value', numpy.float32)])
        u_data = numpy.zeros((2,), u_dtype)
        u_data['name'][0] = 'surfaceCurrentSpeedUncertainty'
        u_data['name'][1] = 'surfaceCurrentDirectionUncertainty'
        u_data['value'][0] = -1.0
        u_data['value'][1] = -1.0
        uncertainty_data = self.feature_instance.create_dataset('uncertainty', (2,), u_dtype)
        uncertainty_data[...] = u_data

    def add_metadata(self):
        """Add metadata.

        Based on grid properties, s111 metadata, and input metadata.
        """

        metadata_xml_reference = numpy.string_('MD_{}.XML'.format(os.path.splitext(self.filename)[0]))

        # Add carrier metadata
        self.h5_file.attrs.create('depthTypeIndex', self.input_metadata.DEPTH_TYPE_INDEX, dtype=numpy.int32)
        self.h5_file.attrs.create('metadata', metadata_xml_reference, dtype=h5py.special_dtype(vlen=str))
        self.h5_file.attrs.create('horizontalDatumValue', self.input_metadata.HORIZONTAL_DATUM_VALUE, dtype=numpy.int32)
        self.h5_file.attrs.create('geographicIdentifier', self.input_metadata.region, dtype=h5py.special_dtype(vlen=str))
        self.h5_file.attrs.create('productSpecification', self.input_metadata.PRODUCT_SPECIFICATION, dtype=h5py.special_dtype(vlen=str))
        self.h5_file.attrs.create('horizontalDatumReference', self.input_metadata.HORIZONTAL_DATUM_REFERENCE, dtype=h5py.special_dtype(vlen=str))

        # Add feature container metadata
        self.feature.attrs.create('methodCurrentsProduct', self.input_metadata.product, dtype=h5py.special_dtype(vlen=str))
        self.feature.attrs.create('dataCodingFormat', self.data_coding_format, dtype=numpy.int32)
        self.feature.attrs.create('commonPointRule', self.input_metadata.COMMON_POINT_RULE, dtype=numpy.int32)
        self.feature.attrs.create('typeOfCurrentData', self.input_metadata.current_datatype, dtype=numpy.int32)
        self.feature.attrs.create('horizontalPositionUncertainty', -1.0, dtype=numpy.float32)
        self.feature.attrs.create('verticalUncertainty', -1.0, dtype=numpy.float32)
        self.feature.attrs.create('timeUncertainty', -1.0, dtype=numpy.float32)

        if self.data_coding_format == 2:

            # Width between first two cells, grid spacing is uniform
            cellsize_x = self.model_index.var_x[1] - self.model_index.var_x[0]
            cellsize_y = self.model_index.var_y[1] - self.model_index.var_y[0]

            if self.subgrid_index is not None:
                if self.subgrid_index < 0 or self.subgrid_index >= self.model_index.dim_subgrid.size:
                    raise Exception(
                        'Subgrid index [{}] out of model index subgrid dimension range [0-{}]'.format(
                            self.subgrid_index,
                            self.model_index.dim_subgrid.size - 1))

                num_points_lon = 1 + self.model_index.var_subgrid_x_max[self.subgrid_index] - self.model_index.var_subgrid_x_min[self.subgrid_index]
                num_points_lat = 1 + self.model_index.var_subgrid_y_max[self.subgrid_index] - self.model_index.var_subgrid_y_min[self.subgrid_index]
                min_lon = self.model_index.var_x[self.model_index.var_subgrid_x_min[self.subgrid_index]]
                min_lat = self.model_index.var_y[self.model_index.var_subgrid_y_min[self.subgrid_index]]
                max_lon = self.model_index.var_x[self.model_index.var_subgrid_x_max[self.subgrid_index]]
                max_lat = self.model_index.var_y[self.model_index.var_subgrid_y_max[self.subgrid_index]]
            else:
                num_points_lon = self.model_index.dim_x.size
                num_points_lat = self.model_index.dim_y.size
                min_lon = numpy.nanmin(self.model_index.var_x)
                min_lat = numpy.nanmin(self.model_index.var_y)
                max_lon = numpy.nanmax(self.model_index.var_x)
                max_lat = numpy.nanmax(self.model_index.var_y)

            # X/Y coordinates are located at the center of each grid cell
            min_lon = numpy.round(min_lon, 7)
            max_lon = numpy.round(max_lon, 7)
            min_lat = numpy.round(min_lat, 7)
            max_lat = numpy.round(max_lat, 7)

            # Add carrier metadata
            self.h5_file.attrs.create('westBoundLongitude', min_lon, dtype=numpy.float32)
            self.h5_file.attrs.create('eastBoundLongitude', max_lon, dtype=numpy.float32)
            self.h5_file.attrs.create('southBoundLatitude', min_lat, dtype=numpy.float32)
            self.h5_file.attrs.create('northBoundLatitude', max_lat, dtype=numpy.float32)

            # Add feature container metadata
            self.feature.attrs.create('sequencingRule.scanDirection', self.input_metadata.SEQUENCING_RULE_SCAN_DIRECTION, dtype=h5py.special_dtype(vlen=str))
            self.feature.attrs.create('sequencingRule.type', self.input_metadata.SEQUENCING_RULE_TYPE, dtype=numpy.int32)
            self.feature.attrs.create('dimension', 2, dtype=numpy.int32)

            # Add feature instance metadata
            self.feature_instance.attrs.create('startSequence', self.input_metadata.START_SEQUENCE, dtype=h5py.special_dtype(vlen=str))
            self.feature_instance.attrs.create('gridOriginLongitude', min_lon, dtype=numpy.float32)
            self.feature_instance.attrs.create('gridOriginLatitude', min_lat, dtype=numpy.float32)
            self.feature_instance.attrs.create('gridSpacingLongitudinal', cellsize_x, dtype=numpy.float32)
            self.feature_instance.attrs.create('gridSpacingLatitudinal', cellsize_y, dtype=numpy.float32)
            self.feature_instance.attrs.create('numPointsLongitudinal', num_points_lon, dtype=numpy.int32)
            self.feature_instance.attrs.create('numPointsLatitudinal', num_points_lat, dtype=numpy.int32)

            # Add feature instance metadata
            self.feature_instance.attrs.create('westBoundLongitude', min_lon, dtype=numpy.float32)
            self.feature_instance.attrs.create('eastBoundLongitude', max_lon, dtype=numpy.float32)
            self.feature_instance.attrs.create('southBoundLatitude', min_lat, dtype=numpy.float32)
            self.feature_instance.attrs.create('northBoundLatitude', max_lat, dtype=numpy.float32)

    def add_feature_instance_group_data(self, datetime_value, speed, direction, cycletime, target_depth):
        """Add data to the S111 file.
        
        As data is added, new groups will be created and relevant attributes updated.

        Args:
            datetime_value: ``datetime.datetime`` instance representing valid time of
                the nowcast/forecast data being added.
            speed: ``numpy.ma.masked_array`` representing current speed.
            direction: ``numpy.ma.masked_array`` representing current.
            cycletime: ``datetime.datetime`` representing model cycle time.
            target_depth: Target depth below the sea surface, in meters, at which
                surface currents are valid. Default target depth is 4.5 meters.
                Must be greater than or equal to 0. For areas shallower than the
                target depth, half the water column height is used instead.
        """
        # Create a list of all feature instance objects and groups
        feature_instance_objs = []
        self.feature_instance.visit(feature_instance_objs.append)
        self.feature_instance_groups = [obj for obj in feature_instance_objs if
                                        isinstance(self.feature_instance[obj], h5py.Group)]

        # Convert time value to string
        time_str = datetime_value.strftime('%Y%m%dT%H%M%SZ')

        # Format speed and get attributes
        min_speed = numpy.nanmin(speed)
        max_speed = numpy.nanmax(speed)
        min_speed = numpy.round(min_speed, decimals=2)
        max_speed = numpy.round(max_speed, decimals=2)

        # Create feature instance groups
        if len(self.feature_instance_groups) == 0:
            self.feature.attrs.create('numInstances', len(self.feature_instance), dtype=numpy.int32)
            feature_group = self.feature_instance.create_group('Group_001')

            # Time attributes updated once
            issuance_time = cycletime.strftime('%H%M%SZ')
            issuance_date = cycletime.strftime('%Y%m%d')
            self.h5_file.attrs.create('issueTime', numpy.string_(issuance_time), dtype=h5py.special_dtype(vlen=str))
            self.h5_file.attrs.create('issueDate', numpy.string_(issuance_date), dtype=h5py.special_dtype(vlen=str))
            self.feature_instance.attrs.create('dateTimeOfFirstRecord', numpy.string_(time_str), dtype=h5py.special_dtype(vlen=str))
            self.feature_instance.attrs.create('dateTimeOfLastRecord', numpy.string_(time_str), dtype=h5py.special_dtype(vlen=str))
            self.feature_instance.attrs.create('timeRecordInterval', 0, dtype=numpy.int32)

            # Add initial speed attributes
            self.feature.attrs.create('minDatasetCurrentSpeed', min_speed, dtype=numpy.float32)
            self.feature.attrs.create('maxDatasetCurrentSpeed', max_speed, dtype=numpy.float32)

            feature_instance_date = int(cycletime.strftime('%Y%m'))

            if feature_instance_date < 199406:
                epoch = 'TRANSIT'
            elif feature_instance_date < 199706:
                epoch = 'G730'
            elif feature_instance_date < 200201:
                epoch = 'G873'
            elif feature_instance_date < 201202:
                epoch = 'G1150'
            elif feature_instance_date < 201310:
                epoch = 'G1674'
            elif feature_instance_date >= 201310:
                epoch = 'G1762'

            self.h5_file.attrs.create('epoch', epoch, dtype=h5py.special_dtype(vlen=str))

        else:
            # num_groups is 0-based
            num_groups = len(self.feature_instance_groups)
            add_group = num_groups + 1
            feature_group = self.feature_instance.create_group('Group_{:03d}'.format(add_group))
            self.feature_instance.attrs.modify('dateTimeOfLastRecord', numpy.string_(time_str))

            # Update speed attributes each time data is added
            prior_min_speed = self.feature.attrs['minDatasetCurrentSpeed']
            prior_max_speed = self.feature.attrs['maxDatasetCurrentSpeed']
            if min_speed < prior_min_speed:
                self.feature.attrs.modify('minDatasetCurrentSpeed', min_speed)
            if max_speed > prior_max_speed:
                self.feature.attrs.modify('maxDatasetCurrentSpeed', max_speed)

        # Add time string to feature instance group compound dataset
        feature_group.attrs.create('timePoint', numpy.string_(time_str), None, h5py.special_dtype(vlen=str))

        # Add speed and direction data to feature instance group compound dataset
        values_dtype = numpy.dtype([('surfaceCurrentSpeed', numpy.float32), ('surfaceCurrentDirection', numpy.float32)])

        # Check if numpy array is masked
        # If numpy array is masked remove nan values
        if numpy.ma.is_masked(speed):
            speed = speed.filled(FILLVALUE)
            direction = direction.filled(FILLVALUE)

        # Format speed/direction
        speed = numpy.round(speed, decimals=2)
        direction = numpy.round(direction, decimals=1)

        # Add speed/direction data
        values = numpy.zeros(speed.shape, dtype=values_dtype)
        values['surfaceCurrentSpeed'] = speed
        values['surfaceCurrentDirection'] = direction
        values_dset = feature_group.create_dataset('values', speed.shape, dtype=values_dtype, chunks=True, compression='gzip', compression_opts=9)
        values_dset[...] = values

        # Update depth attribute
        current_depth = (-abs(target_depth)) + 0

        self.h5_file.attrs.create('surfaceCurrentDepth', current_depth, dtype=numpy.float32)

        # Update chunking attributes
        chunking_str = ','.join(str(x) for x in values_dset.chunks)

        self.groupF_dset.attrs.create('chunking', chunking_str, dtype=h5py.special_dtype(vlen=str))
        self.feature_instance.attrs.create('instanceChunking', numpy.string_(chunking_str))

    def add_positioning(self, longitude, latitude):
        """Add positioning group and data to the S111 file.

        Args:
            longitude: ``numpy.ma.masked_array`` representing longitude.
            latitude: ``numpy.ma.masked_array`` representing latitude.
        """

        # Add longitude/latitude positioning
        geometry_dtype = numpy.dtype([('longitude', numpy.float32), ('latitude', numpy.float32)])
        dim = len(longitude)

        # Create positioning group
        feature_positioning = self.feature_instance.create_group('Positioning')

        # Create lon/lat compound dataset
        geometry = numpy.zeros((dim,), dtype=geometry_dtype)
        geometry['longitude'] = longitude
        geometry['latitude'] = latitude
        geometry_dset = feature_positioning.create_dataset('geometryValues', (dim,), dtype=geometry_dtype,
                                                           chunks=True, compression='gzip', compression_opts=9)
        geometry_dset[...] = geometry

        # X/Y coordinates are located at the center of each grid cell
        min_lon = numpy.nanmin(longitude)
        min_lat = numpy.nanmin(latitude)
        max_lon = numpy.nanmax(longitude)
        max_lat = numpy.nanmax(latitude)

        min_lon = numpy.round(min_lon, 7)
        max_lon = numpy.round(max_lon, 7)
        min_lat = numpy.round(min_lat, 7)
        max_lat = numpy.round(max_lat, 7)

        # Update carrier metadata
        self.h5_file.attrs.create('westBoundLongitude', min_lon, dtype=numpy.float32)
        self.h5_file.attrs.create('eastBoundLongitude', max_lon, dtype=numpy.float32)
        self.h5_file.attrs.create('southBoundLatitude', min_lat, dtype=numpy.float32)
        self.h5_file.attrs.create('northBoundLatitude', max_lat, dtype=numpy.float32)
        # Update feature container metadata
        self.feature.attrs.create('dimension', 1, dtype=numpy.int32)
        # Update feature instance metadata
        self.feature_instance.attrs.create('westBoundLongitude', min_lon, dtype=numpy.float32)
        self.feature_instance.attrs.create('eastBoundLongitude', max_lon, dtype=numpy.float32)
        self.feature_instance.attrs.create('southBoundLatitude', min_lat, dtype=numpy.float32)
        self.feature_instance.attrs.create('northBoundLatitude', max_lat, dtype=numpy.float32)

    def add_model_metadata(self):
        """Model specific metadata"""

        # Update feature container metadata
        self.feature.attrs.create('interpolationType', self.input_metadata.INTERPOLATION_TYPE, dtype=numpy.int32)

        # Update attributes after all the value groups have been added, 0-based
        num_feature_instance_groups = len(self.feature_instance_groups) + 1
        self.feature_instance.attrs.create('numGRP', num_feature_instance_groups, dtype=numpy.int32)
        self.feature_instance.attrs.create('numberOfTimes', num_feature_instance_groups, dtype=numpy.int32)

        if self.h5_file.__contains__('/SurfaceCurrent/SurfaceCurrent.01/Group_003'):
            first_time = datetime.datetime.strptime(
                (self.h5_file['/SurfaceCurrent/SurfaceCurrent.01/Group_001'].attrs['timePoint']), '%Y%m%dT%H%M%SZ')
            second_time = datetime.datetime.strptime(
                (self.h5_file['/SurfaceCurrent/SurfaceCurrent.01/Group_002'].attrs['timePoint']), '%Y%m%dT%H%M%SZ')

            time_interval_secs = (second_time - first_time).total_seconds()
            self.feature_instance.attrs.modify('timeRecordInterval', time_interval_secs)

        if self.data_coding_format == 3:
            nodes = self.h5_file['/SurfaceCurrent/SurfaceCurrent.01/Group_001/values']
            num_nodes = nodes.shape[0]
            self.feature_instance.attrs.create('numberOfNodes', num_nodes, dtype=numpy.int32)

    def add_time_series_metadata(self, datetime_values):
        """Time series specific metadata

        Args: datetime_values: List of datetime objects
        """

        num_feature_instance_groups = len(self.feature_instance_groups) + 1
        self.feature_instance.attrs.create('numGRP', num_feature_instance_groups, dtype=numpy.int32)

        last_time_str = datetime_values[-1].strftime('%Y%m%dT%H%M%SZ')

        # Overwrite last date time record
        self.feature_instance.attrs.modify('dateTimeOfLastRecord', numpy.string_(last_time_str))

        interval = datetime_values[1] - datetime_values[0]
        time_interval = interval.total_seconds()
        num_times = len(datetime_values)

        self.feature_instance.attrs.create('timeRecordInterval', time_interval, dtype=numpy.int32)
        self.feature_instance.attrs.create('numberOfTimes', num_times, dtype=numpy.int32)
        self.feature_instance.attrs.create('numberOfStations', num_feature_instance_groups, dtype=numpy.int32)


class S111Metadata:
    """Contains s111 metadata to pass to S111File.

    PRODUCT_SPECIFICATION: The product specification used to create this dataset.
    HORIZONTAL_DATUM_REFERENCE: Reference to the register from which the horizontal datum value is taken.
    HORIZONTAL_DATUM_VALUE: Horizontal Datum of the entire dataset.
    DEPTH_TYPE_INDEX: 1:Layer average, 2:Sea surface, 3:Vertical datum, 4:Sea bottom.
    INTERPOLATION_TYPE: Interpolation method recommended for evaluation of the S100_GridCoverage.
    COMMON_POINT_RULE: The procedure used for evaluating geometric objects that overlap or lie fall on boundaries.
    DIMENSION: The dimension of the feature instance.
    SEQUENCING_RULE_TYPE: Method to assign values from the sequence of values to the grid coordinates (e.g. "linear").
    SEQUENCING_RULE_SCAN_DIRECTION: AxisNames, comma-separated (e.g. "longitude,latitude").
    START_SEQUENCE: Starting location of the scan.

    """
    PRODUCT_SPECIFICATION = numpy.string_('INT.IHO.S-111.1.0.0')
    HORIZONTAL_DATUM_REFERENCE = numpy.string_('EPSG')
    HORIZONTAL_DATUM_VALUE = 4326
    DEPTH_TYPE_INDEX = 2
    INTERPOLATION_TYPE = 10
    COMMON_POINT_RULE = 3
    SEQUENCING_RULE_TYPE = 1
    SEQUENCING_RULE_SCAN_DIRECTION = numpy.string_('longitude,latitude')
    START_SEQUENCE = numpy.string_('0,0')

    def __init__(self, region, product, current_datatype, producer_code, station_id=None, model_system=None):
        """Initializes S111Metadata object.

        Args:
            region: Geographic identifier.
            product: Description of current meter type, forecast method or
                model, etc.
            current_datatype: 1: Historical observation (O)
                              2: Real-time observation (R)
                              3: Astronomical prediction (A)
                              4: Analysis or hybrid method (Y)
                              5: Hydrodynamic model hindcast (M)
                              6: Hydrodynamic model forecast (F)
            producer_code: Two-character hydrographic office producer code.
            station_id: (Optional, default None) Station identifier.
            model_system:(Optional, default None) Ocean model system
                identifier (e.g. "cbofs").
        """
        self.region = region
        self.product = product
        self.current_datatype = current_datatype
        self.producer_code = producer_code
        self.station_id = station_id
        self.model_system = model_system


class S111TimeSeries:
    """Contains prediction time series data to pass to S111File. """

    def __init__(self, longitude, latitude, speed, direction, datetime_values):
        """Initializes S111TimeSeries object.

        Args:
            latitude: 1d 'numpy.ndarray' containing latitudes.
            longitude: 1d 'numpy.ndarray' containing longitudes.
            speed: 1d 'numpy.ndarray' containing speed.
            direction: 1d 'numpy.ndarray' containing direction.
            datetime_values: List containing datetime objects.
        """
        self.longitude = longitude
        self.latitude = latitude
        self.speed = speed
        self.direction = direction
        self.datetime_values = datetime_values


def model_to_s111(model_index_file, model_files, s111_path_prefix, cycletime, input_metadata, data_coding_format, target_depth):
    """Convert NetCDF hydrodynamic model to S111 format.

    If the supplied model index NetCDF contains information identifying
    subgrids, one S111 file will be generated for each subgrid. Otherwise,
    a single S111 file will be created for the entire domain.

    Args:
        model_index_file:  Instance of ``ModelIndexFile`` (or a subclass)
            containing pre-calculated grid and mask information.
        model_files: List of ``ModelFile`` (or subclasses thereof) instances
            identifying NetCDF model files to be converted. Files should be
            provided in ascending chronological order, as this order will be
            maintained when appending subsequent nowcasts/forecasts to each
            other in individual S111 files.
        s111_path_prefix: Path prefix for desired output location for generated
            S-111 files. If specified path is a directory, file(s) will be
            output to specified directory with autogenerated names. Otherwise,
            generated file(s) will be placed at specified file path, but with a
            filename suffix appended based on the properties of the target
            output grid and the model file.
        cycletime: ``datetime.datetime`` instance representing target cycle
            time of model forecast(s) being processed.
        input_metadata: ``S111Metadata`` instance describing metadata for
            geographic identifier and description of current meter type,
            forecast method, or model identifier.
        data_coding_format: 1:Time series at fixed stations, 2:Regularly gridded arrays,
            3:Ungeorectified gridded arrays, 4:Time series for one moving platform.
        target_depth: The water current at a specified target depth below the sea
            surface in meters.

    Returns:
        List of paths to HDF5 files created.
    """
    # Path format/prefix for output S111 files. Forecast initialization (reference).
    if os.path.isdir(s111_path_prefix):
        if not s111_path_prefix.endswith('/'):
            s111_path_prefix += '/'
        file_issuance = cycletime.strftime('%Y%m%dT%HZ')
        s111_path_prefix += (
            'S111{}_{}_{}_TYP{}'.format(input_metadata.producer_code, file_issuance,
                                        str.upper(input_metadata.model_system), data_coding_format))

    # model_to_s111 requires a target.
    if target_depth is None:
        target_depth = DEFAULT_TARGET_DEPTH

    s111_file_paths = []

    if data_coding_format == 2:

        try:
            model_index_file.open()
            if model_index_file.dim_subgrid is not None and model_index_file.var_subgrid_id is not None:
                # Output to subgrids
                stack = contextlib.ExitStack()
                s111_files = []
                for i in range(model_index_file.dim_subgrid.size):
                    if model_index_file.var_subgrid_name is not None:
                        filename = '{}_{}.h5'.format(s111_path_prefix,
                                                     model_index_file.var_subgrid_name[i])
                    else:
                        filename = '{}_FID_{}.h5'.format(s111_path_prefix,
                                                         model_index_file.var_subgrid_id[i])

                    s111_file = S111File(filename, input_metadata, data_coding_format,
                                         model_index_file, subgrid_index=i, clobber=True)

                    s111_file_paths.append(s111_file.path)
                    stack.enter_context(s111_file)
                    s111_files.append(s111_file)

            else:
                # Output entire domain
                s111_file = S111File('{}.h5'.format(s111_path_prefix), input_metadata, data_coding_format,
                                     model_index_file, clobber=True)
                s111_file_paths.append(s111_file.path)

            for model_file in model_files:
                try:
                    model_file.open()
                    for time_index in range(len(model_file.datetime_values)):
                        # Call model method and convert and interpolate u/v to regular grid
                        # The water current at a specified target depth below the sea surface in meters the default
                        # target depth is 4.5 meters, target interpolation depth must be greater or equal to 0.

                        reg_grid_u, reg_grid_v = model_file.uv_to_regular_grid(model_index_file, time_index, target_depth)

                        reg_grid_u = numpy.ma.masked_array(reg_grid_u, model_index_file.var_mask.mask)
                        reg_grid_v = numpy.ma.masked_array(reg_grid_v, model_index_file.var_mask.mask)

                        # Convert currents at regular grid points from u/v to speed/direction
                        speed, direction = model.regular_uv_to_speed_direction(reg_grid_u, reg_grid_v)

                        # Apply mask
                        direction = numpy.ma.masked_array(direction, model_index_file.var_mask.mask)
                        speed = numpy.ma.masked_array(speed, model_index_file.var_mask.mask)

                        # If any valid data points fall outside of the scipy griddata convex hull
                        # nan values will be used, if nan values are present
                        # add nan values to the original mask
                        if numpy.isnan(speed).any():

                            nan_mask_speed = numpy.ma.masked_invalid(speed)
                            nan_mask_direction = numpy.ma.masked_invalid(direction)
                            speed_mask = numpy.ma.mask_or(model_index_file.var_mask.mask, nan_mask_speed.mask)
                            direction_mask = numpy.ma.mask_or(model_index_file.var_mask.mask, nan_mask_direction.mask)

                            speed = numpy.ma.masked_array(speed, speed_mask)
                            direction = numpy.ma.masked_array(direction, direction_mask)

                        if model_index_file.dim_subgrid is not None and model_index_file.var_subgrid_id is not None:
                            # Output to subgrids
                                for subgrid_index, s111_file in enumerate(s111_files):
                                    if os.path.isfile(s111_file.path):
                                        x_min = model_index_file.var_subgrid_x_min[subgrid_index]
                                        x_max = model_index_file.var_subgrid_x_max[subgrid_index]
                                        y_min = model_index_file.var_subgrid_y_min[subgrid_index]
                                        y_max = model_index_file.var_subgrid_y_max[subgrid_index]
                                        subgrid_speed = speed[y_min:y_max + 1, x_min:x_max + 1]
                                        subgrid_direction = direction[y_min:y_max + 1, x_min:x_max + 1]
                                        if numpy.ma.count(subgrid_speed) >= 20:
                                            s111_file.add_feature_instance_group_data(
                                                model_file.datetime_values[time_index], subgrid_speed,
                                                subgrid_direction, cycletime, target_depth)
                                            s111_file.add_model_metadata()

                                        else:
                                            s111_file.close()
                                            os.remove('{}'.format(s111_file.path))
                                            s111_file_paths.remove(s111_file.path)
                        else:

                            s111_file.add_feature_instance_group_data(model_file.datetime_values[time_index], speed,
                                                                      direction, cycletime, target_depth)
                            s111_file.add_model_metadata()

                finally:
                    model_file.close()
        finally:
            model_index_file.close()

    else:
        with S111File('{}.h5'.format(s111_path_prefix), input_metadata, data_coding_format, clobber=True) as s111_file:
            s111_file_paths.append(s111_file.path)

            for model_file in model_files:
                try:
                    model_file.open()
                    for time_index in range(len(model_file.datetime_values)):

                        # Get native-grid output with invalid/masked values removed
                        u_compressed, v_compressed, lat_compressed, lon_compressed = model_file.output_native_grid(
                            time_index, target_depth)

                        # Convert currents from u/v to speed/direction
                        speed, direction = model.irregular_uv_to_speed_direction(u_compressed, v_compressed)

                        s111_file.add_feature_instance_group_data(model_file.datetime_values[time_index], speed,
                                                                  direction, cycletime, target_depth)

                finally:
                    model_file.close()

            s111_file.add_positioning(lon_compressed, lat_compressed)
            s111_file.add_model_metadata()

    return s111_file_paths


def time_series_to_s111(input_data, s111_path_prefix, input_metadata, data_coding_format, current_depth):
    """Convert oceanographic time series data to S111 format.

    Current observations and predictions at fixed or moving stations.

    Args:
        input_data: List of ``S111TimeSeries`` (or subclasses thereof) instance
            describing observations or predictions time series data, which includes,
            1d `ndarrays` of latitude, longitude, direction and speed.
        s111_path_prefix: Path prefix for desired output location for generated
            S-111 files. If specified path is a directory, file(s) will be
            output to specified directory with autogenerated names. Otherwise,
            generated file(s) will be placed at specified file path, but with a
            filename suffix appended based on the properties of the target
            output grid and the model file.
        input_metadata: ``S111Metadata`` instance describing metadata for
            geographic identifier and description of current meter type,
            forecast method, or station identifier (e.g. "cb0201").
        data_coding_format: 1:Time series at fixed stations, 2:Regularly gridded arrays,
            3:Ungeorectified gridded arrays, 4:Time series for one moving platform.
        current_depth: The water current at a specified target depth below
            the sea surface in meters.

    """
    timestamp = input_data[0].datetime_values[0]

    if os.path.isdir(s111_path_prefix):
        if not s111_path_prefix.endswith('/'):
            s111_path_prefix += '/'
        file_issuance = datetime.datetime.strftime(timestamp, '%Y%m%dT%H%M%SZ')
        s111_path_prefix += (
            'S111{}_{}_{}_TYP{}'.format(input_metadata.producer_code, file_issuance, input_metadata.region, data_coding_format))

        stations_longitude = []
        stations_latitude = []

        with S111File('{}.h5'.format(s111_path_prefix), input_metadata, data_coding_format, clobber=True) as s111_file:

            if data_coding_format == 1:
                for station in input_data:
                    s111_file.add_feature_instance_group_data(station.datetime_values[0], station.speed,
                                                              station.direction, timestamp, current_depth)

                    stations_longitude.append(station.longitude)
                    stations_latitude.append(station.latitude)

                s111_file.add_positioning(stations_longitude, stations_latitude)

                s111_file.add_time_series_metadata(input_data[0].datetime_values)

            else:
                for obs in input_data:
                    s111_file.add_feature_instance_group_data(obs.datetime_values[0], obs.speed, obs.direction, timestamp, current_depth)

                s111_file.add_positioning(obs.longitude, obs.latitude)

                s111_file.add_time_series_metadata(input_data[0].datetime_values)

