"""Utilities for creation and modification of S-111 compliant HDF-5 files.

S-111 is an IHO standard outlining formats for storing and sending surface
water current data and metadata.
"""
import contextlib
import datetime
import math
import os

import h5py
import netCDF4
import numpy
import numpy.ma as ma

from .model import roms

# Default fill value for NetCDF variables
FILLVALUE = -9999.0

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
            self.add_structure_and_metadata()
        else:
            # File already exists, open in append mode
            self.h5_file = h5py.File(self.path, "r+")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5_file.close()

    def close(self):
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
        self.h5_file.attrs.create('metaFeatures', '', dtype=h5py.special_dtype(vlen=str))
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

        # Create s111 structure, feature group, feature type container, and initial feature instance
        self.groupF = self.h5_file.create_group("Group_F")
        self.feature = self.h5_file.create_group("SurfaceCurrent")
        self.feature_instance = self.feature.create_group("SurfaceCurrent.01")

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
        self.feature.attrs.create('sequenceRule.scanDirection', "", dtype=h5py.special_dtype(vlen=str))

        # Add feature instance metadata
        # String types
        self.feature_instance.attrs.create('startSequence', '', dtype=h5py.special_dtype(vlen=str))
        self.feature_instance.attrs.create('dateTimeOfFirstRecord', '', dtype=h5py.special_dtype(vlen=str))
        self.feature_instance.attrs.create('dateTimeOfLastRecord', '', dtype=h5py.special_dtype(vlen=str))
        # Integer types
        self.feature_instance.attrs.create('numPointsLongitudinal', 0, dtype=numpy.int32)
        self.feature_instance.attrs.create('numPointsLatitudinal', 0, dtype=numpy.int32)
        self.feature_instance.attrs.create('timeRecordInterval', 0, dtype=numpy.int32)
        self.feature_instance.attrs.create('numberOfTimes', 0, dtype=numpy.int32)
        self.feature_instance.attrs.create('numGRP', 0, dtype=numpy.int32)
        self.feature_instance.attrs.create('minGridPointLongitudinal', 0, dtype=numpy.int32)
        self.feature_instance.attrs.create('minGridPointLatitudinal', 0, dtype=numpy.int32)
        # Real types
        self.feature_instance.attrs.create('gridOriginLongitude', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('gridOriginLatitude', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('gridSpacingLongitudinal', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('gridSpacingLatitudinal', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('westBoundLongitude', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('eastBoundLongitude', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('southBoundLatitude', 0, dtype=numpy.float32)
        self.feature_instance.attrs.create('northBoundLatitude', 0, dtype=numpy.float32)

    def update_attributes(self, model_index, ofs_product, ofs_region, subgrid_index=None):
        """Update HDF5 attributes based on grid properties.
        
        Args:
            model_index: `ROMSIndexFile` instance representing model index and
                coefficients file.
            subgrid_index: (Optional, default None) Index of subgrid, if any,
                that this S111File represents. Corresponds with index into
                subgrid dimension of model index file.
            ofs_region: Geographic indentifier metadata.
            ofs_product: Model description and type of forecast metadata.
        """

        # Width between first two cells, grid spacing is uniform
        cellsize_x = model_index.var_x[1] - model_index.var_x[0]
        cellsize_y = model_index.var_y[1] - model_index.var_y[0]

        if subgrid_index is not None:
            if subgrid_index < 0 or subgrid_index >= model_index.dim_subgrid.size:
                raise Exception(
                    "Subgrid index [{}] out of model index subgrid dimension range [0-{}]".format(subgrid_index, model_index.dim_subgrid.size - 1))
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

        num_nodes = num_points_lon * num_points_lat

        # Update carrier metadata
        self.h5_file.attrs.modify('geographicIdentifier', ofs_region)
        self.h5_file.attrs.modify('productSpecification', S111Metadata.productSpecification)
        self.h5_file.attrs.modify('epoch', S111Metadata.epoch)
        self.h5_file.attrs.modify('horizontalDatumReference', S111Metadata.horizontalDatumReference)
        self.h5_file.attrs.modify('horizontalDatumValue', S111Metadata.horizontalDatumValue)
        self.h5_file.attrs.modify('depthTypeIndex', S111Metadata.depthTypeIndex)
        self.h5_file.attrs.modify('surfaceCurrentDepth', - 4.5)
        self.h5_file.attrs.modify('westBoundLongitude', min_lon)
        self.h5_file.attrs.modify('eastBoundLongitude', max_lon)
        self.h5_file.attrs.modify('southBoundLatitude', min_lat)
        self.h5_file.attrs.modify('northBoundLatitude', max_lat)

        # Update feature container metadata
        self.feature.attrs.modify('dataCodingFormat', S111Metadata.dataCodingFormat)
        self.feature.attrs.modify('interpolationType', S111Metadata.interpolationType)
        self.feature.attrs.modify('typeOfCurrentData', S111Metadata.typeOfCurrentData)
        self.feature.attrs.modify('commonPointRule', S111Metadata.commonPointRule)
        self.feature.attrs.modify('dimension', S111Metadata.dimension)
        self.feature.attrs.modify('sequenceRule.type', S111Metadata.sequenceRuleType)
        self.feature.attrs.modify('methodCurrentsProduct', ofs_product)
        scan_direction = numpy.string_('latitude,longitude')
        self.feature.attrs.modify('sequenceRule.scanDirection', S111Metadata.sequenceRuleScanDirection)

        # Update feature instance metadata
        if self.feature.attrs['dataCodingFormat'] == S111Metadata.dataCodingFormat:
            self.feature_instance.attrs.modify('gridOriginLongitude', min_lon)
            self.feature_instance.attrs.modify('gridOriginLatitude', min_lat)
            self.feature_instance.attrs.modify('gridSpacingLongitudinal', cellsize_x)
            self.feature_instance.attrs.modify('gridSpacingLatitudinal', cellsize_y)
            self.feature_instance.attrs.modify('numPointsLongitudinal', num_points_lon)
            self.feature_instance.attrs.modify('numPointsLatitudinal', num_points_lat)
            start_sequence = numpy.string_('0,0')
            self.feature_instance.attrs.modify('startSequence', start_sequence)
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
            dtype = numpy.dtype([("code", h5py.special_dtype(vlen=str)),
                                 ("name", h5py.special_dtype(vlen=str)),
                                 ("uom.name", h5py.special_dtype(vlen=str)),
                                 ("fillValue", h5py.special_dtype(vlen=str)),
                                 ("dataType", h5py.special_dtype(vlen=str)),
                                 ("lower", h5py.special_dtype(vlen=str)),
                                 ("upper", h5py.special_dtype(vlen=str)),
                                 ("closure", h5py.special_dtype(vlen=str))])

            fdata = numpy.zeros((2,), dtype=dtype)
            fdata['code'][0] = "surfaceCurrentSpeed"
            fdata['name'][0] = "Surface current speed"
            fdata['uom.name'][0] = "knots"
            fdata['fillValue'][0] = str(FILLVALUE)
            fdata['dataType'][0] = H5T_CLASS_T[h5py.h5t.FLOAT]
            fdata['lower'][0] = 0.0
            fdata['upper'][0] = ""
            fdata['closure'][0] = "geSemilInterval"

            fdata['code'][1] = "surfaceCurrentDirection"
            fdata['name'][1] = "Surface current direction"
            fdata['uom.name'][1] = "degrees"
            fdata['fillValue'][1] = str(FILLVALUE)
            fdata['dataType'][1] = H5T_CLASS_T[h5py.h5t.FLOAT]
            fdata['lower'][1] = 0.0
            fdata['upper'][1] = 359.9
            fdata['closure'][1] = "closedInterval"

            self.groupF_dset = self.groupF.create_dataset("SurfaceCurrent", (2,), dtype=dtype)
            self.groupF_dset[...] = fdata

            # Add group_f feature code dataset
            fc_data = numpy.zeros((1,), dtype=h5py.special_dtype(vlen=str))
            fc_data[0] = "SurfaceCurrent"
            feature_code = self.groupF.create_dataset("featureCode", (1,), dtype=h5py.special_dtype(vlen=str))
            feature_code[...] = fc_data

            # Add feature container dataset
            axis_names = numpy.zeros((2,), dtype=h5py.special_dtype(vlen=str))
            axis_names[0] = "longitude"
            axis_names[1] = "latitude"
            axis_dset = self.feature.create_dataset("axisNames", (2,), dtype=h5py.special_dtype(vlen=str))
            axis_dset[...] = axis_names

        # Create a list of all feature instance objects and groups
        feature_instance_objs = []
        self.feature_instance.visit(feature_instance_objs.append)
        feature_instance_groups = [obj for obj in feature_instance_objs if
                                   isinstance(self.feature_instance[obj], h5py.Group)]

        # Convert time value to string
        time_str = time_value.strftime("%Y%m%dT%H%M%SZ")

        if len(feature_instance_groups) == 0:
            new_group = self.feature_instance.create_group('Group_001')
            print("Creating", "Group_001", "dataset.")

            # Time attributes updated once
            issuance_time = cycletime.strftime("%H%M%SZ")
            issuance_date = cycletime.strftime("%Y%m%d")
            self.h5_file.attrs.modify('issueTime', issuance_time.encode())
            self.h5_file.attrs.modify('issueDate', issuance_date.encode())
            self.feature.attrs.modify('numInstances', len(self.feature_instance))
            self.feature_instance.attrs.modify('dateTimeOfFirstRecord', time_str.encode())
            self.feature_instance.attrs.modify('dateTimeOfLastRecord', time_str.encode())

        else:
            num_grps = len(feature_instance_groups)
            new_grp = num_grps + 1
            new_group = self.feature_instance.create_group('Group_{}'.format(str(new_grp).zfill(3)))
            print("Creating", "Group_{}".format(str(new_grp).zfill(3)), "dataset.")
            self.feature_instance.attrs.modify('dateTimeOfLastRecord', time_str.encode())

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
        new_group.attrs.create('timePoint', time_str.encode())
        self.feature_instance.attrs.modify('numberOfTimes', len(self.feature_instance))
        self.feature_instance.attrs.modify('numGRP', len(self.feature_instance))

        # Write data to empty feature instance group
        values_dtype = numpy.dtype([("SurfaceCurrentSpeed", numpy.float32),
                                    ("SurfaceCurrentDirection", numpy.float32)])

        values = numpy.zeros(speed.shape, dtype=values_dtype)
        values['SurfaceCurrentSpeed'] = speed
        values['SurfaceCurrentDirection'] = direction
        values_dset = new_group.create_dataset('values', speed.shape, dtype=values_dtype, chunks=True,
                                               compression="gzip", compression_opts=9)
        values_dset[...] = values

        # Update group_f attributes
        self.groupF_dset.attrs.create("chunking", str(values_dset.chunks), dtype=h5py.special_dtype(vlen=str))

        if len(self.feature_instance) == 2:
            # Time record interval is the same through out the forecast.
            first_time = datetime.datetime.strptime((self.feature_instance.attrs['dateTimeOfFirstRecord']),
                                                    "%Y%m%dT%H%M%SZ")
            last_time = datetime.datetime.strptime((self.feature_instance.attrs['dateTimeOfLastRecord']),
                                                   "%Y%m%dT%H%M%SZ")
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

        self.h5_file.flush()


class S111Metadata():
    """Contains s111 metadata to pass to s111File.

    product_specification: The product specification used to create this dataset.
    epoch: Code denoting the epoch of the geodetic datum used by the CRS.
    horizontal_datum_reference: Reference to the register from which the horizontal datum value is taken.
    horizontal_datum_value: Horizontal Datum of the entire dataset.
    dataCodingFormat: Time series at fixed stations, regularly gridded arrays, ungeorectified gridded arrays,
                      moving platform
    depthTypeIndex:
    typeOfCurrentData: Historical, real-time, astronomical , hybrid, hindcast, forecast.
    metaFeatures: Name of metafeatures file.
    metadata: Name of XML metadata file.
    interpolationType: Interpolation method recommended for evaluation of the S100_GridCoverage.
    commonPointRule: The procedure used for evaluating geometric objects that overlap or lie fall on boundaries.
    dimension: The dimension of the feature instance
    sequenceRule_type: Method to assign values from the sequence of values to the grid coordinates (e.g. "linear")
    sequenceRule_scan_direction: axisNames, comma-separated (e.g. "latitude,longitude")

    """
    productSpecification = numpy.string_('INT.IHO.S-111.1.0.0')
    epoch = numpy.string_('G1762')
    horizontalDatumReference = numpy.string_('EPSG')
    horizontalDatumValue = 4326
    dataCodingFormat = 2
    depthTypeIndex = 2
    typeOfCurrentData = 6
    metaFeatures = None
    metadata = None
    interpolationType = 10
    commonPointRule = 3
    dimension = 2
    sequenceRuleType = 1
    sequenceRuleScanDirection = numpy.string_('latitude,longitude')


def roms_to_s111(roms_index_path, roms_output_paths, s111_path_prefix, cycletime, ofs_model, ofs_product, ofs_region):
    """Convert ROMS model output to regular grid in S111 format.

    If the supplied ROMS index NetCDF contains information identifying
    subgrids, one S111 file will be generated for each subgrid. Otherwise, a
    single S111 file will be created for the entire domain.

    Note: Only a single time per ROMS file is currently supported. If a ROMS
    NetCDF includes more than one time/forecast, only the first will be
    extracted.

    Args:
        roms_index_path: Path to ROMS index NetCDF file containing
            precalculated grid and interpolation information.
        roms_output_paths: List of paths to one or more NetCDF output files
            from a ROMS-based modeling system. Files should be provided in
            ascending chronological order, as this order will be maintained
            when appending subsequent nowcasts/forecasts to each other in
            individual S111 files.
        s111_path_prefix: Path prefix for desired output location for generated
            S-111 files. If specified path is a directory, file(s) will be
            output to specified directory with autogenerated names. Otherwise,
            generated file(s) will be placed at specified file path, but with a
            filename suffix appended based on the properties of the target
            output grid and the ROMS file.
        cycletime: `datetime.datetime` instance representing target cycle time
            of model forecast(s) being processed.
        ofs_model: Model identifier (e.g. "cbofs").
        ofs_region: Geographic indentifier metadata.
        ofs_product: Model description and type of forecast metadata.
    """
    # Path format/prefix for output S111 files. Forecast initialization (reference).
    if os.path.isdir(s111_path_prefix):
        if not s111_path_prefix.endswith("/"):
            s111_path_prefix += "/"
        file_issuance = cycletime.strftime("%Y%m%dT%HZ")
        s111_path_prefix += ("S111US_{}_{}_TYP2".format(file_issuance, str.upper(ofs_model)))
    with roms.ROMSIndexFile(roms_index_path) as roms_index:
        if roms_index.dim_subgrid is not None and roms_index.var_subgrid_id is not None:
            # Output to subgrids
            with contextlib.ExitStack() as stack:
                s111_files = []
                for i in range(roms_index.dim_subgrid.size):
                    s111_file = S111File("{}_SUBGRID_{}.h5".format(s111_path_prefix, roms_index.var_subgrid_id[i]),
                                         clobber=True)
                    stack.enter_context(s111_file)
                    s111_file.update_attributes(roms_index, ofs_product, ofs_region, i)
                    s111_files.append(s111_file)

                for roms_output_path in roms_output_paths:
                    with roms.ROMSOutputFile(roms_output_path) as roms_file:
                        # Convert gregorian timestamp to datetime timestamp
                        time_val = netCDF4.num2date(roms_file.nc_file.variables['ocean_time'][:],
                                                    roms_file.nc_file.variables['ocean_time'].units)[0]
                        reg_grid_u, reg_grid_v = roms_file.uv_to_regular_grid(roms_index)
                        # Convert currents at regular grid points from u/v to speed
                        # and direction
                        direction, speed = roms.uv_to_speed_direction(reg_grid_u, reg_grid_v)
                        direction = ma.masked_array(direction, roms_index.var_xi1.mask)
                        speed = ma.masked_array(speed, roms_index.var_xi1.mask)

                        for subgrid_index, s111_file in enumerate(s111_files):
                            if os.path.isfile(s111_file.path):
                                x_min = roms_index.var_subgrid_x_min[subgrid_index]
                                x_max = roms_index.var_subgrid_x_max[subgrid_index]
                                y_min = roms_index.var_subgrid_y_min[subgrid_index]
                                y_max = roms_index.var_subgrid_y_max[subgrid_index]
                                subgrid_speed = speed[y_min:y_max + 1, x_min:x_max + 1]
                                subgrid_direction = direction[y_min:y_max + 1, x_min:x_max + 1]
                                if ma.count(subgrid_speed) >= 20:
                                    s111_file.add_data(time_val, subgrid_speed, subgrid_direction, cycletime)
                                else:
                                    s111_file.close()
                                    os.remove("{}".format(s111_file.path))
        else:
            # Output to default grid (no subgrids)
            with S111File("{}.h5".format(s111_path_prefix), clobber=True) as s111_file:
                s111_file.update_attributes(roms_index, ofs_product, ofs_region)
                for roms_output_path in roms_output_paths:
                    with roms.ROMSOutputFile(roms_output_path) as roms_file:
                        # Convert gregorian timestamp to datetime timestamp
                        time_val = netCDF4.num2date(roms_file.nc_file.variables['ocean_time'][:],
                                                    roms_file.nc_file.variables['ocean_time'].units)[0]

                        # Call roms method and convert and interpolate u/v to regular grid
                        reg_grid_u, reg_grid_v = roms_file.uv_to_regular_grid(roms_index)

                        # Convert currents at regular grid points from u/v to speed
                        # and direction
                        direction, speed = roms.uv_to_speed_direction(reg_grid_u, reg_grid_v)
                        direction = ma.masked_array(direction, roms_index.var_xi1.mask)
                        speed = ma.masked_array(speed, roms_index.var_xi1.mask)

                        s111_file.add_data(time_val, speed, direction, cycletime)
