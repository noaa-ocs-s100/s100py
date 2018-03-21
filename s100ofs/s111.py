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
FILLVALUE=-9999.0

class S111File:
    """Create and manage S-111 files.

    This class implements the context manager pattern and should thus be used
    similar to the following:

        with S111File("myfile.h5") as f:
            ...

    Attributes:
        path: Path (relative or absolute) to the file.
    """
    def __init__(self, path):
        """Initializes S111File object and opens h5 file at specified path.

        If `path` has an extension other than ".h5", it is replaced with
        ".h5"

        Args:
            path: Path of target hdf5 file. Must end in ".h5", otherwise its
        extension will be replaced with ".h5".
        """
        filepath, file_extension = os.path.splitext(path)
        self.path = filepath + ".h5"
        if os.path.exists(self.path):
            # File already exists, open in append mode
            self.h5_file = h5py.File(self.path, "r+")
        else:
            # File doesn't exist, open in create (write) mode and add metadata
            self.h5_file = h5py.File(self.path, "w")
            self.add_metadata()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5_file.flush() # is this required??
        self.h5_file.close()

    def close(self):
        self.h5_file.close()

    def add_metadata(self):
        """Add empty global metadata attributes to HDF5 file."""
        
        # Integer types
        self.h5_file.attrs.create('horizDatumValue', 0 , dtype=numpy.int32) 
        self.h5_file.attrs.create('timeRecordInterval', 0 , dtype=numpy.int32)
        self.h5_file.attrs.create('numberOfTimes', 0 , dtype=numpy.int32)
        self.h5_file.attrs.create('numberOfStations', 0 , dtype=numpy.int32)
        self.h5_file.attrs.create('numPointsLongitudinal', 0 , dtype=numpy.int32)
        self.h5_file.attrs.create('numPointsLatitudinal', 0 , dtype=numpy.int32)
        self.h5_file.attrs.create('numberOfNodes', 0 , dtype=numpy.int32)
       
        # Real types
        self.h5_file.attrs.create('surfaceCurrentDepth', 0 , dtype=numpy.float32)
        self.h5_file.attrs.create('gridOriginLongitude', 0 , dtype=numpy.float32)
        self.h5_file.attrs.create('gridOriginLatitude', 0 , dtype=numpy.float32)
        self.h5_file.attrs.create('gridSpacingLongitudinal', 0 , dtype=numpy.float32)
        self.h5_file.attrs.create('gridSpacingLatitudinal', 0 , dtype=numpy.float32)
        self.h5_file.attrs.create('minGridPointLongitudinal', 0 , dtype=numpy.float32)
        self.h5_file.attrs.create('minGridPointLatitudinal', 0 , dtype=numpy.float32)
        self.h5_file.attrs.create('gridLandMaskValue', 0 , dtype=numpy.float32)
        self.h5_file.attrs.create('speedUncertainty', -1.0 , dtype=numpy.float32)
        self.h5_file.attrs.create('directionUncertainty', -1.0 , dtype=numpy.float32)
        self.h5_file.attrs.create('positionUncertainty', -1.0 , dtype=numpy.float32)
        self.h5_file.attrs.create('verticalUncertainty', -1.0 , dtype=numpy.float32)
        self.h5_file.attrs.create('timeUncertainty', -1.0 , dtype=numpy.float32)
        self.h5_file.attrs.create('minDatasetCurrentSpeed', 0 , dtype=numpy.float32)
        self.h5_file.attrs.create('maxDatasetCurrentSpeed', 0 , dtype=numpy.float32)

        # String types
        dt = h5py.special_dtype(vlen=str)
        self.h5_file.attrs.create('productSpecification', 'S-111_v1.11.0' , dtype=dt)
        self.h5_file.attrs.create('dateTimeOfIssue', '' , dtype=dt)
        self.h5_file.attrs.create('nameRegion', '' , dtype=dt)
        self.h5_file.attrs.create('nameSubregion', '' , dtype=dt)
        self.h5_file.attrs.create('horizDatumReference', 'EPSG' , dtype=dt)
        self.h5_file.attrs.create('dateTimeOfFirstRecord', '' , dtype=dt)
        self.h5_file.attrs.create('dateTimeOfLastRecord', '' , dtype=dt)
        self.h5_file.attrs.create('methodCurrentsProduct', '' , dtype=dt)

        # Enumeration types
        self.h5_file.attrs.create('typeOfCurrentData', 0 , dtype=numpy.int32)
        self.h5_file.attrs.create('dataCodingFormat', 0 , dtype=numpy.int32)
        self.h5_file.attrs.create('depthTypeIndex', 0 , dtype=numpy.int32)
        self.h5_file.attrs.create('verticalDatum', 0 , dtype=numpy.int32)

    def update_attributes(self, model_index, subgrid_index=None):
        """Update HDF5 attributes based on grid properties.
        
        *********************
        TODO: MOVE REGION/SUBREGION/CURRENTPRODUCT VALUES OUT OF THIS FUNCTION
        *********************

        Args:
            model_index: `ROMSIndexFile` instance representing model index and
                coefficients file.
            subgrid_index: (Optional, default None) Index of subgrid, if any,
                that this S111File represents. Corresponds with index into
                subgrid dimension of model index file.
        """
        # Grid spacing is uniform, so just look at width between first 2 cells
        cellsize_x = model_index.var_x[1] - model_index.var_x[0]
        cellsize_y = model_index.var_y[1] - model_index.var_y[0]

        if subgrid_index is not None:
            if subgrid_index < 0 or subgrid_index >= model_index.dim_subgrid.size:
                raise Exception("Subgrid index [{}] out of model index subgrid dimension range [0-{}]".format(subgrid_index, model_index.dim_subgrid.size-1))
            num_points_lon = 1 + model_index.subgrid_x_max[subgrid_index] - model_index.subgrid_x_min[subgrid_index]
            num_points_lat = 1 + model_index.subgrid_y_max[subgrid_index] - model_index.subgrid_y_min[subgrid_index]
            #***** need to subtract 0.5*cellsize to get actual min_lon/min_lat ???
            # i.e. does min lon/lat represent center, or bottom-left corner of origin pixel???
            min_lon = model_index.var_x[model_index.subgrid_x_min[subgrid_index]]
            min_lat = model_index.var_y[model_index.subgrid_y_min[subgrid_index]]
        else:
            num_points_lon = model_index.dim_x.size
            num_points_lat = model_index.dim_y.size
            #***** need to subtract 0.5*cellsize to get actual min_lon/min_lat ???
            # i.e. does min lon/lat represent center, or bottom-left corner of origin pixel???
            min_lon = numpy.nanmin(model_index.var_x)
            min_lat = numpy.nanmin(model_index.var_y)

        num_nodes = num_points_lon * num_points_lat
                        
        self.h5_file.attrs.modify('gridSpacingLongitudinal', cellsize_x) 
        self.h5_file.attrs.modify('gridSpacingLatitudinal', cellsize_y) 
        self.h5_file.attrs.modify('horizDatumValue', 4326) 
        self.h5_file.attrs.modify('numPointsLongitudinal', num_points_lon)
        self.h5_file.attrs.modify('numPointsLatitudinal', num_points_lat)
        self.h5_file.attrs.modify('minGridPointLongitudinal', min_lon)
        self.h5_file.attrs.modify('minGridPointLatitudinal', min_lat)
        self.h5_file.attrs.modify('numberOfNodes', num_nodes)
        self.h5_file.attrs.modify('surfaceCurrentDepth', 4.5)
        self.h5_file.attrs.modify('gridOriginLongitude', min_lon)
        self.h5_file.attrs.modify('gridOriginLatitude', min_lat)
        self.h5_file.attrs.modify('gridLandMaskValue', FILLVALUE)
        self.h5_file.attrs.modify('dataCodingFormat', 2)
        self.h5_file.attrs.modify('depthTypeIndex', 2)
        self.h5_file.attrs.modify('typeOfCurrentData', 6)
        
        region = numpy.string_("US_East_Coast")
        subRegion = numpy.string_('Chesapeake_Bay')
        methodCurrentProduct = numpy.string_('ROMS_Hydrodynamic_Model')
        
        self.h5_file.attrs.modify('nameRegion', region)
        self.h5_file.attrs.modify('nameSubregion', subRegion)
        self.h5_file.attrs.modify('methodCurrentsProduct', methodCurrentProduct)

    def add_data(self, time_value, reg_grid_speed, reg_grid_direction):
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
        """
        time_str = time_value.strftime("%Y%m%dT%H%M%SZ")

        if len(self.h5_file.items()) == 0:
            num_groups = 1
            new_group_name = 'Group_001' 
            print("Creating", new_group_name ,"dataset.")
            new_group = self.h5_file.create_group(new_group_name)
            
            group_title = 'Regular Grid at DateTime 1' 
            new_group.attrs.create('Title', group_title.encode())

            self.h5_file.attrs.modify('dateTimeOfIssue', time_str.encode())                                 
            self.h5_file.attrs.modify('dateTimeOfFirstRecord', time_str.encode())
            self.h5_file.attrs.modify('dateTimeOfLastRecord', time_str.encode())
                                  
            # Create the compound datatype for Group F attributes
            DIM0 = 2
            DATASET = "Attributes"

            dtype = numpy.dtype([("0", h5py.special_dtype(vlen=str)), 
                              ("1", h5py.special_dtype(vlen=str)),
                              ("2", h5py.special_dtype(vlen=str)),
                              ("3", h5py.special_dtype(vlen=str)),
                              ("4", h5py.special_dtype(vlen=str)),
                              ("5", h5py.special_dtype(vlen=str))])
            
            fdata = numpy.zeros((DIM0,), dtype=dtype)
            fdata['0'][0] = ("surfaceCurrentSpeed")
            fdata['1'][0] = ("Surface current speed")
            fdata['2'][0] = ("knots")
            fdata['3'][0] = ("-9999.0")
            fdata['4'][0] = ("96,56")
            fdata['5'][0] = ("H5T_FLOAT")
            fdata['0'][1] = ("surfaceCurrentDirection")
            fdata['1'][1] = ("Surface current direction")
            fdata['2'][1] = ("degrees")
            fdata['3'][1] = ("-9999.0")
            fdata['4'][1] = ("96,56")
            fdata['5'][1] = ("H5T_FLOAT")

            groupF = self.h5_file.create_group("Group F")
            dset = groupF.create_dataset(DATASET,(DIM0,), dtype = dtype)
            dset[...] = fdata

        else:
            grps = self.h5_file.items()
            num_groups = len(grps)
    
            new_group_name = 'Group {}'.format(num_groups)
            print("Creating", new_group_name, "dataset.")
            new_group = self.h5_file.create_group(new_group_name)
    
            group_title = 'Regular Grid at DateTime {}'.format(num_groups)
            new_group.attrs.create('Title', group_title.encode())
    
            self.h5_file.attrs.modify('dateTimeOfLastRecord', time_str.encode())
            firstTime = datetime.datetime.strptime((self.h5_file.attrs['dateTimeOfFirstRecord']),"%Y%m%dT%H%M%SZ")
            lastTime = datetime.datetime.strptime((self.h5_file.attrs['dateTimeOfLastRecord']),"%Y%m%dT%H%M%SZ")
            interval = lastTime - firstTime
            timeInterval=  interval.total_seconds()
            self.h5_file.attrs.modify('timeRecordInterval', timeInterval)
         
        # Note: make sure numpy.nanmin works with masked arrays       
        min_speed = numpy.nanmin(reg_grid_speed)
        max_speed = numpy.nanmax(reg_grid_speed)
        min_speed = numpy.round(min_speed,2)
        max_speed = numpy.round(max_speed,2)

        directions = directions.filled(FILLVALUE)
        speeds = speeds.filled(FILLVALUE)
        
        # Write data to empty HDF5 datasets 
        new_group.create_dataset('surfaceCurrentDirection', (directions.shape[0], directions.shape[1]), dtype=numpy.float32, data=directions,  chunks=True, compression="gzip", compression_opts=9, fillvalue=FILLVALUE)
        new_group.create_dataset('surfaceCurrentSpeed', (speeds.shape[0], speeds.shape[1]), dtype=numpy.float32, data=speeds,  chunks=True, compression="gzip", compression_opts=9, fillvalue=FILLVALUE)

        # Update attributes from datasets added
        new_group.attrs.create('DateTime', time_str.encode())
        
        if len(self.h5_file.items()) == 2:
            self.h5_file.attrs.modify('minDatasetCurrentSpeed', min_speed)
            self.h5_file.attrs.modify('maxDatasetCurrentSpeed', max_speed)
            self.h5_file.attrs.modify('numberOfTimes', 1)
        else:
            numberOfTimes = num_groups 
            self.h5_file.attrs.modify('numberOfTimes', numberOfTimes)
            prior_min_speed = self.h5_file.attrs['minDatasetCurrentSpeed']
            prior_max_speed = self.h5_file.attrs['maxDatasetCurrentSpeed']
            if min_speed < prior_min_speed:
                self.h5_file.attrs.modify('minDatasetCurrentSpeed', min_speed)
            if max_speed > prior_max_speed:
                self.h5_file.attrs.modify('maxDatasetCurrentSpeed', max_speed)

def romsToS111(roms_index_file, roms_output_files, s111_path_prefix):
    """Convert ROMS model output to regular grid in S111 format.

    Note: Only a single time per ROMS file is currently supported. If a ROMS
    NetCDF includes more than one time/forecast, only the first will be
    extracted.

    *** TODO: Get name of model (cbofs) from index or output file ***

    Args:
        roms_index_file: `ROMSIndexFile` instance containing precalculated grid
            and interpolation information.
        roms_output_files: List of `ROMSOutputFile` instances pointing to one
            or more NetCDF output files from a ROMS-based modeling system.
            Files should be provided in ascending chronological order, as this
            order will be maintained when appending subsequent
            nowcasts/forecasts to each other in individual S111 files.
        s111_path_prefix: Path prefix for desired output location for generated
            S-111 files. If specified path ends in "/", file(s) will be output
            to specified directory with autogenerated names. Otherwise,
            generated file(s) will be placed at specified file path, but with a
            filename suffix appended based on the properties of the target
            output grid and the ROMS file.
    """
    # Open index file
    # Determine output files to be created (one for whole domain or one per subgrid)
    # Open/initialize each s111 output file
    # For each ROMS file:
    #   interpolate to full grid, convert uv to spd/dir
    #   subset if configured
    #   output to s111 file(s)
    if s111_path_prefix.endswith("/"):
        s111_path_prefix += "cbofs"
    with roms.ROMSIndexFile(roms_index_file) as roms_index:
        if roms_index.dim_subgrid is not None and roms_index.var_subgrid_mask is not None:
            # Output to subgrids
            with contextlib.ExitStack() as stack:
                s111_files = []
                for i in range(roms_index.dim_subgrid.size):
                    s111_file = S111File("{}_subgrid_{}".format(s111_path_prefix, roms_index.var_subgrid_id[i]))
                    stack.enter_context(s111_file)
                    s111_file.update_attributes(roms_index, i)
                    s111_files.append(s111_file)
                roms_files = []
                times = []
                for roms_output_file in roms_output_files:
                    roms_file = roms.ROMSOutputFile(roms_output_file)
                    stack.enter_context(roms_file)
                    roms_files.append(roms_file)
                    # Convert gregorian timestamp to datetime timestamp
                    time_val = netCDF4.num2date(roms_file.nc_file.variables['ocean_time'][:], model_output.nc_file.variables['ocean_time'].units)[0]
                    times.append(time_val)
                
                for i, roms_file in enumerate(roms_files):
                    reg_grid_u, reg_grid_v = roms_file.uvToRegularGrid(roms_index)
                    # Convert currents at regular grid points from u/v to speed
                    # and direction
                    directions, speeds = roms.convertUVToSpeedDirection(reg_grid_u, reg_grid_v)
                    directions = ma.masked_array(directions, roms_index.var_xi1.mask)  
                    speeds = ma.masked_array(speeds, roms_index.var_xi1.mask)
                    
                    for subgrid_index, s111_file in enumerate(s111_files):
                        x_min = roms_index.var_subgrid_x_min[subgrid_index]
                        x_max = roms_index.var_subgrid_x_max[subgrid_index]
                        y_min = roms_index.var_subgrid_y_min[subgrid_index]
                        y_max = roms_index.var_subgrid_y_max[subgrid_index]
                        subgrid_speed = speeds[y_min:y_max+1, x_min:x_max+1]
                        subgrid_direction = directions[y_min:y_max+1, x_min:x_max+1]
                        s111_file.add_data(times[i], subgrid_speed, subgrid_direction)
        else:
            # Output to default grid (no subgrids)
                        


    def add_output(self, model_output_file, model_index_file):
        """Add (append) model output to the S111 file

        Args:
            model_output_file: Path to model output file to be appended.
            model_index_file: Path to model index file required to perform
                regular grid interpolation.
        """
        # Open model output netCDF
        with roms.ROMSOutputFile(model_output_file) as model_output:
            with roms.ROMSIndexFile(model_index_file) as model_index:
                reg_grid_u,reg_grid_v = model_output.uvToRegularGrid(model_index)
                
                        s111_file
        else:
            # Output to default grid (no subgrids)
                        


    def add_output(self, model_output_file, model_index_file):
        """Add (append) model output to the S111 file

        Args:
            model_output_file: Path to model output file to be appended.
            model_index_file: Path to model index file required to perform
                regular grid interpolation.
        """
        # Open model output netCDF
        with roms.ROMSOutputFile(model_output_file) as model_output:
            with roms.ROMSIndexFile(model_index_file) as model_index:
                reg_grid_u,reg_grid_v = model_output.uvToRegularGrid(model_index)
                
                # Create HDF5 groups and datasets
                self.create_group(model_output, model_index, reg_grid_u, reg_grid_v)
                
                # Update HDF5 attributes
                self.update_attributes(model_index, subgrid_id)
                
                print("Data sucessfully added")

