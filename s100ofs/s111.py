"""Utilities for creation and modification of S-111 compliant HDF-5 files.

S-111 is an IHO standard outlining formats for storing and sending surface
water current data and metadata.
"""
from contextlib import ContextDecorator
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
                self.update_attributes(model_index)
                
                print("Data sucessfully added")


    def create_group(self, model_output, model_index, reg_grid_u, reg_grid_v):
        """Create inital HDF5 Group with Speeds and Directions.

        For every additional NetCDF file Create an HDF5 Group containing Speeds
        and Direction Datasets. Update HDF5 time attributes.

        Args:
            model_output: `ROMSOutputFile` representing the source model output
                file.
            model_index: `ROMSIndexFile` representing the model index and
                coefficients file.
            reg_grid_u: `numpy.ma.masked_array` representing U data after
                interpolating to a regular grid.
            reg_grid_v: `numpy.ma.masked_array` representing V data after
                interpolating to a regular grid.
        """
        # Convert gregorian timestamp to datetime timestamp
        time = netCDF4.num2date(model_output.nc_file.variables['ocean_time'][:], model_output.nc_file.variables['ocean_time'].units)
        numberOfTimes = time.shape[0]
        time_val = time[0]
        str_time_val = time_val.strftime("%Y%m%dT%H%M%SZ")
        

        if len(self.h5_file.items()) == 0:  
            new_group_name = 'Group_001' 
            print("Creating", new_group_name ,"dataset.")
            new_group = self.h5_file.create_group(new_group_name)
            
            group_title = 'Regular Grid at DateTime 1' 
            new_group.attrs.create('Title', group_title.encode())

            self.h5_file.attrs.modify('dateTimeOfIssue', str_time_val.encode())                                 
            self.h5_file.attrs.modify('dateTimeOfFirstRecord', str_time_val.encode())
            self.h5_file.attrs.modify('dateTimeOfLastRecord', str_time_val.encode())
                                  
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
    
            new_group_name = 'Group' + ' ' + str(num_groups)
            print("Creating", new_group_name, "dataset.")
            new_group = self.h5_file.create_group(new_group_name)
    
            group_title = 'Regular Grid at DateTime ' + str(num_groups)
            new_group.attrs.create('Title', group_title.encode())
    
            self.h5_file.attrs.modify('dateTimeOfLastRecord', str_time_val.encode())
            firstTime = datetime.datetime.strptime((self.h5_file.attrs['dateTimeOfFirstRecord']),"%Y%m%dT%H%M%SZ")
            lastTime = datetime.datetime.strptime((self.h5_file.attrs['dateTimeOfLastRecord']),"%Y%m%dT%H%M%SZ")
            interval = lastTime - firstTime
            timeInterval=  interval.total_seconds()
            self.h5_file.attrs.modify('timeRecordInterval', timeInterval)
                
        # Convert currents at regular grid points from u/v to speed and
        # direction
        directions, speeds = roms.convertUVToSpeedDirection(reg_grid_u, reg_grid_v, model_index)    

        min_speed = numpy.nanmin(speeds)
        max_speed = numpy.nanmax(speeds)
        min_speed = numpy.round(min_speed,2)
        max_speed = numpy.round(max_speed,2)
        directions = ma.masked_array(directions, model_index.var_xi1.mask)  
        speeds = ma.masked_array(speeds, model_index.var_xi1.mask)
        directions = directions.filled(FILLVALUE)
        speeds = speeds.filled(FILLVALUE)
        
        # Write data to empty HDF5 datasets 
        new_group.create_dataset('surfaceCurrentDirection', (directions.shape[0], directions.shape[1]), dtype=numpy.float32, data=directions,  chunks=True, compression="gzip", compression_opts=9, fillvalue=FILLVALUE)
        new_group.create_dataset('surfaceCurrentSpeed', (speeds.shape[0], speeds.shape[1]), dtype=numpy.float32, data=speeds,  chunks=True, compression="gzip", compression_opts=9, fillvalue=FILLVALUE)

        # Update attributes from datasets added
        new_group.attrs.create('DateTime', str_time_val.encode())
        
        if len(self.h5_file.items()) == 2:
            self.h5_file.attrs.modify('minDatasetCurrentSpeed', min_speed)
            self.h5_file.attrs.modify('maxDatasetCurrentSpeed', max_speed)
            self.h5_file.attrs.modify('numberOfTimes', numberOfTimes)
        else:
            numberOfTimes = num_groups 
            self.h5_file.attrs.modify('numberOfTimes', numberOfTimes)
            prior_min_speed = self.h5_file.attrs['minDatasetCurrentSpeed']
            prior_max_speed = self.h5_file.attrs['maxDatasetCurrentSpeed']
            if min_speed < prior_min_speed:
                self.h5_file.attrs.modify('minDatasetCurrentSpeed', min_speed)
            if max_speed > prior_max_speed:
                self.h5_file.attrs.modify('maxDatasetCurrentSpeed', max_speed)


    def update_attributes(self, model_index):
        """Update HDF5 attributes based on grid properties.
        
        *********************
        TODO: MOVE REGION/SUBREGION/CURRENTPRODUCT VALUES OUT OF THIS FUNCTION
        *********************

        Args:
            model_index: `ROMSIndexFile` instance representing model index and
                coefficients file.
        """
        num_nodes = model_index.dim_x.size * model_index.dim_y.size
        min_lon = numpy.nanmin(model_index.var_x)
        min_lat = numpy.nanmin(model_index.var_y)
        grid_spacing_lon = model_index.var_x[1] - model_index.var_x[0]
        grid_spacing_lat = model_index.var_y[1] - model_index.var_y[0]
                        
        self.h5_file.attrs.modify('gridSpacingLongitudinal',grid_spacing_lon) 
        self.h5_file.attrs.modify('gridSpacingLatitudinal',grid_spacing_lat) 
        self.h5_file.attrs.modify('horizDatumValue', 4326) 
        self.h5_file.attrs.modify('numPointsLongitudinal', model_index.dim_x.size)
        self.h5_file.attrs.modify('numPointsLatitudinal', model_index.dim_y.size)
        self.h5_file.attrs.modify('minGridPointLongitudinal', min_lon)
        self.h5_file.attrs.modify('minGridPointLatitudinal', min_lat)
        self.h5_file.attrs.modify('numberOfNodes', num_nodes)
        self.h5_file.attrs.modify('surfaceCurrentDepth', 4.5)
        self.h5_file.attrs.modify('gridOriginLongitude', min_lon)
        self.h5_file.attrs.modify('gridOriginLatitude', min_lat)
        self.h5_file.attrs.modify('gridLandMaskValue', FILLVALUE )
        self.h5_file.attrs.modify('dataCodingFormat', 2)
        self.h5_file.attrs.modify('depthTypeIndex', 2)
        self.h5_file.attrs.modify('typeOfCurrentData', 6)
        
        region = numpy.string_("US_East_Coast")
        subRegion = numpy.string_('Cheaspeake_Bay')
        methodCurrentProduct = numpy.string_('ROMS_Hydrodynamic_Model')
        
        self.h5_file.attrs.modify('nameRegion', region)
        self.h5_file.attrs.modify('nameSubregion', subRegion)
        self.h5_file.attrs.modify('methodCurrentsProduct', methodCurrentProduct)

