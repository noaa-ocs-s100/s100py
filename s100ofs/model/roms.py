"""Utility classes and methods for working with various ocean models."""
import datetime
import math
import os

import netCDF4
import numpy
import numpy.ma as ma

# Conversion factor for meters/sec to knots
MS2KNOTS = 1.943844

class RegularGridUV:
    """Encapsulate regular grid U, V, and related elements/variables.

    Attributes:
        xgrid: 
        ygrid: 
        grid_x: 
        grid_y: 
        ugrid: 
        vgrid: 
        min_lon: 
        min_lat: 
        ncindex_xi1: 
    """
    def __init__(self, xgrid, ygrid, grid_x, grid_y, ugrid, vgrid, min_lon, min_lat, ncindex_xi1):
        """Initialize RegularGrid object."""
        self.xgrid = xgrid
        self.ygrid = ygrid
        self.grid_x = grid_x
        self.grid_y = grid_y
        self.ugrid = ugrid
        self.vgrid = vgrid
        self.min_lon = min_lon
        self.min_lat = min_lat
        self.ncindex_xi1 = ncindex_xi1

class ROMSIndexFile:
    """Store information about an index file used during interpolation.
    
    Attributes:
        
    """
    def __init__(self, path):
        """Initialize ROMSIndexFile object and open file at specified path.

        Args:
            path: Path of target NetCDF file.
        """
        self.path = path
        if os.path.exists(self.path):
            self.nc_file = netCDF4.Dataset(self.path, "r", format="NETCDF4")
        else:
            # File doesn't exist, raise error
            raise(Exception("NetCDF file does not exist: {}".format(self.path)))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.nc_file.close()

    def close(self):
        self.nc_file.close()

class ROMSOutputFile:
    """Read/process data from a ROMS model output file.

    This class implements the context manager pattern and should thus be used
    similar to the following:

        with ROMSOutputFile("cbofs.nc") as f:
            ...

    Attributes:
        path: Path (relative or absolute) to the file.
    """
    def __init__(self, path):
        """Initialize ROMSOutputFile object and open file at specified path.

        Args:
            path: Path of target NetCDF file.
        """
        self.path = path
        if os.path.exists(self.path):
            self.nc_file = netCDF4.Dataset(self.path, "r", format="NETCDF3_Classic")
        else:
            # File doesn't exist, raise error
            raise(Exception("NetCDF file does not exist: {}".format(self.path)))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.nc_file.close()

    def close(self):
        self.nc_file.close()

    def uvToRegularGrid(self, model_index):
        """Interpolate u/v to regular grid"""
        # Extract the variables from the NetCDF
        ang_rho = self.nc_file.variables['angle'][1:,1:]
        lat_rho = self.nc_file.variables['lat_rho'][1:,1:]
        lon_rho = self.nc_file.variables['lon_rho'][1:,1:]
        u = self.nc_file.variables['u'][0,-1,1:,:]
        v = self.nc_file.variables['v'][0,-1,:,1:]
        mask_u = self.nc_file.variables['mask_u'][1:,:]
        mask_v = self.nc_file.variables['mask_v'][:,1:]
        mask_rho = self.nc_file.variables['mask_rho'][1:,1:]

        #Call masked arrays function and return masked arrays
        water_u,water_v,water_ang_rho,water_lat_rho,water_lon_rho = maskLand(u, v, ang_rho, lat_rho, lon_rho, mask_u, mask_v, mask_rho)
        
        #Call average to rho function u and v scalar values to rho
        u_rho, v_rho = averageUVToRho(water_u, water_v)
       
        #Call rotate function and return rotated u and v vectors
        rot_urho, rot_vrho = rotateUV2D(u_rho, v_rho, water_ang_rho)

        #Call create regular grid function 
        return interpolateUVToRegularGrid(water_lat_rho, water_lon_rho, rot_urho, rot_vrho, model_index)

def convertUVToSpeedDirection(directions, speeds, ugrid, vgrid):
    """Convert u and v averaged and rotated vectors to current speed and direction at regular grid point"""
    for x in range(ugrid.shape[0]):
        for y in range(ugrid.shape[1]):
            v_ms = vgrid[x,y]
            u_ms = ugrid[x,y]

            #Convert from meters per second to knots
            v_knot = v_ms * MS2KNOTS
            u_knot = u_ms * MS2KNOTS

            currentSpeed = math.sqrt(math.pow(u_knot, 2) + math.pow(v_knot, 2))
            currentDirectionRadians = math.atan2(v_knot, u_knot)
            currentDirectionDegrees = math.degrees(currentDirectionRadians)
            currentDirectionNorth = 90.0 - currentDirectionDegrees

            #The direction must always be positive.
            if currentDirectionNorth < 0.0:
                currentDirectionNorth += 360.0

            directions[x,y] = currentDirectionNorth
            speeds[x,y] = currentSpeed

    return directions, speeds

def rotateUV2D(u_rho, v_rho, water_ang_rho):
    """Rotate vectors by geometric angle"""
    angsin = numpy.sin(water_ang_rho)
    angcos = numpy.cos(water_ang_rho)
    rot_urho = u_rho*angcos - v_rho*angsin
    rot_vrho = u_rho*angsin + v_rho*angcos

    return (rot_urho, rot_vrho)
        
def interpolateUVToRegularGrid(water_lat_rho, water_lon_rho, rot_urho, rot_vrho, model_index):
    """Create a regular grid using masked latitude and longitude variables.
       Interpolate averaged, rotated, u/v variables to the regular grid.
       Using the NetCDF Index and Coeffiecient File to obtain weighted coefficents and 
       ROMs indicies to be used for Inverse Distance Weighting Interpolation.

    Args:
        water_lat_rho: 
        water_lon_rho: 
        rot_urho: 
        rot_vrho: 
        model_index: `ROMSIndexFile` from which index values/coefficients will
            be extracted to perform interpolation.
    """
    # Create regular grid using masked lat/long
    min_lat = numpy.nanmin(water_lat_rho)
    max_lat = numpy.nanmax(water_lat_rho)
    min_lon = numpy.nanmin(water_lon_rho)
    max_lon = numpy.nanmax(water_lon_rho)

    grid_x = numpy.linspace(min_lon, max_lon, 445)
    grid_y = numpy.linspace(min_lat, max_lat, 761)

    xgrid, ygrid = numpy.meshgrid(grid_x, grid_y)

    # Read the Index and Coefficient NetCDF variables          
    ncindex_xi1 = model_index.nc_file.variables['xi1'][:,:]
    ncindex_eta1 = model_index.nc_file.variables['eta1'][:,:]
    ncindex_w1 = model_index.nc_file.variables['w1'][:,:]
    ncindex_xi2 = model_index.nc_file.variables['xi2'][:,:]
    ncindex_eta2 = model_index.nc_file.variables['eta2'][:,:]
    ncindex_w2 = model_index.nc_file.variables['w2'][:,:]
    ncindex_xi3 = model_index.nc_file.variables['xi3'][:,:]
    ncindex_eta3 = model_index.nc_file.variables['eta3'][:,:]
    ncindex_w3 = model_index.nc_file.variables['w3'][:,:]
    ncindex_xi4 = model_index.nc_file.variables['xi4'][:,:]
    ncindex_eta4 = model_index.nc_file.variables['eta4'][:,:]
    ncindex_w4 = model_index.nc_file.variables['w4'][:,:]
    ncindex_wsum = model_index.nc_file.variables['wsum'][:,:]

    # Create masked empty regular grid for variable u
    ugrid = numpy.ma.empty(shape=[ncindex_xi1.shape[0],ncindex_xi1.shape[1]]) 
        
    # For each regular grid cell, read the corresponding xi1/eta1/w1,
    # xi2/eta2/w2, xi3/eta3/w3, and xi4/eta4/w4 values
    for y in range(ncindex_xi1.shape[0]):
        for x in range(ncindex_xi1.shape[1]):
            if not ncindex_xi1.mask[y,x]:
                xi1 = ncindex_xi1.data[y,x]
                eta1 = ncindex_eta1.data[y,x]
                xi2 = ncindex_xi2.data[y,x]
                eta2 = ncindex_eta2.data[y,x]
                xi3 = ncindex_xi3.data[y,x]
                eta3 = ncindex_eta3.data[y,x]
                xi4 = ncindex_xi4.data[y,x]
                eta4 = ncindex_eta4.data[y,x]
                u1 = rot_urho.data[eta1,xi1]
                u2 = rot_urho.data[eta2,xi2]
                u3 = rot_urho.data[eta3,xi3]
                u4 = rot_urho.data[eta4,xi4]
                # Use Inverse Distance Weigting algorithm to interpolate u to the regular grid
                ugrid[y,x] = ((((ncindex_w1.data[y,x]) * u1) + ((ncindex_w2.data[y,x]) * u2) + ((ncindex_w3[y,x]) * u3) + ((ncindex_w4[y,x])  * u4)) / ncindex_wsum[y,x])
    
    # Create masked empty regular grid for variable v        
    vgrid = numpy.ma.empty(shape=[ncindex_xi1.shape[0],ncindex_xi1.shape[1]])

    # For each regular grid cell, read the corresponding xi1/eta1/w1,
    # xi2/eta2/w2, xi3/eta3/w3, and xi4/eta4/w4 values
    for y in range(ncindex_xi1.shape[0]):
        for x in range(ncindex_xi1.shape[1]):
            if not ncindex_xi1.mask[y,x]:
                xi1 = ncindex_xi1.data[y,x]
                eta1 = ncindex_eta1.data[y,x]
                xi2 = ncindex_xi2.data[y,x]
                eta2 = ncindex_eta2.data[y,x]
                xi3 = ncindex_xi3.data[y,x]
                eta3 = ncindex_eta3.data[y,x]
                xi4 = ncindex_xi4.data[y,x]
                eta4 = ncindex_eta4.data[y,x]
                v1 = rot_vrho.data[eta1,xi1]
                v2 = rot_vrho.data[eta2,xi2]
                v3 = rot_vrho.data[eta3,xi3]
                v4 = rot_vrho.data[eta4,xi4]
                # Use Inverse Distance Weigting algorithm to interpolate v to
                # the regular grid
                vgrid[y,x] = ((((ncindex_w1.data[y,x]) * v1) + ((ncindex_w2.data[y,x]) * v2) + ((ncindex_w3[y,x]) * v3) + ((ncindex_w4[y,x])  * v4)) / ncindex_wsum[y,x])
                        

    return RegularGridUV(xgrid, ygrid, grid_x, grid_y, ugrid, vgrid, min_lon, min_lat, ncindex_xi1)

def averageUVToRho(water_u, water_v):
    """Average u and v scalars to rho"""
    # Average u values
    u_rho = numpy.ndarray([water_u.shape[0],water_u.shape[1]])
    
    for eta in range(water_u.shape[0]):
        for xi in range(water_u.shape[1]-1):
            u_rho[eta,xi] = (water_u[eta,xi] + water_u[eta,xi+1])/2
        u_rho[water_u.shape[0]-1,xi] = water_u[water_u.shape[0]-1,xi]  
        
    # Average v values
    v_rho = numpy.ndarray([water_v.shape[0],water_v.shape[1]])
    
    for xi in range(water_v.shape[1]):
        for eta in range(water_v.shape[0]-1):
            v_rho[eta,xi] = (water_v[eta,xi] + water_v[eta+1,xi])/2
        v_rho[eta,water_v.shape[1]-1] = water_v[eta,water_v.shape[1]-1]   

    return u_rho, v_rho

def maskLand(u, v, ang_rho, lat_rho, lon_rho, mask_u, mask_v, mask_rho):
    """Create masked arrays for specified variables to mask land values"""
    water_u = ma.masked_array(u, numpy.logical_not(mask_u))
    water_v = ma.masked_array(v, numpy.logical_not(mask_v))
    # u/v masked values need to be set to 0 for averaging 
    water_u = water_u.filled(0)
    water_v = water_v.filled(0)
    water_ang_rho = ma.masked_array(ang_rho, numpy.logical_not(mask_rho))
    water_lat_rho = ma.masked_array(lat_rho, numpy.logical_not(mask_rho))
    water_lon_rho = ma.masked_array(lon_rho, numpy.logical_not(mask_rho))  
           
    return water_u,water_v,water_ang_rho,water_lat_rho,water_lon_rho

