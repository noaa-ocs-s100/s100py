#!/usr/bin/env python
# -*- coding: utf-8 -*-
#******************************************************************************
#
#******************************************************************************
import argparse
import h5py
import numpy
import netCDF4
import math
import datetime
import numpy.ma as ma
from scipy.interpolate import griddata

#******************************************************************************   

def convertVectors(directions, speeds, ugrid, vgrid):

            ms2Knots = 1.943844
            
            for x in range(ugrid.shape[0]):
                for y in range(ugrid.shape[1]):
                
                    v_ms = vgrid[x,y]
                    u_ms = ugrid[x,y]
    
                    #Convert from metres per second to knots
                    v_knot = v_ms * ms2Knots
                    u_knot = u_ms * ms2Knots
    
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
            
#******************************************************************************   

def createGroup(hdf_file, grid_file, ugrid, vgrid, Xv, Yv):

    ocean_time = grid_file.variables['ocean_time']
    #convert gregorian timestamp to datetime timestamp
    times = netCDF4.num2date(ocean_time,ocean_time.units)

 
    numberOfTimes = times.shape[0]

    if hdf_file.items() == []:
        
        for index in range(0, numberOfTimes):

            newGroupName = 'Group_' + str(index + 1)
            print("Creating", newGroupName, "dataset.")
            newGroup = hdf_file.create_group(newGroupName)
            
            groupTitle = 'Regular Grid at DateTime ' + str(index + 1)
            newGroup.attrs.create('Title', groupTitle.encode())

            timeVal = times[index]

            strVal = timeVal.strftime("%Y%m%dT%H%M%SZ")
            newGroup.attrs.create('DateTime', strVal.encode())
            
            #Create dataset containers for speed and direction , with input from function convertVectors
            directions = numpy.empty((ugrid.shape[0],ugrid.shape[1]), dtype=numpy.float64)
            speeds = numpy.empty((ugrid.shape[0],ugrid.shape[1]), dtype=numpy.float64)

            #Call function to convert u and v to current speed and direction
            directions, speeds = convertVectors(directions, speeds, ugrid, vgrid)
    
            #Write data to empty HDF5 datasets 
            direction_dataset = newGroup.create_dataset('Direction', (ugrid.shape[0],ugrid.shape[1]), dtype=numpy.float64, data=directions)
            speed_dataset = newGroup.create_dataset('Speed', (ugrid.shape[0],ugrid.shape[1]), dtype=numpy.float64, data=speeds)
            x_dataset = newGroup.create_dataset('X', (ugrid.shape[0],ugrid.shape[1]), dtype=numpy.float64, data=Xv)
            y_dataset = newGroup.create_dataset('Y', (ugrid.shape[0],ugrid.shape[1]), dtype=numpy.float64, data=Yv)

            #Add CF attributes and geographic coordinates
            direction_dataset.attrs['units'] = 'degrees'
            direction_dataset.attrs['long_name'] = 'surface_current_direction'
            vlen = h5py.special_dtype (vlen = str)
            direction_dataset.attrs.create ('coordinates', data = ['Y', 'X'],
            dtype=vlen)
            
            speed_dataset.attrs['units'] = 'knots'
            speed_dataset.attrs['long_name'] = 'surface_current_speed'
            vlen = h5py.special_dtype (vlen = str)
            speed_dataset.attrs.create ('coordinates', data = ['Y', 'X'],
            dtype=vlen)
            
            y_dataset.attrs["long_name"] = "latitude"
            y_dataset.attrs["units"] = "degrees_north"
            y_dataset.attrs["standard_name"] = "latitude"
            
            x_dataset.attrs["long_name"] = "longitude"
            x_dataset.attrs["units"] = "degrees_east" 
            x_dataset.attrs["standard_name"]= "longitude"
            
            

#******************************************************************************   

def createRegGrid(water_u,water_lat_u, water_lon_u, water_v, water_lat_v, water_lon_v):

    ugridX = numpy.linspace(min(water_lon_u), max(water_lon_u), 200)
    ugridY = numpy.linspace(min(water_lat_u), max(water_lat_u), 200)
    Xu, Yu = numpy.meshgrid(ugridX, ugridY)

    ucoords = numpy.column_stack((water_lon_u, water_lat_u))
    ugrid = griddata(ucoords, water_u, (Xu, Yu), method='nearest', fill_value=numpy.nan)
    #mask = numpy.all(numpy.isnan(ugrid), axis=1)
    #ugrid[~mask]
    #ugrid = ma.masked_where((ugrid < 0) | numpy.isnan(ugrid), ugrid)
    #ugrid = numpy.nan_to_num(ugrid)

    vgridX = numpy.linspace(min(water_lon_v), max(water_lon_v), 200)
    vgridY = numpy.linspace(min(water_lat_v), max(water_lat_v), 200)
    Xv, Yv = numpy.meshgrid(vgridX, vgridY)

    vcoords = numpy.column_stack((water_lon_v, water_lat_v))
    vgrid = griddata(vcoords, water_v, (Xv, Yv), method='nearest', fill_value=numpy.nan)
    #vgrid = ma.masked_where((vgrid < 0) | numpy.isnan(vgrid), vgrid)
    #vgrid = numpy.nan_to_num(vgrid)

    return (ugrid, vgrid, Xv, Yv)

#******************************************************************************   

def maskedArray(grid_file, rot_u, lat_u, lon_u, rot_v, lat_v, lon_v):

    mask_rho = grid_file.variables['mask_rho'][1:,1:].flatten()
    mask_u = grid_file.variables['mask_u'][1:,:].flatten()
    mask_v = grid_file.variables['mask_v'][:,1:].flatten()

    combine_mask= numpy.logical_not(numpy.logical_and(mask_u,mask_v))
   
    water_u = ma.compressed(ma.masked_array(rot_u,combine_mask))
    water_lat_u = ma.compressed(ma.masked_array(lat_u,combine_mask))
    water_lon_u = ma.compressed(ma.masked_array(lon_u,combine_mask))

    water_v = ma.compressed(ma.masked_array(rot_v,combine_mask))
    water_lat_v = ma.compressed(ma.masked_array(lat_v,combine_mask))
    water_lon_v = ma.compressed(ma.masked_array(lon_v,combine_mask))

    return (water_u, water_lat_u, water_lon_u, water_v, water_lat_v, water_lon_v)

#****************************************************************************** 

def rot2d(u, v, ang_rho):
    """rotate vectors by geometric angle"""

    angsin = numpy.sin(ang_rho)
    angcos = numpy.cos(ang_rho)
    rot_u = u*angcos - v*angsin
    rot_v = u*angsin + v*angcos

    return (rot_u, rot_v)
    

#******************************************************************************        
def create_command_line():
    """Create and initialize the command line parser.
    
    :returns: The command line parser.
    """

    parser = argparse.ArgumentParser(description='Add S-111 regular grid Dataset')

    parser.add_argument('-g', '--grid-file', help='The netcdf file containing the regular grid data.', required=True)
    parser.add_argument("inOutFile", nargs=1)

    return parser

#******************************************************************************        
def main():
    
    #Create the command line parser.
    parser = create_command_line()

    #Parse the command line.
    results = parser.parse_args()
    
    #open the HDF5 file.
    with h5py.File(results.inOutFile[0], "r+") as hdf_file:

        #Open the grid file.
        with netCDF4.Dataset(results.grid_file, "r", format="NETCDF3 Classic") as grid_file:
           
            #Grab the variables that we need
            #Flatten data to a 1D array

            ang_rho = grid_file.variables['angle'][1:,1:].flatten()
            lat_rho = grid_file.variables['lat_rho'][1:,1:].flatten()
            lon_rho = grid_file.variables['lon_rho'][1:,1:].flatten()
            lat_u = grid_file.variables['lat_u'][1:,:].flatten()
            lon_u = grid_file.variables['lon_u'][1:,:].flatten()
            lat_v = grid_file.variables['lat_v'][:,1:].flatten()
            lon_v = grid_file.variables['lon_v'][:,1:].flatten()
            u = grid_file.variables['u'][0,-1,1:,:].flatten()
            v = grid_file.variables['v'][0,-1,:,1:].flatten() 

            #Call rotate function and return rotated u and v vectors
            rot_u, rot_v = rot2d(u, v, ang_rho)

            #Call masked arrays function and return masked arrays
            water_u, water_lat_u, water_lon_u, water_v, water_lat_v, water_lon_v = maskedArray(grid_file,rot_u, lat_u, lon_u, rot_v, lat_v, lon_v)
            
            #Create regular grid function using masked u and v arrays
            ugrid, vgrid, Xv, Yv = createRegGrid(water_u,water_lat_u, water_lon_u, water_v, water_lat_v, water_lon_v)

            #Create HDF5 groups and datasets
            createGroup(hdf_file, grid_file, ugrid, vgrid, Xv, Yv)


            print("Data sucessfully added")

            #Flush any edits out.
            hdf_file.flush()
            #hdf_file.close()


if __name__ == "__main__":
    main()
