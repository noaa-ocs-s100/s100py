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

#******************************************************************************   

def updateAttributes(hdf_file, ugrid, gridX, gridY, minLon, minLat):
    """Update HDF5 attributes"""
    
    nodes = ugrid.flatten()
    gridSpacingLon = gridX[1]- gridX[0]
    gridSpacingLat = gridY[1]- gridY[0]
                    
    hdf_file.attrs.modify('gridSpacingLongitudinal',gridSpacingLon) 
    hdf_file.attrs.modify('gridSpacingLatitudinal',gridSpacingLat) 
    hdf_file.attrs.modify('gridSpacingLatitudinal',gridSpacingLat)
    hdf_file.attrs.modify('horizDatumValue', 4326) 
    hdf_file.attrs.modify('numPointsLongitudinal', ugrid.shape[0])
    hdf_file.attrs.modify('numPointsLatitudinal', ugrid.shape[1])
    hdf_file.attrs.modify('minGridPointLongitudinal', minLon)
    hdf_file.attrs.modify('minGridPointLatitudinal', minLat)
    hdf_file.attrs.modify('numberOfNodes', nodes.shape[0])
    hdf_file.attrs.modify('surfaceCurrentDepth', 2)
    hdf_file.attrs.modify('gridOriginLongitude', minLon)
    hdf_file.attrs.modify('gridOriginLatitude', minLat)
    hdf_file.attrs.modify('nameRegion', 'US_East_Coast')
    hdf_file.attrs.modify('nameSubregion', 'Cheaspeake_Bay')
    hdf_file.attrs.modify('methodCurrentsProduct', 'ROMS_Hydrodynamic_Model')
    hdf_file.attrs.modify('gridLandMaskValue', -9999.0 )
    hdf_file.attrs.modify('dataCodingFormat', 2)
    hdf_file.attrs.modify('depthTypeIndex', 2)

#******************************************************************************   

def convertVectors(directions, speeds, ugrid, vgrid):
    """Convert u and v averaged and rotated vectors to current speed and direction at regular grid point"""
    
    ms2Knots = 1.943844
    
    for x in range(ugrid.shape[0]):
        for y in range(ugrid.shape[1]):
        
            v_ms = vgrid[x,y]
            u_ms = ugrid[x,y]

            #Convert from meters per second to knots
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

def createGroup(hdf_file, netcdf_file, ugrid, vgrid, xgrid, ygrid, ncindex_xi1):
    """Create an inital HDF5 Group containing Speeds, Directions, and XY Datasets.
       For every additional NetCDF file Create an HDF5 Group containing Speeds and Direction Datasets.
       Update HDF5 time attributes.
    """
    #Read in time variable 
    ocean_time = netcdf_file.variables['ocean_time']
    
    #Convert gregorian timestamp to datetime timestamp
    times = netCDF4.num2date(ocean_time,ocean_time.units)

    numberOfTimes = times.shape[0]

    if hdf_file.items() == []:
        
        for index in range(0, numberOfTimes):
    
            newGroupName = 'Group' + ' ' + str(index + 1)
            print("Creating", newGroupName ,"dataset.")
            newGroup = hdf_file.create_group(newGroupName)
            
            groupTitle = 'Regular Grid at DateTime ' + str(index + 1)
            newGroup.attrs.create('Title', groupTitle.encode())
    
            timeVal = times[index]
    
            strVal = timeVal.strftime("%Y%m%dT%H%M%SZ")
            newGroup.attrs.create('DateTime', strVal.encode())                    
            hdf_file.attrs.modify('dateTimeOfIssue', strVal)                                 
            hdf_file.attrs.modify('dateTimeOfFirstRecord', strVal)
            hdf_file.attrs.modify('dateTimeOfLastRecord', strVal)
            hdf_file.attrs.modify('typeOfCurrentData', 6)
                                  
            #Create dataset containers for speed and direction , with input from function convertVectors
            directions = numpy.empty((vgrid.shape[0],vgrid.shape[1]), dtype=numpy.float32)
            speeds = numpy.empty((ugrid.shape[0],ugrid.shape[1]), dtype=numpy.float32)
            
            #Calculate Current Speed and Direction at Regular Grid Point
            directions, speeds = convertVectors(directions, speeds, ugrid, vgrid)    

            minSpeed = numpy.nanmin(speeds)
            maxSpeed = numpy.nanmax(speeds)
            minSpeed = numpy.round(minSpeed,2)
            maxSpeed = numpy.round(maxSpeed,2)
            directions = ma.masked_array(directions, ncindex_xi1.mask)  
            speeds = ma.masked_array(speeds, ncindex_xi1.mask)
            directions = directions.filled(-9999.0)
            speeds = speeds.filled(-9999.0)
            xgrid = ma.masked_array(xgrid, ncindex_xi1.mask)
            ygrid = ma.masked_array(ygrid, ncindex_xi1.mask)
            xgrid = xgrid.filled(-9999.0)
            ygrid = ygrid.filled(-9999.0)                    
            
            #Write data to empty HDF5 datasets 
            newGroup.create_dataset('surfaceCurrentDirection', (vgrid.shape[0],vgrid.shape[1]), dtype=numpy.float32, data=directions,  chunks=True, compression="gzip", fillvalue=-9999.0, compression_opts=9)
            newGroup.create_dataset('surfaceCurrentSpeed', (ugrid.shape[0],ugrid.shape[1]), dtype=numpy.float32, data=speeds,  chunks=True, compression="gzip", fillvalue=-9999.0, compression_opts=9)

            hdf_file.attrs.modify('minDatasetCurrentSpeed', minSpeed)
            hdf_file.attrs.modify('maxDatasetCurrentSpeed', maxSpeed)
            hdf_file.attrs.modify('numberOfTimes', numberOfTimes)
                    
    else:
            grps = hdf_file.items()
            numGrps = len(grps)
    
            newGroupName = 'Group' + ' ' + str(numGrps + 1)
            print("Creating", newGroupName, "dataset.")
            newGroup = hdf_file.create_group(newGroupName)
    
            groupTitle = 'Regular Grid at DateTime ' + str(numGrps + 1)
            newGroup.attrs.create('Title', groupTitle.encode())
    
            timeVal = times[0]
    
            strVal = timeVal.strftime("%Y%m%dT%H%M%SZ")
            newGroup.attrs.create('DateTime', strVal.encode())
    
            hdf_file.attrs.modify('dateTimeOfLastRecord', strVal)
            firstTime = datetime.datetime.strptime((hdf_file.attrs['dateTimeOfFirstRecord']),"%Y%m%dT%H%M%SZ")
            lastTime = datetime.datetime.strptime((hdf_file.attrs['dateTimeOfLastRecord']),"%Y%m%dT%H%M%SZ")
            interval = lastTime - firstTime
            timeInterval=  interval.total_seconds()
            hdf_file.attrs.modify('timeRecordInterval', timeInterval)
                
            #Create dataset containers for speed and direction , with input from function convertVectors
            directions = numpy.empty((vgrid.shape[0],vgrid.shape[1]), dtype=numpy.float32)
            speeds = numpy.empty((ugrid.shape[0],ugrid.shape[1]), dtype=numpy.float32)
            
            directions, speeds = convertVectors(directions, speeds, ugrid, vgrid) 
            
            minSpeed = numpy.nanmin(speeds)
            maxSpeed = numpy.nanmax(speeds)
            minSpeed = numpy.round(minSpeed,2)
            maxSpeed = numpy.round(maxSpeed,2)
            directions = ma.masked_array(directions, ncindex_xi1.mask)  
            speeds = ma.masked_array(speeds, ncindex_xi1.mask)
            directions = directions.filled(-9999.0)
            speeds = speeds.filled(-9999.0)
            xgrid = ma.masked_array(xgrid, ncindex_xi1.mask)
            ygrid = ma.masked_array(ygrid, ncindex_xi1.mask)
            xgrid = xgrid.filled(-9999.0)
            ygrid = ygrid.filled(-9999.0)

            #Write data to empty HDF5 datasets 
            newGroup.create_dataset('surfaceCurrentDirection', (vgrid.shape[0],vgrid.shape[1]), dtype=numpy.float32, data=directions,chunks=True, compression="gzip", fillvalue=-9999.0, compression_opts=9)
            newGroup.create_dataset('surfaceCurrentSpeed', (ugrid.shape[0],ugrid.shape[1]), dtype=numpy.float32, data=speeds, chunks=True, compression="gzip", fillvalue=-9999.0, compression_opts=9)

            numberOfTimes = numGrps + 1
            hdf_file.attrs.modify('numberOfTimes', numberOfTimes)
            prior_minSpeed = hdf_file.attrs['minDatasetCurrentSpeed']
            prior_maxSpeed = hdf_file.attrs['maxDatasetCurrentSpeed']
            if minSpeed < prior_minSpeed:
                hdf_file.attrs.modify('minDatasetCurrentSpeed', minSpeed)
            if maxSpeed > prior_maxSpeed:
                hdf_file.attrs.modify('maxDatasetCurrentSpeed', maxSpeed)
            

#******************************************************************************   

def interpolate2RegGrid(water_lat_rho,water_lon_rho,rot_urho, rot_vrho, index_file):
    """Create a regular grid using masked latitude and longitude variables.
       Interpolate averaged, rotated, u/v variables to the regular grid.
       Using the NetCDF Index and Coeffiecient File to obtain weighted coefficents and 
       ROMs indicies to be used for Inverse Distance Weighting Interpolation.
    """
    #Create regular grid using masked lat/long
    minLat = numpy.nanmin(water_lat_rho)
    maxLat = numpy.nanmax(water_lat_rho)
    minLon = numpy.nanmin(water_lon_rho)
    maxLon = numpy.nanmax(water_lon_rho)

    gridX = numpy.linspace(minLon, maxLon, 445)
    gridY = numpy.linspace(minLat, maxLat, 761)

    xgrid, ygrid = numpy.meshgrid(gridX, gridY)

    #Read the Index and Coefficient NetCDF variables          
    ncindex_xi1 = index_file.variables['xi1'][:,:]
    ncindex_eta1 = index_file.variables['eta1'][:,:]
    ncindex_w1 = index_file.variables['w1'][:,:]
    ncindex_xi2 = index_file.variables['xi2'][:,:]
    ncindex_eta2 = index_file.variables['eta2'][:,:]
    ncindex_w2 = index_file.variables['w2'][:,:]
    ncindex_xi3 = index_file.variables['xi3'][:,:]
    ncindex_eta3 = index_file.variables['eta3'][:,:]
    ncindex_w3 = index_file.variables['w3'][:,:]
    ncindex_xi4 = index_file.variables['xi4'][:,:]
    ncindex_eta4 = index_file.variables['eta4'][:,:]
    ncindex_w4 = index_file.variables['w4'][:,:]
    ncindex_wsum = index_file.variables['wsum'][:,:]

    #Create masked empty regular grid for variable u
    ugrid = numpy.ma.empty(shape=[ncindex_xi1.shape[0],ncindex_xi1.shape[1]]) 
        
    #For each regular grid cell, read the corresponding xi1/eta1/w1, xi2/eta2/w2, xi3/eta3/w3, and xi4/eta4/w4 values
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
                #Use Inverse Distance Weigting algorithm to interpolate u to the regular grid
                ugrid[y,x] = ((((ncindex_w1.data[y,x]) * u1) + ((ncindex_w2.data[y,x]) * u2) + ((ncindex_w3[y,x]) * u3) + ((ncindex_w4[y,x])  * u4)) / ncindex_wsum[y,x])
    
    #Create masked empty regular grid for variable v        
    vgrid = numpy.ma.empty(shape=[ncindex_xi1.shape[0],ncindex_xi1.shape[1]])

    #For each regular grid cell, read the corresponding xi1/eta1/w1, xi2/eta2/w2, xi3/eta3/w3, and xi4/eta4/w4 values
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
                #Use Inverse Distance Weigting algorithm to interpolate v to the regular grid
                vgrid[y,x] = ((((ncindex_w1.data[y,x]) * v1) + ((ncindex_w2.data[y,x]) * v2) + ((ncindex_w3[y,x]) * v3) + ((ncindex_w4[y,x])  * v4)) / ncindex_wsum[y,x])
                        

    return xgrid, ygrid, gridX, gridY, ugrid, vgrid, minLon, minLat, ncindex_xi1

#****************************************************************************** 

def rot2d(u_rho, v_rho, water_ang_rho):
    """Rotate vectors by geometric angle"""

    angsin = numpy.sin(water_ang_rho)
    angcos = numpy.cos(water_ang_rho)
    rot_urho = u_rho*angcos - v_rho*angsin
    rot_vrho = u_rho*angsin + v_rho*angcos

    return (rot_urho, rot_vrho)

#****************************************************************************** 

def avg2rho(water_u, water_v):
    """Averge u and v scalars to rho"""

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
    
#******************************************************************************     

def maskingLand(u, v, ang_rho, lat_rho, lon_rho, mask_u, mask_v, mask_rho):
    """Create masked arrays for specified variables to mask land values"""

    water_u = ma.masked_array(u, numpy.logical_not(mask_u))
    water_v = ma.masked_array(v, numpy.logical_not(mask_v))
    #u/v masked values need to be set to 0 for averaging 
    water_u = water_u.filled(0)
    water_v = water_v.filled(0)
    water_ang_rho = ma.masked_array(ang_rho, numpy.logical_not(mask_rho))
    water_lat_rho = ma.masked_array(lat_rho, numpy.logical_not(mask_rho))
    water_lon_rho = ma.masked_array(lon_rho, numpy.logical_not(mask_rho))  
           
    return water_u,water_v,water_ang_rho,water_lat_rho,water_lon_rho


#******************************************************************************  

def createCommandLine():
    """Create and initialize the command line parser.
       returns: The command line parser.
    """

    parser = argparse.ArgumentParser(description='Create S-111 Regular Grid Dataset')

    parser.add_argument('-f', '--netcdf_file', help='Add the ocean model netcdf file.', required=True)
    parser.add_argument('-i', '--index_file', help='Add index and coefficient netcdf file.', required=True)
    parser.add_argument("hdf_file", nargs=1)

    return parser

#******************************************************************************        
def main():
    
    #Create the command line parser
    parser = createCommandLine()

    #Parse the command line
    results = parser.parse_args()
    
    #Open the HDF5 file
    with h5py.File(results.hdf_file[0], "r+") as hdf_file:
        
        #Open the index and coefficient file
        with netCDF4.Dataset(results.index_file, "r", format="NETCDF4") as index_file:

            #Open the NetCDF file
            with netCDF4.Dataset(results.netcdf_file, "r", format="NETCDF3 Classic") as netcdf_file:
               
                #Extract the variables from the NetCDF
                ang_rho = netcdf_file.variables['angle'][1:,1:]
                lat_rho = netcdf_file.variables['lat_rho'][1:,1:]
                lon_rho = netcdf_file.variables['lon_rho'][1:,1:]
                u = netcdf_file.variables['u'][0,-1,1:,:]
                v = netcdf_file.variables['v'][0,-1,:,1:]
                mask_u = netcdf_file.variables['mask_u'][1:,:]
                mask_v = netcdf_file.variables['mask_v'][:,1:]
                mask_rho = netcdf_file.variables['mask_rho'][1:,1:]

                #Call masked arrays function and return masked arrays
                water_u,water_v,water_ang_rho,water_lat_rho,water_lon_rho = maskingLand(u, v, ang_rho, lat_rho, lon_rho, mask_u, mask_v, mask_rho)
                
                #Call average to rho function u and v scalar values to rho
                u_rho, v_rho = avg2rho(water_u, water_v)
               
                #Call rotate function and return rotated u and v vectors
                rot_urho, rot_vrho = rot2d(u_rho, v_rho, water_ang_rho)

                #Call create regular grid function 
                xgrid, ygrid, gridX, gridY, ugrid, vgrid, minLon, minLat, ncindex_xi1 = interpolate2RegGrid(water_lat_rho,water_lon_rho,rot_urho, rot_vrho, index_file)

                #Create HDF5 groups and datasets
                createGroup(hdf_file, netcdf_file, ugrid, vgrid, xgrid, ygrid, ncindex_xi1)

                #Update HDF5 attributes
                updateAttributes(hdf_file, ugrid, gridX, gridY, minLon, minLat)

                print("Data sucessfully added")

            #Close HDF file and flush edits
            hdf_file.flush()


if __name__ == "__main__":
    main()
