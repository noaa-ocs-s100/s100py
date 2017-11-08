#!/usr/bin/env python
# -*- coding: utf-8 -*-
#******************************************************************************
#
#******************************************************************************
import argparse
import h5py
import numpy
import pytz
import netCDF4
import math

ms2Knots = 1.943844

#******************************************************************************        
def create_xy_group(hdf_file, latc, lonc):
    """ Create the XY group containing the position information.

    :param hdf_file: The S-111 HDF file.
    :param latc: A list of latitude values.
    :param lonc: A list of longitude values.
    :returns: A tuple containing minimum x, minimum y, maximum x, maximum y values from the given lists.
    """

    numberOfLat = latc.shape[0]
    numberOfLon = lonc.shape[0]

    xCoordinates =  numpy.empty((1, numberOfLat), dtype=numpy.float64)
    yCoordinates = numpy.empty((1, numberOfLon), dtype=numpy.float64)
    minX = minY = maxX = maxY = None

    for index in range(0, numberOfLat):
        latitude = latc[index]
        longitude = lonc[index]

        #Keep track of the data extents so we can update the metadata.
        if minX == None:
            minX = maxX = longitude
            minY = maxY = latitude
        else:
            minX = min(minX, longitude)
            maxX = max(maxX, longitude)
            minY = min(minY, latitude)
            maxY = max(maxY, latitude)

        xCoordinates[0][index] = longitude
        yCoordinates[0][index] = latitude

    #Add the 'Group XY' to store the position information.
    groupName = 'Group XY'
    print("Creating", groupName, "dataset.")
    xy_group = hdf_file.create_group(groupName)

    #Add the x and y datasets to the xy group.
    xy_group.create_dataset('X', (1, numberOfLat), dtype=numpy.float64, data=xCoordinates)
    xy_group.create_dataset('Y', (1, numberOfLon), dtype=numpy.float64, data=yCoordinates)       

    return (minX, minY, maxX, maxY)
  
#******************************************************************************        
def create_direction_speed(newGroup, ur, vr):
    """ Create the speed and direction datasets.

    :param group: The HDF group to add the speed and direction datasets to.
    :param ua: List of velocity values along the x axis in metres per second.
    :param va: List of velocity values along the y axis in metres per second.
    :returns: A tuple containing the minimum and maximum speed values added.
    """

    min_speed = None
    max_speed = None

    numberOfVaValues = len(vr)

    directions = numpy.empty((1, numberOfVaValues), dtype=numpy.float64)
    speeds = numpy.empty((1, numberOfVaValues), dtype=numpy.float64)
    
    for index in range(0, numberOfVaValues):

        v_ms = vr[index]
        u_ms = ur[index]

        #Convert from metres per second to knots
        v_knot = v_ms * ms2Knots
        u_knot = u_ms * ms2Knots

        windSpeed = math.sqrt(math.pow(u_knot, 2) + math.pow(v_knot, 2))
        windDirectionRadians = math.atan2(v_knot, u_knot)
        windDirectionDegrees = math.degrees(windDirectionRadians)
        windDirectionNorth = 90.0 - windDirectionDegrees

        #The direction must always be positive.
        if windDirectionNorth < 0.0:
            windDirectionNorth += 360.0

        directions[0][index] = windDirectionNorth
        speeds[0][index] = windSpeed

        if min_speed == None:
            min_speed = max_speed = windSpeed
        else:
            min_speed = min(min_speed, windSpeed)
            max_speed = max(max_speed, windSpeed)

    #Create the datasets.
    direction_dataset = group.create_dataset('Direction', (1, numberOfVaValues), dtype=numpy.float64, data=directions)
    speed_dataset = group.create_dataset('Speed', (1, numberOfVaValues), dtype=numpy.float64, data=speeds)

    return min_speed, max_speed


#******************************************************************************        
def create_data_groups(hdf_file, times, ur, vr):
    """Create the data groups in the S-111 file. (One group for each time value)

    :param hdf_file: The S-111 HDF file.
    :param times: The list of time values from the source data.
    :param ua: List of velocity values along the x axis in metres per second. (An array of values per time)
    :param va: List of velocity values along the y axis in metres per second. (An array of values per time)
    :returns: A tuple containing the minimum time, maximum time, time interval, minimum speed, and maximum speed of the source data.
    """

    numberOfTimes = len(times) #times.shape[0]

    interval = None
    minTime = maxTime = None
    minSpeed = maxSpeed = None

    for index in range(0, numberOfTimes):

        newGroupName = 'Group ' + str(index + 1)
        print("Creating", newGroupName, "dataset.")
        newGroup = hdf_file.create_group(newGroupName)
        
        groupTitle = 'Irregular Grid at DateTime ' + str(index + 1)
        newGroup.attrs.create('Title', groupTitle.encode())

        #Store the start time.
        timeVal = times[index]

        #Keep track of the min/max time so we can update the metadata
        if minTime == None:
            minTime = maxTime = timeVal
        else:
            minTime = min(minTime, timeVal)
            maxTime = max(maxTime, timeVal)

        strVal = timeVal.strftime("%Y%m%dT%H%M%SZ")
        newGroup.attrs.create('DateTime', strVal.encode())

        groupMinSpeed, groupMaxSpeed = create_direction_speed(newGroup, ur, vr)

        #Keep track of the min/max speed so we can update the metadata
        if minSpeed == None:
            minSpeed = groupMinSpeed
            maxSpeed = groupMaxSpeed
        else:
            minSpeed = min(minSpeed, groupMinSpeed)
            maxSpeed = max(maxSpeed, groupMaxSpeed)
        
    #Figure out what the interval is between the times (use only the first)
    if numberOfTimes > 1:

        strVal = times[0].tostring().decode()
        firstTimeVal = iso8601.parse_date(strVal)
        firstTimeVal = firstTimeVal.astimezone(pytz.utc)

        strVal = times[1].tostring().decode()
        secondTimeVal = iso8601.parse_date(strVal)
        secondTimeVal = secondTimeVal.astimezone(pytz.utc)

        interval = secondTimeVal - firstTimeVal

    return (minTime, maxTime, interval, minSpeed, maxSpeed)


#******************************************************************************        
def update_metadata(hdf_file, numberOfTimes, numberOfValues, minTime, maxTime, interval, minX, minY, maxX, maxY, minSpeed, maxSpeed):
    """Update the S-111 file's metadata.

    :param hdf_file: The S-111 HDF file.
    :param numberOfTimes: The number of times in the source data.
    :param numberOfValues: The number of values per record in the source data.
    :param minTime: The minimum temporal extents of the source data.
    :param maxTime: The maximum temporal extents of the source data.
    :param interval: The time interval between records of the source data.
    :param minX: The minimum x coordinate of the source data.
    :param minY: The minimum y coordinate of the source data.
    :param maxX: The maximum x coordinate of the source data.
    :param maxY: The maximum y coordinate of the source data.
    :param minSpeed: The minimum surface speed of the source data.
    :param maxSpeed: The maximum surface speed of the source data.
    """

    #Set the correct coding format.
    hdf_file.attrs.create('dataCodingFormat', 3, dtype=numpy.int64)

    #Set the number of times.
    hdf_file.attrs.create('numberOfTimes', numberOfTimes, dtype=numpy.int64)

    #Set the number of nodes.
    hdf_file.attrs.create('numberOfNodes', numberOfValues, dtype=numpy.int64)
    
    #Set the time interval (if we have one)
    if interval != None:
        intervalInSeconds = interval.total_seconds()
        hdf_file.attrs.create('timeRecordInterval', intervalInSeconds, dtype=numpy.int64)

    #Update the temporal extents in the metadata.
    strVal = minTime.strftime("%Y%m%dT%H%M%SZ")
    hdf_file.attrs.create('dateTimeOfFirstRecord', strVal.encode())
    strVal = maxTime.strftime("%Y%m%dT%H%M%SZ")
    hdf_file.attrs.create('dateTimeOfLastRecord', strVal.encode())

    #Update the surface speed values.
    if 'minDatasetCurrentSpeed' in hdf_file.attrs:
        minSpeed = min(minSpeed, hdf_file.attrs['minDatasetCurrentSpeed'])

    if 'maxDatasetCurrentSpeed' in hdf_file.attrs:
        maxSpeed = max(maxSpeed, hdf_file.attrs['maxDatasetCurrentSpeed'])

    hdf_file.attrs.create('minDatasetCurrentSpeed', minSpeed)
    hdf_file.attrs.create('maxDatasetCurrentSpeed', maxSpeed)


#******************************************************************************        
def create_command_line():
    """Create and initialize the command line parser.
    
    :returns: The command line parser.
    """

    parser = argparse.ArgumentParser(description='Add S-111 irregular grid Dataset')

    parser.add_argument('-g', '--grid-file', help='The netcdf file containing the irregular grid data.', required=True)
    parser.add_argument("inOutFile", nargs=1)

    return parser

#****************************************************************************** 

def rot1d(u, v, ang):
    """
    rotate vectors by geometric angle 

    """
    angsin = numpy.sin(ang)
    angcos = numpy.cos(ang)
    ur = u*angcos - v*angsin
    vr = u*angsin + v*angcos
    return ur, vr 

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
            #Grab the data that we need.
            #Flatten data to a 1D array
            latc = grid_file.variables['lat_rho'][1:,1:].flatten()
            lonc = grid_file.variables['lon_rho'][1:,1:].flatten()
            latu = grid_file.variables['lat_u'][1:,:].flatten()
            lonu = grid_file.variables['lon_u'][1:,:].flatten()
            latv = grid_file.variables['lat_v'][:,1:].flatten()
            lonv = grid_file.variables['lon_v'][:,1:].flatten()
            ua = grid_file.variables['u'][0,-1,1:,:].flatten()
            va = grid_file.variables['v'][0,-1,:,1:].flatten()     
            ang = grid_file.variables['angle'][1:,1:].flatten()
            masku = grid_file.variables['mask_u'][1:,:].flatten()
            maskv = grid_file.variables['mask_v'][:,1:].flatten()
            maskc = grid_file.variables['mask_rho'][1:,1:].flatten()
            timesv = grid_file.variables['ocean_time']
            #convert gregorian timestamp to datetime timestamp
            times = netCDF4.num2date(timesv,timesv.units)

            ur, vr = rot1d(ua,va,ang) 

            #Verify that these arrays are the same size.
            numberOfTimes = len(times) #times.shape[0]
            numberOfVaSeries = va.ndim
            numberOfUaSeries = ua.ndim
            if numberOfTimes != numberOfVaSeries or numberOfTimes != numberOfUaSeries:
                raise Exception('The number of time values does not match the number of speed and distance values.')
        
            #Verify that these arrays are the same size.
            numberOfLat = latc.shape[0]
            numberOfLon = lonc.shape[0]
            numberOfVaValues = va.shape[0]
            numberOfUaValues = ua.shape[0]
            if numberOfLat != numberOfLon:
                raise Exception('The input latitude and longitude array are different sizes.')
            elif numberOfLat != numberOfVaValues or numberOfLat != numberOfUaValues:
                raise Exception('The number of positions does not match the number of speed and distance values.')

            print ("Number of VAvalues", numberOfVaValues)

            #Verify that the input data is in the correct units.
            v_units = grid_file.variables['v']
            u_units = grid_file.variables['u']
            vaUnits = v_units.getncattr('units')
            uaUnits = u_units.getncattr('units')
            if vaUnits != uaUnits and vaUnits != 'meter second-1':
                raise Exception('The input velocity data is stored in an unsupported unit.')

            print("Adding irregular grid dataset")
            print("Number of timestamps in source file:", numberOfTimes)
            print("Number of records for each timestamp:", numberOfLat)

            #Add the 'Group XY' to store the position information.
            minX, minY, maxX, maxY = create_xy_group(hdf_file, latc, lonc)
    
            #Add all of the groups
            minTime, maxTime, interval, minSpeed, maxSpeed = create_data_groups(hdf_file, times, ur, vr)

            #Update the s-111 file's metadata
            update_metadata(hdf_file, numberOfTimes, numberOfVaValues,
                            minTime, maxTime, interval, minX, minY, maxX, maxY,
                            minSpeed, maxSpeed)

            print("Dataset successfully added")

        #Flush any edits out.
        hdf_file.flush()


if __name__ == "__main__":
    main()
