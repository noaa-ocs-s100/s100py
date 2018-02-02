#!/usr/bin/env python
# -*- coding: utf-8 -*-
#******************************************************************************
#
#******************************************************************************

import numpy
import netCDF4
import numpy.ma as ma
import argparse

#******************************************************************************

def indexAndcoefficients(mask_rho, lat_rho, lon_rho, xgrid, ygrid, var00, var01, var02, var03, var04, var05, var06, var07, var08, var09, var10, var11, var12, index_file):
    """Add data to output NetCDF variables. 
       For every regular grid point x0,y0 which falls inside four valid irregular grid points, store (eta1,xi),(eta2,xi2),(eta3,xi3),(eta4,xi4) and (w1,w2,w3,w4).
       If all four triangle areas are greater than 0, the point (x0,y0) is located inside the irregular grid cell enclosed by (eta1,xi),(eta2,xi2),(eta3,xi3),(eta4,xi4).
       If it does not fallen within any irregular grid cell, this locations are specified as missing values (-99999.0). If it does fall within a valid irregular cell, 
       calculate the inverse-distance weight to from the four four (eta,xi) locations to the each regular grid point (x0,y0). 
    """
    
    for y in range(ygrid.shape[0]):
        for x in range(xgrid.shape[1]):
            x0 = xgrid[y,x]
            y0 = ygrid[y,x]
            found_cell = False
            for xi1 in range(lat_rho.shape[1]-1):
                if found_cell:
                    break
                for eta1 in range(lat_rho.shape[0]-1):
                    xi2 = xi1 + 1
                    eta2 = eta1
                    xi3 = xi1 + 1
                    eta3 = eta1 + 1
                    xi4 = xi1
                    eta4 = eta1 + 1
                    if mask_rho[eta1,xi1] == 1 and mask_rho[eta2,xi2] == 1 and mask_rho[eta3,xi3] == 1 and mask_rho[eta4,xi4] == 1:
                        x1 = lon_rho[eta1,xi1]
                        y1 = lat_rho[eta1,xi1]
                        x2 = lon_rho[eta2,xi2]
                        y2 = lat_rho[eta2,xi2]
                        x3 = lon_rho[eta3,xi3]
                        y3 = lat_rho[eta3,xi3]
                        x4 = lon_rho[eta4,xi4]
                        y4 = lat_rho[eta4,xi4]
                        s1=0.5*((x1-x0)*(y2-y1)-(x2-x1)*(y1-y0))     
                        s2=0.5*((x2-x0)*(y3-y2)-(x3-x2)*(y2-y0))
                        s3=0.5*((x3-x0)*(y4-y3)-(x4-x3)*(y3-y0))     
                        s4=0.5*((x4-x0)*(y1-y4)-(x1-x4)*(y4-y0))                                        
                        if all([s > 0 for s in [s1,s2,s3,s4]]):
                            #Inverse-distance weight with the power of 1 - (1/numpy.sqrt((x1-x0)**2+(y1-y0)**2))
                            #Inverse-distance weight with the power of 2 - (1/((x1-x0)**2+(y1-y0)**2))
                            w1 = (1/((x1-x0)**2+(y1-y0)**2))
                            w2 = (1/((x2-x0)**2+(y2-y0)**2))
                            w3 = (1/((x3-x0)**2+(y3-y0)**2))
                            w4 = (1/((x4-x0)**2+(y4-y0)**2))
                            wsum = w1 + w2 + w3 +w4
                            var00[y,x] = xi1
                            var01[y,x] = eta1
                            var02[y,x] = w1
                            var03[y,x] = xi2
                            var04[y,x] = eta2
                            var05[y,x] = w2
                            var06[y,x] = xi3
                            var07[y,x] = eta3
                            var08[y,x] = w3
                            var09[y,x] = xi4
                            var10[y,x] = eta4
                            var11[y,x] = w4
                            var12[y,x] = wsum
                            found_cell = True
                            break    


#******************************************************************************        

def createNetCDF(mask_rho, lat_rho, lon_rho, xgrid, ygrid, index_file):
    """Add attributes and variables to output NetCDF file"""                              
    
    index_file.createDimension('eta', len(xgrid))
    index_file.createDimension('xi', len(xgrid[1]))

    var00 = index_file.createVariable('xi1', 'i4', ('eta','xi'),fill_value=-9999)
    var01 = index_file.createVariable('eta1', 'i4', ('eta','xi'),fill_value=-9999)
    var02 = index_file.createVariable('w1', 'f4', ('eta','xi'),fill_value=-9999)
    var03 = index_file.createVariable('xi2', 'i4', ('eta','xi'),fill_value=-9999)
    var04 = index_file.createVariable('eta2', 'i4', ('eta','xi'),fill_value=-9999)
    var05 = index_file.createVariable('w2', 'f4', ('eta','xi'),fill_value=-9999)
    var06 = index_file.createVariable('xi3', 'i4', ('eta','xi'),fill_value=-9999) 
    var07 = index_file.createVariable('eta3', 'i4', ('eta','xi'),fill_value=-9999)
    var08 = index_file.createVariable('w3', 'f4', ('eta','xi'),fill_value=-9999)
    var09 = index_file.createVariable('xi4', 'i4', ('eta','xi'),fill_value=-9999)
    var10 = index_file.createVariable('eta4', 'i4', ('eta','xi'),fill_value=-9999)
    var11 = index_file.createVariable('w4', 'f4', ('eta','xi'),fill_value=-9999)
    var12 = index_file.createVariable('wsum', 'f4', ('eta','xi'),fill_value=-9999)

    index_file.format = "netCDF-4" 


    return mask_rho, lat_rho, lon_rho, xgrid, ygrid, var00, var01, var02, var03, var04, var05, var06, var07, var08, var09, var10, var11, var12, index_file


#******************************************************************************        
def create_command_line():
    """Create and initialize the command line parser.
       Outputfile: Designate an output NetCDF file.
    """
    parser = argparse.ArgumentParser(description='Add netCDF File')
    
    parser.add_argument('-f', '--netcdf_file', help='Add the ocean model netcdf file.', required=True)
    parser.add_argument('-i', '--index_file', help='Add index and coefficient netcdf file.', required=True)

    return parser


#******************************************************************************        
def main():
    """ This script scans a set of regular grid points and determines whether a regular grid point falls within the boundaries 
        of four valid adjancent irregular grid points, forming an irregular grid cell. Then calculates the inverse distance from each 
        regular grid point to the pre determined four irregular grid points and stores the four indicies and four weighted coefficients 
        for each valid regular grid point. These indicies and weighted coefficent will be used to interpolate irregular ROMS variables 
        to a regular grid. 
    """

    #Create the command line parser.
    parser = create_command_line()

    #Parse the command line.
    results = parser.parse_args()
    
        #create output netCDF file.
    with netCDF4.Dataset(results.index_file, "w", format="NETCDF4") as index_file: 

        #Open netCDF file.
        with netCDF4.Dataset(results.netcdf_file, "r", format="NETCDF3 Classic") as netcdf_file:


            #Extract Variables
            lat_rho = netcdf_file.variables['lat_rho'][1:,1:]
            lon_rho = netcdf_file.variables['lon_rho'][1:,1:]
            mask_rho = netcdf_file.variables['mask_rho'][1:,1:] 
            
            water_lat_rho = ma.masked_array(lat_rho, numpy.logical_not(mask_rho))
            water_lon_rho = ma.masked_array(lon_rho, numpy.logical_not(mask_rho))  
           
            #get lat/lon extent
            lat_min = numpy.nanmin(water_lat_rho)
            lat_max = numpy.nanmax(water_lat_rho)
            lon_min = numpy.nanmin(water_lon_rho)
            lon_max = numpy.nanmax(water_lon_rho)
            
            #Create regular grid 
            #500 meters
            gridX = numpy.linspace(lon_min, lon_max, 445)
            gridY = numpy.linspace(lat_min, lat_max, 761)

            xgrid, ygrid = numpy.meshgrid(gridX, gridY)

            #Call function create NetCDF
            mask_rho, lat_rho, lon_rho, xgrid, ygrid, var00, var01, var02, var03, var04, var05, var06, var07, var08, var09, var10, var11, var12, index_file = createNetCDF(mask_rho, lat_rho, lon_rho, xgrid, ygrid, index_file)
            
            #Call function create index and coefficients
            indexAndcoefficients(mask_rho, lat_rho, lon_rho, xgrid, ygrid, var00, var01, var02, var03, var04, var05, var06, var07, var08, var09, var10, var11, var12, index_file)

            print("Index and Coefficients Sucessfully Created")


if __name__ == "__main__":
    main()
