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

def indexAndcoefficients(mask_rho, lat_rho, lon_rho, xgrid, ygrid, var00, var01, var02, var03, var04, var05, var06, var07, var08, var09, var10, var11, var12, out_file):

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
                            print "found one!"
                            #Inverse-distance weight with the power of 1
                            #(1/numpy.sqrt((x1-x0)**2+(y1-y0)**2))
                            # Inverse-distance weight with the power of 2
                            w1 = (1/((x1-x0)**2+(y1-y0)**2))
                            w2 = (1/((x2-x0)**2+(y2-y0)**2))
                            w3 = (1/((x3-x0)**2+(y3-y0)**2))
                            w4 = (1/((x4-x0)**2+(y4-y0)**2))
                            wsum = w1 + w2 + w3 +w4
                            var00.xi1[y,x] = xi1
                            var01.eta1[y,x] = eta1
                            var02.w1[y,x] = w1
                            var03.xi2[y,x][y,x] = xi2
                            var04.eta2[y,x][y,x] = eta2
                            var05.w2[y,x][y,x] = w2
                            var06.xi3[y,x][y,x] = xi3
                            var07.eta3[y,x][y,x] = eta3
                            var08.w3[y,x][y,x] = w3
                            var09.xi4[y,x][y,x] = xi4
                            var10.eta4[y,x][y,x] = eta4
                            var11.w4[y,x][y,x] = w4
                            var12.wsum[y,x][y,x] = wsum
                            found_cell = True
                            break
                
#******************************************************************************        

def createNetCDF(mask_rho, lat_rho, lon_rho, xgrid, ygrid, out_file):

    # Run through every regular grid point test if it inside 4 roms/rho points then calculates the distance 
    # between the regular grid point and the four roms/rho points                                      
    
    out_file.createDimension('eta', len(xgrid))
    out_file.createDimension('xi', len(xgrid[1]))

    var00 = out_file.createVariable('xi1', 'i4', ('eta','xi'),fill_value=-9999)
    var01 = out_file.createVariable('eta1', 'i4', ('eta','xi'),fill_value=-9999)
    var02 = out_file.createVariable('w1', 'f4', ('eta','xi'),fill_value=-9999)
    var03 = out_file.createVariable('xi2', 'i4', ('eta','xi'),fill_value=-9999)
    var04 = out_file.createVariable('eta2', 'i4', ('eta','xi'),fill_value=-9999)
    var05 = out_file.createVariable('x2', 'f4', ('eta','xi'),fill_value=-9999)
    var06 = out_file.createVariable('xi3', 'i4', ('eta','xi'),fill_value=-9999) 
    var07 = out_file.createVariable('eta3', 'i4', ('eta','xi'),fill_value=-9999)
    var08 = out_file.createVariable('w3', 'f4', ('eta','xi'),fill_value=-9999)
    var09 = out_file.createVariable('xi4', 'i4', ('eta','xi'),fill_value=-9999)
    var10 = out_file.createVariable('eta4', 'i4', ('eta','xi'),fill_value=-9999)
    var11 = out_file.createVariable('w4', 'f4', ('eta','xi'),fill_value=-9999)
    var12 = out_file.createVariable('wsum', 'f4', ('eta','xi'),fill_value=-9999)
    
    var00.long_name = "roms_xi1_index"
    var00.standard_name = "xi1"

    var01.long_name = "roms_eta1_index"
    var01.standard_name = "eta1"

    var02.long_name = "weight_1"
    var02.standard_name = "w1"

    var03.long_name = "roms_xi2_index"
    var03.standard_name = "xi2"

    var04.long_name = "roms_eta2_index"
    var04.standard_name = "eta2"
    
    var05.long_name = "weight_2"
    var05.standard_name = "w2"

    var06.long_name = "roms_xi3_index"
    var06.standard_name = "xi3"

    var07.long_name = "roms_eta3_index"
    var07.standard_name = "eta3"
    
    var08.long_name = "weight_3"
    var08.standard_name = "w3"
    
    var09.long_name = "roms_xi4_index"
    var09.standard_name = "xi4"

    var10.long_name = "roms_eta4_index"
    var10.standard_name = "eta4"

    var11.long_name = "weight_4"
    var11.standard_name = "w4"
    
    var12.long_name = "weight_summed"
    var12.standard_name = "w_sum"

    out_file.format = "netCDF-4" 
    out_file.conventions = "CF-1.0"


    indexAndcoefficients(mask_rho, lat_rho, lon_rho, xgrid, ygrid, var00, var01, var02, var03, var04, var05, var06, var07, var08, var09, var10, var11, var12, out_file)

#******************************************************************************        
def create_command_line():
    """Create and initialize the command line parser.
    
    :returns: The command line parser.
    """

    parser = argparse.ArgumentParser(description='Add netCDF File')
    
    parser.add_argument('-i', '--nc-file', help='The netcdf file containing data.', required=True)
    parser.add_argument("outputFile", nargs=1)

    return parser

#******************************************************************************        
def main():
    
    #Create the command line parser.
    parser = create_command_line()

    #Parse the command line.
    results = parser.parse_args()
    
        #create output netCDF file.
    with netCDF4.Dataset(results.outputFile[0], "w", format="NETCDF4") as out_file: 

        #Open netCDF file.
        with netCDF4.Dataset(results.nc_file, "r", format="NETCDF3 Classic") as nc_file:


            #Extract Variables
            lat_rho = nc_file.variables['lat_rho'][1:,1:]
            lon_rho = nc_file.variables['lon_rho'][1:,1:]
            mask_rho = nc_file.variables['mask_rho'][1:,1:] 
            
            water_lat_rho = ma.masked_array(lat_rho, numpy.logical_not(mask_rho))
            water_lon_rho = ma.masked_array(lon_rho, numpy.logical_not(mask_rho))  
           
            #get lat/lon extent
            lat_min = numpy.nanmin(water_lat_rho)
            lat_max = numpy.nanmax(water_lat_rho)
            lon_min = numpy.nanmin(water_lon_rho)
            lon_max = numpy.nanmax(water_lon_rho)
            
            #Then Create regular Grid 
            #500 meters
            gridX = numpy.linspace(lon_min, lon_max, 445)
            gridY = numpy.linspace(lat_min, lat_max, 761)

            xgrid, ygrid = numpy.meshgrid(gridX, gridY)

            #Call function create NetCDF
            createNetCDF(mask_rho, lat_rho, lon_rho, xgrid, ygrid, out_file)

            print("Index and Coefficients Sucessfully Created")

    out_file.close()

if __name__ == "__main__":
    main()