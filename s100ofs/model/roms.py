"""
Utility classes and methods for working with ROMS ocean model output.

The Regional Ocean Modeling System (ROMS) is a 3D hydrodynamic modeling
framework which uses an irregular, curvilinear horizontal grid and a sigma
(bathymetry-following) vertical coordinate system. This module provides
functionality allowing ROMS output to be interpolated to a regular, orthogonal
lat/lon horizontal grid at a given depth-below-surface.
"""
import datetime
import math
import os

import netCDF4
import numpy
import numpy.ma as ma

import gdal
import osr
import ogr
from shapely.geometry import Polygon, Point, MultiPolygon, shape

# Conversion factor for meters/sec to knots
MS2KNOTS = 1.943844
_package = "s100ofs"
_dirPath = os.path.dirname(os.path.realpath(_package))


class RegularGrid:
    """Create a regular grid from an irregular grid.

        This class reads data from a ROMS model output file and creates a regular grid
        at specified resolution. Using either the ocean model extent or a shapefile.  

    """
    
    def __init__(self, path):
        """Initialize ROMSOutputFile object and open file at specified path.

        Args:
            path: Path of target NetCDF file.
        
        Raises:
            Exception: Specified NetCDF file does not exist.
        """
        self.path = path
        if os.path.exists(self.path):
            self.nc_file = netCDF4.Dataset(self.path, "r", format="NETCDF3_Classic")
            self.init_handles()
        else:
            # File doesn't exist, raise error
            raise(Exception("NetCDF file does not exist: {}".format(self.path)))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self.nc_file.close()

    def init_handles(self):
        """Initialize handles to NetCDF variables.

        Because there is no U for the first column (xi=0) and no V for the
        first row (eta=0) in the ROMS staggered grid, we must skip the first
        row and column of every variable associated with rho points as well as
        the first row (eta=0) of U values, and the first column (xi=0) of V
        values. This ensures that u, v, and all other rho variables have the
        same dimensions in xi and eta, with u[eta,xi], v[eta,xi], and rho
        variables e.g. ang_rho[eta,xi] all correspond with the same grid cell
        for a given [eta,xi] coordinate.
        """
        self.lat_rho = self.nc_file.variables['lat_rho'][1:,1:]
        self.lon_rho = self.nc_file.variables['lon_rho'][1:,1:]
        self.mask_rho = self.nc_file.variables['mask_rho'][1:,1:]

    def init_grid(self, res):
        """Initialize grid at specified resolution.

        Args:
            res: Resolution of grid. 
            option: Shapefile Extent.

        Create regular grid and return dimensions, using default ocean model
        extent or using an extent derived from shapefile coverage. 
        """
        self.res = res
        print (res)

        num_cells_y, num_cells_x = ofsExtent(self,res)

def ofsExtent(self, res):
    """For any ocean model extent, determine resolution and create regular grid.

    """
    #Determine ocean model extent
    water_lat_rho = ma.masked_array(self.lat_rho, numpy.logical_not(self.mask_rho))
    water_lon_rho = ma.masked_array(self.lon_rho, numpy.logical_not(self.mask_rho))
    minLon = numpy.nanmin(water_lon_rho)
    maxLon = numpy.nanmax(water_lon_rho)
    minLat = numpy.nanmin(water_lat_rho)
    maxLat = numpy.nanmax(water_lat_rho)

    dist_x = abs(maxLon - minLon)
    dist_y = abs(maxLat - minLat)

    #Midpoint constant 
    #The number of meteres per degree of latitude at the equator
    #Length of 1 degree of Longitude = cosine (latitude) * length of degree (meters) at equator 
    midLat = (minLat + (maxLat-minLat))/2
    radius= 6371000 # Radius of earth in meters 
    c = ((radius * math.pi)/180)*math.cos(midLat)
    
    # Target resolution in decimal degrees
    target_dd = res/c
    
    num_cells_x = int(round(dist_x/target_dd))
    cellsize_x =  dist_x/num_cells_x

    num_cells_y = int(round(dist_y/target_dd))
    cellsize_y = dist_y/num_cells_y
    
    #Creates a regular grid 
    step_x = cellsize_x/2
    step_y = cellsize_y/2
    x = minLon + step_x
    y = minLat + step_y
    
    gridX = []
    while x <= maxLon:
        gridX.append(x)
        x += cellsize_x
    
    gridY = []
    while y <= maxLat:
        gridY.append(y)
        y += cellsize_y
    
    print (num_cells_y)
    print (num_cells_x) 

    return num_cells_y, num_cells_x

def shpExtent(self, res, minLon,maxLon, minLat, maxLat):
    """For any shapefile extent, determine resolution and create regular grid.

    """
    #Read in 160k grid 
    _shpfile="grids_160k.shp"
    shp_path = os.path.realpath(_shpfile)

    shp = ogr.Open(shp_path)
    layer = shp.GetLayer()   

    #Read in 160k grid geometry
    #Export to Json to use FID atttribute
    grid160k_list = []

    for feature in layer:
        geom = feature.GetGeometryRef()
        grid160k_list.append(geom.ExportToJson())

    #Create OGR Geometry from ocean model grid extent   
    #Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(self.minLon, self.maxLat)
    ring.AddPoint(self.minLon, self.minLat)
    ring.AddPoint(self.maxLon, self.minLat)
    ring.AddPoint(self.maxLon, self.maxLat)
    ring.AddPoint(self.minLon, self.maxLat)
    # Create polygon
    ofs_poly = ogr.Geometry(ogr.wkbPolygon)
    ofs_poly.AddGeometry(ring)

    #Find the intersection between 160k grid and ocean model grid extent
    intersect = []
    index = []
    #intersection
    for polygon in range (len(grid160k_list)):
        grid_poly = ogr.CreateGeometryFromJson(grid160k_list[polygon])
        intersection = ofs_poly.Intersection(grid_poly)                
        if intersection.ExportToWkt() != "GEOMETRYCOLLECTION EMPTY":
            intersect.append(intersection.ExportToJson()) 
            index.append(polygon)

    #Create a new polygon with the grid_160k cells that cover the ocean model extent    
    multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
    for id in index: 
        multipolygon.AddGeometry(ogr.CreateGeometryFromJson(grid160k_list[id]))

    #Get coverage extent    
    (minX, maxX, minY, maxY) = multipolygon.GetEnvelope()
    
    #Coverage Extent, X distance & Y distance
    dist_x = abs(maxX - minX)
    dist_y = abs(maxY - minY)

    #Midpoint constant 
    #The number of meteres per degree of latitude at the equator
    #Length of 1 degree of Longitude = cosine (latitude) * length of degree (meters) at equator 
    midLat = (minY + (maxY-minY))/2
    radius= 6371000 # Radius of earth in meters 
    c = ((radius * math.pi)/180)*cos(midLat)
    
    # Target resolution in decimal degrees
    target_dd = self.res/c
    
    num_cells_x = int(round(dist_x/target_dd))
    cellsize_x =  dist_x/num_cells_x

    num_cells_y = int(round(dist_y/target_dd))
    cellsize_y = dist_y/num_cells_y
    
    #Creates a regular grid 
    step_x = cellsize_x/2
    step_y = cellsize_y/2
    x = minX + step_x
    y = minY + step_y
    
    gridX = []
    while x <= maxX:
        gridX.append(x)
        x += cellsize_x
    
    gridY = []
    while y <= maxY:
        gridY.append(y)
        y += cellsize_y

    print (num_cells_y)  
    return num_cells_y, num_cells_x
    
def gridSubset(self, num_cells_y, num_cells_x, xgrid, ygrid):
    """Use 160k grids to subset ocean model output.

    """    
    #subset_masks
    subset_mask = numpy.ma.empty(shape=[num_cells_y,num_cells_x],  dtype=int) 
     
    #Run through intersected polygons and test every regular grid point
    #If the regular grid point falls within an intersected polygon assign it 
    #the fid value
    for id in index: 
        grid_160k = shape(json.loads(poly_list[id]))             
        for eta in range (num_cells_y):
            for xi in range (num_cells_x):
                point = Point(xgrid[eta,xi], ygrid[eta,xi])                    
                if point.within(grid_160k) == True:
                   subset_mask[eta,xi] = id      


def gridShoreline(self, num_cells_y, num_cells_x, xgrid, ygrid): 
    """Use shoreline shp to mask out land at highest resolution.

    """   
    _shoreline ="shoreline.shp"
    shr_path = os.path.realpath(_shoreline)

    shp = ogr.Open(shr_path)
    layer = shp.GetLayer() 

    for feature in layer:
        geom = feature.GetGeometryRef()
        shore_poly = geom.ExportToJson()

    #Test every regular grid point, if it falls within the shoreline polygon
    #Assign it a value of 1
    shoreline_mask = numpy.ma.empty(shape=[num_cells_y, num_cells_x],  dtype=int)     
    
    water = shape(json.loads(shore_poly))  
    for eta in range (num_cells_y):
        for xi in range (num_cells_x):
            point = Point(xgrid[eta,xi], ygrid[eta,xi])                    
            if point.within(water) == True:
               shoreline_mask[eta,xi] = 1
            else:
               shoreline_mask[eta,xi] = -9999.0

class ROMSIndexFile(RegularGrid):
    """Store information about an index file used during interpolation.
    
    Attributes:
        path: Path (relative or absolute) of the file..
        nc_file: Handle to `netCDF4.Dataset` instance for the opened NetCDF
            file.
        dim_x: Handle to x dimension.
        dim_y: Handle to y dimension.
        var_x: Handle to x coordinate variable (longitudes).
        var_y: Handle to y coordinate variable (latitudes).
        var_xi1: Handle to xi1 variable (xi index of point 1).
        var_eta1: Handle to eta1 variable (eta index of point 1).
        var_w1: Handle to w1 variable (weight coefficient of point 1).
        var_xi2: Handle to xi2 variable (xi index of point 2).
        var_eta2: Handle to eta2 variable (eta index of point 2).
        var_w2: Handle to w2 variable (weight coefficient of point 2).
        var_xi3: Handle to xi3 variable (xi index of point 3).
        var_eta3: Handle to eta3 variable (eta index of point 3).
        var_w3: Handle to w3 variable (weight coefficient of point 3).
        var_xi4: Handle to xi4 variable (xi index of point 4).
        var_eta4: Handle to eta4 variable (eta index of point 4).
        var_w1: Handle to w4 variable (weight coefficient of point 4).
        var_wsum: Handle to wsum variable (sum of weight coefficients 1-4).
    """
    def __init__(self, path, clobber=False):
        """Initialize ROMSIndexFile object and open file at specified path.

        If file already exists and `clobber==False`, it will be opened in
        read mode. Otherwise, it will be opened in write mode.
        
        Args:
            path: Path of target NetCDF file.
            clobber: (Optional, default False) If True, existing index file at
                specified path, if any, will be deleted and the new file will
                be opened in write mode.
        """
        self.path = path
        if not os.path.exists(self.path) or clobber:
            self.nc_file = netCDF4.Dataset(self.path, "w", format="NETCDF4")
        else:
            self.nc_file = netCDF4.Dataset(self.path, "r", format="NETCDF4")
            self.init_handles()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self.nc_file.close()

    def init_handles(self):
        """Initialize handles to NetCDF dimensions/variables."""
        self.dim_y = self.nc_file.dimensions['y']
        self.dim_x = self.nc_file.dimensions['x']

        self.var_y = self.nc_file.variables['y'][:]
        self.var_x = self.nc_file.variables['x'][:]

        self.var_xi1 = self.nc_file.variables['xi1'][:,:]
        self.var_eta1 = self.nc_file.variables['eta1'][:,:]
        self.var_w1 = self.nc_file.variables['w1'][:,:]
        self.var_xi2 = self.nc_file.variables['xi2'][:,:]
        self.var_eta2 = self.nc_file.variables['eta2'][:,:]
        self.var_w2 = self.nc_file.variables['w2'][:,:]
        self.var_xi3 = self.nc_file.variables['xi3'][:,:]
        self.var_eta3 = self.nc_file.variables['eta3'][:,:]
        self.var_w3 = self.nc_file.variables['w3'][:,:]
        self.var_xi4 = self.nc_file.variables['xi4'][:,:]
        self.var_eta4 = self.nc_file.variables['eta4'][:,:]
        self.var_w4 = self.nc_file.variables['w4'][:,:]
        self.var_wsum = self.nc_file.variables['wsum'][:,:]
    
    def init_nc(self, roms_file):
        """Initialize NetCDF dimensions/variables/attributes.

        Args:
            roms_file: `ROMSOutputFile` instance containing model output used
                to identify properties of original grid.
            num_cells_x: Number of cells in x dimension of regular grid.
            num_cells_y: Number of cells in y dimension of regular grid.
        """
        self.dim_y = self.nc_file.createDimension('y', RegularGrid.num_cells_y)
        self.dim_x = self.nc_file.createDimension('x', RegularGrid.num_cells_x)

        self.var_y = self.nc_file.createVariable('y', 'f4', ('y',), fill_value=-9999)
        self.var_y.long_name = "latitude of regular grid point"
        self.var_y.units = "degree_north"
        self.var_y.standard_name = "latitude"

        self.var_x = self.nc_file.createVariable('x', 'f4', ('x',), fill_value=-9999)
        self.var_x.long_name = "longitude of regular grid point"
        self.var_x.units = "degree_east"
        self.var_x.standard_name = "longitude"

        # Write out lat/lon coordinates to y/x coordinate variables
        water_lat_rho = ma.masked_array(roms_file.var_lat_rho, numpy.logical_not(roms_file.var_mask_rho))
        water_lon_rho = ma.masked_array(roms_file.var_lon_rho, numpy.logical_not(roms_file.var_mask_rho))
        lon_min = numpy.nanmin(water_lon_rho)
        lon_max = numpy.nanmax(water_lon_rho)
        lat_min = numpy.nanmin(water_lat_rho)
        lat_max = numpy.nanmax(water_lat_rho)
        self.var_x[:] = numpy.linspace(lon_min, lon_max, num_cells_x)
        self.var_y[:] = numpy.linspace(lat_min, lat_max, num_cells_y)

        self.var_xi1 = self.nc_file.createVariable('xi1', 'i4', ('y','x'),fill_value=-9999)
        self.var_eta1 = self.nc_file.createVariable('eta1', 'i4', ('y','x'),fill_value=-9999)
        self.var_w1 = self.nc_file.createVariable('w1', 'f4', ('y','x'),fill_value=-9999)
        self.var_xi2 = self.nc_file.createVariable('xi2', 'i4', ('y','x'),fill_value=-9999)
        self.var_eta2 = self.nc_file.createVariable('eta2', 'i4', ('y','x'),fill_value=-9999)
        self.var_w2 = self.nc_file.createVariable('w2', 'f4', ('y','x'),fill_value=-9999)
        self.var_xi3 = self.nc_file.createVariable('xi3', 'i4', ('y','x'),fill_value=-9999)
        self.var_eta3 = self.nc_file.createVariable('eta3', 'i4', ('y','x'),fill_value=-9999)
        self.var_w3 = self.nc_file.createVariable('w3', 'f4', ('y','x'),fill_value=-9999)
        self.var_xi4 = self.nc_file.createVariable('xi4', 'i4', ('y','x'),fill_value=-9999)
        self.var_eta4 = self.nc_file.createVariable('eta4', 'i4', ('y','x'),fill_value=-9999)
        self.var_w4 = self.nc_file.createVariable('w4', 'f4', ('y','x'),fill_value=-9999)
        self.var_wsum = self.nc_file.createVariable('wsum', 'f4', ('y','x'),fill_value=-9999)

        #self.gridOriginLongitude = 
        #self.gridOriginLatitude = 
        #self.gridSpacingLongitude = 
        #self.gridSpacingLatitude = 
        self.nc_file.model = "CBOFS"
        self.nc_file.format = "netCDF-4"

    def compute_indexes_coefficients(self, roms_file):
        """Compute index/coefficient variables and store in NetCDF file.
        
        For every regular grid point x0,y0 which falls inside four valid
        irregular grid points, store (eta1,xi), (eta2,xi2), (eta3,xi3),
        (eta4,xi4) and (w1,w2,w3,w4). The eta and xi values represent the
        ROMS y and x coordinate indices while the w values represent the
        corresponding coefficients for each of the four associated points
        from the original (irregular) grid. The w coefficients are calculated
        using an inverse-distance weight formula.

        These values are stored in the output index file for use when
        regridding actual model output to the regular grid. As long as the
        regular grid remains the same, these values will not change, thus the
        index file only needs to be generated once per model and kept on the
        data processing system in perpetuity.

        To determine which four irregular grid points (rho points), if any, a
        regular grid point falls betwen, the area of four triangles (created by
        connecting the regular grid point with each pair of adjacent irregular
        grid cells) is computed using determinants. If all four triangle areas
        are greater than 0, the regular grid point (x0,y0) is located inside
        the irregular grid cell enclosed by (eta1,xi1), (eta2,xi2), (eta3,xi3),
        (eta4,xi4).

        If a regular grid point does not fall within the area enclosed by any
        four valid irregular grid points, the location is specified as missing
        by assigning the preconfigured fill value (e.g. -9999.0) for each of
        the index and weight coefficient values. Otherwise, the proper eta, xi,
        and w coefficient values are determined and stored accordingly.

        Args:
            roms_file: `ROMSOutputFile` instance containing irregular grid
                structure to be used to compute index/coefficient values.
        """
        for y in range(self.dim_y.size):
            for x in range(self.dim_x.size):
                x0 = self.var_x[x]
                y0 = self.var_y[y]
                found_cell = False
                for xi1 in range(roms_file.num_xi-1):
                    if found_cell:
                        break
                    for eta1 in range(roms_file.num_eta-1):
                        xi2 = xi1 + 1
                        eta2 = eta1
                        xi3 = xi1 + 1
                        eta3 = eta1 + 1
                        xi4 = xi1
                        eta4 = eta1 + 1
                        if (roms_file.var_mask_rho[eta1,xi1] == 1
                                and roms_file.var_mask_rho[eta2,xi2] == 1
                                and roms_file.var_mask_rho[eta3,xi3] == 1
                                and roms_file.var_mask_rho[eta4,xi4] == 1):
                            x1 = roms_file.var_lon_rho[eta1,xi1]
                            y1 = roms_file.var_lat_rho[eta1,xi1]
                            x2 = roms_file.var_lon_rho[eta2,xi2]
                            y2 = roms_file.var_lat_rho[eta2,xi2]
                            x3 = roms_file.var_lon_rho[eta3,xi3]
                            y3 = roms_file.var_lat_rho[eta3,xi3]
                            x4 = roms_file.var_lon_rho[eta4,xi4]
                            y4 = roms_file.var_lat_rho[eta4,xi4]
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
                                self.var_xi1[y,x] = xi1
                                self.var_eta1[y,x] = eta1
                                self.var_w1[y,x] = w1
                                self.var_xi2[y,x] = xi2
                                self.var_eta2[y,x] = eta2
                                self.var_w2[y,x] = w2
                                self.var_xi3[y,x] = xi3
                                self.var_eta3[y,x] = eta3
                                self.var_w3[y,x] = w3
                                self.var_xi4[y,x] = xi4
                                self.var_eta4[y,x] = eta4
                                self.var_w4[y,x] = w4
                                self.var_wsum[y,x] = wsum
                                found_cell = True
                                break


class ROMSOutputFile:
    """Read/process data from a ROMS model output file.

    This class implements the context manager pattern and should thus be used
    similar to the following:

        with ROMSOutputFile("cbofs.nc") as f:
            ...

    Attributes:
        path: Path (relative or absolute) of the file.
        nc_file: Handle to `netCDF4.Dataset` instance for the opened NetCDF
            file.
    """
    def __init__(self, path):
        """Initialize ROMSOutputFile object and open file at specified path.

        Args:
            path: Path of target NetCDF file.
        
        Raises:
            Exception: Specified NetCDF file does not exist.
        """
        self.path = path
        if os.path.exists(self.path):
            self.nc_file = netCDF4.Dataset(self.path, "r", format="NETCDF3_Classic")
            self.init_handles()
        else:
            # File doesn't exist, raise error
            raise(Exception("NetCDF file does not exist: {}".format(self.path)))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self.nc_file.close()
    
    def init_handles(self):
        """Initialize handles to NetCDF variables.

        Because there is no U for the first column (xi=0) and no V for the
        first row (eta=0) in the ROMS staggered grid, we must skip the first
        row and column of every variable associated with rho points as well as
        the first row (eta=0) of U values, and the first column (xi=0) of V
        values. This ensures that u, v, and all other rho variables have the
        same dimensions in xi and eta, with u[eta,xi], v[eta,xi], and rho
        variables e.g. ang_rho[eta,xi] all correspond with the same grid cell
        for a given [eta,xi] coordinate.
        """
        self.var_ang_rho = self.nc_file.variables['angle'][1:,1:]
        self.var_lat_rho = self.nc_file.variables['lat_rho'][1:,1:]
        self.var_lon_rho = self.nc_file.variables['lon_rho'][1:,1:]
        self.var_u = self.nc_file.variables['u'][0,:,1:,:]
        self.var_v = self.nc_file.variables['v'][0,:,:,1:]
        self.var_mask_u = self.nc_file.variables['mask_u'][1:,:]
        self.var_mask_v = self.nc_file.variables['mask_v'][:,1:]
        self.var_mask_rho = self.nc_file.variables['mask_rho'][1:,1:]
        self.var_zeta = self.nc_file.variables['zeta'][0,1:,1:]
        self.var_h = self.nc_file.variables['h'][1:,1:]
        self.var_s_rho = self.nc_file.variables['s_rho'][:]
        self.var_hc = self.nc_file.variables['hc'][:]
        self.var_cs_r = self.nc_file.variables['Cs_r'][:]
        self.num_eta = self.var_h.shape[0]
        self.num_xi = self.var_h.shape[1]
        self.num_sigma = self.var_s_rho.shape[0]

    def uvToRegularGrid(self, model_index):
        """Interpolate u/v to regular grid"""
        # Extract the variables from the NetCDF

        # Call vertical function and return u and v at target depth
        u_depth, v_depth = vertInterp(self.var_u, self.var_v, self.var_s_rho, self.var_zeta, self.var_h, self.var_hc, self.var_cs_r, self.num_eta, self.num_xi, self.num_sigma)

        # Call masked arrays function and return masked arrays
        water_u,water_v,water_ang_rho,water_lat_rho,water_lon_rho = maskLand(u_depth, v_depth, self.var_ang_rho, self.var_lat_rho, self.var_lon_rho, self.var_mask_u, self.var_mask_v, self.var_mask_rho)

        # Call average to rho function u and v scalar values to rho
        u_rho, v_rho = averageUVToRho(water_u, water_v, self.num_eta, self.num_xi)
       
        # Call rotate function and return rotated u and v vectors
        rot_urho, rot_vrho = rotateUV2D(u_rho, v_rho, water_ang_rho)

        # Call create regular grid function 
        return interpolateUVToRegularGrid(water_lat_rho, water_lon_rho, rot_urho, rot_vrho, model_index)

def convertUVToSpeedDirection(reg_grid_u, reg_grid_v, model_index):
    """Convert u and v averaged/rotated vectors to speed/direction.

    Input u/v values are assumed to be in meters/sec. Output speed values will
    be converted to knots and direction in degrees from true north (0-360).

    Args:
        reg_grid_u: `numpy.ma.masked_array` containing u values interpolated to
            the regular grid.
        reg_grid_v: `numpy.ma.masked_array` containing v values interpolated to
            the regular grid.
        model_index: `ROMSIndexFile` instance representing index and
            coefficients file.
    """
    directions = numpy.empty((reg_grid_u.shape[0],reg_grid_u.shape[1]), dtype=numpy.float32)
    speeds = numpy.empty((reg_grid_u.shape[0],reg_grid_u.shape[1]), dtype=numpy.float32)
    for y in range(reg_grid_u.shape[0]):
        for x in range(reg_grid_u.shape[1]):
            u_ms = reg_grid_u[y,x]
            v_ms = reg_grid_v[y,x]

            #Convert from meters per second to knots
            u_knot = u_ms * MS2KNOTS
            v_knot = v_ms * MS2KNOTS

            currentSpeed = math.sqrt(math.pow(u_knot, 2) + math.pow(v_knot, 2))
            currentDirectionRadians = math.atan2(v_knot, u_knot)
            currentDirectionDegrees = math.degrees(currentDirectionRadians)
            currentDirectionNorth = 90.0 - currentDirectionDegrees

            #The direction must always be positive.
            if currentDirectionNorth < 0.0:
                currentDirectionNorth += 360.0

            directions[y,x] = currentDirectionNorth
            speeds[y,x] = currentSpeed

    return directions, speeds

def rotateUV2D(u_rho, v_rho, water_ang_rho):
    """Rotate vectors by geometric angle.

    Args:
        u_rho: `numpy.ma.masked_array` containing u values averaged to rho.
        v_rho: `numpy.ma.masked_array` containing v values averaged to rho.
        water_ang_rho: `numpy.ma.masked_array` containing angle-of-rotation
            values for rho points.
    """
    angsin = numpy.sin(water_ang_rho)
    angcos = numpy.cos(water_ang_rho)
    rot_urho = u_rho*angcos - v_rho*angsin
    rot_vrho = u_rho*angsin + v_rho*angcos

    return (rot_urho, rot_vrho)
        
def interpolateUVToRegularGrid(water_lat_rho, water_lon_rho, rot_urho, rot_vrho, model_index):
    """Create a regular grid using masked latitude and longitude variables.

    Interpolate averaged, rotated, u/v variables to the regular grid by using
    the NetCDF Indexes and Coeffiecients File to obtain weighted coefficents
    and eta/xi dimension indices to be used for Inverse Distance Weighting
    Interpolation.

    Args:
        water_lat_rho: `numpy.ma.masked_array` containing latitude of rho
            points with NoData/land values masked out.
        water_lon_rho: `numpy.ma.masked_array` containing longitude of rho
            points with NoData/land values masked out.
        rot_urho: `numpy.ma.masked_array` containing u values averaged to rho
            points and angle-of-rotation applied, with NoData/land values
            masked out.
        rot_urho: `numpy.ma.masked_array` containing v values averaged to rho
            points and angle-of-rotation applied, with NoData/land values
            masked out.
        model_index: `ROMSIndexFile` from which index values/coefficients will
            be extracted to perform interpolation.
    """
    # Create masked empty regular grid for variable u
    ugrid = numpy.ma.empty(shape=[model_index.dim_y.size,model_index.dim_x.size])
        
    # For each regular grid cell, read the corresponding xi1/eta1/w1,
    # xi2/eta2/w2, xi3/eta3/w3, and xi4/eta4/w4 values
    for y in range(model_index.dim_y.size):
        for x in range(model_index.dim_x.size):
            if not model_index.var_xi1.mask[y,x]:
                # Need to specify '.data[y,x]' or just '[y,x]' sufficient?
                xi1 = model_index.var_xi1[y,x]
                eta1 = model_index.var_eta1[y,x]
                xi2 = model_index.var_xi2[y,x]
                eta2 = model_index.var_eta2[y,x]
                xi3 = model_index.var_xi3[y,x]
                eta3 = model_index.var_eta3[y,x]
                xi4 = model_index.var_xi4[y,x]
                eta4 = model_index.var_eta4[y,x]
                u1 = rot_urho[eta1,xi1]
                u2 = rot_urho[eta2,xi2]
                u3 = rot_urho[eta3,xi3]
                u4 = rot_urho[eta4,xi4]
                # Use Inverse Distance Weigting algorithm to interpolate u to the regular grid
                ugrid[y,x] = ((((model_index.var_w1[y,x]) * u1) + ((model_index.var_w2[y,x]) * u2) + ((model_index.var_w3[y,x]) * u3) + ((model_index.var_w4[y,x])  * u4)) / model_index.var_wsum[y,x])
    
    # Create masked empty regular grid for variable v        
    vgrid = numpy.ma.empty(shape=[model_index.dim_y.size,model_index.dim_x.size])

    # For each regular grid cell, read the corresponding xi1/eta1/w1,
    # xi2/eta2/w2, xi3/eta3/w3, and xi4/eta4/w4 values
    for y in range(model_index.var_xi1.shape[0]):
        for x in range(model_index.var_xi1.shape[1]):
            if not model_index.var_xi1.mask[y,x]:
                xi1 = model_index.var_xi1[y,x]
                eta1 = model_index.var_eta1[y,x]
                xi2 = model_index.var_xi2[y,x]
                eta2 = model_index.var_eta2[y,x]
                xi3 = model_index.var_xi3[y,x]
                eta3 = model_index.var_eta3[y,x]
                xi4 = model_index.var_xi4[y,x]
                eta4 = model_index.var_eta4[y,x]
                v1 = rot_vrho[eta1,xi1]
                v2 = rot_vrho[eta2,xi2]
                v3 = rot_vrho[eta3,xi3]
                v4 = rot_vrho[eta4,xi4]
                # Use Inverse Distance Weigting algorithm to interpolate v to
                # the regular grid
                vgrid[y,x] = ((((model_index.var_w1[y,x]) * v1) + ((model_index.var_w2[y,x]) * v2) + ((model_index.var_w3[y,x]) * v3) + ((model_index.var_w4[y,x])  * v4)) / model_index.var_wsum[y,x])

    return (ugrid, vgrid)

def averageUVToRho(water_u, water_v, num_eta, num_xi):
    """Average u and v scalars to rho.

    Args:
        water_u: 2D numpy array containing U values to be averaged, with all
            NoData/Land-masked point values set to zero.
        water_v: 2D numpy array containing V values to be averaged, with all
            NoData/Land-masked point values set to zero.
    
    Returns:
        A 2-tuple of 2D `numpy.ndarray`s containing u and v values
        (respectively) averaged to rho points.
    """
    # Average u values
    u_rho = numpy.ndarray([num_eta, num_xi])
    
    for eta in range(num_eta):
        for xi in range(num_xi-1):
            u_rho[eta,xi] = (water_u[eta,xi] + water_u[eta,xi+1])/2
        u_rho[num_eta-1,xi] = water_u[num_eta-1,xi]
        
    # Average v values
    v_rho = numpy.ndarray([num_eta, num_xi])
    
    for xi in range(num_xi):
        for eta in range(num_eta-1):
            v_rho[eta,xi] = (water_v[eta,xi] + water_v[eta+1,xi])/2
        v_rho[eta,num_xi-1] = water_v[eta,num_xi-1]

    return u_rho, v_rho

def maskLand(u_depth, v_depth, ang_rho, lat_rho, lon_rho, mask_u, mask_v, mask_rho):
    """Create masked arrays for specified variables to mask land values.

    Args:
        u_depth: `numpy.ndarray` containing u values for entire grid at
            specified depth.
        v_depth: `numpy.ndarray` containing v values for entire grid at
            specified depth.
        ang_rho: `numpy.ndarray` containing angle-of-rotation values at
            rho points.
        lat_rho: `numpy.ndarray` containing latitude values of rho points.
        lon_rho: `numpy.ndarray` containing longitude values of rho points.
        mask_u: `numpy.ndarray` containing mask values for u points.
        mask_v: `numpy.ndarray` containing mask values for v points.
        mask_rho: `numpy.ndarray` containing mask values for rho points.
    """
    water_u = ma.masked_array(u_depth, numpy.logical_not(mask_u))
    water_v = ma.masked_array(v_depth, numpy.logical_not(mask_v))
    # u/v masked values need to be set to 0 for averaging 
    water_u = water_u.filled(0)
    water_v = water_v.filled(0)
    water_ang_rho = ma.masked_array(ang_rho, numpy.logical_not(mask_rho))
    water_lat_rho = ma.masked_array(lat_rho, numpy.logical_not(mask_rho))
    water_lon_rho = ma.masked_array(lon_rho, numpy.logical_not(mask_rho))  
           
    return water_u,water_v,water_ang_rho,water_lat_rho,water_lon_rho

def vertInterp(u, v, s_rho, zeta, h, hc, cs_r, num_eta, num_xi, num_sigma):
    """Vertically interpolate variables to target depth.

    Args:
        u: `numpy.ndarray` containing u values for entire grid.
        v: `numpy.ndarray` containing v values for entire grid.
        s_rho: `numpy.ndarray` s-coordinate at rho points, -1 min 0 max, positive up.
        zeta: `numpy.ndarray` containing MSL free surface at rho points in meters.
        h: `numpy.ndarray` containing bathymetry at rho points.
        hc: `numpy.ndarray` containing s-coordinate parameter, critical depth".
        cs_r: `numpy.ndarray` containing s-coordinate stretching curves at rho points.
    """

    z = numpy.ma.empty(shape=[num_eta,num_xi, 20])

    # Roms vertical transformation equation 1
    for k in range (num_sigma):
        S = ((hc * s_rho[k]) + (h - hc) * cs_r[k])
        z[:,:,k] = S + zeta * (1 + S/h)

    # For areas shallower than 9m the target depth is half the total depth
    # For areas deeper than 9m the target depth is 4.5m from zeta
    total_depth = h + zeta
    target_depth = zeta - numpy.minimum(9,total_depth)/2

    # For every rho point store z level index values and depth values, 
    # above and below the target depth in a masked array
    z_level = numpy.ma.empty(shape=[num_eta,num_xi, 4])

    for eta in range (num_eta):
        for xi in range(num_xi):
          if zeta.mask[eta,xi] != True:
              # Finds the closest values above and below the target depth
              depth1= (z[eta,xi,:])[(z[eta,xi,:]) >= (target_depth[eta,xi])].min()
              depth2 = (z[eta,xi,:])[(z[eta,xi,:]) <= (target_depth[eta,xi])].max()
            # Identifies the z levels and the depths above and below the target
              zmin = min(enumerate(z[eta,xi,:]), key=lambda x: abs(x[1]-depth1))
              zmax = min(enumerate(z[eta,xi,:]), key=lambda x: abs(x[1]-depth2))
              # Store each variable in the z_level masked array
              z_level[eta,xi,0]= zmin[0]# z level 1
              z_level[eta,xi,1]= zmax[0]# z level 2
              z_level[eta,xi,2]= zmin[1]# z level 1 depth value
              z_level[eta,xi,3]= zmax[1]# z level 2 depth value

    u_depth = numpy.ma.empty(shape=[num_eta,num_xi])

    for eta in range (num_eta):
        for xi in range(num_xi):
          if zeta.mask[eta,xi] != True:
              z1 = z_level[eta,xi,2] # z level 1 depth value
              z2 = z_level[eta,xi,3] # z level 2 depth value
              u_zmin = int(z_level[eta,xi,0])# u sigma level corresponding to z level 1
              u_zmax = int(z_level[eta,xi,1])# u sigma level corresponding to z level 2
              u1 = u[u_zmin,eta,xi]# u sigma level 1
              u2 = u[u_zmax,eta,xi]# u sigma level 2
              td1 = target_depth[eta,xi]
              # Using linear interpolation calculate u at target depth
              u_interp= u2 - ((u2-u1)*((z2 - td1)/(z2-z1)))
              # Store inerpolated u values in a masked array
              u_depth[eta,xi] = u_interp
              u_depth = numpy.nan_to_num(u_depth,-9999.0)

    v_depth = numpy.ma.empty(shape=[num_eta,num_xi])

    for eta in range (num_eta):
        for xi in range(num_xi):
          if zeta.mask[eta,xi] != True:
              z1 = z_level[eta,xi,2] # z level 1 value
              z2 = z_level[eta,xi,3] # z level 2 value
              v_zmin = int(z_level[eta,xi,0])# v sigma level corresponding to z level 1
              v_zmax = int(z_level[eta,xi,1])# v sigma level corresponding to z level 2
              v1 = v[v_zmin,eta,xi]# v sigma level 1
              v2 = v[v_zmax,eta,xi]# v sigma level 2
              d1 = target_depth[eta,xi]
              # Using linear interpolation calculate v at target depth
              v_interp= v2 - ((v2-v1)*((z2 - d1)/(z2-z1)))
              # Store inerpolated v values in a masked array
              v_depth[eta,xi] = v_interp
              v_depth = numpy.nan_to_num(v_depth,-9999.0)

    return u_depth, v_depth