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
import sys

import gdal
import json
import netCDF4
import numpy
import numpy.ma as ma
import osr
import ogr
from shapely.geometry import Polygon, Point, MultiPolygon, shape

_package = "s100ofs"
_dirPath = os.path.dirname(os.path.realpath(_package))

# Conversion factor for meters/sec to knots
MS2KNOTS = 1.943844

# Approximate radius of the Earth spheroid, in meters.
EARTH_RADIUS_METERS=6371000

# Default fill value for NetCDF variables
FILLVALUE=-9999.0

class RegularGrid:
    """Encapsulate information describing a regular lat-lon grid.

    Attributes:
        x_min: X-coordinate of the left edge of the bottom-left grid cell.
        y_min: Y-coordinate of the bottom edge of the bottom-left grid cell.
        x_max: X-coordinate of the right edge of the top-right grid cell.
        y_max: Y-coordinate of the top edge of the top-right grid cell.
        cellsize_x: Width of each regular grid cell in the x-direction.
        cellsize_y: Height of each regular grid cell in the y-direction.
        x_coords: List containing calculated x-coordinate of each grid column.
            The coordinate identifies the center location of each grid cell
            (for example, the x-coordinate of the first grid cell is equal to
            origin_x + cellsize_x/2.
        y_coords: List containing calculated y-coordinate of each grid row. The
            coordinate identifies the center location of each grid cell (for
            example, the y-coordinate of the first grid cell is equal to
            origin_y + cellsize_y/2.
    """
    def __init__(self, x_min, y_min, x_max, y_max, cellsize_x, cellsize_y):
        """Initialize RegularGrid object.

        Args:
            x_min: X-coordinate of the left edge of the bottom-left grid cell.
            y_min: Y-coordinate of the bottom edge of the bottom-left grid
                cell.
            x_max: X-coordinate of the right edge of the top-right grid cell.
            y_max: Y-coordinate of the top edge of the top-right grid cell.
            cellsize_x: Width of each regular grid cell in the x-direction.
            cellsize_y: Height of each regular grid cell in the y-direction.
        """
        self.x_min = x_min
        self.y_min = y_min
        self.x_max = x_max
        self.y_max = y_max
        self.cellsize_x = cellsize_x
        self.cellsize_y = cellsize_y
        self.calc_gridpoints()

    def calc_gridpoints(self):
        """Calculate and store x/y coordinates of grid points."""
        # Cell positions are calculated at center (midpoint) of cell
        x = self.x_min + self.cellsize_x/2
        y = self.y_min + self.cellsize_y/2
        
        self.x_coords = []
        while x <= self.x_max:
            self.x_coords.append(x)
            x += self.cellsize_x
        
        self.y_coords = []
        while y <= self.y_max:
            self.y_coords.append(y)
            y += self.cellsize_y

    @staticmethod
    def calc_cellsizes(lon_min, lat_min, lon_max, lat_max, target_cellsize_meters):
        """Calculate actual x/y cell sizes from an extent and target cell size.
        
        Given a lat/lon extent and target cell size in meters, calculate actual
        x and y cell sizes (in decimal degrees) that approximate the target
        cell size while ensuring the extent is divided into a whole number of
        cells in the x and y directions.

        Because degrees-per-meter varies by latitude, the midpoint of the given
        extent is used to calculate the conversion factor, as it should
        represent the average value for the whole grid.

        Args:
            lon_min: Minimum longitude of extent.
            lat_min: Minimum latitude of extent.
            lon_max: Maximum longitude of extent.
            lat_max: Maximum latitude of extent.
            target_cellsize_meters: Target cell size, in meters. Actual
                calculated cell sizes will be approximations of this.

        Returns:
            A 2-tuple in the form (cellsize_x, cellsize_y)
        """
        grid_width = abs(lon_max - lon_min)
        grid_height = abs(lat_max - lat_min)

        # Calculate number of meters per degree of longitude at the center of
        # the grid. This value can then be used to convert a cellsize in meters
        # to one in decimal degrees. The center of the grid is used as a good
        # approximation, since this value will vary by latitude.
        lat_mid = lat_min + grid_height/2
        meters_per_degree = ((EARTH_RADIUS_METERS * math.pi)/180)*math.cos(lat_mid)
        
        # Target cell size in decimal degrees
        target_cellsize_dd = target_cellsize_meters/meters_per_degree
        
        num_cells_x = int(round(grid_width/target_cellsize_dd))
        cellsize_x =  grid_width/num_cells_x

        num_cells_y = int(round(grid_height/target_cellsize_dd))
        cellsize_y = grid_height/num_cells_y
        
        return cellsize_x, cellsize_y
    
class ROMSIndexFile:
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

    DIMNAME_X = 'x'
    DIMNAME_Y = 'y'
    DIMNAME_SUBGRID = 'subgrid'

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
        self.dim_y = self.nc_file.dimensions[self.DIMNAME_Y]
        self.dim_x = self.nc_file.dimensions[self.DIMNAME_X]

        try:
            self.var_shoreline_mask = self.nc_file.variables['shoreline_mask'][:,:]
        except KeyError:
            pass

        try:
            self.dim_subgrid = self.nc_file.dimensions[self.DIMNAME_SUBGRID]
        except KeyError:
            pass

        try:
            self.var_subgrid_mask = self.nc_file.variables['subgrid_mask'][:,:,:]
        except KeyError:
            pass

        self.var_y = self.nc_file.variables[self.DIMNAME_Y][:]
        self.var_x = self.nc_file.variables[self.DIMNAME_X][:]

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
    
    def create_dims_coordvars(self, num_y, num_x):
        """Create NetCDF dimensions and coordinate variables.

        Args:
            num_y: Number of cells in x dimension.
            num_x: Number of cells in y dimension.
        """
        self.dim_y = self.nc_file.createDimension(self.DIMNAME_Y, num_y)
        self.dim_x = self.nc_file.createDimension(self.DIMNAME_X, num_x)
        
        # Create coordinate variables with same name as dimensions
        self.var_y = self.nc_file.createVariable(self.DIMNAME_Y, 'f4', (self.DIMNAME_Y,), fill_value=FILLVALUE)
        self.var_y.long_name = "latitude of regular grid point"
        self.var_y.units = "degree_north"
        self.var_y.standard_name = "latitude"

        self.var_x = self.nc_file.createVariable(self.DIMNAME_X, 'f4', (self.DIMNAME_X,), fill_value=FILLVALUE)
        self.var_x.long_name = "longitude of regular grid point"
        self.var_x.units = "degree_east"
        self.var_x.standard_name = "longitude"

    def create_index_coefficient_vars(self):
        """Create index/coefficient NetCDF variables."""
        self.var_xi1 = self.nc_file.createVariable('xi1', 'i4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_eta1 = self.nc_file.createVariable('eta1', 'i4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_w1 = self.nc_file.createVariable('w1', 'f4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_xi2 = self.nc_file.createVariable('xi2', 'i4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_eta2 = self.nc_file.createVariable('eta2', 'i4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_w2 = self.nc_file.createVariable('w2', 'f4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_xi3 = self.nc_file.createVariable('xi3', 'i4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_eta3 = self.nc_file.createVariable('eta3', 'i4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_w3 = self.nc_file.createVariable('w3', 'f4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_xi4 = self.nc_file.createVariable('xi4', 'i4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_eta4 = self.nc_file.createVariable('eta4', 'i4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_w4 = self.nc_file.createVariable('w4', 'f4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)
        self.var_wsum = self.nc_file.createVariable('wsum', 'f4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)

    def create_shoreline_mask_var(self):
        """Create shoreline mask NetCDF variable."""
        self.var_shoreline_mask = self.nc_file.createVariable('shoreline_mask', 'i4', (self.DIMNAME_Y,self.DIMNAME_X),fill_value=FILLVALUE)

    def create_subgrid_dims_vars(self, num_subgrids):
        """Create subgrid-related NetCDF dimensions/variables.

        Args:
            num_subgrids: Number of subgrids.
        """
        self.dim_subgrid = self.nc_file.createDimension(self.DIMNAME_SUBGRID, num_subgrids)
        self.var_subgrid_id = self.nc_file.createVariable('subgrid_id', 'i4', (self.DIMNAME_SUBGRID,),fill_value=FILLVALUE)
        self.var_subgrid_mask = self.nc_file.createVariable('subgrid_mask', 'i4', (self.DIMNAME_Y,self.DIMNAME_X,self.DIMNAME_SUBGRID),fill_value=FILLVALUE)
    
    def init_nc(self, roms_file, target_cellsize_meters, shoreline_shp=None, subset_grid_shp=None):
        """Initialize NetCDF dimensions/variables/attributes.

        Args:
            roms_file: `ROMSOutputFile` instance containing model output used
                to identify properties of original grid.
            target_cellsize_meters: Target cell size of grid cells, in meters.
                Actual calculated x/y grid cell sizes will vary slightly from
                this value, since the regular grid uses lat/lon coordinates
                (thus a cell's width/height in meters will vary by latitude),
                and since it will be adjusted in order to fit a whole number of
                grid cells in the x and y directions within the calculated grid
                extent.
            shoreline_shp: (Optional, default None) Path to a polygon shapefile
                containing features identifying land areas. If specified,
                a shoreline mask variable will be created/populated.
            subset_grid_shp: (Optional, default None) Path to a polygon
                shapefile containing orthogonal rectangles identifying areas
                to be used to subset (chop up) the full regular grid into
                tiles. Shapefile is assumed to be in the WGS84 projection. If
                None, the index file will be created assuming no subsets are
                desired and the extent of the model will be used instead.
        """
        # Calculate extent of valid (water) points
        water_lat_rho = ma.masked_array(roms_file.var_lat_rho, numpy.logical_not(roms_file.var_mask_rho))
        water_lon_rho = ma.masked_array(roms_file.var_lon_rho, numpy.logical_not(roms_file.var_mask_rho))
        lon_min = numpy.nanmin(water_lon_rho)
        lon_max = numpy.nanmax(water_lon_rho)
        lat_min = numpy.nanmin(water_lat_rho)
        lat_max = numpy.nanmax(water_lat_rho)
        
        # Populate grid x/y coordinate variables and subset-related variables
        # (if applicable)
        if subset_grid_shp is None:
            reg_grid = self.init_xy(lon_min, lat_min, lon_max, lat_max, target_cellsize_meters)
            self.gridOriginLongitude = reg_grid.x_min
            self.gridOriginLatitude = reg_grid.y_min
        else:
            reg_grid = self.init_xy_with_subsets(lon_min, lat_min, lon_max, lat_max, target_cellsize_meters, subset_grid_shp)
            self.gridOriginLongitude = reg_grid.x_min
            self.gridOriginLatitude = reg_grid.y_min

        # Create NetCDF variables
        self.create_index_coefficient_vars()

        if shoreline_shp is not None:
            land = self.init_shoreline_mask(reg_grid, shoreline_shp)

        self.nc_file.model = "CBOFS"
        self.nc_file.format = "netCDF-4"

        print (len(reg_grid.y_coords),len(reg_grid.x_coords))

        # Calculate the indexes/coefficients - can take many hours
        self.compute_indexes_coefficients(roms_file, land)

    def init_xy(self, lon_min, lat_min, lon_max, lat_max, target_cellsize_meters):
        """Create & initialize x/y dimensions/coordinate vars.
        
        Args:
            lon_min: Minimum longitude of domain.
            lat_min: Minimum latitude of domain.
            lon_max: Maximum longitude of domain.
            lat_max: Maximum latitude of domain.
            target_cellsize_meters: Target cell size, in meters. Actual
                calculated cell sizes will be approximations of this.
        """
        # Calculate actual x/y cell sizes
        cellsize_x, cellsize_y, num_cells_x, num_cells_y = RegularGrid.calc_cellsizes(lon_min, lat_min, lon_max, lat_max, target_cellsize_meters)
        
        # Build a regular grid using calculated cell sizes and given extent
        reg_grid = RegularGrid(lon_min, lat_min, lon_max, lat_max, cellsize_x, cellsize_y)
        
        # Create NetCDF dimensions & coordinate variables using dimension sizes
        # from regular grid
        self.create_dims_coordvars(len(reg_grid.y_coords), len(reg_grid.x_coords))
        
        # Populate NetCDF coordinate variables using regular grid coordinates
        self.var_x[:] = reg_grid.x_coords[:]
        self.var_y[:] = reg_grid.y_coords[:]
        self.gridSpacingLongitude = cellsize_x
        self.gridSpacingLatitude = cellsize_y

        return reg_grid

    def init_xy_with_subsets(self, lon_min, lat_min, lon_max, lat_max, target_cellsize_meters, subset_grid_shp):
        """Create & initialize x/y dimensions/coordinate vars and subset vars.

        Args:
            lon_min: Minimum longitude of domain.
            lat_min: Minimum latitude of domain.
            lon_max: Maximum longitude of domain.
            lat_max: Maximum latitude of domain.
            target_cellsize_meters: Target cell size, in meters. Actual
                calculated cell sizes will be approximations of this.
            subset_grid_shp: Path to subset grid polygon shapefile used to
                define subgrid domains.

        Raises: Exception when given subset grid shapefile does not exist or
            does not include any grid polygons intersecting with given extent.

        Returns: Instance of `RegularGrid` representing the extended generated
            grid whose extent matches the union of all intersecting subset grid
            polygons.
        """

        #shp = ogr.Open(shp_path)
        shp = ogr.Open(subset_grid_shp)
        layer = shp.GetLayer()

        # Create OGR Geometry from ocean model grid extent
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(lon_min, lat_max)
        ring.AddPoint(lon_min, lat_min)
        ring.AddPoint(lon_max, lat_min)
        ring.AddPoint(lon_max, lat_max)
        ring.AddPoint(lon_min, lat_max)
        # Create polygon
        ofs_poly = ogr.Geometry(ogr.wkbPolygon)
        ofs_poly.AddGeometry(ring)

        # Find the intersection between 160k grid and ocean model grid extent
        subset_polys = {}
        fids = []
        fid = 0
        for feature in layer:
            geom = feature.GetGeometryRef()
            intersection = ofs_poly.Intersection(geom)
            if intersection.ExportToWkt() != "GEOMETRYCOLLECTION EMPTY":
                subset_polys[fid] = geom.ExportToJson()
                fids.append(fid)
            fid += 1

        if len(fids) == 0:
            raise Exception("Given subset grid shapefile contains no polygons that intersect with model domain; cannot proceed.")

        # Use a single subset polygon to calculate x/y cell sizes. This ensures
        # that cells do not fall on the border between two grid polygons.
        singlepolygon = ogr.Geometry(ogr.wkbMultiPolygon)
        singlepolygon.AddGeometry(ogr.CreateGeometryFromJson(subset_polys[fids[0]]))
        sp_x_min, sp_x_max, sp_y_min, sp_y_max = singlepolygon.GetEnvelope()

        cellsize_x, cellsize_y = RegularGrid.calc_cellsizes(sp_x_min, sp_y_min, sp_x_max, sp_y_max, target_cellsize_meters)
        # Combine identified subset grid polygons into single multipolygon to
        # calculate full extent of all combined subset grids
        multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
        for fid in fids:
            multipolygon.AddGeometry(ogr.CreateGeometryFromJson(subset_polys[fid]))

        (x_min, x_max, y_min, y_max) = multipolygon.GetEnvelope()
        full_reg_grid = RegularGrid(x_min, y_min, x_max, y_max, cellsize_x, cellsize_y)

        # Create NetCDF dimensions & coordinate variables using dimension sizes
        # from regular grid
        self.create_dims_coordvars(len(full_reg_grid.y_coords), len(full_reg_grid.x_coords))
        # Populate NetCDF coordinate variables using regular grid coordinates
        self.var_x[:] = full_reg_grid.x_coords[:]
        self.var_y[:] = full_reg_grid.y_coords[:]
        self.gridSpacingLongitude = full_reg_grid.cellsize_x
        self.gridSpacingLatitude = full_reg_grid.cellsize_y

        # Create subgrid dimension/variables
        self.create_subgrid_dims_vars(len(subset_polys))

        # Populate subgrid mask variable
        for subgrid_index, fid in enumerate(fids):
            # Convert OGR geometry to shapely geometry
            subset_poly_shape = shape(json.loads(subset_polys[fid]))
            for eta in range(len(self.var_y)):
                for xi in range(len(self.var_x)):
                    point = Point(self.var_x[xi], self.var_y[eta])
                    if point.within(subset_poly_shape):
                        self.var_subgrid_mask[eta,xi,subgrid_index] = 1

        return full_reg_grid

    def init_shoreline_mask(self, reg_grid, shoreline_shp):
        """Use shoreline shapefile to mask out land at highest resolution.
        
        Args:
            reg_grid: `RegularGrid` instance describing the regular grid for
                which the shoreline mask will be created.
            shoreline_shp: Path to a polygon shapefile containing features
                identifying land areas.

        Returns:
            2D numpy array, matching the dimensions of the given RegularGrid,
            containing a value of 1 for water areas and a value of 255 for land
            areas.
        """

        #shp = ogr.Open(shr_path)
        shp = ogr.Open(shoreline_shp)
        layer = shp.GetLayer()

        # Create OGR Geometry from ocean model grid extent
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(reg_grid.x_min, reg_grid.y_max)
        ring.AddPoint(reg_grid.x_min, reg_grid.y_min)
        ring.AddPoint(reg_grid.x_max, reg_grid.y_min)
        ring.AddPoint(reg_grid.x_max, reg_grid.y_max)
        ring.AddPoint(reg_grid.x_min, reg_grid.y_max)
        # Create polygon
        ofs_poly = ogr.Geometry(ogr.wkbPolygon)
        ofs_poly.AddGeometry(ring)

        # Find the intersection between shoreline shapefile and ocean model grid extent
        for feature in layer:
            geom = feature.GetGeometryRef()
            intersection = ofs_poly.Intersection(geom)

        # Create a memory layer to rasterize from.
        driver = ogr.GetDriverByName('Memory')
        memds=driver.CreateDataSource('tmpmemds')
        lyr=memds.CreateLayer('land', geom_type=ogr.wkbPolygon)
        feat = ogr.Feature(lyr.GetLayerDefn())
        feat.SetGeometry(intersection)
        lyr.CreateFeature(feat)

        # Create raster
        pixelWidth = reg_grid.cellsize_x
        pixelHeight = reg_grid.cellsize_y
        cols = len(reg_grid.y_coords)
        rows = len(reg_grid.x_coords)
        target_ds = gdal.GetDriverByName('GTiff').Create("land.tif", rows, cols, 1, gdal.GDT_Float32)
        target_ds.SetGeoTransform((reg_grid.x_min, pixelWidth, 0, reg_grid.y_min, 0,  pixelHeight))
        band = target_ds.GetRasterBand(1)
        NoData_value = 1
        band.SetNoDataValue(1)
        band.FlushCache()

        gdal.RasterizeLayer(target_ds, [1], lyr, None, None)
        target_dsSRS = osr.SpatialReference()
        target_dsSRS.ImportFromEPSG(4326)
        target_ds.SetProjection(target_dsSRS.ExportToWkt())
        target_ds = None

        # Store as numpy array, land = 255, water = 1
        land = gdal.Open("land.tif").ReadAsArray()

        return land

    def compute_indexes_coefficients(self, roms_file, land):
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
            land: 2D numpy array, matching the dimensions of this index file,
                containing a value of 1 for water areas and a value of 255 for
                land areas.
        """
        for y in range(self.dim_y.size):
            for x in range(self.dim_x.size):
                x0 = self.var_x[x]
                y0 = self.var_y[y]
                found_cell = False
                if land[y,x] != 1:
                    continue
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
              u_depth = numpy.nan_to_num(u_depth,FILLVALUE)

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
              v_depth = numpy.nan_to_num(v_depth,FILLVALUE)

    return u_depth, v_depth
