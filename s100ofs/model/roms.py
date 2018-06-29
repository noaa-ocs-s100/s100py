"""
Utility classes and methods for working with ROMS ocean model output.

The Regional Ocean Modeling System (ROMS) is a 3D hydrodynamic modeling
framework which uses an irregular, curvilinear horizontal grid and a sigma
(bathymetry-following) vertical coordinate system. This module provides
functionality allowing ROMS output to be interpolated to a regular, orthogonal
lat/lon horizontal grid at a given depth-below-surface.
"""

import math
import os

import gdal
import json
import netCDF4
import numpy
import numpy.ma as ma
import osr
import ogr
from shapely.geometry import Polygon, Point, MultiPolygon, shape
from scipy import interpolate
import datetime

# Conversion factor for meters/sec to knots
MS2KNOTS = 1.943844

# Approximate radius of the Earth spheroid, in meters.
EARTH_RADIUS_METERS = 6371000

# Default fill value for NetCDF variables
FILLVALUE = -9999.0


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
        lat_mid_radians = math.radians(lat_mid)
        meters_per_degree = ((EARTH_RADIUS_METERS * math.pi)/180)*math.cos(lat_mid_radians)
        
        # Target cell size in decimal degrees
        target_cellsize_dd = target_cellsize_meters/meters_per_degree
        
        num_cells_x = int(round(grid_width/target_cellsize_dd))
        cellsize_x = grid_width/num_cells_x

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
        var_mask: Handle to master_mask variable.
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
            self.dim_subgrid = self.nc_file.dimensions[self.DIMNAME_SUBGRID]
        except KeyError:
            self.dim_subgrid = None

        try:
            self.var_subgrid_id = self.nc_file.variables['subgrid_id'][:]
        except KeyError:
            self.var_subgrid_id = None

        try:
            self.var_subgrid_x_min = self.nc_file.variables['subgrid_x_min'][:]
        except KeyError:
            self.var_subgrid_x_min = None
        try:
            self.var_subgrid_x_max = self.nc_file.variables['subgrid_x_max'][:]
        except KeyError:
            self.var_subgrid_x_max = None
        try:
            self.var_subgrid_y_min = self.nc_file.variables['subgrid_y_min'][:]
        except KeyError:
            self.var_subgrid_y_min = None
        try:
            self.var_subgrid_y_max = self.nc_file.variables['subgrid_y_max'][:]
        except KeyError:
            self.var_subgrid_y_max = None

        self.var_y = self.nc_file.variables[self.DIMNAME_Y][:]
        self.var_x = self.nc_file.variables[self.DIMNAME_X][:]
        self.var_mask = self.nc_file.variables['mask'][:, :]

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
        self.var_y.long_name = "latitude of regular grid point at cell center"
        self.var_y.units = "degree_north"
        self.var_y.standard_name = "latitude"

        self.var_x = self.nc_file.createVariable(self.DIMNAME_X, 'f4', (self.DIMNAME_X,), fill_value=FILLVALUE)
        self.var_x.long_name = "longitude of regular grid point at cell center"
        self.var_x.units = "degree_east"
        self.var_x.standard_name = "longitude"

        self.var_mask = self.nc_file.createVariable('mask', 'i4', (self.DIMNAME_Y, self.DIMNAME_X), fill_value=FILLVALUE)
        self.var_mask.long_name = "regular grid point mask"
        self.var_mask.flag_values = 1, 255;
        self.var_mask.flag_meanings = "land, water";


    def create_subgrid_dims_vars(self, num_subgrids):
        """Create subgrid-related NetCDF dimensions/variables.

        Args:
            num_subgrids: Number of subgrids.
        """
        self.dim_subgrid = self.nc_file.createDimension(self.DIMNAME_SUBGRID, num_subgrids)
        self.var_subgrid_id = self.nc_file.createVariable('subgrid_id', 'i4', (self.DIMNAME_SUBGRID,), fill_value=FILLVALUE)
        self.var_subgrid_x_min = self.nc_file.createVariable('subgrid_x_min', 'i4', (self.DIMNAME_SUBGRID,))
        self.var_subgrid_x_max = self.nc_file.createVariable('subgrid_x_max', 'i4', (self.DIMNAME_SUBGRID,))
        self.var_subgrid_y_min = self.nc_file.createVariable('subgrid_y_min', 'i4', (self.DIMNAME_SUBGRID,))
        self.var_subgrid_y_max = self.nc_file.createVariable('subgrid_y_max', 'i4', (self.DIMNAME_SUBGRID,))
    
    def init_nc(self, roms_file, target_cellsize_meters, ofs_model, shoreline_shp=None, subset_grid_shp=None):
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
            ofs_model: The target model identifier.
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
        else:
            reg_grid = self.init_xy_with_subsets(lon_min, lat_min, lon_max, lat_max, target_cellsize_meters, subset_grid_shp)

        self.nc_file.gridOriginLongitude = reg_grid.x_min
        self.nc_file.gridOriginLatitude = reg_grid.y_min

        land_mask = None
        if shoreline_shp is not None:
            land_mask = self.init_shoreline_mask(reg_grid, shoreline_shp)

        self.nc_file.model = str.upper(ofs_model)
        self.nc_file.format = "netCDF-4"

        print ("Full grid dimensions (y,x): ({},{})".format(len(reg_grid.y_coords), len(reg_grid.x_coords)))

        # Calculate the mask
        self.compute_mask(roms_file, reg_grid, land_mask)

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
        cellsize_x, cellsize_y = RegularGrid.calc_cellsizes(lon_min, lat_min, lon_max, lat_max, target_cellsize_meters)
        
        # Build a regular grid using calculated cell sizes and given extent
        reg_grid = RegularGrid(lon_min, lat_min, lon_max, lat_max, cellsize_x, cellsize_y)
        
        # Create NetCDF dimensions & coordinate variables using dimension sizes
        # from regular grid
        self.create_dims_coordvars(len(reg_grid.y_coords), len(reg_grid.x_coords))
        
        # Populate NetCDF coordinate variables using regular grid coordinates
        self.var_x[:] = reg_grid.x_coords[:]
        self.var_y[:] = reg_grid.y_coords[:]
        self.nc_file.gridSpacingLongitude = cellsize_x
        self.nc_file.gridSpacingLatitude = cellsize_y

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

        # Get the EPSG value from the import shapefile and transform to WGS84
        spatialRef = layer.GetSpatialRef()
        shpSRS = spatialRef.GetAttrValue("AUTHORITY", 1)
        source = osr.SpatialReference()
        source.ImportFromEPSG(int(shpSRS))
        target = osr.SpatialReference()
        target.ImportFromEPSG(4326)
        transform = osr.CoordinateTransformation(source, target)
        ofs_poly.Transform(transform)

        # Find the intersection between grid polygon and ocean model grid extent
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
        self.nc_file.gridSpacingLongitude = full_reg_grid.cellsize_x
        self.nc_file.gridSpacingLatitude = full_reg_grid.cellsize_y

        # Create subgrid dimension/variables
        self.create_subgrid_dims_vars(len(subset_polys))

        print ("start", datetime.datetime.now().time())
        # Calculate subgrid mask ranges, populate subgrid ID
        for subgrid_index, fid in enumerate(fids):
            self.var_subgrid_id[subgrid_index] = fid
            # Start with extreme values for min/max, then narrow down
            subgrid_x_min = len(full_reg_grid.x_coords)
            subgrid_x_max = 0
            subgrid_y_min = len(full_reg_grid.y_coords)
            subgrid_y_max = 0
            # Convert OGR geometry to shapely geometry
            subset_poly_shape = shape(json.loads(subset_polys[fid]))
            for eta in range(len(self.var_y)):
                for xi in range(len(self.var_x)):
                    point = Point(self.var_x[xi], self.var_y[eta])
                    if point.within(subset_poly_shape):
                        if eta < subgrid_y_min:
                            subgrid_y_min = eta
                        if eta > subgrid_y_max:
                            subgrid_y_max = eta
                        if xi < subgrid_x_min:
                            subgrid_x_min = xi
                        if xi > subgrid_x_max:
                            subgrid_x_max = xi
            if subgrid_x_min >= subgrid_x_max or subgrid_y_min >= subgrid_y_max:
                raise Exception("Error calculating subgrid index ranges for subgrid [{} - fid {}]:x_min: {}, x_max: {}, y_min: {}, y_max: {}".format(subgrid_index, fid, subgrid_x_min, subgrid_x_max, subgrid_y_min, subgrid_y_max))
            self.var_subgrid_x_min[subgrid_index] = subgrid_x_min
            self.var_subgrid_x_max[subgrid_index] = subgrid_x_max
            self.var_subgrid_y_min[subgrid_index] = subgrid_y_min
            self.var_subgrid_y_max[subgrid_index] = subgrid_y_max
        print ("end", datetime.datetime.now().time())

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

        # Get the EPSG value from the import shapefile and transform to WGS84
        spatialRef = layer.GetSpatialRef()
        shpSRS = spatialRef.GetAttrValue("AUTHORITY", 1)
        source = osr.SpatialReference()
        source.ImportFromEPSG(int(shpSRS))
        target = osr.SpatialReference()
        target.ImportFromEPSG(4326)
        transform = osr.CoordinateTransformation(source, target)
        ofs_poly.Transform(transform)

        # Find the intersection between shoreline shapefile and ocean model grid extent
        for feature in layer:
            geom = feature.GetGeometryRef()
            intersection = ofs_poly.Intersection(geom)

        # Create a memory layer to rasterize from.
        driver = ogr.GetDriverByName('Memory')
        memds = driver.CreateDataSource('tmpmemds')
        lyr = memds.CreateLayer('land_mask', geom_type=ogr.wkbPolygon)
        feat = ogr.Feature(lyr.GetLayerDefn())
        feat.SetGeometry(intersection)
        lyr.CreateFeature(feat)

        # Create raster
        pixel_width = reg_grid.cellsize_x
        pixel_height = reg_grid.cellsize_y
        cols = len(reg_grid.y_coords)
        rows = len(reg_grid.x_coords)
        target_ds = gdal.GetDriverByName('GTiff').Create("land_mask.tif", rows, cols, 1, gdal.GDT_Float32)
        target_ds.SetGeoTransform((reg_grid.x_min, pixel_width, 0, reg_grid.y_min, 0,  pixel_height))
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(1)
        band.FlushCache()

        gdal.RasterizeLayer(target_ds, [1], lyr, None, None)
        target_dsSRS = osr.SpatialReference()
        target_dsSRS.ImportFromEPSG(4326)
        target_ds.SetProjection(target_dsSRS.ExportToWkt())
        target_ds = None

        # Store as numpy array, land = 255, water = 1
        land_mask = gdal.Open("land_mask.tif").ReadAsArray()

        return land_mask

    def compute_mask(self, roms_file, reg_grid, land_mask):
        """Create model domain mask.

        For every irregular grid point create a polygon from four valid
        grid points, searching counter clockwise (eta1,xi1), (eta2,xi2), (eta3,xi3),
        (eta4,xi4). Rasterize the polygon to create grid domain mask.


        Args:
            roms_file: `ROMSOutputFile` instance containing irregular grid
                structure and variables.
            reg_grid: `RegularGrid` instance describing the regular grid for
                which the mask will be created.
            land_mask: 2D numpy array, matching the dimensions of this index file,
                containing a value of 1 for water areas and a value of 255 for
                land_mask areas.
        """
        # Create shapefile with OGR
        driver = ogr.GetDriverByName('Esri Shapefile')
        ds = driver.CreateDataSource('grid_cell_mask.shp')
        layer = ds.CreateLayer('', None, ogr.wkbMultiPolygon)
        layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))

        # Add spatial reference to polygon
        spatialRef = osr.SpatialReference()
        spatialRef.ImportFromEPSG(4326)
        spatialRef.MorphToESRI()
        file = open('grid_cell_mask.prj', 'w')
        file.write(spatialRef.ExportToWkt())
        file.close()

        # Create shapefile containing polygons for each irregular grid cell
        # using four valid irregular grid points, searching counter
        # clockwise(eta1, xi), (eta2, xi2), (eta3, xi3),(eta4, xi4)
        print ("start", datetime.datetime.now().time())
        for xi1 in range(roms_file.num_xi - 1):
            for eta1 in range(roms_file.num_eta - 1):
                xi2 = xi1 + 1
                eta2 = eta1
                xi3 = xi1 + 1
                eta3 = eta1 + 1
                xi4 = xi1
                eta4 = eta1 + 1

                # Search valid points
                valid_points = []
                for (eta, xi) in ((eta1, xi1), (eta2, xi2), (eta3, xi3), (eta4, xi4)):
                    if roms_file.var_mask_rho[eta, xi] == 1:
                        valid_points.append((eta, xi))
                if len(valid_points) < 3:
                    continue
                ring = ogr.Geometry(ogr.wkbLinearRing)
                for (eta, xi) in valid_points:
                    ring.AddPoint(roms_file.var_lon_rho[eta, xi], roms_file.var_lat_rho[eta, xi])
                (eta, xi) = valid_points[0]
                ring.AddPoint(roms_file.var_lon_rho[eta, xi], roms_file.var_lat_rho[eta, xi])

                # Create polygon
                geom = ogr.Geometry(ogr.wkbPolygon)
                geom.AddGeometry(ring)
                feat = ogr.Feature(layer.GetLayerDefn())
                feat.SetField('id', xi1)
                feat.SetGeometry(geom)
                layer.CreateFeature(feat)
        print ("end", datetime.datetime.now().time())

        # Rasterize grid cell polygons
        pixel_width = reg_grid.cellsize_x
        pixel_height = reg_grid.cellsize_y
        cols = len(reg_grid.y_coords)
        rows = len(reg_grid.x_coords)
        target_ds = gdal.GetDriverByName('GTiff').Create("grid_cell_mask.tif", rows, cols, 1, gdal.GDT_Byte)
        target_ds.SetGeoTransform((reg_grid.x_min, pixel_width, 0, reg_grid.y_min, 0, pixel_height))
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(1)
        band.FlushCache()

        gdal.RasterizeLayer(target_ds, [1], layer, None, None)
        target_dsSRS = osr.SpatialReference()
        target_dsSRS.ImportFromEPSG(4326)
        target_ds.SetProjection(target_dsSRS.ExportToWkt())

        # Save and close everything
        target_ds = None
        feat = geom = None
        ds = layer = feat = geom = None

        # Store as numpy array, # value of 1.0 for invalid areas and
        # a value of 255 for valid areas.
        grid_cell_mask = gdal.Open("grid_cell_mask.tif").ReadAsArray()

        # Use land mask and grid cell mask to create master mask and
        # write to index file
        print ("start_xy", datetime.datetime.now().time())
        for y in range(self.dim_y.size):
            for x in range(self.dim_x.size):
                if land_mask[y, x] != 1:
                    continue
                if grid_cell_mask[y, x] != 1:
                    self.var_mask[y, x] = 1
                else:
                    self.var_mask[y, x] = FILLVALUE
        print ("end_xy", datetime.datetime.now().time())


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
        self.var_ang_rho = self.nc_file.variables['angle'][1:, 1:]
        self.var_lat_rho = self.nc_file.variables['lat_rho'][1:, 1:]
        self.var_lon_rho = self.nc_file.variables['lon_rho'][1:, 1:]
        self.var_u = self.nc_file.variables['u'][0, :, 1:, :]
        self.var_v = self.nc_file.variables['v'][0, :, :, 1:]
        self.var_mask_u = self.nc_file.variables['mask_u'][1:, :]
        self.var_mask_v = self.nc_file.variables['mask_v'][:, 1:]
        self.var_mask_rho = self.nc_file.variables['mask_rho'][1:, 1:]
        self.var_zeta = self.nc_file.variables['zeta'][0, 1:, 1:]
        self.var_h = self.nc_file.variables['h'][1:, 1:]
        self.var_s_rho = self.nc_file.variables['s_rho'][:]
        self.var_hc = self.nc_file.variables['hc'][:]
        self.var_cs_r = self.nc_file.variables['Cs_r'][:]
        self.var_vtransform = self.nc_file.variables['Vtransform'][:]
        self.num_eta = self.var_h.shape[0]
        self.num_xi = self.var_h.shape[1]
        self.num_sigma = self.var_s_rho.shape[0]

    def uv_to_regular_grid(self, model_index, target_depth):
        """Interpolate u/v to regular grid"""
        # Extract the variables from the NetCDF

        # Call vertical function and return u and v at target depth
        u_target_depth, v_target_depth = vertical_interpolation(self.var_u, self.var_v, self.var_s_rho,
                                                                self.var_mask_rho, self.var_mask_u, self.var_mask_v,
                                                                self.var_zeta, self.var_h, self.var_hc, self.var_cs_r,
                                                                self.var_vtransform, self.num_eta, self.num_xi,
                                                                self.num_sigma, target_depth)
        # Call masked arrays function and return masked arrays
        water_u, water_v, water_ang_rho, water_lat_rho, water_lon_rho = mask_land(u_target_depth, v_target_depth,
                                                                                  self.var_ang_rho, self.var_lat_rho,
                                                                                  self.var_lon_rho, self.var_mask_u,
                                                                                  self.var_mask_v, self.var_mask_rho)

        # Call average to rho function u and v scalar values to rho
        u_rho, v_rho = average_uv2rho(water_u, water_v)
        
        # Call rotate function and return rotated u and v vectors
        rot_u_rho, rot_v_rho = rotate_uv2d(u_rho, v_rho, water_ang_rho)
        
        # Call create regular grid function 
        return interpolate_uv_to_regular_grid(rot_u_rho, rot_v_rho, water_lat_rho, water_lon_rho, model_index)


def uv_to_speed_direction(reg_grid_u, reg_grid_v):
    """Convert u and v averaged/rotated vectors to speed/direction.

    Input u/v values are assumed to be in meters/sec. Output speed values will
    be converted to knots and direction in degrees from true north (0-360).

    Args:
        reg_grid_u: `numpy.ma.masked_array` containing u values interpolated to
            the regular grid.
        reg_grid_v: `numpy.ma.masked_array` containing v values interpolated to
            the regular grid.
    """
    direction = numpy.empty((reg_grid_u.shape[0], reg_grid_u.shape[1]), dtype=numpy.float32)
    speed = numpy.empty((reg_grid_u.shape[0], reg_grid_u.shape[1]), dtype=numpy.float32)
    for y in range(reg_grid_u.shape[0]):
        for x in range(reg_grid_u.shape[1]):
            try:
                u_ms = reg_grid_u[y, x]
                v_ms = reg_grid_v[y, x]

                # Convert from meters per second to knots
                u_knot = u_ms * MS2KNOTS
                v_knot = v_ms * MS2KNOTS

                current_speed = math.sqrt(math.pow(u_knot, 2) + math.pow(v_knot, 2))
                current_direction_radians = math.atan2(v_knot, u_knot)
                current_direction_degrees = math.degrees(current_direction_radians)
                current_direction_north = 90.0 - current_direction_degrees
            except OverflowError as e:
                print("OverflowError covering uv to speed/dir at y,x: {},{}".format(y, x))
                print("reg_grid_u.mask[y,x]: {}".format(reg_grid_u.mask[y, x]))
                print("reg_grid_v.mask[y,x]: {}".format(reg_grid_v.mask[y, x]))
                print("u_ms: {}, v_ms: {}".format(u_ms, v_ms))
                print("u_knot: {}, v_knot: {}".format(u_knot, v_knot))
                raise e

            # The direction must always be positive.
            if current_direction_north < 0.0:
                current_direction_north += 360.0

            direction[y, x] = current_direction_north
            speed[y, x] = current_speed

    return direction, speed


def rotate_uv2d(u_rho, v_rho, water_ang_rho):
    """Rotate vectors by geometric angle.

    Args:
        u_rho: `numpy.ma.masked_array` containing u values averaged to rho.
        v_rho: `numpy.ma.masked_array` containing v values averaged to rho.
        water_ang_rho: `numpy.ma.masked_array` containing angle-of-rotation
            values for rho points.
    """
    ang_sin = numpy.sin(water_ang_rho)
    ang_cos = numpy.cos(water_ang_rho)
    rot_u_rho = u_rho*ang_cos - v_rho*ang_sin
    rot_v_rho = u_rho*ang_sin + v_rho*ang_cos

    return rot_u_rho, rot_v_rho


def interpolate_uv_to_regular_grid(rot_u_rho, rot_v_rho, water_lat_rho, water_lon_rho, model_index):
    """Create a regular grid using masked latitude and longitude variables.

    Interpolate averaged, rotated, u/v variables to the regular grid using
    Inverse Distance Weighted Interpolation.

    Args:
        rot_u_rho: `numpy.ma.masked_array` containing u values averaged to rho
            points and angle-of-rotation applied, with NoData/land values
            masked out.
        rot_v_rho: `numpy.ma.masked_array` containing v values averaged to rho
            points and angle-of-rotation applied, with NoData/land values
            masked out.
        water_lat_rho: `numpy.ndarray` containing masked latitude values of rho
            points.
        water_lon_rho: `numpy.ndarray` containing masked longitude values of rho
            points.
        model_index: `ROMSIndexFile` from which index values/coefficients will
            be extracted to perform interpolation.
    """
    # Flatten and compress variables
    rot_u_rho = ma.compressed(rot_u_rho)
    rot_v_rho = ma.compressed(rot_v_rho)
    water_lat_rho = ma.compressed(water_lat_rho)
    water_lon_rho = ma.compressed(water_lon_rho)

    # Create an ogr object containing irregular points eta,xi,u,v and write to memory
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS("WGS84")
    ds = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
    layer = ds.CreateLayer("irregular_points", srs=srs, geom_type=ogr.wkbPoint)
    layer.CreateField(ogr.FieldDefn("u", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("v", ogr.OFTReal))

    for i in range(len(rot_v_rho)):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(water_lon_rho[i], water_lat_rho[i])
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(point)
        feature.SetField("u", rot_u_rho[i])
        feature.SetField("v", rot_v_rho[i])
        layer.CreateFeature(feature)

    # Using ogr object run gdal grid to interpolate irregular grid values to regular grid values
    dst_u = gdal.Grid('u.tif', ds, format='MEM', width=model_index.dim_x.size, height=model_index.dim_y.size,
                      algorithm="invdist:power=2.0:smoothing=0.0:radius1=0.04:radius2=0.04:angle=0.0:max_points=0:min_points=2:nodata=0.0",
                      zfield="u")
    dst_v = gdal.Grid('v.tif', ds, format='MEM', width=model_index.dim_x.size, height=model_index.dim_y.size,
                      algorithm="invdist:power=2.0:smoothing=0.0:radius1=0.04:radius2=0.04:angle=0.0:max_points=0:min_points=2:nodata=0.0",
                      zfield="v")

    reg_grid_u = dst_u.ReadAsArray()
    reg_grid_v = dst_v.ReadAsArray()

    return reg_grid_u, reg_grid_v


def average_uv2rho(water_u, water_v):
    """Average u and v scalars to rho.

    U values at rho [eta,xi] are calculated by averaging u[eta,xi] with
    u[eta,xi+1], except for the final xi, where u at rho [eta,xi_max] is
    simply set to u[eta,xi_max-1].
    
    V values at rho [eta,xi] are calculated by averaging v[eta,xi] with
    u[eta+1,xi], except for the final eta, where v at rho [eta_max,xi] is
    simply set to v[eta_max-1,xi].
    
    Args:
        water_u: 2D numpy array containing U values to be averaged, with all
            NoData/Land-masked point values set to zero.
        water_v: 2D numpy array containing V values to be averaged, with all
            NoData/Land-masked point values set to zero.
    
    Returns:
        A 2-tuple of 2D `numpy.ndarray`s containing u and v values
        (respectively) averaged to rho points.
    """
    num_eta = water_u.shape[0]
    num_xi = water_u.shape[1]

    # Average u values
    u_rho = numpy.ndarray([num_eta, num_xi])
    
    for eta in range(num_eta):
        for xi in range(num_xi-1):
            u_rho[eta, xi] = (water_u[eta, xi] + water_u[eta, xi+1])/2
        u_rho[eta, num_xi-1] = water_u[eta, num_xi-1]

    # Average v values
    v_rho = numpy.ndarray([num_eta, num_xi])
    
    for xi in range(num_xi):
        for eta in range(num_eta-1):
            v_rho[eta, xi] = (water_v[eta, xi] + water_v[eta+1, xi])/2
        v_rho[num_eta-1, xi] = water_v[num_eta-1, xi]

    return u_rho, v_rho


def mask_land(u_target_depth, v_target_depth, ang_rho, lat_rho, lon_rho, mask_u, mask_v, mask_rho):
    """Create masked arrays for specified variables to mask land values.

    Args:
        u_target_depth: `numpy.ndarray` containing u values for entire grid at
            specified depth.
        v_target_depth: `numpy.ndarray` containing v values for entire grid at
            specified depth.
        ang_rho: `numpy.ndarray` containing angle-of-rotation values at
            rho points.
        lat_rho: `numpy.ndarray` containing latitude values of rho points.
        lon_rho: `numpy.ndarray` containing longitude values of rho points.
        mask_u: `numpy.ndarray` containing mask values for u points.
        mask_v: `numpy.ndarray` containing mask values for v points.
        mask_rho: `numpy.ndarray` containing mask values for rho points.
    """
    water_u = ma.masked_array(u_target_depth, numpy.logical_not(mask_u))
    water_v = ma.masked_array(v_target_depth, numpy.logical_not(mask_v))
    # u/v masked values need to be set to 0 for averaging 
    water_u = water_u.filled(0)
    water_v = water_v.filled(0)
    water_ang_rho = ma.masked_array(ang_rho, numpy.logical_not(mask_rho))
    water_lat_rho = ma.masked_array(lat_rho, numpy.logical_not(mask_rho))
    water_lon_rho = ma.masked_array(lon_rho, numpy.logical_not(mask_rho))

    return water_u, water_v, water_ang_rho, water_lat_rho, water_lon_rho


def vertical_interpolation(u, v, s_rho, mask_rho, mask_u, mask_v, zeta, h, hc, cs_r, vtransform, num_eta, num_xi, num_sigma, target_depth):
    """Vertically interpolate variables to target depth.

    Args:
        u: `numpy.ndarray` containing u values for entire grid.
        v: `numpy.ndarray` containing v values for entire grid.
        s_rho: `numpy.ndarray` s-coordinate at rho points, -1 min 0 max, positive up.
        zeta: `numpy.ndarray` containing MSL free surface at rho points in meters.
        h: `numpy.ndarray` containing bathymetry at rho points.
        hc: `numpy.ndarray` containing s-coordinate parameter, critical depth".
        cs_r: `numpy.ndarray` containing s-coordinate stretching curves at rho points.
        mask_u: `numpy.ndarray` containing mask values for u points.
        mask_v: `numpy.ndarray` containing mask values for v points.
        mask_rho: `numpy.ndarray` containing mask values for rho points.
        vtransform: vertical terrain-following transformation equation.
        num_eta: eta dimensions.
        num_xi: xi dimensions.
        num_sigma: sigma dimensions.
        target_depth: The water current at a specified target depth below the sea
            surface in meters, default target depth is 4.5 meters, target interpolation
            depth must be greater or equal to 0.
    """
    zeta = ma.masked_array(zeta, numpy.logical_not(mask_rho))
    total_depth = h + zeta
    z = numpy.ma.empty(shape=[num_sigma, num_eta, num_xi])

    if vtransform == 1:
        # Roms vertical transformation equation 1
        for k in range(num_sigma):
            s = ((hc * s_rho[k]) + (h - hc) * cs_r[k])
            z[k, :, :] = s + zeta * (1 + s/h)
    else:
        # Roms vertical transformation equation 2 GOMOFS Only
        for k in range(num_sigma):
            s = ((hc * s_rho[k]) + (h*cs_r[k]))/(h + hc)
            z[k, :, :] = zeta + ([zeta + h]*s)

    if target_depth < 0:
        raise Exception("Target depth must be positive")
    if target_depth > numpy.nanmax(total_depth):
        raise Exception("Target depth exceeds total depth")

    # For areas shallower than the target depth, depth is half the total depth
    interp_depth = zeta - numpy.minimum(target_depth*2, total_depth)/2

    u_target_depth = numpy.ma.empty(shape=[num_eta, num_xi])
    v_target_depth = numpy.ma.empty(shape=[num_eta, num_xi])
    # Perform vertical linear interpolation on u/v values to target depth
    for eta in range(num_eta):
        for xi in range(num_xi):
            if not zeta.mask[eta, xi]:
                if mask_u[eta, xi] != 0:
                    u_interp_depth = interpolate.interp1d(z[:, eta, xi], u[:, eta, xi], fill_value='extrapolate')
                    u_target_depth[eta, xi] = u_interp_depth(interp_depth.data[eta, xi])

                if mask_v[eta, xi] != 0:
                    v_interp_depth = interpolate.interp1d(z[:, eta, xi], v[:, eta, xi], fill_value='extrapolate')
                    v_target_depth[eta, xi] = v_interp_depth(interp_depth.data[eta, xi])

    return u_target_depth, v_target_depth
