"""
Utility classes and methods for working with numerical ocean models
and deriving water currents.


This module provides functionality allowing ocean models
to be interpolated to a regular, orthogonal lat/lon horizontal grid.

"""

import math
import os

import gdal
import json
import netCDF4
import numpy
import osr
import ogr
from shapely.geometry import shape
from scipy.interpolate import griddata

from s100.regulargrid import RegularGrid

# Conversion factor for meters/sec to knots
MS2KNOTS = 1.943844

# Default fill value for NetCDF variables
FILLVALUE = -9999.0

# Default module for horizontal interpolation
INTERP_METHOD_SCIPY = "scipy"

# Alternative module for horizontal interpolation
INTERP_METHOD_GDAL = "gdal"


class ModelIndexFile:
    """Store a regular grid numpy.array and mask for model of interest in a NetCDF file.
        Created once per model, resolution, model domain extent, and subset model extents.
        These values are stored in the output index file for use when regridding native model
        output to the regular grid. As long as the regular grid remains the same, these values
        will not change, thus the index file only needs to be generated once per model
        and kept on the data processing system in perpetuity.

        Option: Store subset information to index file to be used to subset model
        into multiple output files.
    
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

    def __init__(self, path):
        """Initialize ModelIndexFile object and open file at specified path.

        If file already exists and `clobber==False`, it will be opened in
        read mode. Otherwise, it will be opened in write mode.
        
        Args:
            path: Path of target NetCDF file.
        """
        self.path = path
        self.nc_file = None
        self.dim_y = None
        self.dim_x = None
        self.dim_subgrid = None
        self.var_subgrid_id = None
        self.var_subgrid_name = None
        self.var_subgrid_x_min = None
        self.var_subgrid_x_max = None
        self.var_subgrid_y_min = None
        self.var_subgrid_y_max = None
        self.var_y = None
        self.var_x = None
        self.var_mask = None

    def open(self, clobber=False):
        """Open netCDF file.

        Args:
            clobber: (Optional, default False) If True, existing netCDF file at
                specified path, if any, will be deleted and the new file will
                be opened in write mode.
        """
        # nc_file: Handle to `netCDF4.Dataset` instance for the opened NetCDF
        if not os.path.exists(self.path) or clobber:
            self.nc_file = netCDF4.Dataset(self.path, "w", format="NETCDF4")
        else:
            self.nc_file = netCDF4.Dataset(self.path, "r", format="NETCDF4")
            self.init_handles()

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
            self.var_subgrid_name = self.nc_file.variables['subgrid_name'][:]
        except KeyError:
            self.var_subgrid_name = None
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

    def create_dims_coord_vars(self, num_y, num_x):
        """Create index file NetCDF dimensions and coordinate variables.

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
        self.var_mask.flag_values = -9999.0, 1
        self.var_mask.flag_meanings = "land, water"

    def create_subgrid_dims_vars(self, num_subgrids, subset_grid_field_name=None):
        """Create subgrid-related NetCDF dimensions/variables.

        Args:
            num_subgrids: Number of subgrids.
            subset_grid_field_name: (Optional, default None) Shapefile field name
                to be stored in the index file.
        """
        self.dim_subgrid = self.nc_file.createDimension(self.DIMNAME_SUBGRID, num_subgrids)
        self.var_subgrid_id = self.nc_file.createVariable('subgrid_id', 'i4', (self.DIMNAME_SUBGRID,), fill_value=FILLVALUE)
        self.var_subgrid_x_min = self.nc_file.createVariable('subgrid_x_min', 'i4', (self.DIMNAME_SUBGRID,), fill_value=FILLVALUE)
        self.var_subgrid_x_max = self.nc_file.createVariable('subgrid_x_max', 'i4', (self.DIMNAME_SUBGRID,), fill_value=FILLVALUE)
        self.var_subgrid_y_min = self.nc_file.createVariable('subgrid_y_min', 'i4', (self.DIMNAME_SUBGRID,), fill_value=FILLVALUE)
        self.var_subgrid_y_max = self.nc_file.createVariable('subgrid_y_max', 'i4', (self.DIMNAME_SUBGRID,), fill_value=FILLVALUE)

        if subset_grid_field_name is not None:
            self.var_subgrid_name = self.nc_file.createVariable('subgrid_name', 'S30', (self.DIMNAME_SUBGRID,), fill_value=FILLVALUE)

    def init_nc(self, model_file, target_cellsize_meters, ofs_model, shoreline_shp=None, subset_grid_shp=None,
                subset_grid_field_name=None):
        """Initialize NetCDF dimensions/variables/attributes.

        Args:
            model_file: `ModelFile` instance containing model output used
                to identify properties of original grid.
            target_cellsize_meters: Target cell size of grid cells, in meters.
                Actual calculated x/y grid cell sizes will vary slightly from
                this value, since the regular grid uses lat/lon coordinates
                (thus a cell's width/height in meters will vary by latitude),
                and since it will be adjusted in order to fit a whole number of
                grid cells in the x and y directions within the calculated grid
                extent.
            ofs_model: The target model identifier.
            shoreline_shp: (Optional, default None) Path to a polygon shapefile
                containing features identifying land areas. If specified,
                a shoreline mask variable will be created/populated.
            subset_grid_shp: (Optional, default None) Path to a polygon
                shapefile containing orthogonal rectangles identifying areas
                to be used to subset (chop up) the full regular grid into
                tiles. Shapefile is assumed to be in the WGS84 projection. If
                None, the index file will be created assuming no subsets are
                desired and the extent of the model will be used instead.
            subset_grid_field_name: (Optional, default None) Shapefile
                field name to be stored in the index file.

        """
        # Calculate extent of valid (water) points
        (lon_min, lon_max, lat_min, lat_max) = model_file.get_valid_extent()

        # Populate grid x/y coordinate variables and subset-related variables
        # (if applicable)
        if subset_grid_shp is None:
            reg_grid = self.init_xy(lon_min, lat_min, lon_max, lat_max, target_cellsize_meters)
        elif subset_grid_field_name is None:
            reg_grid = self.init_xy_with_subsets(lon_min, lat_min, lon_max, lat_max, target_cellsize_meters,
                                                 subset_grid_shp)
        else:
            reg_grid = self.init_xy_with_subsets(lon_min, lat_min, lon_max, lat_max, target_cellsize_meters,
                                                 subset_grid_shp, subset_grid_field_name)

        self.nc_file.gridOriginLongitude = reg_grid.x_min
        self.nc_file.gridOriginLatitude = reg_grid.y_min

        land_mask = None
        if shoreline_shp is not None:
            land_mask = self.init_shoreline_mask(reg_grid, shoreline_shp)

        self.nc_file.model = str.upper(ofs_model)
        self.nc_file.format = "netCDF-4"

        print ("Full grid dimensions (y,x): ({},{})".format(len(reg_grid.y_coords), len(reg_grid.x_coords)))

        # Calculate the mask
        grid_cell_mask = self.compute_grid_mask(model_file, reg_grid)
        self.write_mask(land_mask, grid_cell_mask)

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
        self.create_dims_coord_vars(len(reg_grid.y_coords), len(reg_grid.x_coords))
        
        # Populate NetCDF coordinate variables using regular grid coordinates
        self.var_x[:] = reg_grid.x_coords[:]
        self.var_y[:] = reg_grid.y_coords[:]
        self.nc_file.gridSpacingLongitude = cellsize_x
        self.nc_file.gridSpacingLatitude = cellsize_y

        return reg_grid

    def init_xy_with_subsets(self, lon_min, lat_min, lon_max, lat_max, target_cellsize_meters, subset_grid_shp, subset_grid_field_name=None):
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
            subset_grid_field_name: Optional, default None) Shapefile
                field name to be stored in the index file.

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
        spatial_ref = layer.GetSpatialRef()
        shp_srs = spatial_ref.GetAttrValue("AUTHORITY", 1)
        source = osr.SpatialReference()
        source.ImportFromEPSG(int(shp_srs))
        target = osr.SpatialReference()
        target.ImportFromEPSG(4326)
        transform = osr.CoordinateTransformation(source, target)
        ofs_poly.Transform(transform)

        # Find the intersection between grid polygon and ocean model grid extent
        subset_polys = {}
        fids = []
        fields = {}
        fid = 0
        for feature in layer:
            geom = feature.GetGeometryRef()
            intersection = ofs_poly.Intersection(geom)
            if intersection.ExportToWkt() != "GEOMETRYCOLLECTION EMPTY":
                subset_polys[fid] = geom.ExportToJson()
                if subset_grid_field_name is not None:
                    field_name = feature.GetField(str(subset_grid_field_name))
                    fields.update({fid: field_name})
                    fids.append(fid)
                else:
                    fids.append(fid)
            fid += 1

        if len(fids) == 0:
            raise Exception("Given subset grid shapefile contains no polygons that intersect with model domain; cannot proceed.")

        # Use a single subset polygon to calculate x/y cell sizes. This ensures
        # that cells do not fall on the border between two grid polygons.
        single_polygon = ogr.Geometry(ogr.wkbMultiPolygon)
        single_polygon.AddGeometry(ogr.CreateGeometryFromJson(subset_polys[fids[0]]))
        sp_x_min, sp_x_max, sp_y_min, sp_y_max = single_polygon.GetEnvelope()

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
        self.create_dims_coord_vars(len(full_reg_grid.y_coords), len(full_reg_grid.x_coords))
        # Populate NetCDF coordinate variables using regular grid coordinates
        self.var_x[:] = full_reg_grid.x_coords[:]
        self.var_y[:] = full_reg_grid.y_coords[:]
        self.nc_file.gridSpacingLongitude = full_reg_grid.cellsize_x
        self.nc_file.gridSpacingLatitude = full_reg_grid.cellsize_y

        # Create subgrid dimension/variables
        self.create_subgrid_dims_vars(len(subset_polys), subset_grid_field_name)
        # Calculate subgrid mask ranges, populate subgrid ID
        for subgrid_index, fid in enumerate(fids):
            self.var_subgrid_id[subgrid_index] = fid
            if subset_grid_field_name is not None:
                self.var_subgrid_name[subgrid_index] = fields[fid]

            # Convert OGR geometry to shapely geometry
            subset_poly_shape = shape(json.loads(subset_polys[fid]))
            min_x_coord = subset_poly_shape.bounds[0]
            max_x_coord = subset_poly_shape.bounds[2]
            min_y_coord = subset_poly_shape.bounds[1]
            max_y_coord = subset_poly_shape.bounds[3]

            subgrid_x_min = None
            subgrid_x_max = None
            subgrid_y_min = None
            subgrid_y_max = None

            for i, x in enumerate(self.var_x):
                if x >= min_x_coord:
                    subgrid_x_min = i
                    break
            count_x = round(((max_x_coord - min_x_coord) / full_reg_grid.cellsize_x))

            for i, y in enumerate(self.var_y):
                if y >= min_y_coord:
                    subgrid_y_min = i
                    break
            count_y = round(((max_y_coord - min_y_coord) / full_reg_grid.cellsize_y))

            subgrid_x_max = subgrid_x_min + count_x - 1
            subgrid_y_max = subgrid_y_min + count_y - 1

            self.var_subgrid_x_min[subgrid_index] = subgrid_x_min
            self.var_subgrid_x_max[subgrid_index] = subgrid_x_max
            self.var_subgrid_y_min[subgrid_index] = subgrid_y_min
            self.var_subgrid_y_max[subgrid_index] = subgrid_y_max

        return full_reg_grid

    @staticmethod
    def init_shoreline_mask(reg_grid, shoreline_shp):
        """Create a shoreline mask for the region of interest at the
           target resolution.
        
        Args:
            reg_grid: `RegularGrid` instance describing the regular grid for
                which the shoreline mask will be created.
            shoreline_shp: Path to a polygon shapefile containing features
                identifying land areas.

        Returns:
            2D numpy array, matching the dimensions of the given RegularGrid,
            containing a value of 0 for water areas and a value of 255 for land
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
        spatial_ref = layer.GetSpatialRef()
        shp_srs = spatial_ref.GetAttrValue("AUTHORITY", 1)
        source = osr.SpatialReference()
        source.ImportFromEPSG(int(shp_srs))
        target = osr.SpatialReference()
        target.ImportFromEPSG(4326)
        transform = osr.CoordinateTransformation(source, target)
        ofs_poly.Transform(transform)

        # Find the intersection between shoreline shapefile and ocean
        # model grid extent
        for feature in layer:
            geom = feature.GetGeometryRef()
            intersection = ofs_poly.Intersection(geom)

        # Create a regional ogr shoreline polygon layer
        # write to memory
        driver = ogr.GetDriverByName('Memory')
        mem_dset = driver.CreateDataSource('mem_dst')
        dest_srs = ogr.osr.SpatialReference()
        dest_srs.ImportFromEPSG(4326)
        lyr = mem_dset.CreateLayer('land_mask', dest_srs, geom_type=ogr.wkbPolygon)
        feat = ogr.Feature(lyr.GetLayerDefn())
        feat.SetGeometry(intersection)
        lyr.CreateFeature(feat)

        # Rasterize the shoreline polygon layer and write to memory
        pixel_width = reg_grid.cellsize_x
        pixel_height = reg_grid.cellsize_y
        cols = len(reg_grid.y_coords)
        rows = len(reg_grid.x_coords)
        target_dset = gdal.GetDriverByName('MEM').Create("land_mask.tif", rows, cols, 1, gdal.GDT_Byte)
        target_dset.SetGeoTransform((reg_grid.x_min, pixel_width, 0, reg_grid.y_min, 0,  pixel_height))
        target_dset_srs = osr.SpatialReference()
        target_dset_srs.ImportFromEPSG(4326)
        target_dset.SetProjection(target_dset_srs.ExportToWkt())
        band = target_dset.GetRasterBand(1)
        band.SetNoDataValue(FILLVALUE)
        band.FlushCache()

        gdal.RasterizeLayer(target_dset, [1], lyr, burn_values=[2])

        # Store as a numpy array, land = 2, invalid areas = 0
        target_band = target_dset.GetRasterBand(1)
        land_mask = target_band.ReadAsArray(pixel_width, pixel_height, rows, cols).astype(numpy.int)

        return land_mask

    def compute_grid_mask(self, model_file, reg_grid):
        """Create regular grid model domain mask.

        Args:
            model_file: 'ModelFile' instance containing model irregular grid
                structure and variables.
            reg_grid: `RegularGrid` instance describing the regular grid for
                which the mask will be created.
        """
        pass

    @staticmethod
    def rasterize_mask(reg_grid, layer):
        """Rasterize model domain mask from irregular grid cell polygons.

                Args:
                    reg_grid: `RegularGrid` instance describing the regular grid for
                        which the mask will be created.
                    layer: Derived from modeling framework module(i.e. roms, fvcom, etc.)
                        contains irregular or unstructured grid cell polygons.
                """
        # Rasterize the grid cell polygon layer and write to memory
        pixel_width = reg_grid.cellsize_x
        pixel_height = reg_grid.cellsize_y
        cols = len(reg_grid.y_coords)
        rows = len(reg_grid.x_coords)
        target_dset = gdal.GetDriverByName('MEM').Create("grid_cell_mask.tif", rows, cols, 1, gdal.GDT_Byte)
        target_dset.SetGeoTransform((reg_grid.x_min, pixel_width, 0, reg_grid.y_min, 0, pixel_height))
        target_dset_srs = osr.SpatialReference()
        target_dset_srs.ImportFromEPSG(4326)
        target_dset.SetProjection(target_dset_srs.ExportToWkt())
        band = target_dset.GetRasterBand(1)
        band.SetNoDataValue(FILLVALUE)
        band.FlushCache()

        gdal.RasterizeLayer(target_dset, [1], layer, burn_values=[1])

        # Store as numpy array, valid areas = 1, invalid areas = 0
        target_band = target_dset.GetRasterBand(1)
        grid_cell_mask = target_band.ReadAsArray(pixel_width, pixel_height, rows, cols).astype(numpy.int)

        return grid_cell_mask

    def write_mask(self, land_mask, grid_cell_mask):
        """Write master mask to index file.

        Args:
            grid_cell_mask: 2D numpy array, matching the index file dimensions,
                containing regular grid model domain mask values, a value of 1
                for valid water, a value of 0 for invalid areas.
            land_mask: 2D numpy array, matching the index file dimensions,
                containing regular grid masked values, a value of 0 for invalid
                areas and a value of 2 for land.
        """
        # Use land mask and grid cell mask to create a master mask
        # Value of 1 for valid areas, FILLVALUE for invalid areas

        # Write mask to index file
        for y in range(self.dim_y.size):
            for x in range(self.dim_x.size):
                if land_mask[y, x] != 0:
                    continue
                if grid_cell_mask[y, x] == 1:
                    self.var_mask[y, x] = 1
                else:
                    self.var_mask[y, x] = FILLVALUE


class ModelFile:
    """Read/process data from a numerical ocean model file.

    This is a parent class that should be inherited from. This
    class opens a NetCDF model file, reads variables, gets model
    domain extent and converts model variables u/v to a regular
    grid.

    """
    def __init__(self, path):
        """Initialize model file object and opens file at specified path.

        Args:
            path: Path of target NetCDF file.
        
        Raises:
            Exception: Specified NetCDF file does not exist.
        """
        self.path = path
        self.nc_file = None

    def open(self):
        if os.path.exists(self.path):
            self.nc_file = netCDF4.Dataset(self.path, "r", format="NETCDF3_Classic")
            self.init_handles()
        else:
            # File doesn't exist, raise error
            raise(Exception("NetCDF file does not exist: {}".format(self.path)))

    def close(self):
        self.nc_file.close()

    def get_valid_extent(self):
        pass

    def init_handles(self):
        pass

    def uv_to_regular_grid(self, model_index, target_depth, interp=None):
        """Execute functions to process model variables to a regular grid, model specific
           functions and interpolation method derived from modeling framework module of
           interest(i.e. roms, fvcom, etc.).
        """
        pass


def uv_to_speed_direction(reg_grid_u, reg_grid_v):
    """Convert u and v vectors to speed/direction.

    Input u/v values are assumed to be in meters/sec. Output speed values will
    be converted to knots and direction in degrees from true north (0-360).

    Args:
        reg_grid_u: `numpy.ma.masked_array` containing u values interpolated to
            the regular grid.
        reg_grid_v: `numpy.ma.masked_array` containing v values interpolated to
            the regular grid.

    Returns:
        Two tuple of speed and direction
    """
    direction = numpy.ma.empty((reg_grid_u.shape[0], reg_grid_u.shape[1]), dtype=numpy.float32)
    speed = numpy.ma.empty((reg_grid_u.shape[0], reg_grid_u.shape[1]), dtype=numpy.float32)
    for y in range(reg_grid_u.shape[0]):
        for x in range(reg_grid_u.shape[1]):
            if reg_grid_u.mask[y, x]:
                direction[y, x] = numpy.nan
                speed[y, x] = numpy.nan
                continue

            u_ms = reg_grid_u[y, x]
            v_ms = reg_grid_v[y, x]

            # Convert from meters per second to knots
            u_knot = u_ms * MS2KNOTS
            v_knot = v_ms * MS2KNOTS

            current_speed = math.sqrt(math.pow(u_knot, 2) + math.pow(v_knot, 2))
            current_direction_radians = math.atan2(v_knot, u_knot)
            current_direction_degrees = math.degrees(current_direction_radians)
            current_direction_north = 90.0 - current_direction_degrees

            # The direction must always be positive.
            if current_direction_north < 0.0:
                current_direction_north += 360.0

            direction[y, x] = current_direction_north
            speed[y, x] = current_speed

    return direction, speed


def scipy_interpolate_uv_to_regular_grid(u, v, lat, lon, model_index):
    """Create a regular grid using masked latitude and longitude variables.

       Interpolate u/v variables to the regular grid using Linear Interpolation.

       Args:
           u: `numpy.ma.masked_array` containing u values with NoData/land values
               masked out.
           v: `numpy.ma.masked_array` containing v values with NoData/land values
               masked out.
           lat: `numpy.ndarray` containing masked latitude values of rho
               or centroid points.
           lon: `numpy.ndarray` containing masked longitude values of rho
               or centroid points.
           model_index: `ModelIndexFile` from which index values will
               be extracted to perform interpolation.
       """

    # Using scipy to interpolate irregular spaced u/v points to a regular grid
    x, y = numpy.meshgrid(model_index.var_x, model_index.var_y)
    coords = numpy.column_stack((lon, lat))
    reg_grid_u = griddata(coords, u, (x, y), method='linear')
    reg_grid_v = griddata(coords, v, (x, y), method='linear')

    return reg_grid_u, reg_grid_v


def gdal_interpolate_uv_to_regular_grid(u, v, lat, lon, model_index):
    """Create a regular grid using masked latitude and longitude variables.

    Interpolate u/v variables to the regular grid using Linear Interpolation.

    Args:
        u: `numpy.ma.masked_array` containing u values with NoData/land values
           masked out.
        v: `numpy.ma.masked_array` containing v values with NoData/land values
           masked out.
        lat: `numpy.ndarray` containing masked latitude values of rho
           or centroid points.
        lon: `numpy.ndarray` containing masked longitude values of rho
           or centroid points.
        model_index: `ModelIndexFile` from which index values will
           be extracted to perform interpolation.
    """
    # Create an ogr object containing irregularly spaced points for u,v,lat,lon
    # and write to memory
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS("WGS84")
    ds = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Float32)
    layer = ds.CreateLayer("irregular_points", srs=srs, geom_type=ogr.wkbPoint)
    layer.CreateField(ogr.FieldDefn("u", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("v", ogr.OFTReal))

    for i in range(len(v)):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(lon[i], lat[i])
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(point)
        feature.SetField("u", u[i])
        feature.SetField("v", v[i])
        layer.CreateFeature(feature)

    # Input ogr object to gdal grid and interpolate irregularly spaced u/v
    # to a regular grid
    dst_u = gdal.Grid('u.tif', ds, format='MEM', width=model_index.dim_x.size, height=model_index.dim_y.size,
                      algorithm="linear:nodata=0.0", zfield="u")
    dst_v = gdal.Grid('v.tif', ds, format='MEM', width=model_index.dim_x.size, height=model_index.dim_y.size,
                      algorithm="linear:nodata=0.0", zfield="v")

    reg_grid_u = dst_u.ReadAsArray()
    reg_grid_v = dst_v.ReadAsArray()

    return reg_grid_u, reg_grid_v

