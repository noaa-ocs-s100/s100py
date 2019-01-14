"""
Utility classes and methods for working with Global HYCOM ocean model forecast guidance.

The Hybrid Coordinate Ocean Model (HYCOM) is a three-dimensional, primitive
equation, free surface ocean,  numerical ocean model with curvilinear orthogonal
and sigma (terrain-following)coordinates.

"""

import numpy
import numpy.ma as ma
import datetime
import netCDF4
import osr
import ogr
from scipy import interpolate

from s100py.model import model

# Default fill value for NetCDF variables
FILLVALUE = -9999.0

# Default module for horizontal interpolation
INTERP_METHOD_SCIPY = "scipy"

# Alternative module for horizontal interpolation
INTERP_METHOD_GDAL = "gdal"


class HYCOMIndexFile(model.ModelIndexFile):
    """Store a regular grid mask based on HYCOM model grid properties in a NetCDF file."""

    def __init__(self, path):
        super().__init__(path)

    def compute_grid_mask(self, model_file, reg_grid):
        """Create model domain mask and write to index file.

        Args:
            model_file: `HYCOMOutputFile` instance containing irregular grid
                structure and variables.
            reg_grid: `RegularGrid` instance describing the regular grid for
                which the mask will be created.
        """
        # Create shapefile with OGR
        driver = ogr.GetDriverByName('Esri Shapefile')
        dset = driver.CreateDataSource('grid_cell_mask.shp')
        dset_srs = ogr.osr.SpatialReference()
        dset_srs.ImportFromEPSG(4326)
        layer = dset.CreateLayer('', dset_srs, ogr.wkbMultiPolygon)
        layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))

        # Add spatial reference to polygon
        spatial_ref = osr.SpatialReference()
        spatial_ref.ImportFromEPSG(4326)
        spatial_ref.MorphToESRI()
        mask_file = open('grid_cell_mask.prj', 'w')
        mask_file.write(spatial_ref.ExportToWkt())
        mask_file.close()

        # Create shapefile containing polygons for each irregular grid cell
        # searching counter clockwise(eta1, xi), (eta2, xi2), (eta3, xi3),(eta4, xi4)
        for nx1 in range(model_file.num_x-1):
            for ny1 in range(model_file.num_y-1):
                if nx1 == 0:
                    continue
                if ny1 == 0:
                    continue
                if ny1 == model_file.num_y-1:
                    continue
                if nx1 == model_file.num_x-1:
                    continue
                nx2 = nx1 + 1
                ny2 = ny1
                nx3 = nx1 + 1
                ny3 = ny1 + 1
                nx4 = nx1
                ny4 = ny1 + 1

                # Search valid points
                valid_points = []
                for (ny, nx) in ((ny1, nx1), (ny2, nx2), (ny3, nx3), (ny4, nx4)):
                    if model_file.var_u.mask[0, 0, ny, nx]:
                        continue
                    valid_points.append((ny, nx))
                if len(valid_points) < 3:
                    continue
                ring = ogr.Geometry(ogr.wkbLinearRing)
                for (ny, nx) in valid_points:
                    ring.AddPoint(model_file.var_lon[ny, nx], model_file.var_lat[ny, nx])
                (ny, nx) = valid_points[0]
                ring.AddPoint(model_file.var_lon[ny, nx], model_file.var_lat[ny, nx])

                # Create polygon
                geom = ogr.Geometry(ogr.wkbPolygon)
                geom.AddGeometry(ring)
                feat = ogr.Feature(layer.GetLayerDefn())
                feat.SetField('id', 1)
                feat.SetGeometry(geom)
                layer.CreateFeature(feat)

        return self.rasterize_mask(reg_grid, layer)


class HYCOMFile(model.ModelFile):
    """Read/process data from a HYCOM model NetCDF file.

    Attributes:
        path: Path (relative or absolute) of the file.
    """

    def __init__(self, path):
        """Initialize HYCOM file object and open file at specified path.

        Args:
            path: Path of target NetCDF file.

        """
        super().__init__(path)
        self.var_lat = None
        self.var_lon = None
        self.var_u = None
        self.var_v = None
        self.var_date = None
        self.var_x = None
        self.var_y = None
        self.var_depth = None
        self.var_time = None
        self.datetime_values = None
        self.var_mask = None
        self.num_x = None
        self.num_y = None

    def close(self):
        super().close()
        self.release_resources()

    def release_resources(self):
        self.var_lat = None
        self.var_lon = None
        self.var_u = None
        self.var_v = None
        self.var_date = None
        self.var_x = None
        self.var_y = None
        self.var_depth = None
        self.var_time = None
        self.datetime_values = None
        self.var_mask = None
        self.num_x = None
        self.num_y = None

    def get_valid_extent(self):
        """Masked model domain extent."""
        lon_min = numpy.nanmin(self.var_lon)
        lon_max = numpy.nanmax(self.var_lon)
        lat_min = numpy.nanmin(self.var_lat)
        lat_max = numpy.nanmax(self.var_lat)

        return lon_min, lon_max, lat_min, lat_max

    def init_handles(self):
        """Initialize handles to NetCDF variables."""
        self.var_lat = self.nc_file.variables['Latitude'][:, :]
        self.var_lon = self.nc_file.variables['Longitude'][:, :]
        self.var_lat = self.var_lat.astype(numpy.float64)
        self.var_lon = self.var_lon.astype(numpy.float64)
        self.var_u = self.nc_file.variables['u'][:, :, :, :]
        self.var_v = self.nc_file.variables['v'][:, :, :, :]
        self.var_date = self.nc_file.variables['Date'][:]
        self.var_x = self.nc_file.variables['X'][:]
        self.var_y = self.nc_file.variables['Y'][:]
        self.var_depth = self.nc_file.variables['Depth'][:]
        self.var_time = self.nc_file.variables['MT'][:]
        self.num_x = len(self.var_x)
        self.num_y = len(self.var_y)
        # Use the surface layer in the u variable to define a land mask
        self.var_mask = self.nc_file.variables['u'][0, 0, :, :]

        # Convert timestamps to datetime objects and store in a list
        # Rounding to the nearest hour
        self.datetime_values = []
        for time_index in range(len(self.var_time)):
            self.var_time = netCDF4.num2date(self.nc_file.variables['MT'][:], self.nc_file.variables['MT'].units)[time_index]
            print(self.var_time)
            if self.var_time.minute >= 30:
                # round up
                adjusted_time = datetime.datetime(self.var_time.year, self.var_time.month, self.var_time.day, self.var_time.hour, 0, 0) + datetime.timedelta(hours=1)
            elif self.var_time.minute < 30:
                # round down
                adjusted_time = datetime.datetime(self.var_time.year, self.var_time.month, self.var_time.day, self.var_time.hour, 0, 0)

            self.datetime_values.append(adjusted_time)

    def uv_to_regular_grid(self, model_index, time_index, target_depth, interp=INTERP_METHOD_SCIPY):
        """Call grid processing functions and interpolate u/v to a regular grid"""

        u_target_depth, v_target_depth = vertical_interpolation(self.var_u, self.var_v, self.var_depth, self.num_x, self.num_y, time_index, target_depth)

        u_compressed, v_compressed, lat_compressed, lon_compressed = compress_variables(u_target_depth, v_target_depth, self.var_lat, self.var_lon, self.var_mask)

        # Scipy interpolation is default method, change method parameter to change interpolation method
        if interp == model.INTERP_METHOD_SCIPY:
            return model.scipy_interpolate_uv_to_regular_grid(u_compressed, v_compressed, lat_compressed, lon_compressed, model_index)
        elif interp == model.INTERP_METHOD_GDAL:
            return model.gdal_interpolate_uv_to_regular_grid(u_compressed, v_compressed, lat_compressed, lon_compressed, model_index)


def compress_variables(u_target_depth, v_target_depth, lat, lon, mask):
    """Compress masked variables for interpolation.

    Args:
        u_target_depth: `numpy.ma.masked_array` containing u values at target depth.
        v_target_depth: `numpy.ma.masked_array` containing v values at target depth.
        lat: `numpy.ma.masked_array` containing latitude values.
        lon: `numpy.ma.masked_array` containing longitude values.
        mask: `numpy.ma.masked_array` containing mask values.
    """
    water_lat_rho = ma.masked_array(lat, mask.mask)
    water_lon_rho = ma.masked_array(lon, mask.mask)
    water_u = ma.masked_array(u_target_depth, mask.mask)
    water_v = ma.masked_array(v_target_depth, mask.mask)

    u_compressed = ma.compressed(water_u)
    v_compressed = ma.compressed(water_v)
    lat_compressed = ma.compressed(water_lat_rho)
    lon_compressed = ma.compressed(water_lon_rho)

    return u_compressed, v_compressed, lat_compressed, lon_compressed


def vertical_interpolation(u, v, depth, num_x, num_y, time_index, target_depth):
    """Vertically interpolate variables to target depth.

    Args:
        u: `numpy.ndarray` containing u values for entire grid.
        v: `numpy.ndarray` containing v values for entire grid.
        depth: `numpy.1darray` containing depth in meters, positive down.
        time_index: Single forecast time index value.
        num_x: X dimensions
        num_y: Y dimensions
        target_depth: The water current at a specified target depth below the sea
            surface in meters, default target depth is 4.5 meters, target interpolation
            depth must be greater or equal to 0.
    """
    if target_depth < 0:
        raise Exception("Target depth must be positive")
    if target_depth > numpy.nanmax(depth):
        raise Exception("Target depth exceeds total depth")

    u_target_depth = numpy.ma.empty(shape=[num_y, num_x])
    v_target_depth = numpy.ma.empty(shape=[num_y, num_x])

    # Determine the total depth at each valid u/v location
    for ny in range(num_y):
        for nx in range(num_x):
            deepest_valid_depth_layer_index = None
            for i in range(len(depth)):
                if u.mask[time_index, i, ny, nx] or v.mask[time_index, i, ny, nx]:
                    break
                deepest_valid_depth_layer_index = i
            if deepest_valid_depth_layer_index is None:
                continue

            # For areas shallower than the target depth, depth is half the total depth
            interp_depth = numpy.minimum(target_depth * 2, numpy.max(depth[deepest_valid_depth_layer_index])) / 2

            # Perform vertical linear interpolation on u/v values to target depth
            u_interp_depth = interpolate.interp1d(depth[:deepest_valid_depth_layer_index], u[time_index, :deepest_valid_depth_layer_index, ny, nx], fill_value='extrapolate')
            u_target_depth[ny, nx] = u_interp_depth(interp_depth)

            v_interp_depth = interpolate.interp1d(depth[:deepest_valid_depth_layer_index], v[time_index, :deepest_valid_depth_layer_index, ny, nx], fill_value='extrapolate')
            v_target_depth[ny, nx] = v_interp_depth(interp_depth)

    return u_target_depth, v_target_depth
