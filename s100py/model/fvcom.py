"""
Utility classes and methods for working with FVCOM ocean model forecast guidance.

The Finite-Volume, primitive equation Community Ocean Modeling System (FVCOM)
is a unstructured-grid, finite-volume, free-surface, 3D primitive equation
coastal ocean circulation model. The horizontal grid is comprised of unstructured
triangular cells and a generalized(bathymetry-following) vertical coordinate system.
This module provides functionality allowing FVCOM output to be interpolated to
a regular, orthogonal lat/lon horizontal grid at a given depth-below-surface.
"""

import numpy
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


class FVCOMIndexFile(model.ModelIndexFile):
    """Store a regular grid mask based on FVCOM model grid properties in a NetCDF file."""
    def __init__(self, path):
        super().__init__(path)

    def init_nc(self, model_file, target_cellsize_meters, ofs_model, model_type, shoreline_shp=None, subset_grid_shp=None, subset_grid_field_name=None):
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
            model_type: The target model type.
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
        super().init_nc(model_file, target_cellsize_meters, ofs_model, model_type, shoreline_shp, subset_grid_shp, subset_grid_field_name)
        # Determine vertical coordinate type
        (vertical_coordinates) = model_file.get_vertical_coordinate_type()
        self.nc_file.modelVerticalCoordinates = vertical_coordinates

        # For FVCOM models horizontally interpolate sigma values to the centroid
        (siglay_centroid, latc, lonc, num_nele, num_siglay) = model_file.sigma_to_centroid(vertical_coordinates)

        self.nc_file.createDimension("nele", num_nele)
        self.nc_file.createDimension("siglay", num_siglay)

        var_lonc = self.nc_file.createVariable('lon_centroid', 'f4', 'nele', fill_value=-9999)
        var_lonc.long_name = "longitude of unstructured centroid point"
        var_lonc.units = "degree_east"
        var_lonc.standard_name = "longitude"

        var_latc = self.nc_file.createVariable('lat_centroid', 'f4', 'nele', fill_value=-9999)
        var_latc.long_name = "latitude of unstructured centroid point"
        var_latc.units = "degree_north"
        var_latc.standard_name = "latitude"

        var_siglay_centroid = self.nc_file.createVariable('siglay_centroid', 'f4', ("siglay", "nele",), fill_value=FILLVALUE)
        var_siglay_centroid.long_name = "sigma layer at unstructured centroid point"
        var_siglay_centroid.coordinates = "lat_centroid lon_centroid"

        var_siglay_centroid[:, :] = siglay_centroid
        var_latc[:] = latc
        var_lonc[:] = lonc

    def compute_grid_mask(self, model_file, reg_grid):
        """Create model domain mask and write to index file.

        For every centroid create a polygon from three valid node
        points. Rasterize the polygon to create a grid domain mask.

        Args:
            model_file: `FVCOMFile` instance containing irregular grid
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

        # Create shapefile containing polygons for each unstructured
        # triangle using three valid nodes surrounding each centroid
        for node in range(0, model_file.var_nv.shape[1]):
            p1 = model_file.var_nv[0][node] - 1
            p2 = model_file.var_nv[1][node] - 1
            p3 = model_file.var_nv[2][node] - 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(model_file.var_lon_nodal[p1], model_file.var_lat_nodal[p1])
            ring.AddPoint(model_file.var_lon_nodal[p2], model_file.var_lat_nodal[p2])
            ring.AddPoint(model_file.var_lon_nodal[p3], model_file.var_lat_nodal[p3])
            ring.AddPoint(model_file.var_lon_nodal[p1], model_file.var_lat_nodal[p1])
            # Create polygon
            geom = ogr.Geometry(ogr.wkbPolygon)
            geom.AddGeometry(ring)
            feat = ogr.Feature(layer.GetLayerDefn())
            feat.SetField('id', 1)
            feat.SetGeometry(geom)
            layer.CreateFeature(feat)

        return self.rasterize_mask(reg_grid, layer)


class FVCOMFile(model.ModelFile):
    """Read/process data from a FVCOM model NetCDF file.

    Attributes:
        path: Path (relative or absolute) of the file.
    """
    def __init__(self, path, lon_offset=-360):
        """Initialize FVCOM file object and open file at specified path.

        Args:
            path: Path of target NetCDF file.

        """
        super().__init__(path)
        self.lon_offset = lon_offset
        self.var_lat_nodal = None
        self.var_lon_nodal = None
        self.var_lat_centroid = None
        self.var_lon_centroid = None
        self.var_u = None
        self.var_v = None
        self.var_zeta = None
        self.var_h = None
        self.var_nv = None
        self.var_wet_cells = None
        self.var_time = None
        self.var_siglay = None
        self.var_siglev = None
        self.num_siglay = None
        self.num_siglev = None
        self.num_nele = None
        self.num_node = None
        self.num_times = None
        self.datetime_values = None

    def close(self):
        super().close()
        self.release_resources()

    def release_resources(self):
        self.var_lat_nodal = None
        self.var_lon_nodal = None
        self.var_lat_centroid = None
        self.var_lon_centroid = None
        self.var_u = None
        self.var_v = None
        self.var_zeta = None
        self.var_h = None
        self.var_nv = None
        self.var_wet_cells = None
        self.var_time = None
        self.var_siglay = None
        self.var_siglev = None
        self.num_siglay = None
        self.num_siglev = None
        self.num_nele = None
        self.num_node = None
        self.num_times = None
        self.datetime_values = None

    def get_valid_extent(self):
        """Masked model domain extent."""
        lon_min = min(self.var_lon_centroid)
        lon_max = max(self.var_lon_centroid)
        lat_min = min(self.var_lat_centroid)
        lat_max = max(self.var_lat_centroid)

        return lon_min, lon_max, lat_min, lat_max

    def init_handles(self):
        """Initialize handles to NetCDF variables."""
        self.var_lat_nodal = self.nc_file.variables['lat'][:]
        self.var_lon_nodal = self.nc_file.variables['lon'][:] + self.lon_offset
        self.var_lat_nodal = self.var_lat_nodal.astype(numpy.float64)
        self.var_lon_nodal = self.var_lon_nodal.astype(numpy.float64)
        self.var_lat_centroid = self.nc_file.variables['latc'][:]
        self.var_lon_centroid = self.nc_file.variables['lonc'][:] + self.lon_offset
        self.var_lat_centroid = self.var_lat_centroid.astype(numpy.float64)
        self.var_lon_centroid = self.var_lon_centroid.astype(numpy.float64)
        self.var_u = self.nc_file.variables['u'][:, :, :]
        self.var_v = self.nc_file.variables['v'][:, :, :]
        self.var_zeta = self.nc_file.variables['zeta'][:, :]
        self.var_siglay = self.nc_file.variables['siglay'][:, :]
        self.var_siglev = self.nc_file.variables['siglev'][:, :]
        self.var_h = self.nc_file.variables['h'][:]
        self.var_nv = self.nc_file.variables['nv'][:, :]
        self.var_wet_cells = self.nc_file.variables['wet_cells'][:, :]
        self.num_node = len(self.var_h)
        self.num_nele = len(self.var_lat_centroid)
        self.num_siglay = self.var_siglay.shape[0]
        self.num_siglev = self.var_siglev.shape[0]
        self.num_times = self.var_u.shape[0]

        # Convert timestamp to datetime object
        self.datetime_values = []
        for time_index in range(self.num_times):
            self.var_time = netCDF4.num2date(self.nc_file.variables['time'][:], self.nc_file.variables['time'].units)[time_index]

            if self.var_time.minute >= 30:
                # round up
                adjusted_time = datetime.datetime(self.var_time.year, self.var_time.month, self.var_time.day, self.var_time.hour, 0, 0) + datetime.timedelta(hours=1)
            elif self.var_time.minute < 30:
                # round down
                adjusted_time = datetime.datetime(self.var_time.year, self.var_time.month, self.var_time.day, self.var_time.hour, 0, 0)

            self.datetime_values.append(adjusted_time)

        # Determine if sigma values are positive up or down from netCDF metadata
        if self.nc_file.variables['siglay'].positive == 'down':
            self.var_siglay = self.var_siglay * -1

    def get_vertical_coordinate_type(self):
        """Determine FVCOM-based OFS vertical sigma coordinate type"""

        siglay_values = self.var_siglay[:, 0]
        vertical_coordinates = "UNIFORM"
        for i in range(self.var_siglay.shape[1]):
            for s in range(self.var_siglay.shape[0]):
                if self.var_siglay[s, i] != siglay_values[s]:
                    vertical_coordinates = "GENERALIZED"
                    break

        return vertical_coordinates

    def sigma_to_centroid(self, vertical_coordinates):
        """Horizontally interpolate FVCOM-based OFS vertical sigma coordinates to centroid

        Args:
            vertical_coordinates: model vertical coordinate type.
        """

        siglay_centroid = numpy.ma.empty(shape=[self.num_siglay, self.num_nele])

        # Sigma vertical coordinate type
        # Generalized or Uniform
        if vertical_coordinates == "GENERALIZED":
            for i in range(self.num_siglay):
                coords = numpy.column_stack((self.var_lon_nodal, self.var_lat_nodal))
                siglay_centroid[i, :] = interpolate.griddata(coords, self.var_siglay[i, :], (self.var_lon_centroid, self.var_lat_centroid), method='linear')
        else:
            for k in range(self.num_nele):
                    siglay_centroid[:, k] = self.var_siglay[:, 0]

        return siglay_centroid, self.var_lat_centroid, self.var_lon_centroid, self.num_nele, self.num_siglay

    def uv_to_regular_grid(self, model_index, time_index, target_depth, interp=model.INTERP_METHOD_SCIPY):
        """Call grid processing functions and interpolate u/v to a regular grid"""

        h_centroid, zeta_centroid = node_to_centroid(self.var_zeta, self.var_h, self.var_lon_nodal, self.var_lat_nodal,
                                                     self.var_lon_centroid, self.var_lat_centroid, time_index)

        u_target_depth, v_target_depth = vertical_interpolation(self.var_u, self.var_v, h_centroid, zeta_centroid,
                                                                model_index, self.num_nele, self.num_siglay, time_index,
                                                                target_depth)

        # Scipy interpolation is default method, change method parameter to change interpolation method
        if interp == model.INTERP_METHOD_SCIPY:
            return model.scipy_interpolate_uv_to_regular_grid(u_target_depth, v_target_depth, self.var_lat_centroid,
                                                              self.var_lon_centroid, model_index)
        elif interp == model.INTERP_METHOD_GDAL:
            return model.gdal_interpolate_uv_to_regular_grid(u_target_depth, v_target_depth, self.var_lat_centroid,
                                                             self.var_lon_centroid, model_index)


def vertical_interpolation(u, v, h, zeta, model_index, num_nele, num_siglay, time_index, target_depth):
    """Vertically interpolate variables to target depth.

    Args:
        u: `numpy.ndarray` containing u values for entire grid.
        v: `numpy.ndarray` containing v values for entire grid.
        h: `numpy.ndarray` containing bathymetry at centroid points.
        zeta: `numpy.ndarray` containing MSL free surface at centroid points in meters.
        model_index: `ModelIndexFile` instance representing model index file containing siglay.
        num_nele: Number of elements(centroid).
        num_siglay: Number of sigma layers.
        time_index: Single forecast time index value.
        target_depth: The water current at a specified target depth below the sea
            surface in meters, default target depth is 4.5 meters, target interpolation
            depth must be greater or equal to 0.
    """
    true_depth = zeta + h

    sigma_depth_layers = numpy.ma.empty(shape=[num_siglay, num_nele])
    siglay = model_index.nc_file.variables['siglay_centroid'][:, :]

    for k in range(num_siglay):
        sigma_depth_layers[k, :] = siglay[k, :] * true_depth

    if target_depth < 0:
        raise Exception("Target depth must be positive")
    if target_depth > numpy.nanmax(true_depth):
        raise Exception("Target depth exceeds total depth")

    # For areas shallower than the target depth, depth is half the total depth
    interp_depth = zeta - numpy.minimum(target_depth * 2, true_depth) / 2

    u_target_depth = numpy.ma.empty(shape=[num_nele])
    v_target_depth = numpy.ma.empty(shape=[num_nele])

    # Perform vertical linear interpolation on u/v values to target depth
    for nele in range(num_nele):
        u_interp_depth = interpolate.interp1d(sigma_depth_layers[:, nele], u[time_index, :, nele], fill_value='extrapolate')
        u_target_depth[nele] = u_interp_depth(interp_depth.data[nele])

        v_interp_depth = interpolate.interp1d(sigma_depth_layers[:, nele], v[time_index, :, nele], fill_value='extrapolate')
        v_target_depth[nele] = v_interp_depth(interp_depth.data[nele])

    return u_target_depth, v_target_depth


def node_to_centroid(zeta, h, lon_node, lat_node, lon_centroid, lat_centroid, time_index):
    """Horizontally interpolate variables at nodes to centroids(elements).

    Args:
        zeta: `numpy.ndarray` containing MSL free surface at nodal points in meters.
        h: `numpy.ndarray` containing bathymetry at nodal points.
        lon_node: `numpy.ndarray` containing nodal longitude.
        lat_node: `numpy.ndarray` containing nodal latitude.
        lon_centroid: `numpy.ndarray` containing centroid(elements) longitude.
        lat_centroid: `numpy.ndarray` containing centroid(elements) latitude.
        time_index: Single time index value.
    """
    coords = numpy.column_stack((lon_node, lat_node))
    h_centroid = interpolate.griddata(coords, h, (lon_centroid, lat_centroid), method='linear')
    zeta_centroid = interpolate.griddata(coords, zeta[time_index, :], (lon_centroid, lat_centroid), method='linear')

    return h_centroid, zeta_centroid

