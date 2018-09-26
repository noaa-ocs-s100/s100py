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
import math
from scipy import interpolate


from s100.model import model

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
        for node in range(0,model_file.var_nv.shape[1]):
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
            feat.SetField('id', node)
            feat.SetGeometry(geom)
            layer.CreateFeature(feat)

        return self.rasterize_mask(reg_grid, layer)

class FVCOMFile(model.ModelFile):
    """Read/process data from a FVCOM model NetCDF file.

    Attributes:
        path: Path (relative or absolute) of the file.
    """
    def __init__(self, path, lon_offset = -360):
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
        self.var_siglay = None
        self.num_siglay = None
        self.var_siglev = None
        self.num_siglev = None
        self.num_nele = None
        self.num_node = None
        self.wet_cells = None
        self.time_val = None

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
        self.var_siglay = None
        self.num_siglay = None
        self.var_siglev = None
        self.num_siglev = None
        self.num_nele = None
        self.num_node = None
        self.wet_cells = None
        self.time_val = None

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
        self.var_u = self.nc_file.variables['u'][0, :, :]
        self.var_v = self.nc_file.variables['v'][0, :, :]
        self.var_zeta = self.nc_file.variables['zeta'][0,:]
        self.var_siglay = self.nc_file.variables['siglay'][:,:]
        self.var_siglev = self.nc_file.variables['siglev'][:,:]
        self.var_h = self.nc_file.variables['h'][:]
        self.var_nv = self.nc_file.variables['nv'][:,:]
        self.wet_cells = self.nc_file.variables['wet_cells'][0,:]
        self.num_node = self.var_h.shape[0]
        self.num_nele = self.var_u.shape[1]
        self.num_siglay = self.var_siglay.shape[0]
        self.num_siglev = self.var_siglev.shape[0]

        # Convert timestamp to datetime object
        self.time_val = netCDF4.num2date(self.nc_file.variables['time'][:],self.nc_file.variables['time'].units)[0]

        if self.time_val.minute >= 30:
            # round up
            self.time_val = datetime.datetime(self.time_val.year, self.time_val.month, self.time_val.day,self.time_val.hour, 0, 0) + datetime.timedelta(hours=1)
        elif self.time_val.minute < 30:
            # round down
            self.time_val = datetime.datetime(self.time_val.year, self.time_val.month, self.time_val.day,self.time_val.hour, 0, 0)

    def get_vertical_coordinate_type(self):
        """Determine FVCOM-based OFS vertical sigma coordinate type"""

        siglay_values = self.var_siglay[:,0]
        vertical_coordinates = "uniform"
        for i in range(self.var_siglay.shape[1]):
            for s in range(self.var_siglay.shape[0]):
                if self.var_siglay[s,i] != siglay_values[s]:
                    vertical_coordinates = "generalized"
                    break

        return vertical_coordinates

    def uv_to_regular_grid(self, model_index, target_depth, interp=model.INTERP_METHOD_SCIPY):
        """Call grid processing functions and interpolate averaged, rotated u/v to a regular grid"""

        h_centroid, zeta_centroid, siglay_centroid = node_to_element(model_index, self.var_zeta, self.var_h,
                                                                              self.var_siglay, self.var_lon_nodal,
                                                                              self.var_lat_nodal, self.var_lon_centroid,
                                                                              self.var_lat_centroid, self.num_nele, self.num_siglay)

        u_target_depth, v_target_depth = vertical_interpolation(self.var_u, self.var_v, h_centroid, zeta_centroid,
                                                                siglay_centroid, self.num_nele, self.num_siglay,
                                                                target_depth)

        # Scipy interpolation is default method, change method parameter to change interpolation method
        if interp == model.INTERP_METHOD_SCIPY:
            return model.scipy_interpolate_uv_to_regular_grid(u_target_depth, v_target_depth, self.var_lat_centroid,
                                                              self.var_lon_centroid, model_index)
        elif interp == model.INTERP_METHOD_GDAL:
            return model.gdal_interpolate_uv_to_regular_grid(u_target_depth, v_target_depth, self.var_lat_centroid,
                                                             self.var_lon_centroid, model_index)

def vertical_interpolation(u, v, h, zeta, siglay, num_nele, num_siglay, target_depth):
    """Vertically interpolate variables to target depth.

    Args:
        u: `numpy.ndarray` containing u values for entire grid.
        v: `numpy.ndarray` containing v values for entire grid.
        zeta: `numpy.ndarray` containing MSL free surface at centroid points in meters.
        h: `numpy.ndarray` containing bathymetry at centroid points.
        siglay: `numpy.ndarray` containing sigma layers (positive up) at centroid points.
        num_nele: number of elements(centroid).
        num_siglay: number of sigma layers.
        target_depth: The water current at a specified target depth below the sea
            surface in meters, default target depth is 4.5 meters, target interpolation
            depth must be greater or equal to 0.
    """
    true_depth = zeta + h

    sigma_depth_layers = numpy.ma.empty(shape=[num_siglay, num_nele])

    for k in range(num_siglay):
        sigma_depth_layers[k, :] = siglay[k,:] * true_depth

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
        u_interp_depth = interpolate.interp1d(sigma_depth_layers[:, nele], u[:, nele],fill_value='extrapolate')
        u_target_depth[nele] = u_interp_depth(interp_depth.data[nele])

        v_interp_depth = interpolate.interp1d(sigma_depth_layers[:, nele], v[:, nele], fill_value='extrapolate')
        v_target_depth[nele] = v_interp_depth(interp_depth.data[nele])

    return u_target_depth, v_target_depth

def node_to_element(model_index, zeta, h, siglay, lon_node, lat_node, lon_centroid, lat_centroid, num_nele, num_siglay):
    """Horizontally interpolate variables at nodes to elements(centroids).

    Args:
        model_index: `ModelIndexFile` instance representing model index file.
        zeta: `numpy.ndarray` containing MSL free surface at nodal points in meters.
        h: `numpy.ndarray` containing bathymetry at nodal points.
        siglay: `3d numpy.ndarray` containing sigma layers (positive up).
        lon_node: `numpy.ndarray` containing nodal longitude.
        lat_node: `numpy.ndarray` containing nodal latitude.
        lon_centroid: `numpy.ndarray` containing centroid(elements) longitude.
        lat_centroid: `numpy.ndarray` containing centroid(elements) latitude.
        num_nele: number of elements(centroids)
        num_siglay: number of sigma layers

    """
    print ("start", datetime.datetime.now())

    siglay_centroid = numpy.ma.empty(shape=[num_siglay, num_nele])

    # Sigma vertical coordinate type
    if model_index.nc_file.modelVerticalCoordinates == "generalized":
        for i in range (num_siglay):
            coords = numpy.column_stack((lon_node, lat_node))
            siglay_centroid[i,:] = interpolate.griddata(coords, siglay[i,:],(lon_centroid, lat_centroid), method='linear')
    else:
        # uniform vertical coordinates
        for k in range(num_nele):
                siglay_centroid[:, k] = siglay[:, 0]

    coords = numpy.column_stack((lon_node, lat_node))
    h_centroid = interpolate.griddata(coords, h, (lon_centroid, lat_centroid), method='linear')
    zeta_centroid = interpolate.griddata(coords, zeta, (lon_centroid, lat_centroid), method='linear')


    return h_centroid, zeta_centroid, siglay_centroid

