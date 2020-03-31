""" Functions to create S102 data from other sources

"""
import logging
import sys
import os

import numpy
from osgeo import gdal, osr

from ..s1xx import s1xx_sequence
from .api import DEPTH, UNCERTAINTY, S102File, S102Exception


# @todo create a friendly name mapping to s102 nested location, then add s102 functions for "to dictionary" and "from dictionary" to api
# that would make these functions easily invertable

def _get_S102File(output_file):
    """ Small helper function to convert the output_file parameter into a S102File, currently accepting file path as string or S102File instance.
    Could propbably accept h5py.File or other things in the future"""
    if isinstance(output_file, S102File):
        data_file = output_file
    else:  # try everything else -- pathlib, str, tempfile, io.BytesIO
        try:
            data_file = S102File(output_file)
        except TypeError as typeerr:
            msg = "Failed to create S102File using {}".format(str(output_file))
            logging.error(msg)
            raise type(typeerr)(msg).with_traceback(sys.exc_info()[2])

    return data_file


def create_s102(output_file, overwrite=True) -> S102File:
    """ Creates or updates an S102File object.
    Default values are set for any data that don't have options or are mandatory to be filled in the S102 spec.

    Parameters
    ----------
    output_file
        Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
    overwrite
        If updating an existing file then set this option to False in order to retain data (not sure this is needed).

    Returns
    -------
    S102File
        The object created or updated by this function.


    """
    data_file = _get_S102File(output_file)
    # @fixme @todo -- I think this will overwrite no matter what, need to look into that
    data_file.create_empty_metadata()  # init the root with a fully filled out empty metadata set
    root = data_file.root
    bathy_cov_dset = root.feature_information.bathymetry_coverage_dataset
    bathy_depth_info = bathy_cov_dset.append_new_item()  # bathy_cov_dset.append(bathy_cov_dset.metadata_type())
    bathy_depth_info.initialize_properties(True, overwrite=overwrite)
    bathy_depth_info.code = DEPTH
    bathy_depth_info.name = DEPTH
    # these are auto-filled by the api
    # bathy_depth_info.unit_of_measure="metres"
    # bathy_depth_info.fill_value=1000000.0
    # bathy_depth_info.datatype=H5T_NATIVE_FLOAT
    # bathy_depth_info.lower = -12000
    # bathy_depth_info.upper = 12000
    # bathy_depth_info.closure = "closedInterval"

    bathy_uncertainty_info = bathy_cov_dset.append_new_item()
    bathy_uncertainty_info.initialize_properties(True, overwrite=overwrite)
    bathy_uncertainty_info.code = UNCERTAINTY
    bathy_uncertainty_info.name = UNCERTAINTY

    # I'm not sure what to put here, yet
    tracking_cov = root.feature_information.tracking_list_coverage

    track_info = tracking_cov.append_new_item()  # append(tracking_cov.metadata_type())
    track_info.initialize_properties(True, overwrite=overwrite)
    track_info.code = "X"
    track_info.name = "X"
    track_info.unit_of_measure = "N/A"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True, overwrite=overwrite)
    track_info.code = "Y"
    track_info.name = "Y"
    track_info.unit_of_measure = "N/A"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True, overwrite=overwrite)
    track_info.code = "originalValue"
    track_info.name = "Original Value"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True, overwrite=overwrite)
    track_info.code = "trackCode"
    track_info.name = "Track Code"
    track_info.unit_of_measure = "N/A"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True, overwrite=overwrite)
    track_info.code = "listSeries"
    track_info.name = "List Series"
    track_info.unit_of_measure = "N/A"

    root.bathymetry_coverage.axis_names = numpy.array(["longitude", "latitude"])  # row major order means X/longitude first
    root.bathymetry_coverage.common_point_rule = 1  # average
    # root.bathymetry_coverage.data_coding_format = 2  # default
    # root.bathymetry_coverage.dimension = 2  # default value
    root.bathymetry_coverage.sequencing_rule_scan_direction = "Longitude, Latitude"
    root.bathymetry_coverage.interpolation_type = 1  # nearest neighbor
    root.bathymetry_coverage.num_instances = 1  # how many Bathycoverages
    root.bathymetry_coverage.sequencing_rule_type = 1  # linear
    del root.bathymetry_coverage.time_uncertainty

    return data_file


def from_arrays(depth_grid: s1xx_sequence, uncert_grid: s1xx_sequence, output_file, nodata_value=None,
                flip_x: bool = False, flip_y: bool = False, overwrite: bool = True) -> S102File:  # num_array, or list of lists accepted
    """  Creates or updates an S102File object based on numpy array/h5py datasets.
    Calls :any:`create_s102` then fills in the HDF5 datasets with the supplied depth_grid and uncert_grid.
    Fills the number of points areas and any other appropriate places in the HDF5 file per the S102 spec.

    Raises an S102Exception if the shapes of the depth and uncertainty (if not None) grids are not equal.

    Parameters
    ----------
    depth_grid
    uncert_grid
        The uncertainty dataset to embed in the object.
        If None then a numpy.zeros array will be created in the appropriate shape to be stored in the file.
    output_file
        Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
    nodata_value
        Value used to denote an empty cell in the grid.  Used in
    flip_x
        boolean if the data should be mirrored on x coordinate (i.e. the original grid is right to left)
        Flips are done here so we can implement a chunked read/write to save memory
    flip_y
        boolean if the data should be mirrored on y coordinate (i.e. the original grid is top to bottom)
        Flips are done here so we can implement a chunked read/write to save memory
    overwrite
        If updating an existing file then set this option to False in order to retain data (not sure this is needed).

    Returns
    -------
    S102File
        The object created or updated by this function.

    """
    # @todo -- Add logic that if the grids are gdal raster bands then read in blocks and use h5py slicing to write in blocks.  Slower but saves resources
    data_file = create_s102(output_file)
    root = data_file.root
    try:
        bathy_01 = root.bathymetry_coverage.bathymetry_coverage[0]
    except IndexError:
        bathy_01 = root.bathymetry_coverage.bathymetry_coverage.append_new_item()
    bathy_01.initialize_properties(recursively_create_children=True, overwrite=overwrite)

    del bathy_01.grid_spacing_vertical
    del bathy_01.grid_origin_vertical
    del bathy_01.number_of_times
    del bathy_01.time_record_interval
    del bathy_01.date_time_of_last_record
    del bathy_01.date_time_of_first_record
    bathy_01.num_grp = 1

    try:
        bathy_group_object = bathy_01.bathymetry_group[0]
    except IndexError:
        bathy_group_object = bathy_01.bathymetry_group.append_new_item()
    # bathy_group_object.initialize_properties()  # Not creating everything as I'm not sure if the grid attributes should be there

    print("fix here -- row/column order?")
    nx, ny = depth_grid.shape
    if uncert_grid is None:
        uncert_grid = numpy.zeros(depth_grid.shape, dtype=numpy.float32)
    if depth_grid.shape != uncert_grid.shape:
        raise S102Exception("Depth and Uncertainty grids have different shapes")

    bathy_01.num_points_latitudinal = ny
    bathy_01.num_points_longitudinal = nx
    bathy_01.start_sequence = "0,0"
    del bathy_01.num_points_vertical
    del bathy_01.vertical_extent_maximum_z
    del bathy_01.vertical_extent_minimum_z

    bathy_group_object.extent_create()
    bathy_group_object.extent.initialize_properties(True, overwrite=overwrite)
    bathy_group_object.extent.low.coord_values[0:2] = [0, 0]
    bathy_group_object.extent.high.coord_values[0:2] = [nx, ny]

    depth_max = depth_grid[depth_grid != nodata_value].max()
    depth_min = depth_grid[depth_grid != nodata_value].min()
    bathy_group_object.maximum_depth = depth_max
    bathy_group_object.minimum_depth = depth_min

    uncertainty_max = uncert_grid[uncert_grid != nodata_value].max()
    uncertainty_min = uncert_grid[uncert_grid != nodata_value].min()
    bathy_group_object.minimum_uncertainty = uncertainty_min
    bathy_group_object.maximum_uncertainty = uncertainty_max
    bathy_group_object.dimension = 2

    bathy_group_object.origin_create()
    bathy_group_object.origin.initialize_properties(True, overwrite=overwrite)
    bathy_group_object.origin.dimension = 2

    bathy_group_object.values_create()
    grid = bathy_group_object.values
    # @todo -- need to make sure nodata values are correct, especially if converting something other than bag which is supposed to have the same nodata value
    # @todo -- Add logic that if the grids are gdal raster bands then read in blocks and use h5py slicing to write in blocks.  Slower but saves resources
    if flip_x:
        depth_grid = numpy.fliplr(depth_grid)
        uncert_grid = numpy.fliplr(uncert_grid)
    if flip_y:
        depth_grid = numpy.flipud(depth_grid)
        uncert_grid = numpy.flipud(uncert_grid)

    grid.depth = depth_grid
    grid.uncertainty = uncert_grid

    return data_file


def from_arrays_with_metadata(depth_grid: s1xx_sequence, uncert_grid: s1xx_sequence, metadata: dict, output_file, nodata_value=None,
                              overwrite: bool = True) -> S102File:  # raw arrays and metadata accepted
    """ Fills or creates an :any:`S102File` from the given arguments.

    Parameters
    ----------
    depth_grid
        a numpy or hdf5 dataset object of the rectangular grid of depths, lower left corner is the first point
    uncert_grid
        a numpy or hdf5 dataset object of the rectangular grid of uncertainties, lower left corner is the first point
    metadata
        a dictionary of metadata describing the grids passed in,
        metadata should have the following key/value pairs:
            - "origin"; tuple of the position (x,y) or (lon, lat) for the lower left corner resolutions are positive
            - "res": tuple of the resolution (cell size) of each grid cell (x, y).
                If a resolution is negative then the grid will be flipped in that dimension and the origin adjusted accordingly.
            - "horizontalDatumReference": See :any:`S102Root` horizontal_datum_reference, ex: "EPSG".
                "EPSG" is the default value.
            - "horizontalDatumValue":  The value for the horizontal data such as the EPSG code ex: 32611
            - "epoch":
            - "geographicIdentifier": Location of the data, ex: "Long Beach, CA, USA".
                An empty string ("") is the default.
            - "issueDate":
            - "metadataFile": File name for the associated discovery metatadata (xml)
    output_file
        Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
    nodata_value
        the "no data" value used in the grids
    overwrite
        if the output_file was an existing S102File then keep any attributes that might have
    Returns
    -------
    S102File

    """
    # @todo - add logic to see if the coordinate system is lower right, if not then need to mirror the arrays or add flags to do that in from_arrays
    res_x, res_y = metadata["res"]
    flip_x = True if res_x < 0 else False
    flip_y = True if res_y < 0 else False

    # @todo @fixme
    print("Need to determine if projected coords or not - assuming UTM right now")
    nx, ny = depth_grid.shape
    corner_x, corner_y = metadata['origin']

    opposite_corner_x = corner_x + res_x * nx
    opposite_corner_y = corner_y + res_y * ny

    minx = min((corner_x, opposite_corner_x))
    maxx = max((corner_x, opposite_corner_x))
    miny = min((corner_y, opposite_corner_y))
    maxy = max((corner_y, opposite_corner_y))

    data_file = from_arrays(depth_grid, uncert_grid, output_file, nodata_value=nodata_value, overwrite=overwrite, flip_x=flip_x, flip_y=flip_y)

    # now add the additional metadata
    root = data_file.root
    bathy_01 = root.bathymetry_coverage.bathymetry_coverage[0]
    bathy_group_object = bathy_01.bathymetry_group[0]

    root.east_bound_longitude = minx
    root.west_bound_longitude = maxx
    root.south_bound_latitude = miny
    root.north_bound_latitude = maxy
    bathy_01.east_bound_longitude = minx
    bathy_01.west_bound_longitude = maxx
    bathy_01.south_bound_latitude = miny
    bathy_01.north_bound_latitude = maxy
    bathy_01.grid_origin_latitude = miny

    bathy_01.grid_origin_longitude = minx
    bathy_01.grid_origin_latitude = miny
    bathy_01.grid_spacing_longitudinal = abs(res_x)  # we adjust for negative resolution in the from_arrays
    bathy_01.grid_spacing_latitudinal = abs(res_y)

    # @todo  @FIXME
    print("need to determine if this is degrees/metres and set axisNames accordingly")
    # bathy_group_object.axis_names = numpy.array(["longitude", "latitude"])  # row major order means X/longitude first

    bathy_group_object.origin.coordinate = numpy.array([minx, miny])

    # these names are taken from the S100/S102 attribute names
    # but are hard coded here to allow the S102 spec to change but not affect any tools built on these utility functions
    if "horizontalDatumReference" in metadata or overwrite:
        root.horizontal_datum_reference = metadata.get("horizontalDatumReference", "EPSG")
    if "horizontalDatumValue" in metadata or overwrite:
        root.horizontal_datum_value = metadata.get("horizontalDatumValue", 0)
    if "epoch" in metadata or overwrite:
        root.epoch = metadata.get("epoch", "")  # e.g. "G1762"  this is the 2013-10-16 WGS84 used by CRS
    if "geographicIdentifier" in metadata or overwrite:
        root.geographic_identifier = metadata.get("geographicIdentifier", "")
    if "issueDate" in metadata or overwrite:
        root.issue_date = metadata.get('issueDate', "")  # datetime.date.today().isoformat()
    if "metadataFile" in metadata or overwrite:
        root.metadata = metadata.get('metadataFile', "")  # datetime.date.today().isoformat()

    data_file.write()
    return data_file


def from_gdal(input_raster, output_file, metadata: dict = {}) -> S102File:  # gdal instance or filename accepted
    """ Fills or creates an :any:`S102File` from the given arguments.

    Parameters
    ----------
    input_raster
        Either a path to a raster file that GDAL can open or a gdal.Dataset object.
    output_file
        Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
    metadata
        A dictionary of metadata describing the grids passed in.
        All the metadata used in :any:`from_from_arrays_with_metadata` can be specified and
        would override the values that would have been populated based on the GDAL data.

        horizontalDatumReference, horizontalDatumValue, origin, res will be determined from GDAL if not otherwise specified.

    Returns
    -------
    S102File

    """
    if isinstance(input_raster, gdal.Dataset):
        dataset = input_raster
    else:
        dataset = gdal.Open(input_raster)

    # @todo @fixme -- transform the coordinate system to a WGS84.  Strictly this may not end up being square, so how do we handle
    #  transform = osr.CoordinateTransformation( src_srs, tgt_srs)
    print("Glen, how do you want to transform coordinates?")
    if "horizontalDatumReference" not in metadata or "horizontalDatumValue" not in metadata:
        metadata["horizontalDatumReference"] = "EPSG"
        epsg = osr.SpatialReference(dataset.GetProjection()).GetAttrValue("AUTHORITY", 1)
        metadata["horizontalDatumValue"] = epsg

    raster_band = dataset.GetRasterBand(1)
    depth_nodata_value = raster_band.GetNoDataValue()
    uncertainty_band = dataset.GetRasterBand(2)

    nx, ny = raster_band.XSize, raster_band.YSize
    ulx, dxx, dxy, uly, dyx, dyy = dataset.GetGeoTransform()
    if dxy != 0.0 or dyx != 0.0:
        raise S102Exception("raster is not north up but is rotated, this is not handled at this time")

    # not used now as from_arrays moved to origin
    lrx = ulx + (nx * dxx) + (ny * dyx)
    lry = uly + (nx * dxy) + (ny * dyy)

    if "origin" not in metadata:  # gdal convention is upper left -- we need to handle that where we write from arrays
        metadata["origin"] = [ulx, uly]
    if "res" not in metadata:
        metadata["res"] = [dxx, dyy]
    s102_data_file = from_arrays_with_metadata(raster_band.ReadAsArray(), uncertainty_band.ReadAsArray(), metadata, output_file,
                                               nodata_value=depth_nodata_value)

    return s102_data_file


def from_bag(bagfile, output_file) -> S102File:
    """
    Parameters
    ----------
    bagfile
    output_file
        Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
    Returns
    -------

    """
    metadata = {"horizontalDatumValue": 32611,
                "horizontalDatumReference": "EPSG",
                "geographicIdentifier": "Long Beach, CA",
                "epoch": "G1762",
                }
    try:
        pth = os.path.splitext(os.path.basename(output_file))[0] + ".xml"
    except Exception:
        pth = "None"
    metadata["metadataFile"] = pth

    s102_data_file = from_gdal(bagfile, output_file, metadata=metadata)

    raise S102Exception("This function is just for show right now!")


# @todo remove this function once we extract any logic for reading epsg etc
def read_bag(input_bag):
    """
    Extract the required information from the BAG using GDAL and Fuse.  Return
    the metadata in a dictionary and the arrays.
    """
    # check to see if VR
    # if VR, resample at specific resolution
    # if not return the stuff.
    metadata = dict()
    # use gdal to get the rasters and horizontal georeferencing
    bagfile = gdal.Open(input_bag)
    raster_band = bagfile.GetRasterBand(1)
    uncertainty_band = bagfile.GetRasterBand(2)
    raster_band.XSize, raster_band.YSize
    bagfile.GetProjection()
    bagfile.GetMetadata()
    epsg = osr.SpatialReference(bagfile.GetProjection()).GetAttrValue("AUTHORITY", 1)
    if epsg is None:
        raise ValueError('No EPSG code discernible from BAG read with GDAL')
    else:
        if epsg in get_valid_epsg():
            metadata['epsg'] = epsg
        else:
            raise ValueError(f'BAG EPSG code is not within those allowd by S-102 spec: {epsg}')
    bagfile.GetMetadataDomainList()  # ['', 'IMAGE_STRUCTURE', 'DERIVED_SUBDATASETS', 'xml:BAG']
    meta_dict = bagfile.GetMetadata_Dict("xml:BAG")
    print(meta_dict)

    # use fuse to read the XML and get the
    #   Date
    #   resolution (x and y)
    #   bounds (x and y, but lat lon or utm?)
    #   vertical datum
    fuse_bag = bag.BAGSurvey("")
    meta_gdal, bag_version = fuse_bag._parse_bag_gdal(input_bag)
    meta_bagxml = fuse_bag._parse_bag_xml(input_bag, bag_version)
    metadata = {**metadata, **meta_bagxml}

    depth_raster_data = raster_band.ReadAsArray()
    metadata['nodata'] = raster_band.GetNoDataValue()
    uncertainty_raster_data = uncertainty_band.ReadAsArray()

    return depth_raster_data, uncertainty_raster_data, metadata
