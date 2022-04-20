""" Functions to create S102 data from other sources

If the utils module is run as __main__ then it will run  :func:`from_bag` or :func:`from_gdal` on the filename given in the command line arguments.

"""
# this also works but doesn't put parenthesis in the html   :any:`from_bag`
# fully qualified would also work   :func:`s100py.s102.utils.from_gdal`

import os
import sys

if getattr(sys, 'frozen', False):
    # in a frozen exe os.pyc is at the root level
    proj_db_path = os.path.join(os.path.dirname(os.__file__), "library\\share\\proj")
    # print("frozen exe", proj_db_path)
    os.environ["PROJ_LIB"] = proj_db_path
    # print(os.listdir(proj_db_path))  # os.listdir(os.path.dirname(os.__file__)))
else:
    pass
    # print("running as script")

import logging
import warnings
import argparse
from xml.etree import ElementTree as et
import tkinter as tk
from tkinter import filedialog, messagebox
import numpy
from osgeo import gdal, osr
import h5py
import re

try:
    from matplotlib import pyplot
    from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
    from matplotlib.colorbar import ColorbarBase
except:
    if not getattr(sys, 'frozen', False):  # we expect the frozen exe to not have matplotlib
        print("matplotlib.pyplot failed to import, plotting will not work")

from s100py.s1xx import s1xx_sequence, S1XXFile
from s100py.s102.api import DEPTH, UNCERTAINTY, S102File, S102Exception

__all__ = ['plot_depth_using_h5py', 'create_s102', 'from_arrays', 'from_arrays_with_metadata',
           'from_gdal', 'from_bag', 'get_valid_epsg']

gco = "{http://www.isotc211.org/2005/gco}"

# @todo create a friendly name mapping to s102 nested location, then add s102 functions for "to dictionary" and "from dictionary" to api
#   that would make these functions easily invertable

r"""
from s100py.s102 import utils
navo_name = r"G:\Data\S102 Data\LA_LB_Area_GEO_reprojected.bag.navo_%d.h5"
utils.plot_depth_using_h5py(navo_name)
bag_names = [r"G:\Data\S102 Data\NBS_US5NYCAH_20200430.bag", r"G:\Data\S102 Data\NBS_US5NYCBH_20200429.bag", r"G:\Data\S102 Data\LA_LB_Area_GEO_reprojected.bag"]
fout = utils.from_bag(bag_name, bag_name+r".noaa.h5")
for bag_name in bag_names:
    fout = utils.from_bag(bag_name, bag_name+r".noaa.h5")
    fout.close()
    utils.plot_depth_using_h5py(bag_name + r".noaa.h5")
"""


def plot_depth_using_h5py(filename, enc_color=False):
    # filename = r"G:\Data\S102 Data\GlenS102Test\102USA15NYCAH200430.H5"
    # h5py.File(r"G:\Data\S102 Data\LA_LB_Area_GEO_reprojected.bag_%d.h5", mode="r", driver="family", memb_size=681574400)
    try:
        f = h5py.File(filename, mode="r", driver="family")
    except OSError as e:
        # see if the member size isn't right and then fall back to standard opening
        # OSError: Unable to open file (Family member size should be 681574400.  But the size from file access property is 2147483647)
        try:
            error_string = str(e)
            m = re.search(r'Family member size should be\s+(\d+)', error_string)
            if m:
                sz = int(m[1])
            f = h5py.File(filename, mode="r", driver="family", memb_size=int(sz))
        except:
            f = h5py.File(filename, mode="r")
    fill_val = 1000000
    try:
        d = f["BathymetryCoverage/BathymetryCoverage.01/Group.001/values"]['depth']
    except KeyError:
        try:
            d = f["BathymetryCoverage/BathymetryCoverage.001/Group.001/values"]['depth']
        except KeyError:
            d = f["SurfaceCurrent/SurfaceCurrent.01/Group_001/values"]['surfaceCurrentSpeed']
            fill_val = -9999
    d[d==fill_val] = numpy.nan

    # ud = numpy.flipud(d)
    if enc_color:
        colors = numpy.array([(255, 255, 255), (201, 237, 252), (167, 217, 251), (130, 202, 255), (97, 180, 255)]) / 255
        # bins = [-1000, -9.1, -5.4, -3.6, -1.8, ]
        bins = [-1000, -30, -22, -15, -12]
        norm = BoundaryNorm(bins, colors.shape[0])
        cmap = ListedColormap(colors)
        pyplot.imshow(d, interpolation='nearest', cmap=cmap, norm=norm)
    else:
        pyplot.imshow(d)
    pyplot.gca().invert_yaxis()
    pyplot.show()
    # cm = LinearSegmentedColormap.from_list("enc", colors, N=4)
    # im = pyplot.imshow(ud, cmap=cm, interpolation='nearest')
    # pyplot.colorbar(im)


def _get_S102File(output_file):
    """ Small helper function to convert the output_file parameter into a S102File, currently accepting file path as string or S102File instance.
    Could propbably accept h5py.File or other things in the future"""
    if isinstance(output_file, S1XXFile):
        data_file = output_file
    else:  # try everything else -- pathlib, str, tempfile, io.BytesIO
        try:
            data_file = S102File(output_file, "w")
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

    root.bathymetry_coverage.axis_names = numpy.array(["Longitude", "Latitude"])  # row major order means X/longitude first
    root.bathymetry_coverage.sequencing_rule_scan_direction = "Longitude, Latitude"
    root.bathymetry_coverage.common_point_rule = 1  # average
    # root.bathymetry_coverage.data_coding_format = 2  # default
    # root.bathymetry_coverage.dimension = 2  # default value
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
        Value used to denote an empty cell in the grid.  Used in finding the min/max and then converted to the S102 fillValue.
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

    # @todo @fixme fix here -- row/column order?
    rows, cols = depth_grid.shape
    if uncert_grid is None:
        uncert_grid = numpy.full(depth_grid.shape, nodata_value, dtype=numpy.float32)
    if depth_grid.shape != uncert_grid.shape:
        raise S102Exception("Depth and Uncertainty grids have different shapes")

    bathy_01.num_points_latitudinal = rows
    bathy_01.num_points_longitudinal = cols
    bathy_01.start_sequence = "0,0"
    del bathy_01.num_points_vertical
    del bathy_01.vertical_extent_maximum_z
    del bathy_01.vertical_extent_minimum_z

    bathy_group_object.extent_create()
    bathy_group_object.extent.initialize_properties(True, overwrite=overwrite)
    bathy_group_object.extent.low.coord_values[0:2] = [0, 0]
    bathy_group_object.extent.high.coord_values[0:2] = [rows, cols]

    depth_max = depth_grid[depth_grid != nodata_value].max()
    depth_min = depth_grid[depth_grid != nodata_value].min()
    bathy_group_object.maximum_depth = depth_max
    bathy_group_object.minimum_depth = depth_min

    try:
        uncertainty_max = uncert_grid[uncert_grid != nodata_value].max()
        uncertainty_min = uncert_grid[uncert_grid != nodata_value].min()
    except ValueError:  # an empty uncertainty array (all values == nodata) will cause this
        uncertainty_max = uncertainty_min = nodata_value

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
    if nodata_value != root.feature_information.bathymetry_coverage_dataset[0].fill_value:
        depth_grid = numpy.copy(depth_grid)
        depth_grid[depth_grid == nodata_value] = root.feature_information.bathymetry_coverage_dataset[0].fill_value
        uncert_grid = numpy.copy(uncert_grid)
        uncert_grid[uncert_grid == nodata_value] = root.feature_information.bathymetry_coverage_dataset[1].fill_value

    grid.depth = depth_grid
    grid.uncertainty = uncert_grid

    return data_file


def from_arrays_with_metadata(depth_grid: s1xx_sequence, uncert_grid: s1xx_sequence, metadata: dict, output_file, nodata_value=None,
                              overwrite: bool = True) -> S102File:  # raw arrays and metadata accepted
    """ Fills or creates an :any:`S102File` from the given arguments.

    Parameters
    ----------
    depth_grid
        a numpy or hdf5 dataset object of the rectangular grid of depths
    uncert_grid
        a numpy or hdf5 dataset object of the rectangular grid of uncertainties, lower left corner is the first point
    metadata
        a dictionary of metadata describing the grids passed in,
        metadata should have the following key/value pairs:
            - "origin": tuple of the position (x,y) or (lon, lat) for the reference corner node.
                Other corners are calulated from this corner using the resolution and size of the data array.
            - "res": tuple of the resolution (cell size) of each grid cell (x, y).
                Lower left corner is the first point of both resolutions are positive.
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

    rows, cols = depth_grid.shape
    corner_x, corner_y = metadata['origin']

    # S-102 is node based, so distance to far corner is res * (n -1)
    opposite_corner_x = corner_x + res_x * (cols - 1)
    opposite_corner_y = corner_y + res_y * (rows - 1)

    minx = min((corner_x, opposite_corner_x))
    maxx = max((corner_x, opposite_corner_x))
    miny = min((corner_y, opposite_corner_y))
    maxy = max((corner_y, opposite_corner_y))

    data_file = from_arrays(depth_grid, uncert_grid, output_file, nodata_value=nodata_value, overwrite=overwrite, flip_x=flip_x, flip_y=flip_y)

    # now add the additional metadata
    root = data_file.root
    bathy_01 = root.bathymetry_coverage.bathymetry_coverage[0]
    bathy_group_object = bathy_01.bathymetry_group[0]

    root.east_bound_longitude = maxx
    root.west_bound_longitude = minx
    root.south_bound_latitude = miny
    root.north_bound_latitude = maxy
    bathy_01.east_bound_longitude = maxx
    bathy_01.west_bound_longitude = minx
    bathy_01.south_bound_latitude = miny
    bathy_01.north_bound_latitude = maxy
    bathy_01.grid_origin_latitude = miny

    bathy_01.grid_origin_longitude = minx
    bathy_01.grid_origin_latitude = miny
    bathy_01.grid_spacing_longitudinal = abs(res_x)  # we adjust for negative resolution in the from_arrays
    bathy_01.grid_spacing_latitudinal = abs(res_y)


    bathy_group_object.origin.coordinate = numpy.array([minx, miny])

    # these names are taken from the S100/S102 attribute names
    # but are hard coded here to allow the S102 spec to change but not affect any tools built on these utility functions
    if "horizontalDatumReference" in metadata or overwrite:
        root.horizontal_datum_reference = metadata.get("horizontalDatumReference", "EPSG")
    if "horizontalDatumValue" in metadata or overwrite:
        source_epsg = int(metadata.get("horizontalDatumValue", 0))
        if source_epsg in get_valid_epsg():
            root.horizontal_datum_value = source_epsg
        else:
            raise ValueError(f'The provided EPSG code {source_epsg} is not within the S102 specified values.')
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(root.horizontal_datum_value)
    if srs.IsProjected():
        axes = ["Easting", "Northing"]  # ["Northing", "Easting"]  # row major instead of
    else:
        axes = ["Longitude", "Latitude"]  # ["Latitude", "Longitude"]  # row major instead of

    bathy_group_object.axis_names = numpy.array(axes)  # row major order means X/longitude first
    root.bathymetry_coverage.axis_names = numpy.array(axes)  # row major order means X/longitude first
    root.bathymetry_coverage.sequencing_rule_scan_direction = ", ".join(axes)

    if "epoch" in metadata or overwrite:
        root.epoch = metadata.get("epoch", "")  # e.g. "G1762"  this is the 2013-10-16 WGS84 used by CRS
    if "geographicIdentifier" in metadata or overwrite:
        root.geographic_identifier = metadata.get("geographicIdentifier", "")
    if "issueDate" in metadata or overwrite:
        root.issue_date = metadata.get('issueDate', "")  # datetime.date.today().isoformat()
    if "metadataFile" in metadata or overwrite:
        root.metadata = metadata.get('metadataFile', "")  # datetime.date.today().isoformat()

    data_file.write()
    data_file.flush()
    
    return data_file


def from_gdal(input_raster, output_file, metadata: dict = None) -> S102File:  # gdal instance or filename accepted
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
    if metadata is None:
        metadata = {}
    else:
        metadata = metadata.copy()

    if isinstance(input_raster, gdal.Dataset):
        dataset = input_raster
    else:
        dataset = gdal.Open(input_raster)

    # @todo @fixme -- transform the coordinate system to a WGS84.  Strictly this may not end up being square, so how do we handle
    #  transform = osr.CoordinateTransformation( src_srs, tgt_srs)
    # Until we have a working datum engine this module should not do datum transformations - GR 20200402
    if "horizontalDatumReference" not in metadata or "horizontalDatumValue" not in metadata:
        metadata["horizontalDatumReference"] = "EPSG"
        sr = osr.SpatialReference(dataset.GetProjection())
        epsg = sr.GetAuthorityCode(None)
        # FIXME: this is likely incorrect. We probably don't want to get the code of the geographic CRS when the CRS is projected
        if epsg is None and sr.IsProjected():
            sr = sr.CloneGeogCS()
        if epsg:
            metadata["horizontalDatumValue"] = int(epsg)
        else:
            if sr.GetAttrValue("GEOGCS") == 'WGS 84':
                metadata["horizontalDatumValue"] = 4326
            #elif sr.GetAttrValue("GEOGCS") == 'North_American_Datum_1983':
            #    metadata["horizontalDatumValue"] = 4269
            else:
                raise S102Exception("Projection not understood, was searching for an EPSG code and found " + osr.SpatialReference(dataset.GetProjection()).ExportToWkt())

    if "epoch" not in metadata:
        # @todo We should be able to pull this from the WKT
        pass

    raster_band = dataset.GetRasterBand(1)
    depth_nodata_value = raster_band.GetNoDataValue()
    uncertainty_band = dataset.GetRasterBand(2)

    ulx, dxx, dxy, uly, dyx, dyy = dataset.GetGeoTransform()
    if dxy != 0.0 or dyx != 0.0:
        raise S102Exception("raster is not north up but is rotated, this is not handled at this time")

    if "origin" not in metadata:
        # shift the gdal geotransform corner point to reference the node (pixel is center) rather than cell (pixel is area)
        metadata["origin"] = [ulx + dxx/2, uly + dyy/2]
    if "res" not in metadata:
        metadata["res"] = [dxx, dyy]
    s102_data_file = from_arrays_with_metadata(raster_band.ReadAsArray(), uncertainty_band.ReadAsArray(), metadata, output_file,
                                               nodata_value=depth_nodata_value)

    return s102_data_file


def from_bag(bagfile, output_file, metadata: dict = None) -> S102File:
    """
    Parameters
    ----------
    bagfile
        Either a path to a raster file that GDAL can open or a gdal.Dataset object.
    output_file
        Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
    metadata
        Supports the metadata options in :any:`from_from_arrays_with_metadata`.
        In addition, 'resample_resolution' can supplied to use a particular resolution using gdal "MODE=RESAMPLED_GRID"
    Returns
    -------

    """
    # @todo update method docstring for possible metadata fields
    if metadata is None:
        metadata = {}
    else:
        metadata = metadata.copy()

    if isinstance(bagfile, gdal.Dataset):
        bag = bagfile
    else:
        bag = gdal.Open(bagfile)
        
    # check for and resample variable resolution BAG if able
    gdal_metadata = bag.GetMetadata()
    if 'HAS_SUPERGRIDS' in gdal_metadata and gdal_metadata['HAS_SUPERGRIDS'] == 'TRUE':
        bag_filename = bag.GetFileList()[0]
        if "resample_resolution" in metadata:
            res = metadata["resample_resolution"]
            bag = None
            bag = gdal.OpenEx(bag_filename, open_options=['MODE=RESAMPLED_GRID', f'RESX={res}',f'RESY={res}'])
        else:
            warnings.warn(f'No resampling resolution provided for variable resolution bag {bag_filename}.  Using overview resolution.', category=RuntimeWarning)

    # populate the issueDate if possible from a simple string search
    if 'issueDate' not in metadata:
        xml_str = bag.GetMetadata('xml:BAG')[0]
        root = et.fromstring(xml_str)
        elem = root.find(".//" + gco + "Date")
        if elem is not None and elem.text:
            metadata['issueDate'] = elem.text

    s102_data_file = from_gdal(bag, output_file, metadata=metadata)
    
    return s102_data_file


def get_valid_epsg() -> list:
    """
    Create and return the list of valid EPSG codes for S-102 version 2.0.
    """
    valid_epsg = [4326, 5041, 5042]
    valid_epsg += list(numpy.arange(32601, 32660 + 1))
    valid_epsg += list(numpy.arange(32701, 32760 + 1))
    return valid_epsg

def browse_files(question):
    # using tkinter since it is built in to python and smaller to distribute than PySide2 or wxPython in an executable
    root = tk.Tk()
    root.withdraw()
    # root.filename = tkFileDialog.askopenfilename(initialdir="/", title="Select file",
    #                                              filetypes=(("jpeg files", "*.jpg"), ("all files", "*.*")))
    file_path = filedialog.askopenfilename(title=question)
    return file_path

def bool_question(question, title="", icon="warning"):
    root = tk.Tk()
    root.withdraw()
    result = messagebox.askquestion(title, question, icon=icon)
    return result == 'yes'


def make_parser():
    parser = argparse.ArgumentParser(description='Convert a georeferenced file to S102')
    parser.add_argument("-?", "--show_help", action="store_true", help="show this help message and exit")
    parser.add_argument("-i", "--input_filename", help="full path to the file to be processed")
    parser.add_argument("-o", "--output_filename", help="output filename, default is same name as input with .h5 appended")
    parser.add_argument("-r", "--res", help="Resolution.  If the input file is a BAG then use attempt to use the given resolution" )
    return parser


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    if args.show_help:
        parser.print_help()
        sys.exit()

    if not args.input_filename:
        path = browse_files('Choose file to convert to S102')
        if path:
            args.input_filename = path

    if args.input_filename:
        # get the output file path from the command line or the user if it wasn't specified
        if args.output_filename:
            output_name = args.output_filename
        else:
            output_name = args.input_filename + ".h5"
            if os.path.exists(output_name):
                if bool_question(f"{output_name} already exists, overwrite?", "Overwrite File"):
                    os.remove(output_name)

        # check if the data is a bag and should be sent to the from_bag function or just raster and send to from_gdal
        ds = gdal.Open(args.input_filename)
        drv = ds.GetDriver()
        if drv.GetDescription() == "BAG":
            from_bag(ds, output_name)
        else:
            from_gdal(ds, output_name)

