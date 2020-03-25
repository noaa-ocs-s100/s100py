import os
from osgeo import gdal, osr
import tempfile
import logging

import numpy

from .s102 import make_s102
from fuse.raw_read.noaa import bag  # this is from NBS (National Bathymetric Source)

def get_valid_epsg() -> list:
    """
    Create and return the list of valid EPSG codes for S-102 version 2.0.
    """
    valid_epsg = [4326, 5041, 5042]
    valid_epsg += list(numpy.arange(32601, 32660 + 1))
    valid_epsg += list(numpy.arange(32701, 32760 + 1))
    return valid_epsg

def bag_to_s102(input_bag, region, output_path=""):
    """
    Read a BAG and write it to S-102 v2.
    """
    elev, uncert, metadata = read_bag(input_bag)
    if not output_path:
        output_path = input_bag + ".s102.h5"
    make_s102(output_path, elev, uncert, metadata)

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

def NAVO_convert_bag(bag_path, output_path, path_to_convertor=".\\BAG_to_S102.exe", buffer=False):
    import subprocess
    cmd = '"' + path_to_convertor + '" "' + bag_path + '" "' + output_path + '"'
    print(cmd)
    if buffer:
        std_out = tempfile.TemporaryFile()  # deleted on exit from function
        std_err = tempfile.TemporaryFile()
    else:
        std_out = None
        std_err = None
    p = subprocess.Popen(cmd, stdout=std_out, stderr=std_err)
    p.wait()
    if buffer:
        std_out.seek(0)
        std_err.seek(0)
        out = std_out.read()
        err = std_err.read()
        print(out)
        print(err)
