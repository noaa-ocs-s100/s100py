import os
from osgeo import gdal, osr
import tempfile
import pprint
import logging

import numpy

from .s102 import S102File, DEPTH, UNCERTAINTY
from fuse.raw_read.noaa import bag  # this is from NBS (National Bathymetric Source)


def sr_bag_to_s102(input_bag, output_path=""):
    use_gdal = True
    use_hdf5 = False
    use_fuse = True

    if not output_path:
        output_path = input_bag + ".s102.h5"
    if use_gdal:
        bagfile = gdal.Open(input_bag)
        nx, ny = bagfile.RasterXSize, bagfile.RasterYSize
        raster_band = bagfile.GetRasterBand(1)
        uncertainty_band = bagfile.GetRasterBand(2)
        raster_band.XSize, raster_band.YSize
        bagfile.GetProjection()
        bagfile.GetMetadata()
        epsg = osr.SpatialReference(bagfile.GetProjection()).GetAttrValue("AUTHORITY", 1)
        bagfile.GetMetadataDomainList()  # ['', 'IMAGE_STRUCTURE', 'DERIVED_SUBDATASETS', 'xml:BAG']
        meta_dict = bagfile.GetMetadata_Dict("xml:BAG")
        print(meta_dict)

    if use_hdf5:
        import h5py
        import xml
        hdf_bag = h5py.File(input_bag)
        meta = hdf_bag["BAG_root/metadata"]
        meta_xml = b"".join(meta)
        et_xml = xml.etree.ElementTree.fromstring(meta_xml)
        logging.debug(meta_xml.decode("UTF-8").replace("\\n", "\n"))

        def dumpxml(e, indent=0):
            for i, c in enumerate(e):
                print(indent * "\t", i, c)
                dumpxml(c, indent + 1)

        dumpxml(et_xml)
        # print(et_xml.findall('.//CharacterString'))

    if use_fuse:
        fuse_bag = bag.BAGSurvey("")
        meta_gdal, bag_version = fuse_bag._parse_bag_gdal(input_bag)
        meta_bagxml = fuse_bag._parse_bag_xml(input_bag, bag_version)
        print(bag_version)
        pprint.pprint(meta_gdal)
        pprint.pprint(meta_bagxml)

    depth_raster_data = raster_band.ReadAsArray()
    depth_nodata_value = raster_band.GetNoDataValue()
    uncertainty_raster_data = uncertainty_band.ReadAsArray()

    sfile = S102File(output_path, "w", driver=None)
    sfile.create_empty_metadata()  # init the root with a fully filled out empty metadata set
    root = sfile.root
    del root.feature_information.feature_code  # Guessing at the right dataset name to keep (keep featureName)
    # root.feature_information.feature_code_remove()
    root.product_specification = "INT.IHO.S-102.2.0.0"
    root.metadata = os.path.splitext(os.path.basename(output_path))[0] + ".xml"
    root.vertical_datum = meta_bagxml['from_vert_datum']

    bathy_cov_dset = root.feature_information.bathymetry_coverage_dataset
    bathy_depth_info = bathy_cov_dset.append_new_item()  # bathy_cov_dset.append(bathy_cov_dset.metadata_type())
    bathy_depth_info.initialize_properties(True)
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
    bathy_uncertainty_info.initialize_properties(True)
    bathy_uncertainty_info.code = UNCERTAINTY
    bathy_uncertainty_info.name = UNCERTAINTY

    # I'm not sure what to put here, yet
    tracking_cov = root.feature_information.tracking_list_coverage

    track_info = tracking_cov.append_new_item()  # append(tracking_cov.metadata_type())
    track_info.initialize_properties(True)
    track_info.code = "X"
    track_info.name = "X"
    track_info.unit_of_measure = "N/A"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True)
    track_info.code = "Y"
    track_info.name = "Y"
    track_info.unit_of_measure = "N/A"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True)
    track_info.code = "originalValue"
    track_info.name = "Original Value"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True)
    track_info.code = "trackCode"
    track_info.name = "Track Code"
    track_info.unit_of_measure = "N/A"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True)
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

    bathy_01 = root.bathymetry_coverage.bathymetry_coverage.append_new_item()
    bathy_01.initialize_properties(True)

    # @todo @fixme
    print("Need to determine if projected coords or not - assuming UTM right now")
    minx = meta_bagxml['lon_min']
    maxx = meta_bagxml['lon_max']
    miny = meta_bagxml['lat_min']
    maxy = meta_bagxml['lat_max']

    minx = min([meta_bagxml['bounds'][0][0], meta_bagxml['bounds'][1][0]])
    maxx = max([meta_bagxml['bounds'][0][0], meta_bagxml['bounds'][1][0]])
    miny = min([meta_bagxml['bounds'][0][1], meta_bagxml['bounds'][1][1]])
    maxy = max([meta_bagxml['bounds'][0][1], meta_bagxml['bounds'][1][1]])

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
    bathy_01.grid_spacing_longitudinal = meta_bagxml["res"][0]
    bathy_01.grid_spacing_latitudinal = meta_bagxml["res"][1]
    del bathy_01.grid_spacing_vertical
    del bathy_01.grid_origin_vertical
    bathy_01.num_grp = 1

    bathy_group_object = bathy_01.bathymetry_group.append_new_item()
    # bathy_group_object.initialize_properties()  # Not creating everything as I'm not sure if the grid attributes shoul dbe thereTrue)
    print(bathy_group_object.get_standard_properties())
    # @todo  @FIXME
    print("need to determine if this is degrees/metres and set axisNames accordingly")
    # bathy_group_object.axis_names = numpy.array(["longitude", "latitude"])  # row major order means X/longitude first
    # use default dimension =2

    print("fix here -- row/column order?")
    nx, ny = meta_bagxml['shape']

    bathy_01.num_points_latitudinal = ny
    bathy_01.num_points_longitudinal = nx
    del bathy_01.num_points_vertical
    bathy_01.start_sequence = "0,0"
    del bathy_01.vertical_extent_maximum_z
    del bathy_01.vertical_extent_minimum_z

    bathy_group_object.extent_create()
    bathy_group_object.extent.initialize_properties(True)
    bathy_group_object.extent.low.coord_values[0:2] = [0, 0]
    bathy_group_object.extent.high.coord_values[0:2] = [nx, ny]

    depth_max = depth_raster_data[depth_raster_data != depth_nodata_value].max()
    depth_min = depth_raster_data[depth_raster_data != depth_nodata_value].min()
    bathy_group_object.maximum_depth = depth_max
    bathy_group_object.minimum_depth = depth_min

    uncertainty_max = uncertainty_raster_data[uncertainty_raster_data != depth_nodata_value].max()
    uncertainty_min = uncertainty_raster_data[uncertainty_raster_data != depth_nodata_value].min()
    bathy_group_object.minimum_uncertainty = uncertainty_min
    bathy_group_object.maximum_uncertainty = uncertainty_max
    bathy_group_object.dimension = 2

    bathy_group_object.origin_create()
    bathy_group_object.origin.initialize_properties(True)
    bathy_group_object.origin.dimension = 2
    bathy_group_object.origin.coordinate = numpy.array([minx, miny])

    bathy_group_object.values_create()
    grid = bathy_group_object.values
    # @todo -- need to make sure nodata values are correct, especially if converting something other than bag which is supposed to have the same nodata value
    grid.depth = depth_raster_data
    grid.uncertainty = uncertainty_raster_data

    print("hard coding datum for now")
    # @todo @fixme hardcoded stuff....
    root.horizontal_datum_reference = "EPSG"
    root.horizontal_datum_value = 32611
    root.epoch = "G1762"  # this is the 2013-10-16 WGS84 used by CRS
    root.geographic_identifier = "Long Beach, CA"
    root.issue_date = meta_bagxml['date_stamp']  # datetime.date.today().isoformat()

    sfile.write()


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
