import pytest

import os
import logging
import tempfile

import numpy
from osgeo import gdal

from s100py import s100, s102

@pytest.fixture(scope="module")
def s102_file():
    tf = tempfile.TemporaryFile(suffix=".h5")
    f = s102.S102File(tf)
    return f

@pytest.fixture(scope="module")
def bagname():
    fname = "F00788_SR_8m.bag"
    yield fname

@pytest.fixture(scope="module")
def tifname():
    # made with:   gdalwarp -of GTiff -srcnodata 1000000.0 -co "COMPRESS=LZW" -dstnodata 9999.0 F00788_SR_8m.bag F00788_SR_8m.tif
    fname = "F00788_SR_8m.tif"
    yield fname

@pytest.fixture(scope="module")
def temp_bagname(bagname):
    fname = bagname+".copy.bag"
    yield fname
    try:
        os.remove(fname)
    except (FileNotFoundError, PermissionError):
        pass

@pytest.fixture(scope="module")
def output_path(bagname):
    out_path = bagname + ".s102_output_test.h5"
    yield out_path
    try:
        os.remove(out_path)
    except (FileNotFoundError, PermissionError):
        pass

@pytest.fixture(scope="module")
def copy_path(bagname):
    out_path = bagname + ".s102_output_test.copy.h5"
    yield out_path
    try:
        os.remove(out_path)
    except (FileNotFoundError, PermissionError):
        pass

def check_s102_data(s102obj):
    assert s102obj.root
    assert s102obj.root.horizontal_datum_reference == "EPSG"
    assert s102obj.root.horizontal_datum_value == 26910
    assert s102obj.root.east_bound_longitude > 523812
    assert s102obj.root.east_bound_longitude < 523813
    b = s102obj.root.bathymetry_coverage.bathymetry_coverage[0]
    assert b.num_points_latitudinal == 179
    assert b.east_bound_longitude == s102obj.root.east_bound_longitude
    group = b.bathymetry_group[0]
    assert group.origin.coordinate[0] == s102obj.root.east_bound_longitude
    assert group.values.depth.shape == (179, 179)
    assert numpy.min(group.values.depth) == pytest.approx(-68.44306, 0.0001)

def test_make_from_gdal(bagname, output_path):
    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass
    new_s102 = s102.from_gdal(bagname, output_path)
    check_s102_data(new_s102)

def test_read_s102(output_path):
    s102_read_test = s102.S102File(output_path, "r")
    check_s102_data(s102_read_test)
    s102_read_test.close()

def test_copy_s102(output_path, copy_path):
    s102_read_test = s102.S102File(output_path, "r")
    check_s102_data(s102_read_test)
    try:
        os.remove(copy_path)
    except FileNotFoundError:
        pass
    s102_copy_root_test = s102.S102File(copy_path, "w")
    s102_copy_root_test.root = s102_read_test.root
    s102_copy_root_test.write()
    check_s102_data(s102_copy_root_test)
    s102_read_test.close()
    s102_copy_root_test.close()


def test_tif_conversion(tifname, temp_bagname):
    # test that nodata is changed to 1000000 and that passing a gdal instance works instead of filename
    try:
        os.remove(temp_bagname)
    except FileNotFoundError:
        pass
    gdal_data = gdal.Open(tifname)
    depth_raster = gdal_data.GetRasterBand(1)
    orig = depth_raster.ReadAsArray()
    rows, cols = orig.shape
    # the corner at zero zero in the test data is empty
    empty_row = 0
    empty_col = 0
    min_data = numpy.min(orig)
    min_row, min_col = numpy.where(orig == min_data)
    ulx, dxx, dxy, uly, dyx, dyy = gdal_data.GetGeoTransform()
    if dyy < 0:
        min_row = rows -1 - min_row
        empty_row = -1
    if dxx < 0:
        min_col = cols - 1 - min_col
        empty_col = -1
    assert depth_raster.GetNoDataValue() == 9999
    new_s102 = s102.from_gdal(gdal_data, temp_bagname)
    empty_corner = new_s102.root.bathymetry_coverage.bathymetry_coverage[0].bathymetry_group[0].values.depth[empty_row, empty_col]
    assert new_s102.root.feature_information.bathymetry_coverage_dataset[0].fill_value == empty_corner
    assert orig[-1][0] == depth_raster.GetNoDataValue()
    assert min_data == new_s102.root.bathymetry_coverage.bathymetry_coverage[0].bathymetry_group[0].values.depth[min_row, min_col]


