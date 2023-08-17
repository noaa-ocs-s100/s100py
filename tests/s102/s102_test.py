import pytest

import os
import pathlib
import logging
import tempfile

import numpy
from osgeo import gdal

from s100py import s100, s102
from s100py.s102 import v2_0
from s100py.s102 import v2_1
from s100py.s102 import v2_2

local_path = pathlib.Path(__file__).parent

# FIXME reinstate the other versions
@pytest.fixture(scope="module", params=[v2_2,]) # v2_1, v2_0])
def s102(request):
    yield request.param

def h5py_string_comp(h5py_val, cmp_str):
    # h5py <3.0 returns a string, >3.0 returns bytes
    return h5py_val in (cmp_str, bytes(cmp_str, "utf-8"))


@pytest.fixture(scope="module")
def s102_file():
    tf = tempfile.TemporaryFile(suffix=".h5")
    f = s102.S102File(tf)
    return f


@pytest.fixture(scope="module")
def bagname():
    fname = str(local_path.joinpath("F00788_SR_8m.bag"))
    yield fname


@pytest.fixture(scope="module")
def tifname():
    # made with:   gdalwarp -of GTiff -srcnodata 1000000.0 -co "COMPRESS=LZW" -dstnodata 9999.0 F00788_SR_8m.bag F00788_SR_8m.tif
    fname = str(local_path.joinpath("F00788_SR_8m.tif"))
    yield fname


@pytest.fixture(scope="module")
def temp_bagname(bagname):
    fname = bagname + ".copy.bag"
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
    # read with lower level h5py access too
    assert h5py_string_comp(s102obj.attrs['horizontalDatumReference'], "EPSG")

    assert s102obj.root.horizontal_datum_value == 32610
    assert s102obj.root.west_bound_longitude > 523816
    assert s102obj.root.west_bound_longitude < 523817
    assert s102obj.root.east_bound_longitude > 525240
    assert s102obj.root.east_bound_longitude < 525241
    b = s102obj.root.bathymetry_coverage.bathymetry_coverage[0]
    assert b.num_points_latitudinal == 179
    assert b.east_bound_longitude == s102obj.root.east_bound_longitude
    group = b.bathymetry_group[0]
    assert group.origin.coordinate[0] == s102obj.root.west_bound_longitude
    assert group.values.depth.shape == (179, 179)
    # depending on if the z is positive up or down the min depth can be 68.4 or 36.1, using min avoids the 1000000 nodata value
    assert numpy.min(group.values.depth) == pytest.approx(-68.44306, 0.0001) or numpy.min(group.values.depth) == pytest.approx(36.18454, 0.0001)


def test_make_from_gdal(bagname, output_path):
    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass
    # the sample data is in NAD83 so does not meet spec - test that it's caught
    pytest.raises(s102.S102Exception, s102.from_gdal, *(bagname, output_path))

    # override the metadata for the datum to WGS84 zone 10N and go from there
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
    new_s102 = s102.from_gdal(bagname, output_path, metadata=metadata)

    check_s102_data(new_s102)


def test_make_from_bag(bagname, output_path):
    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass
    # the sample data is in NAD83 so does not meet spec - test that it's caught
    pytest.raises(s102.S102Exception, s102.from_bag, *(bagname, output_path))

    # override the metadata for the datum to WGS84 zone 10N and go from there
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
    new_s102 = s102.from_bag(bagname, output_path, metadata=metadata)

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
        min_row = rows - 1 - min_row
        empty_row = -1
    if dxx < 0:
        min_col = cols - 1 - min_col
        empty_col = -1
    assert depth_raster.GetNoDataValue() == 9999
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
    new_s102 = s102.from_gdal(gdal_data, temp_bagname, metadata=metadata)
    empty_corner = new_s102.root.bathymetry_coverage.bathymetry_coverage[0].bathymetry_group[0].values.depth[empty_row, empty_col]
    assert new_s102.root.feature_information.bathymetry_coverage_dataset[0].fill_value == empty_corner
    assert orig[-1][0] == depth_raster.GetNoDataValue()
    assert min_data == new_s102.root.bathymetry_coverage.bathymetry_coverage[0].bathymetry_group[0].values.depth[min_row, min_col]


def test_subdivide(output_path):
    # output_path = r"C:\Pydro22_Dev\NOAA\site-packages\Python38\git_repos\s100py\tests\s102\F00788_SR_8m.bag.s102.h5"
    s102_file = s102.S102File(output_path, "r")
    s102_file.subdivide(output_path, 2, 3)


# @TODO iterate the versions using fixtures so all the functions run against all the versions
def test_s102_versions(bagname):
    from s100py.s102 import v2_0, v2_1
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
    bagname = pathlib.Path(bagname)
    f20 = v2_0.api.S102File.from_bag(bagname, bagname.with_suffix(".2_0.h5"), metadata=metadata)
    assert "2.0" in str(f20.root.product_specification)
    f20.close()
    f21 = v2_1.api.S102File.from_bag(bagname, bagname.with_suffix(".2_1.h5"), metadata=metadata)
    assert "2.1" in str(f21.root.product_specification)
    f21.close()
    open20 = s100.open(bagname.with_suffix(".2_0.h5"))
    open21 = s100.open(bagname.with_suffix(".2_1.h5"))
    assert "2.0" in str(open20.root.product_specification)
    assert isinstance(open20, v2_0.api.S102File)
    assert "2.1" in str(open21.root.product_specification)
    assert isinstance(open21, v2_1.api.S102File)
    open20.close()
    open21.close()


def test_s102_version_upgrade(bagname):
    from s100py.s102 import v2_0, v2_1
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
    bagname = pathlib.Path(bagname)
    f20 = v2_0.api.S102File.from_bag(bagname, bagname.with_suffix(".2_0_to_2_1.h5"), metadata=metadata)
    assert "2.0" in str(f20.root.product_specification)
    f20.close()
    f21 = v2_1.api.S102File.upgrade(bagname.with_suffix(".2_0_to_2_1.h5"))
    assert "2.1" in str(f21.root.product_specification)
    f21.close()


# tiffname = r"C:\Data\BlueTopo\RATs\BlueTopo_BC25M26L_20221102.tiff"
# output_path = r"C:\Data\BlueTopo\RATs\BlueTopo_BC25M26L_20221102.h5"
# @TODO reduce the size of the test dataset and add it to the test directory
def test_rat(s102, tifname=r"C:\Data\BlueTopo\RATs\BlueTopo_BC25M26L_20221102b.tiff", output_path=r"C:\Data\BlueTopo\RATs\BlueTopo_BC25M26L_20221102b.h5"):
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
    try:
        os.remove(output_path)
    except (FileNotFoundError, PermissionError):
        pass

    new_outname = str(output_path)+f".{s102.api.EDITION}.h5"
    new_s102_20 = s102.utils.from_gdal(tifname, new_outname, metadata=metadata)
    try:
        os.remove(new_outname)
    except (FileNotFoundError, PermissionError):
        pass

# test_rat(str(local_path.joinpath("F00788_SR_8m.tif")), output_path)
# test_rat(tiffname, output_path)
metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
out_path = r"C:\Data\BlueTopo\RATs\BlueTopo_BC25M26L_20221102b.2_2.h5"
try:
    os.remove(out_path)
except (FileNotFoundError, PermissionError):
    pass
new_s102_20 = v2_2.utils.from_gdal(r"C:\Data\BlueTopo\RATs\BlueTopo_BC25M26L_20221102b.tiff", out_path, metadata=metadata)
