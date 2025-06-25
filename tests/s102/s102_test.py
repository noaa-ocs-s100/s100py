import pytest
import logging

import os
import pathlib
import logging
import tempfile

import numpy
from osgeo import gdal

from s100py import s100
from s100py.s102 import v2_0
from s100py.s102 import v2_1
from s100py.s102 import v2_2
from s100py.s102 import v3_0
latest = v3_0

local_path = pathlib.Path(__file__).parent


def remove_file(pth, quiet=True):
    """Remove a file if it exists"""
    try:
        os.remove(pth)
    except (FileNotFoundError, PermissionError):  #
        if not quiet:
            logging.warning(f"{pth} not found")

@pytest.fixture(scope="module", params=[v3_0, v2_2, v2_1, v2_0])
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
def tif_with_rat_name():
    # made with:   gdalwarp -of GTiff -srcnodata 1000000.0 -co "COMPRESS=LZW" -dstnodata 9999.0 F00788_SR_8m.bag F00788_SR_8m.tif
    fname = str(local_path.joinpath("BlueTopo_BC25M26L_20221102b.tiff"))
    yield fname


@pytest.fixture(scope="module")
def temp_bagname(bagname):
    fname = bagname + ".copy.bag"
    yield fname
    remove_file(fname)

@pytest.fixture(scope="module")
def output_path(bagname):
    out_path = bagname + ".s102_output_test.h5"
    yield out_path
    remove_file(out_path)

@pytest.fixture(scope="module")
def output_with_rat_path(tif_with_rat_name):
    out_path = tif_with_rat_name + ".s102_output_rat_test.h5"
    yield out_path
    remove_file(out_path)

@pytest.fixture(scope="module")
def copy_path(bagname):
    out_path = bagname + ".s102_output_test.copy.h5"
    yield out_path
    remove_file(out_path)

def check_s102_data(s102, s102obj):
    assert s102obj.root
    if s102.api.EDITION <= 2.1:
        assert s102obj.root.horizontal_datum_reference == "EPSG"
        # read with lower level h5py access too
        assert h5py_string_comp(s102obj.attrs['horizontalDatumReference'], "EPSG")
        assert s102obj.root.horizontal_datum_value == 32610
    else:  # v2.2 and later
        assert s102obj.root.horizontal_crs == 32610
    assert s102obj.epsg == 32610

    if s102.api.EDITION == 2.0:
        assert s102obj.root.west_bound_longitude == pytest.approx(523816, 1)
        assert s102obj.root.east_bound_longitude == pytest.approx(525240, 1)
        assert s102obj.root.north_bound_latitude == pytest.approx(5332689, 1)
        assert s102obj.root.south_bound_latitude == pytest.approx(5334113, 1)
    else:
        assert s102obj.root.west_bound_longitude == pytest.approx(-122.679, .001)
        assert s102obj.root.east_bound_longitude == pytest.approx(-122.660, .001)
        assert s102obj.root.north_bound_latitude == pytest.approx(48.159, .001)
        assert s102obj.root.south_bound_latitude == pytest.approx(48.147, .001)

    b = s102obj.root.bathymetry_coverage.bathymetry_coverage[0]

    assert b.west_bound_longitude == pytest.approx(523816, 1)
    assert b.east_bound_longitude == pytest.approx(525240, 1)
    assert b.north_bound_latitude == pytest.approx(5332689, 1)
    assert b.south_bound_latitude == pytest.approx(5334113, 1)
    assert s102obj.grid_origin_latitude == pytest.approx(5332689, 1)
    assert s102obj.grid_origin_longitude == pytest.approx(523816, 1)
    assert s102obj.grid_spacing_latitudinal == 8.0
    assert s102obj.grid_spacing_longitudinal == 8.0

    assert b.num_points_latitudinal == 179
    group = b.bathymetry_group[0]
    assert group.values.depth.shape == (179, 179)
    # depending on if the z is positive up or down the min depth can be 68.4 or 36.1, using min avoids the 1000000 nodata value
    assert numpy.min(group.values.depth) == pytest.approx(-68.44306, 0.0001) or numpy.min(group.values.depth) == pytest.approx(36.18454, 0.0001)
    assert numpy.min(s102obj.depth) == pytest.approx(-68.44306, 0.0001) or numpy.min(s102obj.depth) == pytest.approx(36.18454, 0.0001)
    check_s102_optional_data(s102, s102obj)


def check_s102_rat_data(s102, s102obj):
    if s102.api.EDITION >= 2.2:
        assert s102obj.root
        if hasattr(s102obj.root, 'horizontal_datum_reference'):  # v2.1 and prior
            assert s102obj.root.horizontal_datum_reference == "EPSG"
            # read with lower level h5py access too
            assert h5py_string_comp(s102obj.attrs['horizontalDatumReference'], "EPSG")
            assert s102obj.root.horizontal_datum_value == 32610
        else:  # v2.2 and later
            assert s102obj.root.horizontal_crs == 32615
        assert s102obj.root.west_bound_longitude == pytest.approx(-95.993, .001)
        assert s102obj.root.east_bound_longitude == pytest.approx(-94.806, .001)
        assert s102obj.root.north_bound_latitude == pytest.approx(26.414, .001)
        assert s102obj.root.south_bound_latitude == pytest.approx(25.186, .001)

        b = s102obj.root.bathymetry_coverage.bathymetry_coverage[0]
        if s102.api.EDITION >= 3.0:
            q = s102obj.root.quality_of_bathymetry_coverage.quality_of_bathymetry_coverage[0]
        else:
            q = s102obj.root.quality_of_survey.quality_of_survey[0]
        for instance in (b, q):
            assert instance.west_bound_longitude == pytest.approx(198254, 1)
            assert instance.east_bound_longitude == pytest.approx(319873, 1)
            assert instance.north_bound_latitude == pytest.approx(2922835, 1)
            assert instance.south_bound_latitude == pytest.approx(2788956, 1)

            assert instance.grid_origin_latitude == pytest.approx(2788956, 1)
            assert instance.grid_origin_longitude == pytest.approx(198254, 1)
            assert instance.grid_spacing_latitudinal == pytest.approx(1352.32, .01)
            assert instance.grid_spacing_longitudinal == pytest.approx(1228.48, .01)

        assert b.num_points_latitudinal == 100
        group = b.bathymetry_group[0]
        assert group.values.depth.shape == (100, 100)
        # depending on if the z is positive up or down the min depth can be 68.4 or 36.1, using min avoids the 1000000 nodata value
        assert numpy.min(group.values.depth) == pytest.approx(791.93, 0.0001) or numpy.min(group.values.depth) == pytest.approx(-3541.02, 0.0001)
        qgroup = q.quality_group[0]
        assert numpy.min(qgroup.values) == 0
        assert numpy.max(qgroup.values) == 1188907
        assert numpy.min(qgroup.values[qgroup.values > 0]) == 11134
    check_s102_optional_data(s102, s102obj)


def check_s102_optional_data(s102, s102obj):
    if s102.api.EDITION < 3.0:
        qual_name = "QualityOfSurvey"
    else:
        qual_name = "QualityOfBathymetryCoverage"
        # older apis left an empty QualityOfSurvey object - only checking 3.0 and above
        if qual_name in s102obj.keys():
            assert f"{qual_name}.01" in s102obj[qual_name].keys()  # must have feature instance QualityOfBathymetryCoverage.01
            assert qual_name in s102obj['Group_F'].keys()  # must have the Group_F entry
        else:
            with pytest.raises(KeyError):
                s102obj['Group_F'].attrs[qual_name]


def test_make_from_gdal(s102, bagname, output_path):
    remove_file(output_path)
    # the sample data is in NAD83 so does not meet spec - test that it's caught
    with pytest.raises(s102.S102Exception):
        s102.from_gdal(bagname, output_path)

    # override the metadata for the datum to WGS84 zone 10N and go from there
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
    new_s102 = s102.from_gdal(bagname, output_path, metadata=metadata)

    check_s102_data(s102, new_s102)


def test_make_from_bag(s102, bagname, output_path):
    remove_file(output_path)
    # the sample data is in NAD83 so does not meet spec - test that it's caught
    pytest.raises(s102.S102Exception, s102.from_bag, *(bagname, output_path))

    # override the metadata for the datum to WGS84 zone 10N and go from there
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
    new_s102 = s102.from_bag(bagname, output_path, metadata=metadata)

    check_s102_data(s102, new_s102)


def test_no_uncertainty(s102, bagname, output_path):
    if s102.api.EDITION >= 3.0:
        remove_file(output_path)
        # the sample data is in NAD83 so does not meet spec - test that it's caught
        pytest.raises(s102.S102Exception, s102.from_bag, *(bagname, output_path))

        # override the metadata for the datum to WGS84 zone 10N and go from there
        metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
        new_s102 = s102.from_bag(bagname, output_path, metadata=metadata)
        assert new_s102['Group_F']['BathymetryCoverage'].shape[0] == 2
        del new_s102.root.bathymetry_coverage.bathymetry_coverage[0].bathymetry_group[0].values.uncertainty
        new_s102.write()
        assert new_s102['Group_F']['BathymetryCoverage'].shape[0] == 1


def test_read_s102(s102, output_path):
    s102_read_test = s102.S102File(output_path, "r")
    check_s102_data(s102, s102_read_test)
    s102_read_test.close()


def test_copy_s102(s102, output_path, copy_path):
    s102_read_test = s102.S102File(output_path, "r")
    check_s102_data(s102, s102_read_test)
    remove_file(copy_path)
    s102_copy_root_test = s102.S102File(copy_path, "w")
    s102_copy_root_test.root = s102_read_test.root
    s102_copy_root_test.write()
    check_s102_data(s102, s102_copy_root_test)
    s102_read_test.close()
    s102_copy_root_test.close()


def test_tif_conversion(s102, tifname, temp_bagname):
    # test that nodata is changed to 1000000 and that passing a gdal instance works instead of filename
    remove_file(temp_bagname)
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
    new_s102.to_geotiffs(local_path)


def test_subdivide(s102, output_path):
    # output_path = r"C:\Pydro22_Dev\NOAA\site-packages\Python38\git_repos\s100py\tests\s102\F00788_SR_8m.bag.s102.h5"
    s102_file = s102.S102File(output_path, "r")
    s102_file.subdivide(output_path, 2, 3)


# @TODO iterate the versions using fixtures so all the functions run against all the versions
def test_s102_versions(s102, bagname):
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
    bagname = pathlib.Path(bagname)
    outname = bagname.with_suffix(f".{s102.api.EDITION}.h5")
    f20 = s102.api.S102File.from_bag(bagname, outname, metadata=metadata)
    assert f".{s102.api.EDITION}" in str(f20.root.product_specification)
    f20.close()
    open20 = s100.open(outname)
    assert f".{s102.api.EDITION}" in str(open20.root.product_specification)
    assert isinstance(open20, s102.api.S102File)
    open20.close()


def test_s102_version_upgrade(bagname):
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
    bagname = pathlib.Path(bagname)
    upgrade_name = bagname.with_suffix(".2_0_to_3_0.h5")
    f20 = v2_0.api.S102File.from_bag(bagname, upgrade_name, metadata=metadata)
    assert "2.0" in str(f20.root.product_specification)
    f20.close()
    f21 = v2_1.api.S102File.upgrade(upgrade_name)
    assert "2.1" in str(f21.root.product_specification)
    f21.close()
    f22 = v2_2.api.S102File.upgrade(upgrade_name)
    assert "2.2" in str(f22.root.product_specification)
    f22.close()
    f30 = v3_0.api.S102File.upgrade(upgrade_name)
    assert "3.0" in str(f30.root.product_specification)
    f30.close()


# tiffname = r"C:\Data\BlueTopo\RATs\BlueTopo_BC25M26L_20221102.tiff"
# output_path = r"C:\Data\BlueTopo\RATs\BlueTopo_BC25M26L_20221102.h5"
# @TODO reduce the size of the test dataset and add it to the test directory
def test_rat(s102, tif_with_rat_name, output_with_rat_path):
    # the sample data is in NAD83 so does not meet spec - test that it's caught
    pytest.raises(s102.S102Exception, s102.from_gdal, *(tif_with_rat_name, output_with_rat_path))

    # override the metadata for the datum to WGS84 zone 10N and go from there
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32615}
    # The tiff is in elevation so flip the z
    new_s102 = s102.from_gdal(tif_with_rat_name, output_with_rat_path, metadata=metadata, flip_z=True)
    check_s102_rat_data(s102, new_s102)
    new_s102.to_geotiff(local_path.joinpath("test_rat.tif"))
    upgraded = latest.api.S102File.upgrade(output_with_rat_path)


def test_edition3_changes(s102, bagname, output_path):
    """ Try an illegal vertical datum, like 50, one that is legal for s100 but not s102 and a verticalDatumReference that is not '1'"""
    if s102.api.EDITION >= 3.0:
        metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
        remove_file(output_path)
        new_s102 = s102.from_bag(bagname, output_path, metadata=metadata)
        # make sure the verticalDatum is not '1' fails
        assert new_s102.root.vertical_datum_reference.value == 1
        with pytest.raises(s102.S102Exception):
            new_s102.root.vertical_datum_reference = 2
        # Confirm that the bathy coverages doesn't have a vertical datum reference (only allowed if differs from root)
        with pytest.raises(KeyError):  # @TODO these should raise an attribute error not key error
            new_s102.root.bathymetry_coverage.bathymetry_coverage[0].vertical_datum
        with pytest.raises(KeyError):
            new_s102.root.bathymetry_coverage.bathymetry_coverage[0].vertical_datum_reference
        # root.Metadata (xml file path) is now optional - make sure an empty string is not written
        with pytest.raises(KeyError):
            new_s102.root.metadata
        # @TODO Name changed from QualityOfSurvey to QualityOfBathymetryCoverage
        assert new_s102.root.quality_of_bathymetry_coverage
        orig = new_s102.root.vertical_datum
        with pytest.raises(ValueError):
            new_s102.root.vertical_datum = 47
        assert orig == new_s102.root.vertical_datum


# def test_multiple_vertical_datums(s102, bagname, output_path):
#     raise NotImplementedError("This test is not implemented yet, see the TODOs in the function")
#     if s102.api.EDITION >= 3.0:
#         metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
#         remove_file(output_path)
#         new_s102 = s102.from_bag(bagname, output_path, metadata=metadata)
#         # @TODO test two vertical datums, Bathymetry.01 and BathymetryCoverage.02 should appear and numInstances of the BathymetryCoverage parent should be 2
#         # @TODO Confirm that all BathymetryCoverage instances have the same extents (shape)
#         # @TODO Make two quality of bathymetry that match the two bathymetry coverages
#         # @TODO make sure one of the vertical datums is against the file reference datum and doesn't have the verticalDatum attribute in it
#         # @TODO There is only one QualityOfBathymetryCoverage for all BathymetryCoverages
#         # @TODO Must use a domainExtent.polygon for each BathymetryCoverage.NN but the QualityOfBathymetryCoverage should use the bounding box
#         pass
    