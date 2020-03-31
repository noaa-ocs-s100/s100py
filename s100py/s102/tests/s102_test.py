import pytest

import os
import logging
import tempfile

import numpy

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
    assert s102obj.root.horizontal_datum_value == '26910'
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
    s102_copy_root_test = s102.S102File(copy_path, "w")
    s102_copy_root_test.root = s102_read_test.root
    s102_copy_root_test.write()
    check_s102_data(s102_copy_root_test)
    s102_read_test.close()
    check_s102_data.close()

