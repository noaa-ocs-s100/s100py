import pytest

import os
import logging
import tempfile

from s100py import s100, s102, s1xx, bag_to_s102

@pytest.fixture(scope="module")
def s100_file():
    h, fstr = tempfile.mkstemp(".h5", dir=os.path.split(__file__)[0])
    os.close(h)
    # create a S100 file to do basic tests, each test of a specific spec (102, 111) will need to create it's own.
    f = s100.S100File(fstr)
    yield f
    f.close()
    os.remove(fstr)

def test_create_attrs(s100_file):
    s100_file.root.west_bound_longitude = -100.5
    s100_file.root.east_bound_longitude = 5
    assert s100_file.root.east_bound_longitude == 5
    assert s100_file.root.west_bound_longitude == -100.5

def test_enumeration(s100_file):
    s100_file.root.vertical_datum = "MLLW"
    assert s100_file.root.vertical_datum == s100_file.root.vertical_datum_type["MLLW"]
    assert s100_file.root.vertical_datum == s100_file.root.vertical_datum_type(12)
    assert s100_file.root.vertical_datum == s100.VERTICAL_DATUM["meanLowerLowWater"]

def test_write(s100_file):
    s100_file.root.add_metadata("bogus", "testing")
    s100_file.write()

def test_read(s100_file):
    read_file = s100.S100File(s100_file.filename, "r")
    read_file.read()
    assert read_file.root.east_bound_longitude == s100_file.root.east_bound_longitude
    assert read_file.root.west_bound_longitude == s100_file.root.west_bound_longitude
    assert read_file.root.vertical_datum == s100_file.root.vertical_datum
    assert read_file.root.get_metadata("bogus") == "testing"

def test_initialize_props(s100_file):
    s100_file.root.initialize_properties()
    s100_file.root.initialize_properties(True)
    assert isinstance(s100_file.root.north_bound_latitude, float)
