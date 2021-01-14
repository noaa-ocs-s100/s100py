import pytest

import os
import datetime
import logging
import tempfile

import numpy
import h5py

from s100py import s100


def test_feature_information_conversions():
    i = s100.FeatureInformation()
    assert i._python_datatype() == str
    i.fill_value = 10.0
    assert i.fill_value == "10.0"
    i.lower = 5
    assert i.lower == '5'
    i.upper = "hello"
    assert i.upper == "hello"
    i.datatype = 0
    assert i.datatype == 'H5T_INTEGER'
    assert i.lower == 5
    with pytest.raises(ValueError):
        i.fill_value  # this was set as a float but now changed to an int
    i.fill_value = 10.0  # now this will work since it knows the datatype
    assert i.fill_value == 10
    assert isinstance(i.fill_value, int)
    i.datatype = 'H5T_FLOAT'
    i.fill_value = 10.01
    assert i._attributes['fillValue'] == '10.01'
    i.fill_value = 10.00  # make sure the truncation of insignificant zeros works
    assert i._attributes['fillValue'] == '10'


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
    s100_file.root.issue_date = datetime.datetime(2001, 7, 17, 6, 45)
    assert s100_file.root.issue_date == datetime.date(2001, 7, 17)
    s100_file.root.issue_date = "2001-07-17T06:45:00Z"
    assert s100_file.root.issue_date == datetime.date(2001, 7, 17)
    s100_file.root.issue_time = "2001-07-17T06:45:00Z"
    s100_file.root.issue_time = "2001-07-17T06:45:00+05:30"
    s100_file.root.issue_date = "2001-07-17"
    assert s100_file.root.east_bound_longitude == 5
    assert s100_file.root.west_bound_longitude == -100.5


def test_enumeration(s100_file):
    s100_file.root.vertical_datum = "MLLW"
    assert s100_file.root.vertical_datum == s100_file.root.vertical_datum_type["MLLW"]
    assert s100_file.root.vertical_datum == s100_file.root.vertical_datum_type(12)
    assert s100_file.root.vertical_datum == s100.VERTICAL_DATUM["meanLowerLowWater"]


def test_write(s100_file):
    s100_file.root.add_data("bogus", "testing")
    s100_file.write()


def test_read(s100_file):
    read_file = s100.S100File(s100_file.filename, "r")
    read_file.read()
    assert read_file.root.east_bound_longitude == s100_file.root.east_bound_longitude
    assert read_file.root.west_bound_longitude == s100_file.root.west_bound_longitude
    assert read_file.root.vertical_datum == s100_file.root.vertical_datum
    assert read_file.root.get_data("bogus") == "testing"
    assert read_file.root.issue_date == datetime.date(2001, 7, 17)


def test_initialize_props(s100_file):
    s100_file.root.initialize_properties()
    s100_file.root.initialize_properties(True)
    assert isinstance(s100_file.root.north_bound_latitude, float)


class TestDataset(s100.FeatureInformationDataset):
    @property
    def metadata_name(self) -> str:
        return "FeatInfoTest"


def test_dataset(s100_file):
    td = TestDataset()
    rec1 = td.metadata_type()
    rec1.initialize_properties()
    rec1.code = "unit"
    rec1.name = "Latitude"
    rec2 = td.metadata_type()
    rec2.initialize_properties()
    rec2.code = "test"
    td.append(rec1)
    td.append(rec2)
    s100_file.root.add_data(td.metadata_name, td)
    s100_file.write()

    read_file = s100.S100File(s100_file.filename, "r")
    read_file.read()
    data = read_file["/" + td.metadata_name]
    assert data.shape == (2,)
    assert data.dtype[0].type == numpy.object_  # make sure the strings were stored as utf8 which comes back as an object type
    assert data[0]['code'] == "unit"
    assert data[0]['name'] == "Latitude"
    assert data[1]['code'] == "test"


def test_numpy_string(s100_file):
    numpy_strings = numpy.array(["sLat", "slong"], dtype='S')
    s100_file.root.add_data("ndStrings", numpy_strings)

    numpy_strings2 = numpy.array(["pLat", "pLon"])
    s100_file.root.add_data("narray", numpy_strings2)

    h5py_specials = numpy.ndarray([2], dtype=h5py.special_dtype(vlen=str))
    h5py_specials[0] = "Lat"
    h5py_specials[1] = "Long"
    s100_file.root.add_data("h5pyStrings", h5py_specials)
    s100_file.write()

    read_file = s100.S100File(s100_file.filename, "r")
    assert read_file['/ndStrings'][0] == "sLat"
    assert read_file['/ndStrings'][1] == "slong"
    assert read_file['/narray'][1] == "pLon"
    assert read_file['/h5pyStrings'][0] == "Lat"
