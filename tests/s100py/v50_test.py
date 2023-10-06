import pytest
from s100py.s100.v5_0 import api

def test_vertical_datum():
    assert api.VERTICAL_DATUM["meanLowerLowWater"].value == 12  # old values still work?
    assert api.VERTICAL_DATUM["balticSeaChartDatum2000"].value == 44  # new value exists
    # use the dot notation too
    assert api.VERTICAL_DATUM.balticSeaChartDatum2000.value == 44

def test_interpolation_type():
    assert api.INTERPOLATION_TYPE.nearestneighbor.value == 1
    assert api.INTERPOLATION_TYPE.bilinear.value == 5
    with pytest.raises(KeyError):
        api.INTERPOLATION_TYPE['linear']

