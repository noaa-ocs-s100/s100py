import pytest
from s100py.v5_0 import s100

def test_vertical_datum():
    assert s100.VERTICAL_DATUM["meanLowerLowWater"] == 12  # old values still work?
    assert s100.VERTICAL_DATUM["balticSeaChartDatum2000"] == 44  # new value exists
    # use the dot notation too
    assert s100.VERTICAL_DATUM.balticSeaChartDatum2000 == 44

def test_interpolation_type():
    assert s100.INTERPOLATION_TYPE.nearestneighbor == 1
    assert s100.INTERPOLATION_TYPE.bilinear == 5
    pytest.raises(KeyError, s100.INTERPOLATION_TYPE['linear'])

