"""
Regression tests for GitHub issue #22: the S-102 root (WGS84) bounding box of a projected grid must enclose the grid.

The grid used here mimics the BlueTopo UTM 17N tiles in which the defect was found: 16 m cells near 25.8N, 80.1W,
where grid convergence is about 0.6 deg.  The old two-corner transform left the SE and NW corners ~300 m outside the
box in latitude; test_two_corner_box_misses_corners checks that this grid actually exercises that failure.
"""
import pathlib
import tempfile

import h5py
import numpy
import pytest
from osgeo import osr

from s100py.s102 import v2_1, v2_2, v3_0
from s100py.s102.bounds import geographic_bounds

EPSG = 32617
RES = 16.0
ROWS, COLS = 600, 500  # 9.6 km x 8 km
ORIGIN = (590008.0, 2850008.0)  # centre of the first (lower left) cell, metres

ROOT_NAMES = ("westBoundLongitude", "southBoundLatitude", "eastBoundLongitude", "northBoundLatitude")


def _extent(cell_edges):
    """(minx, miny, maxx, maxy) of the grid: outer cell edges (S-102 3.0) or node extent (S-102 2.1/2.2)."""
    ox, oy = ORIGIN
    if cell_edges:
        return ox - RES / 2, oy - RES / 2, ox + RES * COLS - RES / 2, oy + RES * ROWS - RES / 2
    return ox, oy, ox + RES * (COLS - 1), oy + RES * (ROWS - 1)


def _outline_lonlat(minx, miny, maxx, maxy, epsg=EPSG, n=2001):
    """Densely sampled outline of the projected rectangle, as lon/lat arrays."""
    t = numpy.linspace(0.0, 1.0, n)
    xs = numpy.concatenate([minx + (maxx - minx) * t, numpy.full(n, maxx), maxx - (maxx - minx) * t, numpy.full(n, minx)])
    ys = numpy.concatenate([numpy.full(n, miny), miny + (maxy - miny) * t, numpy.full(n, maxy), maxy - (maxy - miny) * t])
    src = osr.SpatialReference()
    src.ImportFromEPSG(epsg)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    wgs = osr.SpatialReference()
    wgs.ImportFromEPSG(4326)
    wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    pts = numpy.array(osr.CoordinateTransformation(src, wgs).TransformPoints(list(zip(xs, ys))))
    return pts[:, 0], pts[:, 1]


def _assert_encloses(box, lon, lat):
    west, south, east, north = (float(v) for v in box)
    assert west <= east  # this grid does not cross the antimeridian
    assert lon.min() >= west and lon.max() <= east
    assert lat.min() >= south and lat.max() <= north


def _write(api_module, path, cell_edges_version):
    depth = numpy.random.default_rng(22).uniform(1.0, 30.0, (ROWS, COLS)).astype(numpy.float32)
    uncert = numpy.full((ROWS, COLS), 0.5, dtype=numpy.float32)
    if cell_edges_version:
        metadata = {"origin": ORIGIN, "res": (RES, RES), "horizontalCRS": EPSG}
    else:
        metadata = {"origin": ORIGIN, "res": (RES, RES), "horizontalDatumReference": "EPSG", "horizontalDatumValue": EPSG}
    # nodata_value is passed explicitly: the default (None) raises in numpy.isnan inside load_arrays (separate bug)
    f = api_module.S102File.from_arrays_with_metadata(depth, uncert, metadata, str(path), nodata_value=1000000.0)
    f.close()


def _root_box(path, dtype=numpy.float32):
    with h5py.File(path, "r") as h5:
        for name in ROOT_NAMES:
            assert h5.attrs[name].dtype == dtype, name
        return tuple(h5.attrs[name] for name in ROOT_NAMES)


@pytest.fixture
def tmpdir_path():
    with tempfile.TemporaryDirectory() as d:
        yield pathlib.Path(d)


def test_two_corner_box_misses_corners():
    """Guard: the old SW/NE-only transform must fail on this grid, otherwise the tests below prove nothing."""
    minx, miny, maxx, maxy = _extent(cell_edges=True)
    lon, lat = _outline_lonlat(minx, miny, maxx, maxy)
    # n=1 samples only the start of each edge: SW, SE, NE, NW.  The old code used SW (index 0) and NE (index 2).
    corner_lon, corner_lat = _outline_lonlat(minx, miny, maxx, maxy, n=1)
    south, north = corner_lat[0], corner_lat[2]
    assert lat.min() < south or lat.max() > north
    # the miss scales with tile width x tan(convergence): ~57 m for this 8 km wide grid (280-365 m for the wider
    # BlueTopo tiles), far above float32 noise
    assert max(south - lat.min(), lat.max() - north) * 111_000 > 30


# S-102 2.1 is built on S-100 4.0, whose bound attributes are float64; 2.2 and 3.0 (S-100 5.x) use float32
@pytest.mark.parametrize("api_module,cell_edges,dtype", [(v3_0.api, True, numpy.float32), (v2_2.api, False, numpy.float32),
                                                      (v2_1.api, False, numpy.float64)],
                         ids=["v3_0", "v2_2", "v2_1"])
def test_writer_root_bounds_enclose_projected_grid(api_module, cell_edges, dtype, tmpdir_path):
    path = tmpdir_path / "bbox.h5"
    _write(api_module, path, cell_edges)
    box = _root_box(path, dtype)
    lon, lat = _outline_lonlat(*_extent(cell_edges))
    _assert_encloses(box, lon, lat)
    # and it is not grossly oversized: within ~1 m of the true envelope on every side
    west, south, east, north = (float(v) for v in box)
    assert lon.min() - west < 1e-5 and east - lon.max() < 1e-5
    assert lat.min() - south < 1e-5 and north - lat.max() < 1e-5


def test_upgrade_2_2_to_3_0_writes_enclosing_root_bounds(tmpdir_path):
    path = tmpdir_path / "upgrade.h5"
    _write(v2_2.api, path, cell_edges_version=False)
    v22_box = _root_box(path)
    upgraded = v3_0.api.S102File.upgrade(str(path))
    upgraded.close()
    v30_box = _root_box(path)
    # upgrade moves the instance box out to the cell edges, so the root box must grow to match
    assert v30_box != v22_box
    _assert_encloses(v30_box, *_outline_lonlat(*_extent(cell_edges=True)))


def test_geographic_crs_rounds_outward_to_float32():
    # values chosen so that nearest-float32 rounding would move each edge inward
    west, south, east, north = -80.1028115, 25.700000001, -80.0000001, 25.799999999
    box = geographic_bounds(4326, west, south, east, north)
    assert all(isinstance(v, numpy.float32) for v in box)
    assert box[0] <= west and box[1] <= south and box[2] >= east and box[3] >= north
    # outward by at most one float32 ulp
    for got, want in zip(box, (west, south, east, north)):
        assert abs(float(got) - want) <= abs(float(numpy.spacing(numpy.float32(want))))


def test_float64_dtype_is_unrounded():
    box = geographic_bounds(4326, -80.1028115, 25.7, -80.0, 25.8, dtype=numpy.float64)
    assert box == (-80.1028115, 25.7, -80.0, 25.8)
    assert all(isinstance(v, numpy.float64) for v in box)


def test_ups_box_containing_pole():
    """UPS North (EPSG:5041) box around the pole: all longitudes, north bound at 90."""
    west, south, east, north = geographic_bounds(5041, 1990000.0, 1990000.0, 2010000.0, 2010000.0)
    assert float(north) == 90.0
    assert float(west) == -180.0 and float(east) == 180.0
    assert 89.8 < float(south) < 90.0


def test_utm_box_crossing_antimeridian():
    """UTM 60N box straddling 180 deg (Aleutians): returned with west > east, the ISO 19115 crossing convention."""
    west, south, east, north = geographic_bounds(32660, 690000.0, 5760000.0, 720000.0, 5770000.0)
    assert float(west) > float(east)
    assert 179.0 < float(west) <= 180.0 and -180.0 <= float(east) < -179.0
