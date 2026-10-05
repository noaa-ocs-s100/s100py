"""
Geographic (EPSG:4326) bounding box for an S-102 grid extent.

Shared by the S-102 v2.1, v2.2 and v3.0 writers and by the v2.2 -> v3.0 upgrade
(see GitHub issue #22).

Why this exists
---------------
The root ``westBoundLongitude`` / ``eastBoundLongitude`` / ``southBoundLatitude`` /
``northBoundLatitude`` attributes are always geographic, while S-102 grids are often
stored in a projected CRS (UTM or UPS).  The writers previously projected only the
SW and NE corners of the grid to lon/lat.  Grid convergence rotates a projected
rectangle relative to the graticule, so the SE and NW corners (and the interior of a
constant-northing edge, which is curved in lat/lon) fall outside a box built from two
corners.  For BlueTopo UTM 17N tiles (convergence 0.57-0.73 deg) the error was
280-365 m in latitude, leaving valid depth cells outside the box.

S-102 Ed 3.0.0 Table 10-2 requires the box to *encompass* all data, using outer cell
boundaries; S-100 Ed 5.2.1 Table 10c-6 requires the min/max of the data.

Decisions
---------
* ``osr.CoordinateTransformation.TransformBounds`` (GDAL >= 3.4) is used rather than
  hand-rolled corner/edge sampling.  It densifies every edge and also handles the
  two cases a min/max of sampled points gets wrong, both reachable with the CRSs
  S-102 allows: a box containing a pole (UPS, EPSG:5041/5042) and a box crossing
  the antimeridian (UTM zones 1 and 60).  In the antimeridian case it returns
  ``west > east``, which is the ISO 19115 convention for EX_GeographicBoundingBox.
* There is deliberately no fallback for older GDAL: a silently wrong bounding box is
  the defect being fixed, so an old GDAL raises instead.
* Both axis orders are pinned to traditional GIS order (x/easting, y/northing;
  lon, lat) so the result does not depend on the EPSG authority axis order of
  EPSG:4326, which GDAL 3 honours by default.
* From S-100 5.0 (S-102 2.2 and 3.0) the attributes are stored as float32.  Rounding a
  float64 value to the nearest float32 can move an edge inward by up to half a float32
  ulp (about 4e-6 deg at 180 deg, ~0.4 m), so each edge is rounded *outward* to the
  next representable value of the storage dtype.  This applies to geographic grids
  too, which pass through unchanged apart from that rounding.  S-102 2.1 (S-100 4.0)
  stores float64 and passes ``dtype=numpy.float64``, for which the rounding is a no-op.
* The outward-rounding comparison is done in float64.  Under NumPy 2 (NEP 50) a
  comparison between a float32 scalar and a Python float casts the Python float to
  float32 first, so ``numpy.float32(v) < v`` is False even when the float32 is
  smaller - the first version of this module made exactly that mistake, and the
  regression tests caught it.
"""
import numpy

try:
    from osgeo import osr
except ImportError:  # GDAL is optional at import time elsewhere in s100py; fail when actually needed
    osr = None

#: Points added along each edge by TransformBounds (GDAL's own recommended default).
DENSIFY_POINTS = 21


def _round_down(value, dtype):
    """Largest value of ``dtype`` that is <= value (compared in float64, see module notes on NEP 50)."""
    f = dtype(value)
    if numpy.float64(f) > numpy.float64(value):
        f = numpy.nextafter(f, dtype(-numpy.inf))
    return f


def _round_up(value, dtype):
    """Smallest value of ``dtype`` that is >= value (compared in float64, see module notes on NEP 50)."""
    f = dtype(value)
    if numpy.float64(f) < numpy.float64(value):
        f = numpy.nextafter(f, dtype(numpy.inf))
    return f


def geographic_bounds(epsg, minx, miny, maxx, maxy, dtype=numpy.float32):
    """ Compute the EPSG:4326 bounding box that encloses a grid extent.

    Parameters
    ----------
    epsg
        EPSG code of the grid's horizontal CRS (projected or geographic).
    minx, miny, maxx, maxy
        Grid extent in that CRS, x = easting/longitude, y = northing/latitude.
        For S-102 v3.0 these are the outer cell edges (Table 10-2).
    dtype
        Storage type of the root bound attributes: numpy.float32 (S-100 5.x) or numpy.float64 (S-100 4.0).

    Returns
    -------
    tuple of dtype
        (west_lon, south_lat, east_lon, north_lat), each rounded outward to ``dtype``.
        For a projected box crossing the antimeridian, west_lon > east_lon.

    Raises
    ------
    RuntimeError
        If GDAL is not installed, or is older than 3.4 (no TransformBounds) and the CRS is projected.
    """
    if osr is None:
        raise RuntimeError("GDAL (osgeo.osr) is required to compute S-102 geographic bounds")
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg))
    if srs.IsProjected():
        if not hasattr(osr.CoordinateTransformation, "TransformBounds"):
            from osgeo import gdal
            raise RuntimeError("GDAL >= 3.4 is required (CoordinateTransformation.TransformBounds) to compute the "
                               f"geographic bounds of a projected S-102 grid; found GDAL {gdal.__version__}")
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        wgs = osr.SpatialReference()
        wgs.ImportFromEPSG(4326)  # S-102 root bounds are WGS84 geographic
        wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        transform = osr.CoordinateTransformation(srs, wgs)
        west, south, east, north = transform.TransformBounds(minx, miny, maxx, maxy, DENSIFY_POINTS)
    else:
        west, south, east, north = minx, miny, maxx, maxy
    return _round_down(west, dtype), _round_down(south, dtype), _round_up(east, dtype), _round_up(north, dtype)
