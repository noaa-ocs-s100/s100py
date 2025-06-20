import re
import os
import pathlib
from abc import ABC, abstractmethod
from typing import Callable, Iterator, Union, Optional, List, Type
import logging
import datetime
from enum import Enum

import h5py
from osgeo import gdal, osr, ogr
# @todo - consider removing the numpy dependence
import numpy

try:
    from .. import s1xx
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py.s100"

from s100py.s1xx import s1xx_sequence, S1xxObject, S1xxDatasetBase, S1XXFile, S1xxGridsBase, h5py_string_dtype, is_sub_class, h5py_string_comp

EDITION = 5.1
PRODUCT_SPECIFICATION = 'INT.IHO.S-100.5.1'


# @TODO add support for Positioning group (S100 v5.0 10c-9.10.1, 10c-9.10.2 and table 10c-16) for DataCodingFormats 1,3,4,7,8
#   Basically makes a 1D dataset array of the positions and attributes (like Z).
#   For TINs  DCF=7
#   Optionally has triangles dataset NumTriangles x 3 -- array of the 3 node indices for positions
#   Optionally has adjacency dataset NumTriangles x3 -- array of the indices or the triangles that share an edge with each triangle side
#   adjacency[i][0] = triangle adjacent to the edge specified by triangles[i][0] & triangles[i][1]
#   adjacency[i][1] = triangle adjacent to edge triangles[i][1] & triangles[i][2]
#   adjacency[i][2] = triangle adjacent to edge triangles[i][2] & triangles[i][0]
#   Elements for edges without adjacent triangles are filled with the value -1
#   Applicable only for dataCodingFormat = 7 (TIN), but optional even for TIN.


class S100Exception(Exception):
    pass


H5T_CLASS_T = {
    h5py.h5t.NO_CLASS: 'H5T_NO_CLASS',
    h5py.h5t.INTEGER: 'H5T_INTEGER',
    h5py.h5t.FLOAT: 'H5T_NATIVE_FLOAT',
    h5py.h5t.TIME: 'H5T_TIME',
    h5py.h5t.STRING: 'H5T_STRING',
    h5py.h5t.BITFIELD: 'H5T_BITFIELD',
    h5py.h5t.OPAQUE: 'H5T_OPAQUE',
    h5py.h5t.COMPOUND: 'H5T_COMPOUND',
    h5py.h5t.REFERENCE: 'H5T_REFERENCE',
    h5py.h5t.ENUM: 'H5T_ENUM',
    h5py.h5t.VLEN: 'H5T_VLEN',
    h5py.h5t.ARRAY: 'H5T_ARRAY',
    h5py.h5t.NATIVE_INT8: 'H5T_NATIVE_INT8',
    h5py.h5t.NATIVE_UINT8: 'H5T_NATIVE_UINT8',
    h5py.h5t.NATIVE_INT16: 'H5T_NATIVE_INT16',
    h5py.h5t.NATIVE_UINT16: 'H5T_NATIVE_UINT16',
    h5py.h5t.NATIVE_INT32: 'H5T_NATIVE_INT32',
    h5py.h5t.NATIVE_UINT32: 'H5T_NATIVE_UINT32',
    h5py.h5t.NATIVE_INT64: 'H5T_NATIVE_INT64',
    h5py.h5t.NATIVE_UINT64: 'H5T_NATIVE_UINT64',
    h5py.h5t.C_S1: 'H5T_C_S1'
}

# a dictionary of what python type each string representation of an HDF5 type goes to.  e.g. h5py.h5t.INTEGER -> int
_H5T_Types = {
    int: [H5T_CLASS_T[v] for v in (h5py.h5t.INTEGER, h5py.h5t.NATIVE_INT8, h5py.h5t.NATIVE_UINT8, h5py.h5t.NATIVE_INT16, h5py.h5t.NATIVE_UINT16,
                                   h5py.h5t.NATIVE_INT32, h5py.h5t.NATIVE_UINT32, h5py.h5t.NATIVE_INT64, h5py.h5t.NATIVE_UINT64)],
    float: [H5T_CLASS_T[v] for v in (h5py.h5t.FLOAT,)] + ["H5T_FLOAT"],
    str: [H5T_CLASS_T[v] for v in (h5py.h5t.C_S1, h5py.h5t.STRING)],
}


# noinspection PyPep8Naming
class VERTICAL_DATUM(Enum):
    """ Note: while a Vertical Datum can be created with the shorthand aliases, ex: MLWS, the string written and
    returned from the file/S100 object will be the official long name, e.g. "meanLowWaterSprings" etc.
    S100 Part 4a Metadata

    S100 v5 Part 17 Vertical and Sounding Datum
    Added balticSeaChartDatum2000 = 44 to hydrographicZero = 49
    """
    meanLowWaterSprings = 1
    MLWS = 1
    meanLowerLowWaterSprings = 2
    meanSeaLevel = 3
    MSL = 3
    lowestLowWater = 4
    meanLowWater = 5
    MLW = 5
    lowestLowWaterSprings = 6
    approximateMeanLowWaterSprings = 7
    indianSpringLowWater = 8
    lowWaterSprings = 9
    approximateLowestAstronomicalTide = 10
    nearlyLowestLowWater = 11
    meanLowerLowWater = 12
    MLLW = 12
    lowWater = 13
    LW = 13
    approximateMeanLowWater = 14
    approximateMeanLowerLowWater = 15
    meanHighWater = 16
    MHW = 16
    meanHighWaterSprings = 17
    MHWS = 17
    highWater = 18
    approximateMeanSeaLevel = 19
    highWaterSprings = 20
    meanHigherHighWater = 21
    MHHW = 21
    equinoctialSpringLowWater = 22
    lowestAstronomicalTide = 23
    LAT = 23
    localDatum = 24
    internationalGreatLakesDatum1985 = 25
    meanWaterLevel = 26
    lowerLowWaterLargeTide = 27
    higherHighWaterLargeTide = 28
    nearlyHighestHighWater = 29
    highestAstronomicalTide = 30
    HAT = 30
    balticSeaChartDatum2000 = 44
    internationalGreatLakesDatum2020 = 46
    seaFloor = 47
    seaSurface = 48
    hydrographicZero = 49


HORIZONTAL_DATUM_REFERENCE = 'EPSG'
REGULAR = 'Regularly-gridded arrays2'
# Table 10c-4 and Table 10c-23, added values 8 and 9 for Fixed stations - stationwise and Feature oriented regular grid
# v5.0 added camel case names
DATA_CODING_FORMAT = Enum(value="DATA_CODING_FORMAT",
                          names=[
                              ('fixedStations', 1),
                              ('Time series at fixed stations', 1),
                              ('regularGrid', 2),
                              ('Regularly-gridded arrays', 2),
                              ('ungeorectifiedGrid', 3),
                              ('Ungeorectified gridded arrays', 3),
                              ('movingPlatform', 4),
                              ('Moving platform', 4),
                              ('irregularGrid', 5),
                              ('Irregular grid', 5),
                              ('variableCellSize', 6),
                              ('Variable cell size', 6),
                              ('TIN', 7),
                              ('Triangulated irregular network', 7),
                              ('stationwiseFixed', 8),
                              ('Fixed stations stationwise', 8),
                              ('featureOrientedRegularGrid', 9),
                              ('Feature oriented regular grid', 9),
                              # alternate shortcut names that also show up in sphinx, these will be stored with full names including spaces
                              ('TIME', 1),
                              ('REGULAR', 2),
                              ('UNGEORECTIFIED', 3),
                              ('MOVING', 4),
                              ('IRREGULAR', 5),
                              ('VARIABLE', 6),
                              ('STATIONWISE', 8),
                              ('FEATURE_REGULAR', 9)
                          ]
                          )


"""
Sphinx is not interpreting the enum names properly when there are spaces. The correct enum names with spaces are::

  ('Time series at fixed stations', 1),
  ('Regularly-gridded arrays', 2),
  ('Ungeorectified gridded arrays', 3),
  ('Moving platform', 4),
  ('Irregular grid', 5),
  ('Variable cell size', 6),
  ('TIN', 7),
  ('Fixed stations - stationwise', 8),
  ('Feature oriented regular grid', 9),
"""


# noinspection PyPep8Naming
class INTERPOLATION_TYPE(Enum):
    """ S100 v4.0 table 10c-21
    Enumeration S100_CV_InterpolationMethod Codes for interpolation methods between known feature attribute
    values associated with geometric objects in the domain of the discrete coverage
    Extension of ISO 19123
    CV_InterpolationMethod

    Literal nearestneighbor
    Assign the feature attribute value associated with the nearest domain object in the domain of the coverage
    1 Any type of coverage

    Literal linear
    Assign the value computed by a linear function along a line segment connecting two point value pairs, or along a curve with positions are described by values
    of an arc-length parameter
    2 Only segmented curves

    Literal quadratic
    Assign the value computed by a quadratic function of distance along a value segment
    3 Only segmented curves

    Literal cubic
    Assign the value computed by a cubic function of distance along a value segment
    4 Only segmented curves

    Literal bilinear
    Assign a value computed by using a bilinear function of position within the grid cell
    5 Only quadrilateral grids

    Literal biquadratic
    Assign a value computed by using a biquadratic function of position within the grid cell
    6 Only quadrilateral grids

    Literal bicubic
    Assign a value computed by using a bicubic function of position within the grid cell
    7 Only quadrilateral grids

    Literal lostarea
    Assign a value computed by using the lost area method described in ISO 19123
    8 Only Thiessen polygons

    Literal barycentric
    Assign a value computed by using the barycentric method described in ISO 19123
    9 Only TIN

    Literal discrete
    No interpolation method applies to the coverage
    10
    """
    nearestneighbor = 1
    # linear = 2  # these are struck through in the spec
    # quadratic = 3  # these are struck through in the spec
    # cubic = 4  # these are struck through in the spec
    bilinear = 5
    biquadratic = 6
    bicubic = 7
    lostarea = 8
    barycentric = 9
    discrete = 10


# noinspection PyPep8Naming
class COMMON_POINT_RULE(Enum):
    """ S100 v4.0 Table 10c-19
    """
    average = 1
    low = 2
    high = 3
    all = 4
    # start = 5  # these are struck through in the spec
    # end = 6  # these are struck through in the spec


# noinspection PyPep8Naming
class SEQUENCING_RULE_TYPE(Enum):
    """ S100 v4.0 Table 10c-20
    """
    linear = 1
    boustrophedonic = 2
    CantorDiagonal = 3
    spiral = 4
    Morton = 5
    Hilbert = 6


SEQUENCING_RULE_SCAN_DIRECTION = 'longitude,latitude'  # See S100 v4.0 Table 10c-10
START_SEQUENCE = '0,0'

# Table 10c-24
TYPE_OF_HORIZONTAL_CRS = Enum(value="TYPE_OF_HORIZONTAL_CRS",
                                names=[
                                    ('geodeticCRS2D', 1),
                                    ('projectedCRS', 2),
                              ])
# Table 10c-25
VERTICAL_COORDINATE_BASE = Enum(value="VERTICAL_COORDINATE_BASE",
                                names=[
                                    ('seaSurface', 1),
                                    ('verticalDatum', 2),
                                    ('seaBottom', 3),
                              ])

# Table 10c-26
VERTICAL_DATUM_REFERENCE = Enum(value="VERTICAL_DATUM_REFERENCE",
                                names=[
                                    ('s100VerticalDatum', 1),
                                    ('EPSG', 2),
                                ])

# Table 10c-27
PROJECTION_METHOD = Enum(value="PROJECTION_METHODS",
                            names=[
                                ('Mercator', 9805),
                                ('Transverse Mercator', 9807),
                                ('Oblique Mercator', 9815),
                                ('Hotine Oblique Mercator', 9812),
                                ('Lambert Conic Conformal 1SP', 9801),
                                ('Lambert Conic Conformal 2SP', 9802),
                                ('Oblique Stereographic', 9809),
                                ('Polar Stereographic', 9810),
                                ('Krovak Oblique Conic Conformal', 9819),
                                ('American Polyconic', 9818),
                                ('Albers Equal Area', 9822),
                                ('Lambert Azimuthal Equal Area', 9820),
                          ])

# Table 10c-6
HORIZONTAL_CS = Enum(value="HORIZONTAL_CS",
                        names=[
                            ('Lat, Lon (degrees)', 6422),
                            ('Easting, Northing (m)', 4400),
                            ('Northing, Easting (m)', 4500),
                     ])

#Table 10c-6
VERTICAL_CS = Enum(value="VERTICAL_CS",
                        names=[
                            ('Depth', 6498),
                            ('Height', 6499),
                     ])

DATA_OFFSET_CODE = Enum(value="DATA_OFFSET_CODE",
                        names=[
                            ('XMin, YMin', 1),
                            ('XMax, YMax', 2),
                            ('XMax, YMin', 3),
                            ('XMin, YMax', 4),
                            ('Barycenter', 5),
                            ('Lower left', 1),
                            ('Upper right', 2),
                            ('Lower right', 3),
                            ('Upper left', 4),
                            ('Centroid', 5),
                        ])


# Note there is a DirectPosition in the S100 spec which is possibly different than the DirectPosition in the S102 spec
class DirectPosition(S1xxObject):
    """ S102 4.2.1.1.4 of v2.0.0
    """
    __coordinate_hdf_name__ = "coordinate"  #: HDF5 naming
    __dimension_hdf_name__ = "dimension"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def coordinate(self) -> s1xx_sequence:
        return self._attributes[self.__coordinate_hdf_name__]

    @coordinate.setter
    def coordinate(self, val: s1xx_sequence):
        self._attributes[self.__coordinate_hdf_name__] = val

    @property
    def __coordinate_type__(self):
        return numpy.ndarray

    def coordinate_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.coordinate = self.__coordinate_type__([2], numpy.float64)

    @property
    def dimension(self) -> int:
        return self._attributes[self.__dimension_hdf_name__]

    @dimension.setter
    def dimension(self, val: int):
        self._attributes[self.__dimension_hdf_name__] = val

    @property
    def __dimension_type__(self):
        return int

    def dimension_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.dimension = self.__dimension_type__()


class GridCoordinate(S1xxObject):
    """ 4.2.1.1.6 of S102 v2.0.0, references ISO 19123
    """

    __coord_values_hdf_name__ = "coordValues"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def coord_values(self) -> s1xx_sequence:
        """The attribute coordValues has the value class Sequence Integer that shall hold one integer value for each dimension of the grid.
        The ordering of these coordinate values shall be the same as that of the elements of axisNames.
        The value of a single coordinate shall be the number of offsets from the origin of the grid in the direction of a specific axis"""
        return self._attributes[self.__coord_values_hdf_name__]

    @coord_values.setter
    def coord_values(self, val: s1xx_sequence):
        self._attributes[self.__coord_values_hdf_name__] = val

    @property
    def __coord_values_type__(self):
        return numpy.ndarray

    def coord_values_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.coord_values = self.__coord_values_type__([2], numpy.int64)


class GridEnvelope(S1xxObject):
    """ S102  4.2.1.1.5 of v2.0.0, , references ISO 19123
    While I would think that the envelope would describe the real world extents of the grid,
    in the docs it describes the envelope as specifying the row/column offsets for the lower left and upper right
    coordinates using the integer indices (S100 and ISO 19123 sec. 8.3).  The real world coordinates are in the origin and offsetVectors instead.
    So this seems unnecessary since the value seems the same as the size of the matrix held, which can be gotten by reading that instead.

    https://www.fgdc.gov/standards/projects/frameword-data-standard/GI_FrameworkDataStandard_Part3_Elevation.doc/at_download/file&usg=AOvVaw07QEsNy5urachwIO1e4ALU
    """

    __low_hdf_name__ = "low"  #: HDF5 naming
    __high_hdf_name__ = "high"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def low(self) -> S1xxObject:
        return self._attributes[self.__low_hdf_name__]

    @low.setter
    def low(self, val: S1xxObject):
        self._attributes[self.__low_hdf_name__] = val

    @property
    def __low_type__(self):
        return GridCoordinate

    def low_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.low = self.__low_type__()

    @property
    def high(self) -> S1xxObject:
        return self._attributes[self.__high_hdf_name__]

    @high.setter
    def high(self, val: S1xxObject):
        self._attributes[self.__high_hdf_name__] = val

    @property
    def __high_type__(self):
        return GridCoordinate

    def high_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.high = self.__high_type__()


class SequenceRule(S1xxObject):
    """ 4.2.1.1.7 (and .8) of v2.0.0
    CV_SequenceRule specified in ISO 19123
    """

    __type_hdf_name__ = "type"  #: HDF5 naming
    __scan_direction_hdf_name__ = "scanDirection"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def type(self) -> str:
        """From S100 - CV_SequenceRule specified in ISO 19123. Only the values "linear" (for a simple regular cell size grid) and "Morton"
        (for a Quad Tree Grid) shall be used for data that conforms to this standard.
        While S102 further specifies - The default value is "linear". No other options are allowed.

        CodeList types are sets of strings (enumerations if all options are known).
        For SequenceType linear is lowercase while Morton is capitalized.
        """
        return self._attributes[self.__type_hdf_name__]

    @type.setter
    def type(self, val: str):
        self._attributes[self.__type_hdf_name__] = val

    @property
    def __type_type__(self):
        return str

    def type_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.type = self.__type_type__("linear")

    @property
    def scan_direction(self) -> s1xx_sequence:
        """The attribute scanDirection has the value class Sequence<CharacterString> a list of axis names that indicates
        the order in which grid points shall be mapped to position within the sequence of records of feature attribute values.
        The scan direction for all layers in S-102 is "Longitude" and "Latitude" or west to east, then south to north.
        """
        return self._attributes[self.__scan_direction_hdf_name__]

    @scan_direction.setter
    def scan_direction(self, val: s1xx_sequence):
        self._attributes[self.__scan_direction_hdf_name__] = val

    @property
    def __scan_direction_type__(self):
        return str

    def scan_direction_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.scan_direction = self.__scan_direction_type__("Longitude, Latitude")


class Point(S1xxObject):
    """ 4.2.1.1.11 of v2.0.0
    The class GM_Point is taken from ISO 19107 and is the basic data type for a geometric object consisting of one and only one point.
    """
    __position_hdf_name__ = "position"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def position(self) -> DirectPosition:
        """ DirectPosition - see Figure 7-3 in S100 v4.0.0
        """
        return self._attributes[self.__position_hdf_name__]

    @position.setter
    def position(self, val: DirectPosition):
        self._attributes[self.__position_hdf_name__] = val

    @property
    def __position_type__(self):
        return DirectPosition

    def position_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.position = self.__position_type__()


# class GeographicExtent(S1xxObject):
#     """ 4.2.1.1.12 of v2.0.0
#     The class EX_GeographicExtent is a metadata class from ISO 19115.
#     It is a component of the metaclass EX_Extent.
#     The use of EX_GeographicExtent is optional.
#     When used it describes the spatial boundaries of the Tracking List elements within the bounds established by CV_GridEnvelope for the BathymetryCoverage.
#     That is, the tracking list may carry information corresponding only to a portion of the spatial extent covered by the BathymetryCoverage.
#     There is one attribute and one subtype.
#     """
#
#     __extent_type_code_hdf_name__ = "extentTypeCode"  #: HDF5 naming
#
#     @property
#     def __version__(self) -> int:
#         return 1
#
#     @property
#     def extent_type_code(self) -> bool:
#         """ The attribute extentTypeCode is a Boolean value.
#         It is used to indicate whether the bounding polygon/box encompasses an area covered by the data or an area where data is not present.
#         In S-102 it is set to 1.
#         """
#         return self._attributes[self.__extent_type_code_hdf_name__]
#
#     @extent_type_code.setter
#     def extent_type_code(self, val: bool):
#         self._attributes[self.__extent_type_code_hdf_name__] = val
#
#     @property
#     def __extent_type_code_type__(self):
#         return bool
#
#     def extent_type_code_create(self):
#         # noinspection PyAttributeOutsideInit
#         # pylint: disable=attribute-defined-outside-init
#         self.extent_type_code = self.__extent_type_code_type__()


class GeographicBoundingBox(S1xxObject):  # GeographicExtent
    """ S100 Tables 10C-6 and 10c-12
    see also 4.2.1.1.13 of S102 v2.0.0
    The class EX_GeographicBoundingBox is a metadata class from ISO 19115.
    It is a subtype of the abstract class EX_GeographicExtent.
    It defines a bounding box used to indicate the spatial boundaries of the tracking list elements within the
    bounds established by CV_GridEnvelope for the BathymetryCoverage.

    From S100:
    The geographic extent of the grid, as a bounding box
    Ref. domainExtent: EX_GeographicExtent > EX_GeographicBoundingBox
    Either this or the domainExtent dataset must be populated
    The bounds must either all be populated or all omitted
    """

    __west_bound_longitude_hdf_name__ = "westBoundLongitude"  #: HDF5 naming
    __east_bound_longitude_hdf_name__ = "eastBoundLongitude"  #: HDF5 naming
    __south_bound_latitude_hdf_name__ = "southBoundLatitude"  #: HDF5 naming
    __north_bound_latitude_hdf_name__ = "northBoundLatitude"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def west_bound_longitude(self) -> float:
        """Western extent"""
        return self._attributes[self.__west_bound_longitude_hdf_name__]

    @west_bound_longitude.setter
    def west_bound_longitude(self, val: float):
        self._attributes[self.__west_bound_longitude_hdf_name__] = val

    @property
    def __west_bound_longitude_type__(self):
        return numpy.float32

    def west_bound_longitude_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.west_bound_longitude = self.__west_bound_longitude_type__()

    @property
    def east_bound_longitude(self) -> float:
        """Eastern extent"""
        return self._attributes[self.__east_bound_longitude_hdf_name__]

    @east_bound_longitude.setter
    def east_bound_longitude(self, val: float):
        self._attributes[self.__east_bound_longitude_hdf_name__] = val

    @property
    def __east_bound_longitude_type__(self):
        return numpy.float32

    def east_bound_longitude_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.east_bound_longitude = self.__east_bound_longitude_type__()

    @property
    def south_bound_latitude(self) -> float:
        """Southern extent"""
        return self._attributes[self.__south_bound_latitude_hdf_name__]

    @south_bound_latitude.setter
    def south_bound_latitude(self, val: float):
        self._attributes[self.__south_bound_latitude_hdf_name__] = val

    @property
    def __south_bound_latitude_type__(self):
        return numpy.float32

    def south_bound_latitude_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.south_bound_latitude = self.__south_bound_latitude_type__()

    @property
    def north_bound_latitude(self) -> float:
        """Northern extent"""
        return self._attributes[self.__north_bound_latitude_hdf_name__]

    @north_bound_latitude.setter
    def north_bound_latitude(self, val: float):
        self._attributes[self.__north_bound_latitude_hdf_name__] = val

    @property
    def __north_bound_latitude_type__(self):
        return numpy.float32

    def north_bound_latitude_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.north_bound_latitude = self.__north_bound_latitude_type__()


class VertexPoint(S1xxObject):
    """ From Figure 8-21 in S100 v4.0.0

    """

    __geometry_hdf_name__ = "geometry"  #: HDF5 naming
    __value_hdf_name__ = "value"  # HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def geometry(self) -> Point:
        """ Derived from ISO 19107, referenced figure 7-3 and 8-A-5 of S100 v4.0.0"""
        return self._attributes[self.__geometry_hdf_name__]

    @geometry.setter
    def geometry(self, val: Point):
        self._attributes[self.__geometry_hdf_name__] = val

    @property
    def __geometry_type__(self):
        return Point

    def geometry_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.geometry = self.__geometry_type__()

    @property
    def value(self) -> s1xx_sequence:
        """The attribute value has the value class Record which is a sequence of value items that shall assign values to the discrete grid point.
        There are two values in each record in the S102_TrackingListValues class.
        These are the depth and the uncertainty values that were overridden in corresponding grid coverages

        'value' in tracking list should be HDF5 dataset of (depth, uncertainty)
        which are matched to the listSeries which holds the indices of the data locations

        It is an ISO 19103 class Record
        """
        return self._attributes[self.__value_hdf_name__]

    @value.setter
    def value(self, val: s1xx_sequence):
        self._attributes[self.__value_hdf_name__] = val

    @property
    def __value_type__(self):
        return numpy.ndarray

    def value_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.value = self.__value_type__([2, ], numpy.float64)


# FIXME @TODO Add base class (maybe full implementation for many of the datasets) for FeatureInstanceBase
#  Positioning, Group_IDX, Group_TL, cellGeometry, uncertainty, extent, domainExtent.verticalElement, domainExtent.polygon
#  This would be added to the FeatureInstanceBase - See S100 v4.0 table 10c-11 and 10c-12


class FeatureInstanceBase(GeographicBoundingBox):
    """ The feature instance group attributes from table 10c-12 in S100 spec
    """

    __vertical_extent_minimum_z_hdf_name__ = "verticalExtent.minimumZ"
    __vertical_extent_maximum_z_hdf_name__ = "verticalExtent.maximumZ"
    __num_grp_hdf_name__ = "numGRP"
    __number_of_times_hdf_name__ = "numberOfTimes"
    __time_record_interval_hdf_name__ = "timeRecordInterval"
    # @TODO  @FIXME -- first and last records are supposed to be datetime but S100 doc says 'character'  Need to create a datetime handler
    __date_time_of_first_record_hdf_name__ = "dateTimeOfFirstRecord"
    __date_time_of_last_record_hdf_name__ = "dateTimeOfLastRecord"
    __extent_hdf_name__ = "extent"  #: HDF5 naming

    def write(self, hdf5_object):
        super().write(hdf5_object)

    @property
    def vertical_extent_minimum_z(self) -> float:
        """Vertical extent of 3-D grids
        minimumZ, maximumZ: Minimum and maximum values of the grid’s spatial extent
        along the vertical direction. They are encoded as separate attributes"""
        return self._attributes[self.__vertical_extent_minimum_z_hdf_name__]

    @vertical_extent_minimum_z.setter
    def vertical_extent_minimum_z(self, val: float):
        self._attributes[self.__vertical_extent_minimum_z_hdf_name__] = val

    @property
    def __vertical_extent_minimum_z_type__(self) -> Type[float]:
        return float

    def vertical_extent_minimum_z_create(self):
        """ Creates a blank, empty or zero value for vertical_extent_minimum_z"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_extent_minimum_z = self.__vertical_extent_minimum_z_type__()

    @property
    def vertical_extent_maximum_z(self) -> float:
        """Vertical extent of 3-D grids
        minimumZ, maximumZ: Minimum and maximum values of the grid’s spatial extent
        along the vertical direction. They are encoded as separate attributes"""
        return self._attributes[self.__vertical_extent_maximum_z_hdf_name__]

    @vertical_extent_maximum_z.setter
    def vertical_extent_maximum_z(self, val: float):
        self._attributes[self.__vertical_extent_maximum_z_hdf_name__] = val

    @property
    def __vertical_extent_maximum_z_type__(self) -> Type[float]:
        return float

    def vertical_extent_maximum_z_create(self):
        """ Creates a blank, empty or zero value for vertical_extent_maximum_z"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_extent_maximum_z = self.__vertical_extent_maximum_z_type__()

    @property
    def num_grp(self) -> int:
        """The number of data values groups contained in this instance group"""
        return self._attributes[self.__num_grp_hdf_name__]

    @num_grp.setter
    def num_grp(self, val: int):
        self._attributes[self.__num_grp_hdf_name__] = val

    @property
    def __num_grp_type__(self) -> Type[int]:
        return numpy.uint32

    def num_grp_create(self):
        """ Creates a blank, empty or zero value for num_grp"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.num_grp = self.__num_grp_type__()

    @property
    def number_of_times(self) -> int:
        """The total number of time records.
        Time series data only"""
        return self._attributes[self.__number_of_times_hdf_name__]

    @number_of_times.setter
    def number_of_times(self, val: int):
        self._attributes[self.__number_of_times_hdf_name__] = val

    @property
    def __number_of_times_type__(self) -> Type[int]:
        return int

    def number_of_times_create(self):
        """ Creates a blank, empty or zero value for number_of_times"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.number_of_times = self.__number_of_times_type__()

    @property
    def time_record_interval(self) -> int:
        """The interval between time records. Units: Seconds.
        Time series data only"""
        return self._attributes[self.__time_record_interval_hdf_name__]

    @time_record_interval.setter
    def time_record_interval(self, val: int):
        self._attributes[self.__time_record_interval_hdf_name__] = val

    @property
    def __time_record_interval_type__(self) -> Type[int]:
        return int

    def time_record_interval_create(self):
        """ Creates a blank, empty or zero value for time_record_interval"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.time_record_interval = self.__time_record_interval_type__()

    @property
    def date_time_of_first_record(self) -> str:
        """The validity time of the earliest time record. Units: DateTime.
        Time series data only"""
        return self._attributes[self.__date_time_of_first_record_hdf_name__]

    @date_time_of_first_record.setter
    def date_time_of_first_record(self, val: str):
        self._attributes[self.__date_time_of_first_record_hdf_name__] = val

    @property
    def __date_time_of_first_record_type__(self) -> Type[str]:
        return str

    def date_time_of_first_record_create(self):
        """ Creates a blank, empty or zero value for date_time_of_first_record"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.date_time_of_first_record = self.__date_time_of_first_record_type__()

    @property
    def date_time_of_last_record(self) -> str:
        """The validity time of the latest time record. Units: DateTime.
        Time series data only"""
        return self._attributes[self.__date_time_of_last_record_hdf_name__]

    @date_time_of_last_record.setter
    def date_time_of_last_record(self, val: str):
        self._attributes[self.__date_time_of_last_record_hdf_name__] = val

    @property
    def __date_time_of_last_record_type__(self) -> Type[str]:
        return str

    def date_time_of_last_record_create(self):
        """ Creates a blank, empty or zero value for date_time_of_last_record"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.date_time_of_last_record = self.__date_time_of_last_record_type__()

    @property
    def extent(self) -> GridEnvelope:
        """ v5.0 Table 10c-11 Structure fo Feature Instance Groups (optional)

        From 4.2.1.1.1.12,
        The attribute extent has the value class CV_GridEnvelope that shall contain the extent of the spatial domain of the coverage.
        It uses the value class CV_GridEnvelope which provides the grid coordinate values for the diametrically opposed corners of the grid.
        The default is that this value is derived from the bounding box for the data set or tile in a multi tile data set"""
        return self._attributes[self.__extent_hdf_name__]

    @extent.setter
    def extent(self, val: GridEnvelope):
        self._attributes[self.__extent_hdf_name__] = val

    @property
    def __extent_type__(self) -> Type[GridEnvelope]:
        # FIXME is this a numpy array or s1xx list and what is its layout?
        return numpy.ndarray  # GridEnvelope

    @property
    def extent_dtype(self) -> Type[int]:
        return numpy.int64

    def extent_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.extent = self.__extent_type__([], self.extent_dtype)


# noinspection PyUnresolvedReferences
class GridOrigin:
    """ Mixin class for gridOriginLatitude/Longitude/Vertical.
    Used in Data Conding Formats 2,5,6

    S100 v5.0
    8-6.2.8 Grid cell structure S-100 utilizes the same view of grid cell structure as Section 8.2.2 of ISO 19123.
    The grid data in S-100 grid coverages are nominally situated exactly at the grid points defined by the grid coordinates.
    The grid points are therefore the “sample points.” Data values at a sample point represent measurements over a neighbourhood of the sample point.
    This neighbourhood is assumed to extend a half-cell in each dimension.
    The effect is that the sample space corresponding to each grid point is a cell centred at the grid point.

    """
    __grid_origin_longitude_hdf_name__ = "gridOriginLongitude"
    __grid_origin_latitude_hdf_name__ = "gridOriginLatitude"
    __grid_origin_vertical_hdf_name__ = "gridOriginVertical"

    @property
    def grid_origin_longitude(self) -> float:
        """The longitude of the grid origin. Unit: Arc Degrees"""
        return self._attributes[self.__grid_origin_longitude_hdf_name__]

    @grid_origin_longitude.setter
    def grid_origin_longitude(self, val: float):
        self._attributes[self.__grid_origin_longitude_hdf_name__] = val

    @property
    def __grid_origin_longitude_type__(self) -> Type[float]:
        return numpy.float64

    def grid_origin_longitude_create(self):
        """ Creates a blank, empty or zero value for grid_origin_longitude"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_origin_longitude = self.__grid_origin_longitude_type__()

    @property
    def grid_origin_latitude(self) -> float:
        """The latitude of the grid origin. Arc Degrees"""
        return self._attributes[self.__grid_origin_latitude_hdf_name__]

    @grid_origin_latitude.setter
    def grid_origin_latitude(self, val: float):
        self._attributes[self.__grid_origin_latitude_hdf_name__] = val

    @property
    def __grid_origin_latitude_type__(self) -> Type[float]:
        return numpy.float64

    def grid_origin_latitude_create(self):
        """ Creates a blank, empty or zero value for grid_origin_latitude"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_origin_latitude = self.__grid_origin_latitude_type__()

    @property
    def grid_origin_vertical(self) -> float:
        """The grid origin in the vertical dimension. Only for 3-D grids. Units specified by product specifications"""
        return self._attributes[self.__grid_origin_vertical_hdf_name__]

    @grid_origin_vertical.setter
    def grid_origin_vertical(self, val: float):
        self._attributes[self.__grid_origin_vertical_hdf_name__] = val

    @property
    def __grid_origin_vertical_type__(self) -> Type[float]:
        return float

    def grid_origin_vertical_create(self):
        """ Creates a blank, empty or zero value for grid_origin_vertical"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_origin_vertical = self.__grid_origin_vertical_type__()


# noinspection PyUnresolvedReferences
class GridSpacing:
    """Mixin class for gridSpacingLongitudinal/Latitudinal/Vertical.  Probably used with :class:`GridOrigin`
    in Data Conding Formats 2,5,6"""
    __grid_spacing_longitudinal_hdf_name__ = "gridSpacingLongitudinal"
    __grid_spacing_latitudinal_hdf_name__ = "gridSpacingLatitudinal"
    __grid_spacing_vertical_hdf_name__ = "gridSpacingVertical"

    @property
    def grid_spacing_longitudinal(self) -> float:
        """Cell size in the X/longitude dimension. This is the X/longitudinal component of the
        offset vector (8-7.1.4). Units: Arc Degrees"""
        return self._attributes[self.__grid_spacing_longitudinal_hdf_name__]

    @grid_spacing_longitudinal.setter
    def grid_spacing_longitudinal(self, val: float):
        self._attributes[self.__grid_spacing_longitudinal_hdf_name__] = val

    @property
    def __grid_spacing_longitudinal_type__(self) -> Type[float]:
        return numpy.float64

    def grid_spacing_longitudinal_create(self):
        """ Creates a blank, empty or zero value for grid_spacing_longitudinal"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_spacing_longitudinal = self.__grid_spacing_longitudinal_type__()

    @property
    def grid_spacing_latitudinal(self) -> float:
        """Cell size in the Y/latitude dimension. This is the Y/latitudinal component of the offset
        vector (8-7.1.4). Units: Arc Degrees"""
        return self._attributes[self.__grid_spacing_latitudinal_hdf_name__]

    @grid_spacing_latitudinal.setter
    def grid_spacing_latitudinal(self, val: float):
        self._attributes[self.__grid_spacing_latitudinal_hdf_name__] = val

    @property
    def __grid_spacing_latitudinal_type__(self) -> Type[float]:
        return numpy.float64

    def grid_spacing_latitudinal_create(self):
        """ Creates a blank, empty or zero value for grid_spacing_latitudinal"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_spacing_latitudinal = self.__grid_spacing_latitudinal_type__()

    @property
    def grid_spacing_vertical(self) -> float:
        """Cell size in the vertical dimension. Only for 3-D grids. Units specified by product specifications."""
        return self._attributes[self.__grid_spacing_vertical_hdf_name__]

    @grid_spacing_vertical.setter
    def grid_spacing_vertical(self, val: float):
        self._attributes[self.__grid_spacing_vertical_hdf_name__] = val

    @property
    def __grid_spacing_vertical_type__(self) -> Type[float]:
        return float

    def grid_spacing_vertical_create(self):
        """ Creates a blank, empty or zero value for grid_spacing_vertical"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_spacing_vertical = self.__grid_spacing_vertical_type__()


# noinspection PyUnresolvedReferences
class StartSequence:
    """Mixin class for startSequence.  Data Coding Formats 2,5,6 """
    __start_sequence_hdf_name__ = "startSequence"

    @property
    def start_sequence(self) -> str:
        """ Grid coordinates of the grid point to which the first in the sequence of values is to be
        assigned. The choice of a valid point for the start sequence is determined by the
        sequencing rule. Format: n, n… (comma-separated list of grid points, one per
        dimension – For example, 0,0)
        """
        return self._attributes[self.__start_sequence_hdf_name__]

    @start_sequence.setter
    def start_sequence(self, val: str):
        self._attributes[self.__start_sequence_hdf_name__] = val

    @property
    def __start_sequence_type__(self) -> Type[str]:
        return str

    def start_sequence_create(self):
        """ Creates a blank, empty or zero value for start_sequence"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.start_sequence = self.__start_sequence_type__()


class NumberOfStations:
    """ Mixin class for numberOfStations.  Data Coding Formats 4,8
    """
    __number_of_stations_hdf_name__ = "numberOfStations"  #: HDF5 naming

    @property
    def number_of_stations(self) -> int:
        return self._attributes[self.__number_of_stations_hdf_name__]

    @number_of_stations.setter
    def number_of_stations(self, val: int):
        self._attributes[self.__number_of_stations_hdf_name__] = val

    @property
    def __number_of_stations_type__(self) -> Type[int]:
        return int

    def number_of_stations_create(self):
        """ Creates a blank, empty or zero value for number_of_stations"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.number_of_stations = self.__number_of_stations_type__()


class NumberOfNodes:
    """ Mixin class for numberOfNodes.  Data Coding Formats 3, 7 """
    __number_of_nodes_hdf_name__ = "numberOfNodes"  #: HDF5 naming

    @property
    def number_of_nodes(self) -> int:
        return self._attributes[self.__number_of_nodes_hdf_name__]

    @number_of_nodes.setter
    def number_of_nodes(self, val: int):
        self._attributes[self.__number_of_nodes_hdf_name__] = val

    @property
    def __number_of_nodes_type__(self) -> Type[int]:
        return int

    def number_of_nodes_create(self):
        """ Creates a blank, empty or zero value for number_of_nodes"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.number_of_nodes = self.__number_of_nodes_type__()


class GeometryValuesDataset(S1xxGridsBase):
    __longitude_hdf_name__ = "longitude"
    __latitude_hdf_name__ = "latitude"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        """ The plain text name of the dataset"""
        return "geometryValues"

    @property
    def longitude(self) -> s1xx_sequence:
        """Get the data"""
        return self._attributes[self.__longitude_hdf_name__]

    @longitude.setter
    def longitude(self, val: s1xx_sequence):
        """Potential validation or other checks/changes to incoming data"""
        self._attributes[self.__longitude_hdf_name__] = val

    @property
    def __longitude_type__(self) -> s1xx_sequence:
        """S100 Datatype"""
        return numpy.ndarray

    @property
    def longitude_dtype(self) -> Type[float]:
        """S100 Datatype"""
        return numpy.float64

    def longitude_create(self):
        """ Creates a blank, empty or zero value for longitude"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.longitude = self.__longitude_type__([], self.longitude_dtype)

    @property
    def latitude(self) -> s1xx_sequence:
        """Get the data"""
        return self._attributes[self.__latitude_hdf_name__]

    @latitude.setter
    def latitude(self, val: s1xx_sequence):
        """Potential validation or other checks/changes to incoming data"""
        self._attributes[self.__latitude_hdf_name__] = val

    @property
    def __latitude_type__(self) -> s1xx_sequence:
        """S100 Datatype"""
        return numpy.ndarray

    @property
    def latitude_dtype(self) -> Type[float]:
        """S100 Datatype"""
        return numpy.float64

    def latitude_create(self):
        """ Creates a blank, empty or zero value for latitude"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.latitude = self.__latitude_type__([], self.latitude_dtype)

    def get_write_order(self):
        """Specify order of attributes for ordered dict"""
        return [self.__longitude_hdf_name__, self.__latitude_hdf_name__]

    def get_compound_dtype(self):
        return [self.longitude_dtype, self.latitude_dtype]


class PositioningGroup(S1xxObject):
    __geometry_values_hdf_name__ = "geometryValues"
    __triangles_hdf_name__ = "triangles"
    __adjacency_hdf_name__ = "adjacency"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        """ The plain text name of the group

        Returns:
            str: Name of the group
        """
        return "Positioning"

    @property
    def metadata_type(self) -> type:
        return GeometryValuesDataset

    @property
    def geometry_values(self) -> GeometryValuesDataset:
        """Get the data"""
        return self._attributes[self.__geometry_values_hdf_name__]

    @geometry_values.setter
    def geometry_values(self, val: GeometryValuesDataset):
        self._attributes[self.__geometry_values_hdf_name__] = val

    @property
    def __geometry_values_type__(self) -> Type[GeometryValuesDataset]:
        """S100 Datatype"""
        return GeometryValuesDataset

    def geometry_values_create(self):
        """ Creates a blank, empty or zero value for geometry_values"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.geometry_values = self.__geometry_values_type__()

    @property
    def triangles(self) -> s1xx_sequence:
        return self._attributes[self.__triangles_hdf_name__]

    @triangles.setter
    def triangles(self, val: s1xx_sequence):
        self._attributes[self.__triangles_hdf_name__] = val

    @property
    def __triangles_type__(self) -> Type[numpy.array]:
        return numpy.ndarray

    def triangles_create(self):
        """ Creates a blank, empty or zero value for triangles"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.triangles = numpy.array([0, 0, 0], dtype=numpy.int32)

    @property
    def adjacency(self) -> s1xx_sequence:
        return self._attributes[self.__adjacency_hdf_name__]

    @adjacency.setter
    def adjacency(self, val: s1xx_sequence):
        self._attributes[self.__adjacency_hdf_name__] = val

    @property
    def __adjacency_type__(self) -> Type[numpy.array]:
        return numpy.ndarray

    def adjacency_create(self):
        """ Creates a blank, empty or zero value for adjacency"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.adjacency = numpy.array([0, 0, 0], dtype=numpy.int32)


class Positions:
    __positioning_hdf_name__ = "Positioning"

    @property
    def positioning(self) -> S1xxObject:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.__positioning_hdf_name__]

    @positioning.setter
    def positioning(self, val: S1xxObject):
        self._attributes[self.__positioning_hdf_name__] = val

    @property
    def __positioning_type__(self):
        """Defines datatype"""
        return PositioningGroup

    def positioning_create(self):
        """ Creates a blank, empty or zero value for positioning"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.positioning = self.__positioning_type__()


class FeatureInstanceDCF2(StartSequence, GridSpacing, GridOrigin, FeatureInstanceBase):
    """ Data Coding Format 2 is the grid format from table 10c-12 in S100 spec.  Used in S102 for example.
    """

    __num_points_longitudinal_hdf_name__ = "numPointsLongitudinal"
    __num_points_latitudinal_hdf_name__ = "numPointsLatitudinal"
    __num_points_vertical_hdf_name__ = "numPointsVertical"

    @property
    def num_points_longitudinal(self) -> int:
        """Number of grid points in the X/longitude dimension. (iMax)"""
        return self._attributes[self.__num_points_longitudinal_hdf_name__]

    @num_points_longitudinal.setter
    def num_points_longitudinal(self, val: int):
        self._attributes[self.__num_points_longitudinal_hdf_name__] = val

    @property
    def __num_points_longitudinal_type__(self) -> Type[int]:
        return numpy.uint32

    def num_points_longitudinal_create(self):
        """ Creates a blank, empty or zero value for num_points_longitudinal"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.num_points_longitudinal = self.__num_points_longitudinal_type__()

    @property
    def num_points_latitudinal(self) -> int:
        """Number of grid points in the Y/latitude dimension. (jMax)"""
        return self._attributes[self.__num_points_latitudinal_hdf_name__]

    @num_points_latitudinal.setter
    def num_points_latitudinal(self, val: int):
        self._attributes[self.__num_points_latitudinal_hdf_name__] = val

    @property
    def __num_points_latitudinal_type__(self) -> Type[int]:
        return numpy.uint32

    def num_points_latitudinal_create(self):
        """ Creates a blank, empty or zero value for num_points_latitudinal"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.num_points_latitudinal = self.__num_points_latitudinal_type__()

    @property
    def num_points_vertical(self) -> int:
        """Number of grid points in the vertical dimension. (kMax)"""
        return self._attributes[self.__num_points_vertical_hdf_name__]

    @num_points_vertical.setter
    def num_points_vertical(self, val: int):
        self._attributes[self.__num_points_vertical_hdf_name__] = val

    @property
    def __num_points_vertical_type__(self) -> Type[int]:
        return int

    def num_points_vertical_create(self):
        """ Creates a blank, empty or zero value for num_points_vertical"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.num_points_vertical = self.__num_points_vertical_type__()


# @TODO add support for Positioning group (S100 v5.0 10c-9.10.1, 10c-9.10.2 and table 10c-16) for DataCodingFormats 1,3,4,7,8
#   Basically makes a 1D dataset array of the positions and attributes (like Z).
#   For TINs  DCF=7
#   Optionally has triangles dataset NumTriangles x 3 -- array of the 3 node indices for positions
#   Optionally has adjacency dataset NumTriangles x3 -- array of the indices or the triangles that share an edge with each triangle side
#   adjacency[i][0] = triangle adjacent to the edge specified by triangles[i][0] & triangles[i][1]
#   adjacency[i][1] = triangle adjacent to edge triangles[i][1] & triangles[i][2]
#   adjacency[i][2] = triangle adjacent to edge triangles[i][2] & triangles[i][0]
#   Elements for edges without adjacent triangles are filled with the value -1
#   Applicable only for dataCodingFormat = 7 (TIN), but optional even for TIN.


class FeatureInstanceDCF1(Positions, NumberOfStations, FeatureInstanceBase):
    """ Data Coding Format 1 is the Fixed Stations from table 10c-12 in S100 spec.
    """


class FeatureInstanceDCF3(Positions, NumberOfNodes, FeatureInstanceBase):
    """ Data Coding Format 3 is the Ungeorectified grid format from table 10c-12 in S100 spec.
    """


class FeatureInstanceDCF4(Positions, NumberOfStations, FeatureInstanceBase):
    """ Data Coding Format 4 is the Moving Platform format from table 10c-12 in S100 spec.
    """


class FeatureInstanceDCF5(StartSequence, NumberOfNodes, GridSpacing, GridOrigin, FeatureInstanceBase):
    """ Data Coding Format 5 is the Irregular grid format from table 10c-12 in S100 spec.
    """


""" Data Coding Format 6 is the Variable Cell Size grid format from table 10c-12 in S100 spec. """
FeatureInstanceDCF6 = FeatureInstanceDCF5


class FeatureInstanceDCF7(Positions, NumberOfNodes, FeatureInstanceBase):
    """ Data Coding Format 7 is the Triangulated Irregular Network (TIN) format from table 10c-12 in S100 spec.
    """
    __number_of_triangles_hdf_name__ = "numberOfTriangles"  #: HDF5 naming

    @property
    def number_of_triangles(self) -> int:
        return self._attributes[self.__number_of_triangles_hdf_name__]

    @number_of_triangles.setter
    def number_of_triangles(self, val: int):
        self._attributes[self.__number_of_triangles_hdf_name__] = val

    @property
    def __number_of_triangles_type__(self) -> Type[int]:
        return int

    def number_of_triangles_create(self):
        """ Creates a blank, empty or zero value for number_of_triangles"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.number_of_triangles = self.__number_of_triangles_type__()


class FeatureInstanceDCF8(Positions, NumberOfStations, FeatureInstanceBase):
    """ Fixed stations - stationwise from S100 v5.0 Table 10c-12"""


class FeatureInstanceDCF9(FeatureInstanceDCF2):
    pass


class FeatureInformation(S1xxObject):
    """  In S100, table 10c-8.
    In S102, 10.2.1 and table 10.2 and Table 10.1 of v2.0.0

    FeatureInformation (GroupF) is used to describe the FeatureInstances at the root of the file,
    ex: BathymetryCoverage in S102 or SurfaceCurrent in S111.

    The actual data is stored in the top level FeatureInstance element while basic metadata is stored in this FeatureInformation element.

    Note that the data contained in this class are stored in the HDF5 as strings
    but are translated by s100py to appropriate python types (int, float etc)

    The “code” and “datatype” components encode the rangeType attribute of the coverage features in Part 8.

    “lower”, “upper”, and “closure” encode any constraints on attribute values as encoded in the
    feature catalogue (see “S100_FC_SimpleAttribute>constraints” in Part 5 and
    S100_NumericRange in Part 1)
    """
    __code_hdf_name__ = "code"
    __name_hdf_name__ = "name"
    __unit_of_measure_hdf_name__ = "uom.name"
    __fill_value_hdf_name__ = "fillValue"
    __datatype_hdf_name__ = "datatype"
    __lower_hdf_name__ = "lower"
    __upper_hdf_name__ = "upper"
    __closure_hdf_name__ = "closure"

    @property
    def __version__(self) -> int:
        return 1

    def get_write_order(self):
        return [self.__code_hdf_name__,
                self.__name_hdf_name__,
                self.__unit_of_measure_hdf_name__,
                self.__fill_value_hdf_name__,
                self.__datatype_hdf_name__,
                self.__lower_hdf_name__,
                self.__upper_hdf_name__,
                self.__closure_hdf_name__]

    @property
    def code(self) -> str:
        """ Camel case code of attribute as in feature catalogue.
        The “code” and “datatype” components encode the rangeType attribute of the coverage features in Part 8.
        """
        return self._attributes[self.__code_hdf_name__]

    @code.setter
    def code(self, val: str):
        self._attributes[self.__code_hdf_name__] = val

    @property
    def __code_type__(self):
        return str

    def code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.code = self.__code_type__()

    @property
    def name(self) -> str:
        """ Long name as in feature catalogue
        """
        return self._attributes[self.__name_hdf_name__]

    @name.setter
    def name(self, val: str):
        self._attributes[self.__name_hdf_name__] = val

    @property
    def __name_type__(self):
        return str

    def name_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.name = self.__name_type__()

    @property
    def unit_of_measure(self) -> str:
        """ Units of measurement for the dataset.  (uom>name from S-100 feature catalogue)
        """
        return self._attributes[self.__unit_of_measure_hdf_name__]

    @unit_of_measure.setter
    def unit_of_measure(self, val: str):
        self._attributes[self.__unit_of_measure_hdf_name__] = val

    @property
    def __unit_of_measure_type__(self):
        return str

    def unit_of_measure_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.unit_of_measure = self.__unit_of_measure_type__()

    def _python_datatype(self):
        """ Determine what kind of python type best fits the HDF5 type held in the S100 'datatype' attribute.
        The datatype attribute can be set by each S100+ spec individually (S102 and S111 are floats).
        For unhandled types (like H5T.OPAQUE) returns str.

        Returns
        -------
        A Python type: int, str, float are supported currently
        """
        try:
            self.datatype
        except:  # datatype not set yet, so save as a string by default
            val = str
        else:
            if self.datatype in _H5T_Types[int]:
                val = int
            elif self.datatype in _H5T_Types[float]:
                val = float
            else:
                val = str
        return val

    def _convert_from_string_based_on_datatype(self, str_val):
        use_datatype = self._python_datatype()
        if use_datatype is int:
            try:
                val = int(str_val)
            except ValueError as e:
                if str_val == "":
                    val = None
                else:
                    raise e
        elif use_datatype is float:
            try:
                val = float(str_val)
            except ValueError as e:
                if str_val == "":
                    val = None
                else:
                    raise e
        else:
            val = str_val
        return val

    def _convert_to_string_based_on_datatype(self, val):
        use_datatype = self._python_datatype()
        if use_datatype in (int, float) and not isinstance(val, str):
            if val is None or val == "":
                str_val = ""
            elif use_datatype is int:
                str_val = str(int(val))  # this extra conversion gives python a chance to convert scientific notation to standard
            elif use_datatype is float:
                str_val = str(float(val))  # this extra conversion gives python a chance to convert scientific notation to standard
                if str_val[-2:] == ".0":  # remove trailing '.0' so a 12000.0 becomes 12000
                    str_val = str_val[:-2]
        else:
            str_val = str(val)
        return str_val

    @property
    def fill_value(self) -> Union[float, int, str]:
        """ Value denoting missing data.  Fill value (integer or float value, string representation)
        """
        return self._convert_from_string_based_on_datatype(self._attributes[self.__fill_value_hdf_name__])

    @fill_value.setter
    def fill_value(self, val: Union[float, int, str]):
        self._attributes[self.__fill_value_hdf_name__] = self._convert_to_string_based_on_datatype(val)

    @property
    def __fill_value_type__(self):
        return self._python_datatype()

    def fill_value_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.fill_value = self.__fill_value_type__()

    @property
    def datatype(self) -> str:
        return self._attributes[self.__datatype_hdf_name__]

    @datatype.setter
    def datatype(self, val: Union[str, int]):
        """ The “code” and “datatype” components encode the rangeType attribute of the coverage features in Part 8

        Parameters
        ----------
        val
            Either the string name (ex: 'H5T_INTEGER') of the datatype or the h5py constant (ex: h5py.h5t.INTEGER)
        """
        if isinstance(val, int):
            val = H5T_CLASS_T[val]
        self._attributes[self.__datatype_hdf_name__] = val

    @property
    def __datatype_type__(self):
        return str

    def datatype_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.datatype = self.__datatype_type__()

    @property
    def lower(self) -> Union[float, int, str]:
        """ Lower bound on value of attribute """
        return self._convert_from_string_based_on_datatype(self._attributes[self.__lower_hdf_name__])

    @lower.setter
    def lower(self, val: Union[float, int, str]):
        self._attributes[self.__lower_hdf_name__] = self._convert_to_string_based_on_datatype(val)

    @property
    def __lower_type__(self):
        return self._python_datatype()

    def lower_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.lower = self.__lower_type__()

    @property
    def upper(self) -> Union[float, int, str]:
        """ Upper bound on attribute value """
        return self._convert_from_string_based_on_datatype(self._attributes[self.__upper_hdf_name__])

    @upper.setter
    def upper(self, val: Union[float, int, str]):
        self._attributes[self.__upper_hdf_name__] = self._convert_to_string_based_on_datatype(val)

    @property
    def __upper_type__(self):
        return self._python_datatype()

    def upper_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.upper = self.__upper_type__()

    @property
    def closure(self) -> str:
        """ type of closure from S100 Table 1-3 — Interval Types::

        openInterval    The open interval             (a,b)   a < x < b
        geLtInterval    The right half-open interval  [a,b)   a ≤ x < b
        gtLeInterval    The left half-open interval   (a,b]   a < x ≤ b
        closedInterval  The closed interval           [a,b]   a≤ x ≤ b
        gtSemiInterval  The left half-open ray        (a,∞)   a < x
        geSemiInterval  The left closed ray           [a,∞)   a ≤ x
        ltSemiInterval  The right half-open ray       (-∞,a)  x < a
        leSemiInterval  The right closed ray          (-∞,a]  x ≤ a
        """
        return self._attributes[self.__closure_hdf_name__]

    @closure.setter
    def closure(self, val: str):
        self._attributes[self.__closure_hdf_name__] = val

    @property
    def __closure_type__(self):
        return str

    def closure_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.closure = self.__closure_type__()


# FIXME @TODO - this FeatureInformationDataset class in v4.0 writes the compound array to the HDF5 file.
#  Do we need a new Compound type that does this instead (we didn't implement any other compounds in v4.0)
#  The FeatureContainer's new feature_attribute_table needs to be stored similarly.  Also be able to store an array of compounds

class FeatureInformationDataset(S1xxDatasetBase, ABC):
    """ This class comes from S100 -- 10c-9.5 Feature information group.
    This class serves to keep a list of FeatureInformation objects which will be turned into a compound array
    of strings in the HDF5 file.

    The metadata_name property must be overridden.
    The metadata_type will likely be overridden with a specific subclass for the s100+ spec
    """

    @property
    def metadata_type(self) -> Type[FeatureInformation]:
        return FeatureInformation


class CommonPointRule:
    __common_point_rule_hdf_name__ = "commonPointRule"

    @property
    def common_point_rule(self) -> COMMON_POINT_RULE:
        """ The procedure used for evaluating the coverage at a position that falls on the
        boundary or in an area of overlap between geometric objects
        Values from CV_CommonPointRule (Table 10c-19).

        see :data:`~COMMON_POINT_RULE`
        """
        return self._attributes[self.__common_point_rule_hdf_name__]

    @common_point_rule.setter
    def common_point_rule(self, val: Union[int, str, COMMON_POINT_RULE]):
        self.set_enum_attribute(val, self.__common_point_rule_hdf_name__, self.__common_point_rule_type__)

    @property
    def __common_point_rule_type__(self) -> Type[Enum]:
        return COMMON_POINT_RULE

    def common_point_rule_create(self):
        """ Creates a blank, empty or zero value for common_point_rule"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.common_point_rule = self.__common_point_rule_type__["average"]


class FeatureContainer(CommonPointRule, S1xxObject):
    """ Feature Container from S100 v5.0 Table 10c-9
    The feature attribute table was added, primarily for S102 to supply extended information about the origin of the data.
    The Data Coding Format must be 9 to use the feature attribute table.
    The Columns will be "id" followed by an arbitrary number of column names.
    Each "id" is an unsigned integer that matches the values in the values record of a data values group.

    Additional columns
    The number of additional columns depends on the information to be transported.  The number must always be greater than zero (>0).
    Table 10c-8 defines that the name of the "feature container group" must always correspond to the
    camelCase name of an object of the type "Feature Type" of the IHO Geospatial Information (GI) Registry, Data Dictionary Register.
    The names of the other columns are based on the attributes of the selected feature type of the Data Dictionary Register.
    Due to the mapping form of the attributes in the HDF5 format, only simple data types (see clause 10c-7) can be used.
    For the mapping of complex data types, a corresponding notation by means of a dot is necessary.
    The respective levels of the complex data type are separated by a dot up to the simple data type;
    for example, "surveyDateRange.dateEnd" & "surveyDateRange.dateStart".

    From v4.0:  This class comes from S100 in Table 10c-9 – Structure of feature container groups and
    Table 10c-10 – Attributes of feature container groups
    """
    __axis_names_hdf_name__ = "axisNames"
    __coordinate_size_hdf_name__ = "coordinateSize"  #: HDF5 naming
    __data_coding_format_hdf_name__ = "dataCodingFormat"
    __dimension_hdf_name__ = "dimension"
    __horizontal_position_uncertainty_hdf_name__ = "horizontalPositionUncertainty"
    __vertical_uncertainty_hdf_name__ = "verticalUncertainty"
    __time_uncertainty_hdf_name__ = "timeUncertainty"
    __num_instances_hdf_name__ = "numInstances"
    __feature_attribute_table_hdf_name__ = "featureAttributeTable"  #: HDF5 naming

    def __init__(self, *args, **opts):
        super().__init__(*args, **opts)
        self.data_coding_format_create()  # this is defined by the subclass and is constant, so we will automatically set it here

    @property
    def __version__(self) -> int:
        return 1

    @property
    def axis_names(self) -> s1xx_sequence:
        """sequence of character strings

        S100 Spec: Array (1-D): 0..D-1 where D is the value of the dimension attribute
        Axes should be in major-minor order; that is, if storage is to be in row-major order the
        X/longitude axis should be first.
        """
        return self._attributes[self.__axis_names_hdf_name__]

    @axis_names.setter
    def axis_names(self, val: s1xx_sequence):
        self._attributes[self.__axis_names_hdf_name__] = val

    @property
    def __axis_names_type__(self) -> Type[numpy.array]:
        return numpy.array

    def axis_names_create(self):
        """ The attribute axisNames has the value class Sequence<CharacterString> that shall be used to assign names to the grid axis.
        The grid axis names shall be "Latitude" and "Longitude" for unprojected data sets or “Northing” and “Easting” in a projected space.

        Returns
        -------

        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.axis_names = self.__axis_names_type__(["", ""], dtype=h5py_string_dtype)

    @property
    def coordinate_size(self) -> s1xx_sequence:
        return self._attributes[self.__coordinate_size_hdf_name__]

    @coordinate_size.setter
    def coordinate_size(self, val: s1xx_sequence):
        self._attributes[self.__coordinate_size_hdf_name__] = val

    @property
    def __coordinate_size_type__(self) -> Type[s1xx_sequence]:
        return numpy.ndarray

    @property
    def coordinate_size_dtype(self) -> Type[int]:
        return numpy.int64

    def coordinate_size_create(self):
        """ Creates a blank, empty or zero value for coordinate_size"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.coordinate_size = self.__coordinate_size_type__([], self.coordinate_size_dtype)

    @property
    def data_coding_format(self) -> DATA_CODING_FORMAT:
        """ Indication of the type of coverage in instances of this feature. Used to read the
        data (see Table 10c-4) or :data:`~DATA_CODING_FORMAT`
        """
        return self._attributes[self.__data_coding_format_hdf_name__]

    @data_coding_format.setter
    def data_coding_format(self, val: int):
        self.set_enum_attribute(val, self.__data_coding_format_hdf_name__, self.__data_coding_format_type__)

    @property
    def __data_coding_format_type__(self) -> DATA_CODING_FORMAT:
        return DATA_CODING_FORMAT

    @abstractmethod
    def data_coding_format_create(self):
        """ Creates a blank, empty or zero value for data_coding_format"""
        raise NotImplementedError("each s100+ spec implementation must override this  data coding format with the correct default")

    @property
    def dimension(self) -> int:
        """ The dimension of the feature instances
        This is the number of coordinate axes, not the rank of the HDF5 arrays storing
        coordinates or values. For example, a fixed stations dataset with positions in
        latitude and longitude will have dimension=2
        """
        return self._attributes[self.__dimension_hdf_name__]

    @dimension.setter
    def dimension(self, val: int):
        self._attributes[self.__dimension_hdf_name__] = val

    @property
    def __dimension_type__(self) -> Type[int]:
        return numpy.uint8

    def dimension_create(self):
        """ Creates a blank, empty or zero value for dimension"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.dimension = self.__dimension_type__()

    @property
    def horizontal_position_uncertainty(self) -> float:
        """ The uncertainty in horizontal coordinates.
        For example, -1.0 (unknown/inapplicable) or positive value (m)
        """
        return self._attributes[self.__horizontal_position_uncertainty_hdf_name__]

    @horizontal_position_uncertainty.setter
    def horizontal_position_uncertainty(self, val: float):
        self._attributes[self.__horizontal_position_uncertainty_hdf_name__] = val

    @property
    def __horizontal_position_uncertainty_type__(self) -> Type[float]:
        return numpy.float32

    def horizontal_position_uncertainty_create(self):
        """ Creates a blank, empty or zero value for horizontal_position_uncertainty"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_position_uncertainty = self.__horizontal_position_uncertainty_type__()

    @property
    def vertical_uncertainty(self) -> float:
        """ The uncertainty in vertical coordinate(s).
        For example, -1.0 (unknown/inapplicable) or positive value (m)
        """
        return self._attributes[self.__vertical_uncertainty_hdf_name__]

    @vertical_uncertainty.setter
    def vertical_uncertainty(self, val: float):
        self._attributes[self.__vertical_uncertainty_hdf_name__] = val

    @property
    def __vertical_uncertainty_type__(self) -> Type[float]:
        return numpy.float32

    def vertical_uncertainty_create(self):
        """ Creates a blank, empty or zero value for vertical_uncertainty"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_uncertainty = self.__vertical_uncertainty_type__()

    @property
    def time_uncertainty(self) -> float:
        """ Uncertainty in time values.
        For example, -1.0 (unknown/inapplicable) or positive value (s)

        Only for time series data
        """
        return self._attributes[self.__time_uncertainty_hdf_name__]

    @time_uncertainty.setter
    def time_uncertainty(self, val: float):
        self._attributes[self.__time_uncertainty_hdf_name__] = val

    @property
    def __time_uncertainty_type__(self) -> Type[float]:
        return numpy.float32

    def time_uncertainty_create(self):
        """ Creates a blank, empty or zero value for time_uncertainty"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.time_uncertainty = self.__time_uncertainty_type__()

    @property
    def num_instances(self) -> int:
        """ Number of instances of the feature
        (Records in the same time series or moving platform sequence are counted as a
        single instance, not as separate instances)
        """
        return self._attributes[self.__num_instances_hdf_name__]

    @num_instances.setter
    def num_instances(self, val: int):
        self._attributes[self.__num_instances_hdf_name__] = val

    @property
    def __num_instances_type__(self) -> Type[int]:
        return numpy.uint32

    def num_instances_create(self):
        """ Creates a blank, empty or zero value for num_instances"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.num_instances = self.__num_instances_type__()

    @property
    def feature_attribute_table(self) -> s1xx_sequence:
        return self._attributes[self.__feature_attribute_table_hdf_name__]

    @feature_attribute_table.setter
    def feature_attribute_table(self, val: s1xx_sequence):
        self._attributes[self.__feature_attribute_table_hdf_name__] = val

    @property
    def __feature_attribute_table_type__(self) -> Type[numpy.array]:
        return numpy.array

    def feature_attribute_table_create(self):
        """ Creates a blank, empty or zero value for feature_attribute_table"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        # FIXME @TODO -- this has to be an arbitrary number of strings, is this limiting it to two?
        self.feature_attribute_table = self.__feature_attribute_table_type__(["", ""], dtype=h5py_string_dtype)


class SequencingRule:
    """ Mixin class for Sequencing Rule.  At least used in Data Coding Format 2,5,6
    """
    __sequencing_rule_type_hdf_name__ = "sequencingRule.type"
    __sequencing_rule_scan_direction_hdf_name__ = "sequencingRule.scanDirection"

    @property
    def sequencing_rule_type(self) -> SEQUENCING_RULE_TYPE:
        # @todo -- clean up formatting
        """ table 10c-20 of S100

        Enumeration CV_SequenceType Codes that identify the method of ordering grid
        points or value records
        ISO 19123 CV_SequenceType
        Literal linear Sequencing is consecutive along grid lines,
        starting with the first grid axis listed in
        scanDirection

        1 For example, for 2-D
        grids with scan
        direction=(x,y), scanning
        will be in row-major order
        Literal boustrophedonic Variant of linear sequencing in which the
        direction of the scan is reversed on alternating
        grid lines. For grids of dimension > 2, it is also
        reversed on alternating planes

        2
        Literal CantorDiagonal Sequencing in alternating directions along
        parallel diagonals of the grid. For dimension > 2,
        it is repeated in successive planes

        3
        Literal spiral Sequencing in spiral order 4
        S-100 Edition 4.0.0 December 2018
        40 Part 10c – HDF5 Data Format
        Literal Morton Sequencing along a Morton curve 5
        Literal Hilbert Sequencing along a Hilbert curve 6
        Morton

        Returns
        -------

        """
        return self._attributes[self.__sequencing_rule_type_hdf_name__]

    @sequencing_rule_type.setter
    def sequencing_rule_type(self, val: Union[int, str, SEQUENCING_RULE_TYPE]):
        self.set_enum_attribute(val, self.__sequencing_rule_type_hdf_name__, self.__sequencing_rule_type_type__)

    @property
    def __sequencing_rule_type_type__(self) -> Type[Enum]:
        return SEQUENCING_RULE_TYPE

    def sequencing_rule_type_create(self):
        """ Creates a blank, empty or zero value for sequencing_rule_type"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.sequencing_rule_type = self.__sequencing_rule_type_type__["linear"]

    @property
    def sequencing_rule_scan_direction(self) -> str:
        return self._attributes[self.__sequencing_rule_scan_direction_hdf_name__]

    @sequencing_rule_scan_direction.setter
    def sequencing_rule_scan_direction(self, val: str):
        self._attributes[self.__sequencing_rule_scan_direction_hdf_name__] = val

    @property
    def __sequencing_rule_scan_direction_type__(self) -> Type[str]:
        return str

    def sequencing_rule_scan_direction_create(self):
        """ Creates a blank, empty or zero value for sequencing_rule_scan_direction"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.sequencing_rule_scan_direction = self.__sequencing_rule_scan_direction_type__()


class InterpolationType:
    """ Mixin class for Interpolation Type.  At least used in Data Coding Format 2,3,4,5,6,7
    """
    __interpolation_type_hdf_name__ = "interpolationType"

    @property
    def interpolation_type(self) -> Type[int]:
        """ S100 table 10c-21

        Returns
        -------

        """
        return self._attributes[self.__interpolation_type_hdf_name__]

    @interpolation_type.setter
    def interpolation_type(self, val: Union[int, str, INTERPOLATION_TYPE]):
        self.set_enum_attribute(val, self.__interpolation_type_hdf_name__, self.__interpolation_type_type__)

    @property
    def __interpolation_type_type__(self) -> Type[Enum]:
        return INTERPOLATION_TYPE

    def interpolation_type_create(self):
        """ Creates a blank, empty or zero value for interpolation_type"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.interpolation_type = self.__interpolation_type_type__['nearestneighbor']


class DataOffset:
    """ Mixin class for Data Offset from S100 v5.0 Table 10c-10
    Added dataOffsetCode and dataOffsetVector to the featureContainer for DataCodingFormats 2,5,6,9.
    """
    __data_offset_code_hdf_name__ = "dataOffsetCode"  #: HDF5 naming
    __data_offset_vector_hdf_name__ = "dataOffsetVector"  #: HDF5 naming

    @property
    def data_offset_code(self) -> DATA_OFFSET_CODE:
        return self._attributes[self.__data_offset_code_hdf_name__]

    @data_offset_code.setter
    def data_offset_code(self, val: Union[int, str, DATA_OFFSET_CODE]):
        self.set_enum_attribute(val, self.__data_offset_code_hdf_name__, self.__data_offset_code_type__)

    @property
    def __data_offset_code_type__(self) -> Type[DATA_OFFSET_CODE]:
        return DATA_OFFSET_CODE

    def data_offset_code_create(self):
        """ Creates a blank, empty or zero value for data_offset_code
        """
        # Make the enum into a list and take the first value
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_offset_code = list(self.__data_offset_code_type__)[0]

    # FIXME - this is an array not just a float
    @property
    def data_offset_vector(self) -> s1xx_sequence:
        return self._attributes[self.__data_offset_vector_hdf_name__]

    @data_offset_vector.setter
    def data_offset_vector(self, val: s1xx_sequence):
        self._attributes[self.__data_offset_vector_hdf_name__] = val

    @property
    def __data_offset_vector_type__(self) -> Type[numpy.array]:
        return numpy.ndarray

    def data_offset_vector_create(self):
        """ Creates a blank, empty or zero value for data_offset_vector"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_offset_vector =numpy.array([0, 0], dtype=float)


class FeatureContainerDCF1(FeatureContainer):
    """ Container for Data Coding Format 1 """

    def data_coding_format_create(self):
        self.data_coding_format = self.__data_coding_format_type__(1)


class FeatureContainerDCF2(DataOffset, SequencingRule, InterpolationType, FeatureContainer):
    """ Container for Data Coding Format 2 """

    def data_coding_format_create(self):
        self.data_coding_format = self.__data_coding_format_type__(2)


class FeatureContainerDCF3(InterpolationType, FeatureContainer):
    """ Container for Data Coding Format 3 """

    def data_coding_format_create(self):
        self.data_coding_format = self.__data_coding_format_type__(3)


class FeatureContainerDCF4(FeatureContainer):
    """ Container for Data Coding Format 4 """

    def data_coding_format_create(self):
        self.data_coding_format = self.__data_coding_format_type__(4)


class FeatureContainerDCF5(DataOffset, SequencingRule, InterpolationType, FeatureContainer):
    """ Container for Data Coding Format 5 """

    def data_coding_format_create(self):
        self.data_coding_format = self.__data_coding_format_type__(5)


class FeatureContainerDCF6(DataOffset, SequencingRule, InterpolationType, FeatureContainer):
    """ Container for Data Coding Format 6 """

    def data_coding_format_create(self):
        self.data_coding_format = self.__data_coding_format_type__(6)


class FeatureContainerDCF7(InterpolationType, FeatureContainer):
    """ Container for Data Coding Format 7 """

    def data_coding_format_create(self):
        self.data_coding_format = self.__data_coding_format_type__(7)


class FeatureContainerDCF8(FeatureContainer):
    """ Container for Data Coding Format 8 """

    def data_coding_format_create(self):
        self.data_coding_format = self.__data_coding_format_type__(8)

class FeatureContainerDCF9(DataOffset, SequencingRule, FeatureContainer):
    """ Feature Oriented Regular Grid from S100 v5.0 Table 10c-12
    Container for Data Coding Format 9 """

    def data_coding_format_create(self):
        self.data_coding_format = self.__data_coding_format_type__(9)


class GroupFBase(S1xxObject):
    """ From S100 Table 10c-8 – Components of feature information group

    There will also be a :class:`FeatureInformationDataset` holding a list of :class:`FeatureInformation`
    which will be defined by the subclasses of this base class.
    """
    __feature_code_hdf_name__ = "featureCode"

    @property
    def __feature_code_type__(self):
        return numpy.array

    @abstractmethod
    def feature_code_create(self):
        raise NotImplementedError("must overload feature_code_create")

    @property
    def feature_code(self) -> s1xx_sequence:
        """Array (1-d): i=0, F-1.
        Values = codes of feature classes
        (F is the number of feature classes in the application schema.)
        """
        return self._attributes[self.__feature_code_hdf_name__]

    @feature_code.setter
    def feature_code(self, val: s1xx_sequence):
        self._attributes[self.__feature_code_hdf_name__] = val


class S100Root(GeographicBoundingBox):
    """ From table 10c-6 in S100 v4.0 spec.

    The root of the S100 v5.0 schema.  There are restrictions on many of the CRS attributes based on other attributes.
    For example, the horizontal_cs has different value options based on horizontal_crs.

    """
    __feature_information_hdf_name__ = "Group_F"
    __epoch_hdf_name__ = "epoch"
    __geographic_identifier_hdf_name__ = "geographicIdentifier"
    __meta_features_hdf_name__ = "metaFeatures"
    __metadata_hdf_name__ = "metadata"
    __product_specification_hdf_name__ = "productSpecification"
    __issue_time_hdf_name__ = "issueTime"
    __issue_date_hdf_name__ = "issueDate"
    # revised/added in S100 v5.0
    __horizontal_crs_hdf_name__ = "horizontalCRS"  #: HDF5 naming
    __name_of_horizontal_crs_hdf_name__ = "nameOfHorizontalCRS"  #: HDF5 naming
    __type_of_horizontal_crs_hdf_name__ = "typeOfHorizontalCRS"  #: HDF5 naming
    __horizontal_cs_hdf_name__ = "horizontalCS"  #: HDF5 naming
    __horizontal_datum_hdf_name__ = "horizontalDatum"  #: HDF5 naming
    __name_of_horizontal_datum_hdf_name__ = "nameOfHorizontalDatum"  #: HDF5 naming
    __prime_meridian_hdf_name__ = "primeMeridian"  #: HDF5 naming
    __spheriod_hdf_name__ = "spheriod"  #: HDF5 naming
    __projection_method_hdf_name__ = "projectionMethod"  #: HDF5 naming
    __projection_parameter_1_hdf_name__ = "projectionParameter1"  #: HDF5 naming
    __projection_parameter_2_hdf_name__ = "projectionParameter2"  #: HDF5 naming
    __projection_parameter_3_hdf_name__ = "projectionParameter3"  #: HDF5 naming
    __projection_parameter_4_hdf_name__ = "projectionParameter4"  #: HDF5 naming
    __projection_parameter_5_hdf_name__ = "projectionParameter5"  #: HDF5 naming
    __false_easting_hdf_name__ = "falseEasting"  #: HDF5 naming
    __false_northing_hdf_name__ = "falseNorthing"  #: HDF5 naming
    __vertical_cs_hdf_name__ = "verticalCS"  #: HDF5 naming
    __vertical_coordinate_base_hdf_name__ = "verticalCoordinateBase"  #: HDF5 naming
    __vertical_datum_reference_hdf_name__ = "verticalDatumReference"  #: HDF5 naming
    __vertical_datum_hdf_name__ = "verticalDatum"  #: HDF5 naming

    # Removed in v5.0
    # __horizontal_datum_reference_hdf_name__ = "horizontalDatumReference"
    # __horizontal_datum_value_hdf_name__ = "horizontalDatumValue"

    @property
    def __version__(self) -> int:
        return 1

    def write(self, group_object):
        super().write(group_object)
        feat_info = None
        for property_name in self.get_standard_properties():
            if is_sub_class(self.__getattribute__("__" + property_name + "_type__"), GroupFBase):
                feat_info = self.__getattribute__(property_name)

    @property
    def feature_information(self) -> GroupFBase:
        return self._attributes[self.__feature_information_hdf_name__]

    @feature_information.setter
    def feature_information(self, val: GroupFBase):
        self._attributes[self.__feature_information_hdf_name__] = val

    @property
    def __feature_information_type__(self):
        print("Supposed to override feature_information before using")
        return float
        raise NotImplementedError()

    def feature_information_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_information = self.__feature_information_type__()

    @property
    def product_specification(self) -> str:
        return self._attributes[self.__product_specification_hdf_name__]

    @product_specification.setter
    def product_specification(self, val: str):
        self._attributes[self.__product_specification_hdf_name__] = val

    @property
    def __product_specification_type__(self) -> Type[str]:
        return str

    def product_specification_create(self):
        """ Creates a blank, empty or zero value for product_specification"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.product_specification = self.__product_specification_type__()

    @property
    def issue_time(self) -> datetime.time:
        return self._attributes[self.__issue_time_hdf_name__]

    @issue_time.setter
    def issue_time(self, val: Union[datetime.time, datetime.datetime, str]):
        self.set_datetime_attribute(val, self.__issue_time_hdf_name__, self.__issue_time_type__)

    @property
    def __issue_time_type__(self) -> Type[datetime.time]:
        return datetime.time

    def issue_time_create(self):
        """ Creates a blank, empty or zero value for issue_time"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.issue_time = self.__issue_time_type__(0, 0)  # midnight

    # example of how to write a custom format to the HDF5 file -- make sure it will pass the iso parser to be able to read it back from the file
    # @property
    # def __issue_date_repr__(self) -> str:
    #     return self._attributes[self.__issue_date_hdf_name__].strftime("%Y%m%dBarry")

    @property
    def issue_date(self) -> datetime.date:
        return self._attributes[self.__issue_date_hdf_name__]

    @issue_date.setter
    def issue_date(self, val: Union[datetime.date, datetime.datetime, str]):
        self.set_datetime_attribute(val, self.__issue_date_hdf_name__, self.__issue_date_type__)

    @property
    def __issue_date_type__(self) -> Type[datetime.date]:
        return datetime.date

    def issue_date_create(self):
        """ Creates a blank, empty or zero value for issue_date"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.issue_date = self.__issue_date_type__(1970, 1, 1)

    @property
    def epoch(self) -> str:
        return self._attributes[self.__epoch_hdf_name__]

    @epoch.setter
    def epoch(self, val: str):
        self._attributes[self.__epoch_hdf_name__] = val

    @property
    def __epoch_type__(self) -> Type[str]:
        return str

    def epoch_create(self):
        """ Creates a blank, empty or zero value for epoch"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.epoch = self.__epoch_type__()

    @property
    def geographic_identifier(self) -> str:
        return self._attributes[self.__geographic_identifier_hdf_name__]

    @geographic_identifier.setter
    def geographic_identifier(self, val: str):
        self._attributes[self.__geographic_identifier_hdf_name__] = val

    @property
    def __geographic_identifier_type__(self) -> Type[str]:
        return str

    def geographic_identifier_create(self):
        """ Creates a blank, empty or zero value for geographic_identifier"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.geographic_identifier = self.__geographic_identifier_type__()

    @property
    def metadata(self) -> str:
        return self._attributes[self.__metadata_hdf_name__]

    @metadata.setter
    def metadata(self, val: str):
        self._attributes[self.__metadata_hdf_name__] = val

    @property
    def __metadata_type__(self) -> Type[str]:
        return str

    def metadata_create(self):
        """ Creates a blank, empty or zero value for metadata"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.metadata = self.__metadata_type__()

    @property
    def meta_features(self) -> str:
        return self._attributes[self.__meta_features_hdf_name__]

    @meta_features.setter
    def meta_features(self, val: str):
        self._attributes[self.__meta_features_hdf_name__] = val

    @property
    def __meta_features_type__(self) -> Type[str]:
        return str

    def meta_features_create(self):
        """ Creates a blank, empty or zero value for meta_features"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.meta_features = self.__meta_features_type__()

    @property
    def name_of_horizontal_crs(self) -> str:
        return self._attributes[self.__name_of_horizontal_crs_hdf_name__]

    @name_of_horizontal_crs.setter
    def name_of_horizontal_crs(self, val: str):
        self._attributes[self.__name_of_horizontal_crs_hdf_name__] = val

    @property
    def __name_of_horizontal_crs_type__(self) -> Type[str]:
        return str

    def name_of_horizontal_crs_create(self):
        """ Creates a blank, empty or zero value for name_of_horizontal_crs"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.name_of_horizontal_crs = self.__name_of_horizontal_crs_type__()

    @property
    def type_of_horizontal_crs(self) -> TYPE_OF_HORIZONTAL_CRS:
        return self._attributes[self.__type_of_horizontal_crs_hdf_name__]

    @type_of_horizontal_crs.setter
    def type_of_horizontal_crs(self, val: Union[int, str, TYPE_OF_HORIZONTAL_CRS]):
        self.set_enum_attribute(val, self.__type_of_horizontal_crs_hdf_name__, self.__type_of_horizontal_crs_type__)

    @property
    def __type_of_horizontal_crs_type__(self) -> Type[TYPE_OF_HORIZONTAL_CRS]:
        return TYPE_OF_HORIZONTAL_CRS

    def type_of_horizontal_crs_create(self):
        """ Creates a blank, empty or zero value for type_of_horizontal_crs
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.type_of_horizontal_crs = list(self.__type_of_horizontal_crs_type__)[0]

    @property
    def horizontal_crs(self) -> int:
        return self._attributes[self.__horizontal_crs_hdf_name__]

    @horizontal_crs.setter
    def horizontal_crs(self, val: int):
        self._attributes[self.__horizontal_crs_hdf_name__] = val

    @property
    def __horizontal_crs_type__(self) -> Type[int]:
        return numpy.int32

    def horizontal_crs_create(self):
        """ Creates a blank, empty or zero value for horizontal_crs"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_crs = self.__horizontal_crs_type__()

    # Even though the horizontal_cs is an enum, it is stored as an int per the spec in table 10c-6
    @property
    def horizontal_cs(self) -> int:
        return self._attributes[self.__horizontal_cs_hdf_name__]

    @horizontal_cs.setter
    def horizontal_cs(self, val: Union[int, str, HORIZONTAL_CS]):
        # Use an enumeration to control the values, though not officially an enum it essentially is, throws an exception if not in the enum
        self.set_enum_attribute(val, self.__horizontal_cs_hdf_name__, HORIZONTAL_CS)
        # convert the legal values into an integer
        self._attributes[self.__horizontal_cs_hdf_name__] = self._attributes[self.__horizontal_cs_hdf_name__].value

    @property
    def __horizontal_cs_type__(self) -> Type[int]:
        return numpy.int32

    def horizontal_cs_create(self):
        """ Creates a blank, empty or zero value for horizontal_cs
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_cs = list(HORIZONTAL_CS)[0].value  # this is limited to an enumerated list but they specify it as an int

    @property
    def horizontal_datum(self) -> int:
        return self._attributes[self.__horizontal_datum_hdf_name__]

    @horizontal_datum.setter
    def horizontal_datum(self, val: int):
        self._attributes[self.__horizontal_datum_hdf_name__] = val

    @property
    def __horizontal_datum_type__(self) -> Type[int]:
        return numpy.int32

    def horizontal_datum_create(self):
        """ Creates a blank, empty or zero value for horizontal_datum"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_datum = self.__horizontal_datum_type__()

    @property
    def name_of_horizontal_datum(self) -> str:
        return self._attributes[self.__name_of_horizontal_datum_hdf_name__]

    @name_of_horizontal_datum.setter
    def name_of_horizontal_datum(self, val: str):
        self._attributes[self.__name_of_horizontal_datum_hdf_name__] = val

    @property
    def __name_of_horizontal_datum_type__(self) -> Type[str]:
        return str

    def name_of_horizontal_datum_create(self):
        """ Creates a blank, empty or zero value for name_of_horizontal_datum"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.name_of_horizontal_datum = self.__name_of_horizontal_datum_type__()


    @property
    def prime_meridian(self) -> int:
        return self._attributes[self.__prime_meridian_hdf_name__]

    @prime_meridian.setter
    def prime_meridian(self, val: int):
        self._attributes[self.__prime_meridian_hdf_name__] = val

    @property
    def __prime_meridian_type__(self) -> Type[int]:
        return numpy.int32

    def prime_meridian_create(self):
        """ Creates a blank, empty or zero value for prime_meridian"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.prime_meridian = self.__prime_meridian_type__()


    """NOTE: All latitudes and longitudes of the projection parameters must be given in degrees (south and
west negative). Azimuths are given in degrees. For detailed description of the projection method refer to
the EPSG documentation."""

    @property
    def spheriod(self) -> int:
        return self._attributes[self.__spheriod_hdf_name__]

    @spheriod.setter
    def spheriod(self, val: int):
        self._attributes[self.__spheriod_hdf_name__] = val

    @property
    def __spheriod_type__(self) -> Type[int]:
        return numpy.int32

    def spheriod_create(self):
        """ Creates a blank, empty or zero value for spheriod"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.spheriod = self.__spheriod_type__()

    @property
    def projection_method(self) -> PROJECTION_METHOD:
        return self._attributes[self.__projection_method_hdf_name__]

    @projection_method.setter
    def projection_method(self, val: Union[int, str, PROJECTION_METHOD]):
        # Use an enumeration to control the values, though not officially an enum it essentially is
        self.set_enum_attribute(val, self.__projection_method_hdf_name__, PROJECTION_METHOD)
        # convert the legal values into an integer
        self._attributes[self.__projection_method_hdf_name__] = self._attributes[self.__projection_method_hdf_name__].value

    @property
    def __projection_method_type__(self) -> Type[PROJECTION_METHOD]:
        return numpy.int32  # PROJECTION_METHOD

    def projection_method_create(self):
        """ Creates a blank, empty or zero value for projection_method
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_method = list(PROJECTION_METHOD)[0]

    @property
    def projection_parameter_1(self) -> float:
        return self._attributes[self.__projection_parameter_1_hdf_name__]

    @projection_parameter_1.setter
    def projection_parameter_1(self, val: float):
        self._attributes[self.__projection_parameter_1_hdf_name__] = val

    @property
    def __projection_parameter_1_type__(self) -> Type[float]:
        return numpy.float64

    def projection_parameter_1_create(self):
        """ Creates a blank, empty or zero value for projection_parameter_1"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_parameter_1 = self.__projection_parameter_1_type__()

    @property
    def projection_parameter_2(self) -> float:
        return self._attributes[self.__projection_parameter_2_hdf_name__]

    @projection_parameter_2.setter
    def projection_parameter_2(self, val: float):
        self._attributes[self.__projection_parameter_2_hdf_name__] = val

    @property
    def __projection_parameter_2_type__(self) -> Type[float]:
        return numpy.float64

    def projection_parameter_2_create(self):
        """ Creates a blank, empty or zero value for projection_parameter_2"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_parameter_2 = self.__projection_parameter_2_type__()

    @property
    def projection_parameter_3(self) -> float:
        return self._attributes[self.__projection_parameter_3_hdf_name__]

    @projection_parameter_3.setter
    def projection_parameter_3(self, val: float):
        self._attributes[self.__projection_parameter_3_hdf_name__] = val

    @property
    def __projection_parameter_3_type__(self) -> Type[float]:
        return numpy.float64

    def projection_parameter_3_create(self):
        """ Creates a blank, empty or zero value for projection_parameter_3"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_parameter_3 = self.__projection_parameter_3_type__()

    @property
    def projection_parameter_4(self) -> float:
        return self._attributes[self.__projection_parameter_4_hdf_name__]

    @projection_parameter_4.setter
    def projection_parameter_4(self, val: float):
        self._attributes[self.__projection_parameter_4_hdf_name__] = val

    @property
    def __projection_parameter_4_type__(self) -> Type[float]:
        return numpy.float64

    def projection_parameter_4_create(self):
        """ Creates a blank, empty or zero value for projection_parameter_4"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_parameter_4 = self.__projection_parameter_4_type__()

    @property
    def projection_parameter_5(self) -> float:
        return self._attributes[self.__projection_parameter_5_hdf_name__]

    @projection_parameter_5.setter
    def projection_parameter_5(self, val: float):
        self._attributes[self.__projection_parameter_5_hdf_name__] = val

    @property
    def __projection_parameter_5_type__(self) -> Type[float]:
        return numpy.float64

    def projection_parameter_5_create(self):
        """ Creates a blank, empty or zero value for projection_parameter_5"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_parameter_5 = self.__projection_parameter_5_type__()

    @property
    def false_northing(self) -> float:
        return self._attributes[self.__false_northing_hdf_name__]

    @false_northing.setter
    def false_northing(self, val: float):
        self._attributes[self.__false_northing_hdf_name__] = val

    @property
    def __false_northing_type__(self) -> Type[float]:
        return numpy.float64

    def false_northing_create(self):
        """ Creates a blank, empty or zero value for false_northing"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.false_northing = self.__false_northing_type__()

    @property
    def false_easting(self) -> float:
        return self._attributes[self.__false_easting_hdf_name__]

    @false_easting.setter
    def false_easting(self, val: float):
        self._attributes[self.__false_easting_hdf_name__] = val

    @property
    def __false_easting_type__(self) -> Type[float]:
        return numpy.float64

    def false_easting_create(self):
        """ Creates a blank, empty or zero value for false_easting"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.false_easting = self.__false_easting_type__()

    @property
    def vertical_cs(self) -> VERTICAL_CS:
        return self._attributes[self.__vertical_cs_hdf_name__]

    @vertical_cs.setter
    def vertical_cs(self, val: Union[int, str, VERTICAL_CS]):
        # Use an enumeration to control the values, though not officially an enum it essentially is
        self.set_enum_attribute(val, self.__vertical_cs_hdf_name__, VERTICAL_CS)
        # convert the legal values into an integer
        self._attributes[self.__vertical_cs_hdf_name__] = self._attributes[self.__vertical_cs_hdf_name__].value

    @property
    def __vertical_cs_type__(self) -> Type[VERTICAL_CS]:
        return numpy.int32  # VERTICAL_CS

    def vertical_cs_create(self):
        """ Creates a blank, empty or zero value for vertical_cs
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_cs = list(VERTICAL_CS)[0]


    @property
    def vertical_coordinate_base(self) -> VERTICAL_COORDINATE_BASE:
        return self._attributes[self.__vertical_coordinate_base_hdf_name__]

    @vertical_coordinate_base.setter
    def vertical_coordinate_base(self, val: Union[int, str, VERTICAL_COORDINATE_BASE]):
        self.set_enum_attribute(val, self.__vertical_coordinate_base_hdf_name__, self.__vertical_coordinate_base_type__)

    @property
    def __vertical_coordinate_base_type__(self) -> Type[VERTICAL_COORDINATE_BASE]:
        return VERTICAL_COORDINATE_BASE

    def vertical_coordinate_base_create(self):
        """ Creates a blank, empty or zero value for vertical_coordinate_base
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_coordinate_base = list(self.__vertical_coordinate_base_type__)[0]

    @property
    def vertical_datum_reference(self) -> VERTICAL_DATUM_REFERENCE:
        return self._attributes[self.__vertical_datum_reference_hdf_name__]

    @vertical_datum_reference.setter
    def vertical_datum_reference(self, val: Union[int, str, VERTICAL_DATUM_REFERENCE]):
        self.set_enum_attribute(val, self.__vertical_datum_reference_hdf_name__, self.__vertical_datum_reference_type__)

    @property
    def __vertical_datum_reference_type__(self) -> Type[VERTICAL_DATUM_REFERENCE]:
        return VERTICAL_DATUM_REFERENCE

    def vertical_datum_reference_create(self):
        """ Creates a blank, empty or zero value for vertical_datum_reference
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum_reference = list(VERTICAL_DATUM_REFERENCE)[0]

    @property
    def vertical_datum(self) -> Enum:
        val = self._attributes[self.__vertical_datum_hdf_name__]
        return val

    @vertical_datum.setter
    def vertical_datum(self, val: (int, str, VERTICAL_DATUM)):
        # NOTE: When reading from a file h5py gets attributes alphabetically so we can't rely on vertical_datum_reference being set before vertical_datum
        # verticalDatumReference when 1 is from the enumeration but when 2 is an EPSG code that we don't have a check for.
        try:
            self.vertical_datum_reference  # see if the attribute exists
        except KeyError:
            pass
        else:
            if self.vertical_datum_reference == VERTICAL_DATUM_REFERENCE(1):
                try:
                    self.set_enum_attribute(val, self.__vertical_datum_hdf_name__, self.__vertical_datum_restriction__)
                except S100Exception as e:
                    raise S100Exception(f"When vertical_datum_reference is '1' then vertical_datum must be a value given in the enumeration {VERTICAL_DATUM}, the supplied {val} was not found")
                # convert the enumeration back to an integer
                val = self._attributes[self.__vertical_datum_hdf_name__].value
        try:
            self._attributes[self.__vertical_datum_hdf_name__] = int(val)
        except ValueError:
            raise S100Exception(f"vertical_datum must be an integer when vertical_datum_reference is '2' or not set, the supplied {val} was not an integer")

    @property
    def __vertical_datum_restriction__(self):
        return VERTICAL_DATUM

    @property
    def __vertical_datum_type__(self) -> Type[int]:
        return int

    def vertical_datum_create(self):
        """ Creates a blank, empty or zero value for vertical_datum"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum = VERTICAL_DATUM.MLLW.value  # "MLLW"


class S100File(S1XXFile):
    PRODUCT_SPECIFICATION = PRODUCT_SPECIFICATION

    def __init__(self, *args, **kywrds):
        if 'root' not in kywrds:
            kywrds['root'] = S100Root  # inherited classes will specify their own root type
        super().__init__(*args, **kywrds)

    @property
    def epsg(self):
        return self.root.horizontal_crs

    def iter_groups(self):
        """Create a 2-Band GeoTIFF for every speed and direction compound dataset
           within each HDF5 file(s).

        Args:
            input_path: Path to a single S-111 HDF5 file or a directory containing
                one or more.
            output_path: Path to a directory where GeoTIFF file(s) will be
                generated.
        """
        # @TODO this could be a generic method of S100 files but we need to formalize how to find the datasets
        #   It is close and with some hardcoded assumptions it would work,
        #   but we need to know what the feature name(s) are (from GroupF which is doable)
        #   then know the feature infomation instance (like BathymetryCoverage.01 but the formating changes depending on version)
        #   then get the feature group instances (Group_001)
        #   then the value names (depth, uncertainty)
        #   The biggest issue is
        dataname = self.root.feature_information.feature_code[0]
        if isinstance(dataname, bytes):
            dataname = dataname.decode()
        data = self.root.get_s1xx_attr(dataname)
        for instance_key in self[data._hdf5_path].keys():
            if re.match(dataname+r"[._]\d{2,3}", instance_key):
                feature_instance = self["/".join([data._hdf5_path, instance_key])]
                num_grp = feature_instance.attrs['numGRP']

                for idx in range(1, num_grp + 1):
                    for group_key in feature_instance.keys():
                        if re.match(f"Group[._]0*{idx}", group_key):  # Find anything from Group_1 to Group.001
                            group_instance = feature_instance[group_key]
                            yield data, dataname, feature_instance, group_instance

    def to_raster_datasets(self):
        """ Return a GDAL raster dataset with a band for each gridded dataset in the HDF5 file.
        The returned dataset will use the upper left corner and have a negative dyy value in the geotransform like the TIF convention.

        Returns
        -------
        tuple
            (dataset, group_instance) -- a GDAL raster and a h5py group instance
            which could be queried for additional attributes the raster was created from

        """
        if self.epsg < 0:
            raise ValueError("Unable to convert to GeoTIFF.  EPSG code is not valid.")
        for data, dataname, feature_instance, group_instance in self.iter_groups():
            if data.data_coding_format.value not in (2, 9):
                raise S100Exception(f"Unable to convert to GeoTIFF.  Data coding format ({data.data_coding_format.value}) must be regular grid (2 or 9).")
            values = group_instance['values']
            bands = []
            for i, name in enumerate(self['Group_F'][dataname]['code']):
                if isinstance(name, bytes):  # in h5py 3.0+ the string datasets are bytes
                    name = name.decode()
                if name in values.dtype.names:
                    bands.append([name, float(self['Group_F'][dataname]['fillValue'][i])])

            y_dim, x_dim = values[bands[0][0]].shape
            dxx = feature_instance.attrs['gridSpacingLongitudinal']
            dyy = feature_instance.attrs['gridSpacingLatitudinal']
            x_min = feature_instance.attrs['gridOriginLongitude']
            y_min = feature_instance.attrs['gridOriginLatitude']
            # TIFF convention is to list the upper left corner of the upper left pixel
            # S100 convention is to list the center of the bottom left pixel
            # But rather than assuming, check the gridspacing and flip if necessary
            flipx, flipy = False, False
            # Again, TIF convention is to use upper left which yields a negative delta Y.
            # That doesn't seem significant but some software (like QGIS) will not display properly if dyy is positive.
            # QGIS is internally reversing the value but also modifying it in the fourth decimal place making a slight visual distortion.
            if dyy > 0:  # the origin is the bottom corner - so flip the dyy and move the origin to the other side
                # remove one from dimension to stay to inside the grid.  Imagine a 1 cell rester -- the move has to be zero.
                y_min = y_min + dyy * (y_dim - 1)
                dyy *= -1
                flipy = True
            if dxx < 0:  # the origin is the right corner - so flip the dxx and move the origin to the other side
                # remove one from dimension to stay to inside the grid.  Imagine a 1 cell rester -- the move has to be zero.
                x_min = x_min + dxx * (x_dim - 1)
                dxx *= -1
                flipx = True
            # now we know the origin is upper left, so we can adjust from center to upper left of the cell
            geoTransform = [x_min - dxx / 2.0,  # dxx is positive so this moves the origin left by half a cell
                            dxx,
                            0,  # dxy - zero because we are not rotating the raster
                            y_min - dyy / 2.0,  # dyy is negative so this is adding half a cell height
                            0,  # dyx - zero because we are not rotating the raster
                            dyy
                            ]

            srs = osr.SpatialReference()
            srs.ImportFromEPSG(int(self.epsg))

            num_bands = len(bands)
            dataset = gdal.GetDriverByName('MEM').Create('', x_dim, y_dim, num_bands, gdal.GDT_Float32)
            dataset.SetGeoTransform(geoTransform)
            dataset.SetProjection(srs.ExportToWkt())
            for band_num, (band_name, fill_value) in enumerate(bands):
                band_values = values[band_name]
                if flipx:
                    band_values = numpy.fliplr(band_values)
                if flipy:
                    band_values = numpy.flipud(band_values)
                dataset.GetRasterBand(band_num+1).WriteArray(band_values)
                dataset.GetRasterBand(band_num+1).SetDescription(band_name)
                dataset.GetRasterBand(band_num+1).SetNoDataValue(fill_value)
            # Forces geotiffs to use nodata value from first band rather than the last band
            # that contains a different fillvalue
            dataset.GetRasterBand(2).SetNoDataValue(bands[0][1])
            yield dataset, group_instance, flipx, flipy

    def to_geotiffs(self, output_directory: (str, pathlib.Path), creation_options: list=None):
        """ Creates a GeoTIFF file for each regularly gridded dataset in the HDF5 file.
        If only one dataset is found, a single GeoTIFF file will be created named the same as the .h5 but with a .tif extension.
        If multiple datasets are found, multiple GeoTIFF files will be created named the same as the .h5 but with a _{timepoint}.tif extension.
        Supports DCF 2 and 9.

        Parameters
        ----------
        output_directory
            Directory of the location to save the GeoTIFF file(s)
        creation_options
            List of GDAL creation options

        Returns
        -------
        list
            List of filenames created
        """
        filenames = []
        for gdal_dataset, group_instance, flipx, flipy in self.to_raster_datasets():
            split_path = os.path.split(self.filename)
            filename = os.path.splitext(split_path[1])
            try:
                timepoint = group_instance.attrs['timePoint']
            except KeyError:
                datetime_str = ""
            else:
                try:
                    timepoint_str = datetime.datetime.strptime(timepoint, "%Y-%m-%dT%H:%M:%S")
                except ValueError:
                    timepoint_str = datetime.datetime.strptime(timepoint, "%Y%m%dT%H%M%SZ")

                datetime_str = timepoint_str.strftime("_%Y%m%dT%H%M%SZ")

            name = os.path.join(str(output_directory), '{}{}.tif'.format(filename[0], datetime_str))

            if creation_options is None:
                creation_options = []
            gdal.GetDriverByName('GTiff').CreateCopy(name, gdal_dataset, options=creation_options)
            filenames.append(name)
        return filenames

    def to_vector_dataset(self):
        """ Create an osgeo.ogr vector datasource with a layer for each dataset in the HDF5 file.
        Currently only supports 'Ungeorectified gridded arrays' (DCF=3).

        Returns
        -------
        osgeo.ogr.DataSource

        """
        if self.epsg < 0:
            raise ValueError(f"Unable to convert to vector layers.  EPSG code {self.epsg} is not valid.")
        ds = ogr.GetDriverByName("MEMORY").CreateDataSource('memData')
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(self.epsg))
        for data, dataname, feature_instance, group_instance in self.iter_groups():
            dcf = data.data_coding_format.value
            if dcf not in (3, 7):
                raise S100Exception("Unable to convert to Geopackage.  Data coding format must be 3 or 7.")
            values = group_instance['values']
            positions = feature_instance['Positioning']['geometryValues']
            layer_name = "_".join([os.path.split(feature_instance.name)[-1], os.path.split(group_instance.name)[-1]])
            layer = ds.CreateLayer(layer_name, srs, ogr.wkbPoint)
            fields = []
            col_info = []
            for col_name, col_type in values.dtype.descr:
                if isinstance(col_type, (list, tuple)):  # enumerations have something like ('|u1', {'enum': {'Decreasing': 1, ...}})
                    if numpy.dtype(col_type[0]) in (numpy.int32, numpy.int64, numpy.int8, numpy.int16, numpy.int64,
                                                    numpy.uint32, numpy.uint64, numpy.uint8, numpy.uint16, numpy.uint):
                        col_type = ogr.OFTInteger
                        converter = int
                    else:
                        raise S100Exception(f"Unable to convert to GeoPackage.  Data type {col_type} is not supported.")
                elif numpy.dtype(col_type) in (numpy.float32, numpy.float64):
                    col_type = ogr.OFTReal
                    converter = float
                elif numpy.dtype(col_type) in (numpy.int32, numpy.int64, numpy.int8, numpy.int16, numpy.int64):
                    col_type = ogr.OFTInteger
                    converter = int
                else:
                    raise S100Exception(f"Unable to convert to GeoPackage.  Data type {col_type} is not supported.")
                col_info.append([col_name, col_type, converter])
                fields.append(ogr.FieldDefn(col_name, col_type))

            layer.CreateFields(fields)
            layer_defn = layer.GetLayerDefn()
            layer.StartTransaction()
            for pt_data, (longitude, latitude) in zip(values, positions):
                feat = ogr.Feature(layer_defn)
                pt = ogr.Geometry(ogr.wkbPoint)  # or wkbMultiPoint
                pt.AddPoint_2D(float(longitude), float(latitude))  # lat/lon or x/y
                for (col_name, col_type, converter), col_val in zip(col_info, pt_data):
                    feat.SetField(col_name, converter(col_val))
                feat.SetGeometry(pt)
                layer.CreateFeature(feat)
            layer.CommitTransaction()
            if dcf == 7:  # create a triangles layer with TriangleZM to hold the two expected additional values
                if False:
                    # Create test polygon
                    ring = ogr.Geometry(ogr.wkbLinearRing)
                    ring.AddPoint(1179091.1646903288, 712782.8838459781)
                    ring.AddPoint(1161053.0218226474, 667456.2684348812)
                    ring.AddPoint(1214704.933941905, 641092.8288590391)
                    ring.AddPoint(1228580.428455506, 682719.3123998424)
                    ring.AddPoint(1218405.0658121984, 721108.1805541387)
                    ring.AddPoint(1179091.1646903288, 712782.8838459781)
                    poly = ogr.Geometry(ogr.wkbPolygon)
                    poly.AddGeometry(ring)

                    # Create the output Driver
                    outDriver = ogr.GetDriverByName('GeoJSON')

                    # # Create the output GeoJSON
                    # outDataSource = outDriver.CreateDataSource('test.geojson')
                    # outLayer = outDataSource.CreateLayer('test.geojson', geom_type=ogr.wkbPolygon)

                    # Get the output Layer's Feature Definition
                    featureDefn = layer.GetLayerDefn()

                    # create a new feature
                    outFeature = ogr.Feature(featureDefn)

                    # Set new geometry
                    outFeature.SetGeometry(poly)

                    # Add new feature to output Layer
                    layer.CreateFeature(outFeature)

                    # dereference the feature
                    outFeature = None

                    # Save and close DataSources
                    outDataSource = None
                else:
                    layer = ds.CreateLayer(layer_name+"_tin", srs, ogr.wkbPolygonZM)
                    # layer.CreateFields(fields)
                    layer_defn = layer.GetLayerDefn()
                    layer.StartTransaction()
                    triangles = feature_instance['Positioning']['triangles']
                    for nodes in triangles:
                        feat = ogr.Feature(layer_defn)
                        ring = ogr.Geometry(ogr.wkbLinearRing)  # ogr.wkbTriangle)  # or wkbMultiPoint
                        # tri = ogr.Geometry(ogr.wkbTriangle)  # Triangle was coming up with non-standard extension message, so use Polygon
                        for node in list(nodes)+[nodes[0]]:
                            longitude, latitude = positions[node]
                            pt_data = values[node]
                            vals = []
                            for (col_name, col_type, converter), col_val in zip(col_info, pt_data):
                                vals.append(converter(col_val))
                            z, m = vals[:2]  # add the two attributes the vertices have as Z and M values
                            ring.AddPointZM(float(longitude), float(latitude), z, m)  # lat/lon or x/y
                        poly = ogr.Geometry(ogr.wkbPolygonZM)
                        poly.AddGeometry(ring)
                        feat.SetGeometry(poly)
                        # feat.SetGeometry(tri)
                        layer.CreateFeature(feat)
                    layer.CommitTransaction()

        return ds

    def to_geopackage(self, output_path: (str, pathlib.Path)=None):
        """ Create a geopackage file with a layer for each dataset in the HDF5 file.
        Based on the to_vector_dataset method, so only supports ungeorectified gridded data (DCF=3).

        Parameters
        ----------
        output_path
            Full path of the geopackage file, if None then the same name as the HDF5 file will be used with a .gpkg extension

        Returns
        -------
        None

        """
        if output_path is None:
            output_path = os.path.splitext(self.filename)[0] + '.gpkg'
        ds = self.to_vector_dataset()
        fix = ogr.GetDriverByName("GPKG").CopyDataSource(ds, str(output_path))

