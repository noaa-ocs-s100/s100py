from abc import ABC, abstractmethod
from typing import Callable, Iterator, Union, Optional, List, Type
import logging
import datetime
from enum import Enum

import h5py
# @todo - consider removing the numpy dependence
import numpy

try:
    from . import s1xx
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py"

from .s1xx import s1xx_sequence, S1xxAttributesBase, S1xxDatasetBase, S1XXFile, Chunking


class S100Exception(Exception):
    pass


H5T_CLASS_T = {
    h5py.h5t.NO_CLASS: 'H5T_NO_CLASS',
    h5py.h5t.INTEGER: 'H5T_INTEGER',
    h5py.h5t.FLOAT: 'H5T_FLOAT',
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


class VERTICAL_DATUM(Enum):
    """ Note: while a Vertical Datum can be created with the shorthand aliases, ex: MLWS, the string written and
    returned from the file/S100 object will be the official long name, e.g. "meanLowWaterSprings" etc.
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


HORIZONTAL_DATUM_REFERENCE = numpy.string_('EPSG')
DATA_CODING_FORMAT = Enum(value="DATA_CODING_FORMAT",
                          names=[
                              ('Time series at fixed stations', 1),
                              ('Regularly-gridded arrays', 2),
                              ('Ungeorectified gridded arrays', 3),
                              ('Moving platform', 4),
                              ('Irregular grid', 5),
                              ('Variable cell size', 6),
                              ('TIN', 7),
                          ]
                          )


class INTERPOLATION_TYPE(Enum):
    """
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
    linear = 2
    quadratic = 3
    cubic = 4
    bilinear = 5
    biquadratic = 6
    bicubic = 7
    lostarea = 8
    barycentric = 9
    discrete = 10


class COMMON_POINT_RULE(Enum):
    average = 1
    low = 2
    high = 3
    all = 4


class SEQUENCING_RULE_TYPE(Enum):
    linear = 1
    boustrophedonic = 2
    CantorDiagonal = 3
    spiral = 4
    Morton = 5
    Hilbert = 6


SEQUENCING_RULE_SCAN_DIRECTION = numpy.string_('longitude,latitude')
START_SEQUENCE = numpy.string_('0,0')


class DirectPosition(S1xxAttributesBase):
    """ 4.2.1.1.4 of v2.0.0
    """
    coordinate_attribute_name = "coordinate"  #: HDF5 naming
    dimension_attribute_name = "dimension"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def coordinate(self) -> s1xx_sequence:
        return self._attributes[self.coordinate_attribute_name]

    @coordinate.setter
    def coordinate(self, val: s1xx_sequence):
        self._attributes[self.coordinate_attribute_name] = val

    @property
    def coordinate_type(self):
        return numpy.ndarray

    def coordinate_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.coordinate = self.coordinate_type([2], numpy.float64)

    @property
    def dimension(self) -> int:
        return self._attributes[self.dimension_attribute_name]

    @dimension.setter
    def dimension(self, val: int):
        self._attributes[self.dimension_attribute_name] = val

    @property
    def dimension_type(self):
        return int

    def dimension_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.dimension = self.dimension_type()


class GridCoordinate(S1xxAttributesBase):
    """ 4.2.1.1.6 of v2.0.0
    """

    coord_values_attribute_name = "coordValues"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def coord_values(self) -> s1xx_sequence:
        """The attribute coordValues has the value class Sequence Integer that shall hold one integer value for each dimension of the grid.
        The ordering of these coordinate values shall be the same as that of the elements of axisNames.
        The value of a single coordinate shall be the number of offsets from the origin of the grid in the direction of a specific axis"""
        return self._attributes[self.coord_values_attribute_name]

    @coord_values.setter
    def coord_values(self, val: s1xx_sequence):
        self._attributes[self.coord_values_attribute_name] = val

    @property
    def coord_values_type(self):
        return numpy.ndarray

    def coord_values_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.coord_values = self.coord_values_type([2], numpy.int)


class GridEnvelope(S1xxAttributesBase):
    """ 4.2.1.1.5 of v2.0.0
    While I would think that the envelope would describe the real world extents of the grid,
    in the docs it describes the envelope as specifying the row/column offsets for the lower left and upper right
    coordinates using the integer indices (S100 and ISO 19123 sec. 8.3).  The real world coordinates are in the origin and offsetVectors instead.
    So this seems unnecessary since the value seems the same as the size of the matrix held, which can be gotten by reading that instead.

    https://www.fgdc.gov/standards/projects/frameword-data-standard/GI_FrameworkDataStandard_Part3_Elevation.doc/at_download/file&usg=AOvVaw07QEsNy5urachwIO1e4ALU
    """

    low_attribute_name = "low"  #: HDF5 naming
    high_attribute_name = "high"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def low(self) -> S1xxAttributesBase:
        return self._attributes[self.low_attribute_name]

    @low.setter
    def low(self, val: S1xxAttributesBase):
        self._attributes[self.low_attribute_name] = val

    @property
    def low_type(self):
        return GridCoordinate

    def low_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.low = self.low_type()

    @property
    def high(self) -> S1xxAttributesBase:
        return self._attributes[self.high_attribute_name]

    @high.setter
    def high(self, val: S1xxAttributesBase):
        self._attributes[self.high_attribute_name] = val

    @property
    def high_type(self):
        return GridCoordinate

    def high_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.high = self.high_type()


class SequenceRule(S1xxAttributesBase):
    """ 4.2.1.1.7 (and .8) of v2.0.0
    CV_SequenceRule specified in ISO 19123
    """

    type_attribute_name = "type"  #: HDF5 naming
    scan_direction_attribute_name = "scanDirection"  #: HDF5 naming

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
        return self._attributes[self.type_attribute_name]

    @type.setter
    def type(self, val: str):
        self._attributes[self.type_attribute_name] = val

    @property
    def type_type(self):
        return str

    def type_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.type = self.type_type("linear")

    @property
    def scan_direction(self) -> s1xx_sequence:
        """The attribute scanDirection has the value class Sequence<CharacterString> a list of axis names that indicates
        the order in which grid points shall be mapped to position within the sequence of records of feature attribute values.
        The scan direction for all layers in S-102 is "Longitude" and "Latitude" or west to east, then south to north.
        """
        return self._attributes[self.scan_direction_attribute_name]

    @scan_direction.setter
    def scan_direction(self, val: s1xx_sequence):
        self._attributes[self.scan_direction_attribute_name] = val

    @property
    def scan_direction_type(self):
        return str

    def scan_direction_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.scan_direction = self.scan_direction_type("Longitude, Latitude")


class Point(S1xxAttributesBase):
    """ 4.2.1.1.11 of v2.0.0
    The class GM_Point is taken from ISO 19107 and is the basic data type for a geometric object consisting of one and only one point.
    """
    position_attribute_name = "position"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def position(self) -> DirectPosition:
        """ DirectPosition - see Figure 7-3 in S100 v4.0.0
        Returns
        -------

        """
        return self._attributes[self.position_attribute_name]

    @position.setter
    def position(self, val: DirectPosition):
        self._attributes[self.position_attribute_name] = val

    @property
    def position_type(self):
        return DirectPosition

    def position_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.position = self.position_type()


class GeographicExtent(S1xxAttributesBase):
    """ 4.2.1.1.12 of v2.0.0
    The class EX_GeographicExtent is a metadata class from ISO 19115.
    It is a component of the metaclass EX_Extent.
    The use of EX_GeographicExtent is optional.
    When used it describes the spatial boundaries of the Tracking List elements within the bounds established by CV_GridEnvelope for the BathymetryCoverage.
    That is, the tracking list may carry information corresponding only to a portion of the spatial extent covered by the BathymetryCoverage.
    There is one attribute and one subtype.
    """

    extent_type_code_attribute_name = "extentTypeCode"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def extent_type_code(self) -> bool:
        """ The attribute extentTypeCode is a Boolean value.
        It is used to indicate whether the bounding polygon/box encompasses an area covered by the data or an area where data is not present.
        In S-102 it is set to 1.
        """
        return self._attributes[self.extent_type_code_attribute_name]

    @extent_type_code.setter
    def extent_type_code(self, val: bool):
        self._attributes[self.extent_type_code_attribute_name] = val

    @property
    def extent_type_code_type(self):
        return bool

    def extent_type_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.extent_type_code = self.extent_type_code_type()


class GeographicBoundingBox(GeographicExtent):
    """ 4.2.1.1.13 of v2.0.0
    see also S100 Tables 10C-6 and 10c-12
    The class EX_GeographicBoundingBox is a metadata class from ISO 19115.
    It is a subtype of the abstract class EX_GeographicExtent.
    It defines a bounding box used to indicate the spatial boundaries of the tracking list elements within the
    bounds established by CV_GridEnvelope for the BathymetryCoverage.
    """

    west_bound_longitude_attribute_name = "westBoundLongitude"  #: HDF5 naming
    east_bound_longitude_attribute_name = "eastBoundLongitude"  #: HDF5 naming
    south_bound_latitude_attribute_name = "southBoundLatitude"  #: HDF5 naming
    north_bound_latitude_attribute_name = "northBoundLatitude"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def west_bound_longitude(self) -> float:
        return self._attributes[self.west_bound_longitude_attribute_name]

    @west_bound_longitude.setter
    def west_bound_longitude(self, val: float):
        self._attributes[self.west_bound_longitude_attribute_name] = val

    @property
    def west_bound_longitude_type(self):
        return float

    def west_bound_longitude_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.west_bound_longitude = self.west_bound_longitude_type()

    @property
    def east_bound_longitude(self) -> float:
        return self._attributes[self.east_bound_longitude_attribute_name]

    @east_bound_longitude.setter
    def east_bound_longitude(self, val: float):
        self._attributes[self.east_bound_longitude_attribute_name] = val

    @property
    def east_bound_longitude_type(self):
        return float

    def east_bound_longitude_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.east_bound_longitude = self.east_bound_longitude_type()

    @property
    def south_bound_latitude(self) -> float:
        return self._attributes[self.south_bound_latitude_attribute_name]

    @south_bound_latitude.setter
    def south_bound_latitude(self, val: float):
        self._attributes[self.south_bound_latitude_attribute_name] = val

    @property
    def south_bound_latitude_type(self):
        return float

    def south_bound_latitude_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.south_bound_latitude = self.south_bound_latitude_type()

    @property
    def north_bound_latitude(self) -> float:
        return self._attributes[self.north_bound_latitude_attribute_name]

    @north_bound_latitude.setter
    def north_bound_latitude(self, val: float):
        self._attributes[self.north_bound_latitude_attribute_name] = val

    @property
    def north_bound_latitude_type(self):
        return float

    def north_bound_latitude_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.north_bound_latitude = self.north_bound_latitude_type()


class VertexPoint(S1xxAttributesBase):
    """ From Figure 8-21 in S100 v4.0.0

    """

    geometry_attribute_name = "geometry"  #: HDF5 naming
    value_attribute_name = "value"  # HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def geometry(self) -> Point:
        """ Derived from ISO 19107, referenced figure 7-3 and 8-A-5 of S100 v4.0.0
        Returns
        -------

        """
        return self._attributes[self.geometry_attribute_name]

    @geometry.setter
    def geometry(self, val: Point):
        self._attributes[self.geometry_attribute_name] = val

    @property
    def geometry_type(self):
        return Point

    def geometry_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.geometry = self.geometry_type()

    @property
    def value(self) -> s1xx_sequence:
        """The attribute value has the value class Record which is a sequence of value items that shall assign values to the discrete grid point.
        There are two values in each record in the S102_TrackingListValues class.
        These are the depth and the uncertainty values that were overridden in corresponding grid coverages

        'value' in tracking list should be HDF5 dataset of (depth, uncertainty)
        which are matched to the listSeries which holds the indices of the data locations

        It is an ISO 19103 class Record
        """
        return self._attributes[self.value_attribute_name]

    @value.setter
    def value(self, val: s1xx_sequence):
        self._attributes[self.value_attribute_name] = val

    @property
    def value_type(self):
        return numpy.ndarray

    def value_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.value = self.value_type([2, ], numpy.float)


class FeatureInstanceBase(GeographicBoundingBox):
    """ The feature instance group attributes from table 10c-12 in S100 spec
    """

    vertical_extent_minimum_z_attribute_name = "verticalExtent.minimumZ"
    vertical_extent_maximum_z_attribute_name = "verticalExtent.maximumZ"
    num_grp_attribute_name = "numGRP"
    instance_chunking_attribute_name = "instanceChunking"
    number_of_times_attribute_name = "numberOfTimes"
    time_record_interval_attribute_name = "timeRecordInterval"
    # @TODO  @FIXME -- first and last records are supposed to be datetime but S100 doc says 'character'  Need to create a datetime handler
    date_time_of_first_record_attribute_name = "dateTimeOfFirstRecord"
    date_time_of_last_record_attribute_name = "dateTimeOfLastRecord"


    def write(self, hdf5_object):
        super().write(hdf5_object)
        # find any group_NNN objects
        chunking = None
        for _pattern, group_attrib in self.get_standard_list_properties().items():
            group_list = self.__getattribute__(group_attrib)
            for grp in group_list:
                # now we are going to take advantage of the h5py interface to get the chunks attribute from each dataset
                # the S100 spec says things should be written with datasets named 'values'
                # if this does not hold true in the future then we could search for datasets generically here
                try:
                    chunking = hdf5_object[grp._hdf5_path.split("/")[-1] + '/values'].chunks
                except KeyError:
                    pass
        if chunking is not None:
            self.instance_chunking = chunking
            # now that we updated the chunking attribute we need to re-write them (but not the datasets etc)
            self.write_simple_attributes(hdf5_object)

    @property
    def vertical_extent_minimum_z(self) -> float:
        return self._attributes[self.vertical_extent_minimum_z_attribute_name]

    @vertical_extent_minimum_z.setter
    def vertical_extent_minimum_z(self, val: float):
        self._attributes[self.vertical_extent_minimum_z_attribute_name] = val

    @property
    def vertical_extent_minimum_z_type(self) -> Type[float]:
        return float

    def vertical_extent_minimum_z_create(self):
        """ Creates a blank, empty or zero value for vertical_extent_minimum_z"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_extent_minimum_z = self.vertical_extent_minimum_z_type()

    @property
    def vertical_extent_maximum_z(self) -> float:
        return self._attributes[self.vertical_extent_maximum_z_attribute_name]

    @vertical_extent_maximum_z.setter
    def vertical_extent_maximum_z(self, val: float):
        self._attributes[self.vertical_extent_maximum_z_attribute_name] = val

    @property
    def vertical_extent_maximum_z_type(self) -> Type[float]:
        return float

    def vertical_extent_maximum_z_create(self):
        """ Creates a blank, empty or zero value for vertical_extent_maximum_z"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_extent_maximum_z = self.vertical_extent_maximum_z_type()

    @property
    def num_grp(self) -> int:
        return self._attributes[self.num_grp_attribute_name]

    @num_grp.setter
    def num_grp(self, val: int):
        self._attributes[self.num_grp_attribute_name] = val

    @property
    def num_grp_type(self) -> Type[int]:
        return int

    def num_grp_create(self):
        """ Creates a blank, empty or zero value for num_grp"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.num_grp = self.num_grp_type()

    @property
    def instance_chunking(self) -> str:
        return self._attributes[self.instance_chunking_attribute_name]

    @instance_chunking.setter
    def instance_chunking(self, val: Union[str, list, tuple]):
        if isinstance(val, str):
            pass
        else:
            val = ",".join(str(a) for a in val)
        self._attributes[self.instance_chunking_attribute_name] = val

    @property
    def instance_chunking_type(self) -> Type[str]:
        return str

    def instance_chunking_create(self):
        """ Creates a blank, empty or zero value for instance_chunking"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.instance_chunking = self.instance_chunking_type()

    @property
    def number_of_times(self) -> int:
        return self._attributes[self.number_of_times_attribute_name]

    @number_of_times.setter
    def number_of_times(self, val: int):
        self._attributes[self.number_of_times_attribute_name] = val

    @property
    def number_of_times_type(self) -> Type[int]:
        return int

    def number_of_times_create(self):
        """ Creates a blank, empty or zero value for number_of_times"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.number_of_times = self.number_of_times_type()

    @property
    def time_record_interval(self) -> int:
        return self._attributes[self.time_record_interval_attribute_name]

    @time_record_interval.setter
    def time_record_interval(self, val: int):
        self._attributes[self.time_record_interval_attribute_name] = val

    @property
    def time_record_interval_type(self) -> Type[int]:
        return int

    def time_record_interval_create(self):
        """ Creates a blank, empty or zero value for time_record_interval"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.time_record_interval = self.time_record_interval_type()

    @property
    def date_time_of_first_record(self) -> str:
        return self._attributes[self.date_time_of_first_record_attribute_name]

    @date_time_of_first_record.setter
    def date_time_of_first_record(self, val: str):
        self._attributes[self.date_time_of_first_record_attribute_name] = val

    @property
    def date_time_of_first_record_type(self) -> Type[str]:
        return str

    def date_time_of_first_record_create(self):
        """ Creates a blank, empty or zero value for date_time_of_first_record"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.date_time_of_first_record = self.date_time_of_first_record_type()

    @property
    def date_time_of_last_record(self) -> str:
        return self._attributes[self.date_time_of_last_record_attribute_name]

    @date_time_of_last_record.setter
    def date_time_of_last_record(self, val: str):
        self._attributes[self.date_time_of_last_record_attribute_name] = val

    @property
    def date_time_of_last_record_type(self) -> Type[str]:
        return str

    def date_time_of_last_record_create(self):
        """ Creates a blank, empty or zero value for date_time_of_last_record"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.date_time_of_last_record = self.date_time_of_last_record_type()


class FeatureInstanceDCF2(FeatureInstanceBase):
    """ Data Coding Format 2 is the grid format from table 10c-12 in S100 spec.  Used in S102 for example.
    """

    grid_origin_longitude_attribute_name = "gridOriginLongitude"
    grid_origin_latitude_attribute_name = "gridOriginLatitude"
    grid_origin_vertical_attribute_name = "gridOriginVertical"
    grid_spacing_longitudinal_attribute_name = "gridSpacingLongitudinal"
    grid_spacing_latitudinal_attribute_name = "gridSpacingLatitudinal"
    grid_spacing_vertical_attribute_name = "gridSpacingVertical"
    num_points_longitudinal_attribute_name = "numPointsLongitudinal"
    num_points_latitudinal_attribute_name = "numPointsLatitudinal"
    num_points_vertical_attribute_name = "numPointsVertical"
    start_sequence_attribute_name = "startSequence"

    @property
    def grid_origin_longitude(self) -> float:
        return self._attributes[self.grid_origin_longitude_attribute_name]

    @grid_origin_longitude.setter
    def grid_origin_longitude(self, val: float):
        self._attributes[self.grid_origin_longitude_attribute_name] = val

    @property
    def grid_origin_longitude_type(self) -> Type[float]:
        return float

    def grid_origin_longitude_create(self):
        """ Creates a blank, empty or zero value for grid_origin_longitude"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_origin_longitude = self.grid_origin_longitude_type()

    @property
    def grid_origin_latitude(self) -> float:
        return self._attributes[self.grid_origin_latitude_attribute_name]

    @grid_origin_latitude.setter
    def grid_origin_latitude(self, val: float):
        self._attributes[self.grid_origin_latitude_attribute_name] = val

    @property
    def grid_origin_latitude_type(self) -> Type[float]:
        return float

    def grid_origin_latitude_create(self):
        """ Creates a blank, empty or zero value for grid_origin_latitude"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_origin_latitude = self.grid_origin_latitude_type()

    @property
    def grid_origin_vertical(self) -> float:
        return self._attributes[self.grid_origin_vertical_attribute_name]

    @grid_origin_vertical.setter
    def grid_origin_vertical(self, val: float):
        self._attributes[self.grid_origin_vertical_attribute_name] = val

    @property
    def grid_origin_vertical_type(self) -> Type[float]:
        return float

    def grid_origin_vertical_create(self):
        """ Creates a blank, empty or zero value for grid_origin_vertical"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_origin_vertical = self.grid_origin_vertical_type()

    @property
    def grid_spacing_longitudinal(self) -> float:
        return self._attributes[self.grid_spacing_longitudinal_attribute_name]

    @grid_spacing_longitudinal.setter
    def grid_spacing_longitudinal(self, val: float):
        self._attributes[self.grid_spacing_longitudinal_attribute_name] = val

    @property
    def grid_spacing_longitudinal_type(self) -> Type[float]:
        return float

    def grid_spacing_longitudinal_create(self):
        """ Creates a blank, empty or zero value for grid_spacing_longitudinal"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_spacing_longitudinal = self.grid_spacing_longitudinal_type()

    @property
    def grid_spacing_latitudinal(self) -> float:
        return self._attributes[self.grid_spacing_latitudinal_attribute_name]

    @grid_spacing_latitudinal.setter
    def grid_spacing_latitudinal(self, val: float):
        self._attributes[self.grid_spacing_latitudinal_attribute_name] = val

    @property
    def grid_spacing_latitudinal_type(self) -> Type[float]:
        return float

    def grid_spacing_latitudinal_create(self):
        """ Creates a blank, empty or zero value for grid_spacing_latitudinal"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_spacing_latitudinal = self.grid_spacing_latitudinal_type()

    @property
    def grid_spacing_vertical(self) -> float:
        return self._attributes[self.grid_spacing_vertical_attribute_name]

    @grid_spacing_vertical.setter
    def grid_spacing_vertical(self, val: float):
        self._attributes[self.grid_spacing_vertical_attribute_name] = val

    @property
    def grid_spacing_vertical_type(self) -> Type[float]:
        return float

    def grid_spacing_vertical_create(self):
        """ Creates a blank, empty or zero value for grid_spacing_vertical"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.grid_spacing_vertical = self.grid_spacing_vertical_type()

    @property
    def num_points_longitudinal(self) -> int:
        return self._attributes[self.num_points_longitudinal_attribute_name]

    @num_points_longitudinal.setter
    def num_points_longitudinal(self, val: int):
        self._attributes[self.num_points_longitudinal_attribute_name] = val

    @property
    def num_points_longitudinal_type(self) -> Type[int]:
        return int

    def num_points_longitudinal_create(self):
        """ Creates a blank, empty or zero value for num_points_longitudinal"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.num_points_longitudinal = self.num_points_longitudinal_type()

    @property
    def num_points_latitudinal(self) -> int:
        return self._attributes[self.num_points_latitudinal_attribute_name]

    @num_points_latitudinal.setter
    def num_points_latitudinal(self, val: int):
        self._attributes[self.num_points_latitudinal_attribute_name] = val

    @property
    def num_points_latitudinal_type(self) -> Type[int]:
        return int

    def num_points_latitudinal_create(self):
        """ Creates a blank, empty or zero value for num_points_latitudinal"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.num_points_latitudinal = self.num_points_latitudinal_type()

    @property
    def num_points_vertical(self) -> int:
        return self._attributes[self.num_points_vertical_attribute_name]

    @num_points_vertical.setter
    def num_points_vertical(self, val: int):
        self._attributes[self.num_points_vertical_attribute_name] = val

    @property
    def num_points_vertical_type(self) -> Type[int]:
        return int

    def num_points_vertical_create(self):
        """ Creates a blank, empty or zero value for num_points_vertical"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.num_points_vertical = self.num_points_vertical_type()

    @property
    def start_sequence(self) -> str:
        return self._attributes[self.start_sequence_attribute_name]

    @start_sequence.setter
    def start_sequence(self, val: str):
        self._attributes[self.start_sequence_attribute_name] = val

    @property
    def start_sequence_type(self) -> Type[str]:
        return str

    def start_sequence_create(self):
        """ Creates a blank, empty or zero value for start_sequence"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.start_sequence = self.start_sequence_type()


class FeatureInformation(S1xxAttributesBase):
    """  In S100, table 10c-8.
    In S102, 10.2.1 and table 10.2 and Table 10.1 of v2.0.0

    This is used to describe the BathymetryCoverage and TrackingListCoverage within the GroupF feature listing.
    The features described under GroupF have a matching named entry parallel to GroupF (top level).
    The actual data (depths etc) is stored in the top level element while basic metadata is stored in this element.
    """
    code_attribute_name = "code"
    name_attribute_name = "name"
    unit_of_measure_attribute_name = "uom.name"
    fill_value_attribute_name = "fillValue"
    datatype_attribute_name = "datatype"
    lower_attribute_name = "lower"
    upper_attribute_name = "upper"
    closure_attribute_name = "closure"

    @property
    def __version__(self) -> int:
        return 1

    def get_write_order(self):
        return [self.code_attribute_name,
                self.name_attribute_name,
                self.unit_of_measure_attribute_name,
                self.fill_value_attribute_name,
                self.datatype_attribute_name,
                self.lower_attribute_name,
                self.upper_attribute_name,
                self.closure_attribute_name]

    @property
    def code(self) -> str:
        """ The camel case name of the data

        Returns
        -------
        str
            The name of the dataset ("depth" or "uncertainty")
        """
        return self._attributes[self.code_attribute_name]

    @code.setter
    def code(self, val: str):
        self._attributes[self.code_attribute_name] = val

    @property
    def code_type(self):
        return str

    def code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.code = self.code_type()

    @property
    def name(self) -> str:
        """ The plain text name of the data
        Returns
        -------
        str
            Name of the dataset ("depth" or "uncertainty")
        """
        return self._attributes[self.name_attribute_name]

    @name.setter
    def name(self, val: str):
        self._attributes[self.name_attribute_name] = val

    @property
    def name_type(self):
        return str

    def name_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.name = self.name_type()

    @property
    def unit_of_measure(self) -> str:
        """ Units of measurement for the dataset
        Returns
        -------
        str
            "metres"
        """
        return self._attributes[self.unit_of_measure_attribute_name]

    @unit_of_measure.setter
    def unit_of_measure(self, val: str):
        self._attributes[self.unit_of_measure_attribute_name] = val

    @property
    def unit_of_measure_type(self):
        return str

    def unit_of_measure_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.unit_of_measure = self.unit_of_measure_type("metres")

    @property
    def fill_value(self) -> float:
        """ Value denoting missing data
        Returns
        -------
        float
            1000000
        """
        return self._attributes[self.fill_value_attribute_name]

    @fill_value.setter
    def fill_value(self, val: float):
        self._attributes[self.fill_value_attribute_name] = val

    @property
    def fill_value_type(self):
        return float

    def fill_value_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.fill_value = self.fill_value_type(1000000)

    @property
    def datatype(self) -> int:
        """
        Returns
        -------
        string
            H5T_NATIVE_FLOAT
        """
        return self._attributes[self.datatype_attribute_name]

    @datatype.setter
    def datatype(self, val: int):
        self._attributes[self.datatype_attribute_name] = val

    @property
    def datatype_type(self):
        return str

    def datatype_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.datatype = self.datatype_type("H5T_NATIVE_FLOAT")

    @property
    def lower(self) -> float:
        """
        Returns
        -------
        float
            -12000
        """
        return self._attributes[self.lower_attribute_name]

    @lower.setter
    def lower(self, val: float):
        self._attributes[self.lower_attribute_name] = val

    @property
    def lower_type(self):
        return float

    def lower_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.lower = self.lower_type(-12000)

    @property
    def upper(self) -> float:
        """
        Returns
        -------
        float
            12000
        """
        return self._attributes[self.upper_attribute_name]

    @upper.setter
    def upper(self, val: float):
        self._attributes[self.upper_attribute_name] = val

    @property
    def upper_type(self):
        return float

    def upper_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.upper = self.upper_type(12000)

    @property
    def closure(self) -> str:
        """
        Returns
        -------
        str
            closedInterval
        """
        return self._attributes[self.closure_attribute_name]

    @closure.setter
    def closure(self, val: str):
        self._attributes[self.closure_attribute_name] = val

    @property
    def closure_type(self):
        return str

    def closure_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.closure = self.closure_type("closedInterval")


class FeatureInformationDataset(Chunking, S1xxDatasetBase, ABC):
    """ This class comes from S100 -- 10c-9.5 Feature information group.

    metadata_type will likely be overridden with a specific subclass for the s100+ spec
    """

    @property
    def metadata_type(self) -> Type[FeatureInformation]:
        return FeatureInformation

    # chunking_attribute_name = "chunking"  #: HDF5 naming
    # @property
    # def chunking(self) -> str:
    #     return self._attributes[self.chunking_attribute_name]
    #
    # @chunking.setter
    # def chunking(self, val: str):
    #     self._attributes[self.chunking_attribute_name] = val
    #
    # @property
    # def chunking_type(self) -> Type[str]:
    #     return str
    #
    # def chunking_create(self):
    #     # noinspection PyAttributeOutsideInit
    #     # pylint: disable=attribute-defined-outside-init
    #     self.chunking = self.chunking_type()

class FeatureContainer(S1xxAttributesBase):
    axis_names_attribute_name = "axisNames"
    data_coding_format_attribute_name = "dataCodingFormat"
    dimension_attribute_name = "dimension"
    common_point_rule_attribute_name = "commonPointRule"
    horizontal_position_uncertainty_attribute_name = "horizontalPositionUncertainty"
    vertical_uncertainty_attribute_name = "verticalUncertainty"
    time_uncertainty_attribute_name = "timeUncertainty"
    num_instances_attribute_name = "numInstances"
    sequencing_rule_type_attribute_name = "sequencingRule.type"
    sequencing_rule_scan_direction_attribute_name = "sequencingRule.scanDirection"
    interpolation_type_attribute_name = "interpolationType"

    def __init__(self, *args, **opts):
        super().__init__(*args, **opts)
        self.data_coding_format_create()  # this is defined by the subclass and is constant, so we will automatically set it here

    @property
    def axis_names(self) -> s1xx_sequence:
        """sequence of character strings"""
        return self._attributes[self.axis_names_attribute_name]

    @axis_names.setter
    def axis_names(self, val: s1xx_sequence):
        self._attributes[self.axis_names_attribute_name] = val

    @property
    def axis_names_type(self) -> Type[str]:
        return numpy.ndarray

    def axis_names_create(self):
        """ The attribute axisNames has the value class Sequence<CharacterString> that shall be used to assign names to the grid axis.
        The grid axis names shall be "Latitude" and "Longitude" for unprojected data sets or “Northing” and “Easting” in a projected space.

        Returns
        -------

        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.axis_names = self.axis_names_type([2], dtype='S')

    @property
    def data_coding_format(self) -> DATA_CODING_FORMAT:
        return self._attributes[self.data_coding_format_attribute_name]

    @data_coding_format.setter
    def data_coding_format(self, val: int):
        self.set_enum_attribute(val, self.data_coding_format_attribute_name, self.data_coding_format_type)

    @property
    def data_coding_format_type(self) -> DATA_CODING_FORMAT:
        return DATA_CODING_FORMAT

    @abstractmethod
    def data_coding_format_create(self):
        """ Creates a blank, empty or zero value for data_coding_format"""
        raise NotImplementedError("each s100+ spec implementation must override this  data coding format with the correct default")

    @property
    def dimension(self) -> int:
        return self._attributes[self.dimension_attribute_name]

    @dimension.setter
    def dimension(self, val: int):
        self._attributes[self.dimension_attribute_name] = val

    @property
    def dimension_type(self) -> Type[int]:
        return int

    def dimension_create(self):
        """ Creates a blank, empty or zero value for dimension"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.dimension = self.dimension_type()

    @property
    def common_point_rule(self) -> COMMON_POINT_RULE:
        return self._attributes[self.common_point_rule_attribute_name]

    @common_point_rule.setter
    def common_point_rule(self, val: Union[int, str, COMMON_POINT_RULE]):
        self.set_enum_attribute(val, self.common_point_rule_attribute_name, self.common_point_rule_type)

    @property
    def common_point_rule_type(self) -> Type[Enum]:
        return COMMON_POINT_RULE

    def common_point_rule_create(self):
        """ Creates a blank, empty or zero value for common_point_rule"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.common_point_rule = self.common_point_rule_type["average"]

    @property
    def horizontal_position_uncertainty(self) -> float:
        return self._attributes[self.horizontal_position_uncertainty_attribute_name]

    @horizontal_position_uncertainty.setter
    def horizontal_position_uncertainty(self, val: float):
        self._attributes[self.horizontal_position_uncertainty_attribute_name] = val

    @property
    def horizontal_position_uncertainty_type(self) -> Type[float]:
        return float

    def horizontal_position_uncertainty_create(self):
        """ Creates a blank, empty or zero value for horizontal_position_uncertainty"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_position_uncertainty = self.horizontal_position_uncertainty_type()

    @property
    def vertical_uncertainty(self) -> float:
        return self._attributes[self.vertical_uncertainty_attribute_name]

    @vertical_uncertainty.setter
    def vertical_uncertainty(self, val: float):
        self._attributes[self.vertical_uncertainty_attribute_name] = val

    @property
    def vertical_uncertainty_type(self) -> Type[float]:
        return float

    def vertical_uncertainty_create(self):
        """ Creates a blank, empty or zero value for vertical_uncertainty"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_uncertainty = self.vertical_uncertainty_type()

    @property
    def time_uncertainty(self) -> float:
        return self._attributes[self.time_uncertainty_attribute_name]

    @time_uncertainty.setter
    def time_uncertainty(self, val: float):
        self._attributes[self.time_uncertainty_attribute_name] = val

    @property
    def time_uncertainty_type(self) -> Type[float]:
        return float

    def time_uncertainty_create(self):
        """ Creates a blank, empty or zero value for time_uncertainty"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.time_uncertainty = self.time_uncertainty_type()

    @property
    def num_instances(self) -> int:
        return self._attributes[self.num_instances_attribute_name]

    @num_instances.setter
    def num_instances(self, val: int):
        self._attributes[self.num_instances_attribute_name] = val

    @property
    def num_instances_type(self) -> Type[int]:
        return int

    def num_instances_create(self):
        """ Creates a blank, empty or zero value for num_instances"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.num_instances = self.num_instances_type()

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
        return self._attributes[self.sequencing_rule_type_attribute_name]

    @sequencing_rule_type.setter
    def sequencing_rule_type(self, val: Union[int, str, SEQUENCING_RULE_TYPE]):
        self.set_enum_attribute(val, self.sequencing_rule_type_attribute_name, self.sequencing_rule_type_type)

    @property
    def sequencing_rule_type_type(self) -> Type[Enum]:
        return SEQUENCING_RULE_TYPE

    def sequencing_rule_type_create(self):
        """ Creates a blank, empty or zero value for sequencing_rule_type"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.sequencing_rule_type = self.sequencing_rule_type_type["linear"]

    @property
    def sequencing_rule_scan_direction(self) -> str:
        return self._attributes[self.sequencing_rule_scan_direction_attribute_name]

    @sequencing_rule_scan_direction.setter
    def sequencing_rule_scan_direction(self, val: str):
        self._attributes[self.sequencing_rule_scan_direction_attribute_name] = val

    @property
    def sequencing_rule_scan_direction_type(self) -> Type[str]:
        return str

    def sequencing_rule_scan_direction_create(self):
        """ Creates a blank, empty or zero value for sequencing_rule_scan_direction"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.sequencing_rule_scan_direction = self.sequencing_rule_scan_direction_type()

    @property
    def interpolation_type(self) -> Type[int]:
        """ S100 table 10c-21

        Returns
        -------

        """
        return self._attributes[self.interpolation_type_attribute_name]

    @interpolation_type.setter
    def interpolation_type(self, val: Union[int, str, INTERPOLATION_TYPE]):
        self.set_enum_attribute(val, self.interpolation_type_attribute_name, self.interpolation_type_type)

    @property
    def interpolation_type_type(self) -> Type[Enum]:
        return INTERPOLATION_TYPE

    def interpolation_type_create(self):
        """ Creates a blank, empty or zero value for interpolation_type"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.interpolation_type = self.interpolation_type_type['nearestneighbor']

class GroupFBase(S1xxAttributesBase):
    """ Table 10.3 and sect 10.2.2 of v1.0.1
    """
    feature_code_attribute_name = "featureCode"

    @property
    def feature_code_type(self):
        return numpy.array

    @abstractmethod
    def feature_code_create(self):
        raise NotImplementedError("must overload feature_code_create")

    @property
    def feature_code(self) -> s1xx_sequence:
        return self._attributes[self.feature_code_attribute_name]

    @feature_code.setter
    def feature_code(self, val: s1xx_sequence):
        self._attributes[self.feature_code_attribute_name] = val


class S100Root(GeographicBoundingBox):
    """ From table 10c-6 in S100 spec.
    """
    feature_information_attribute_name = "Group_F"
    horizontal_datum_reference_attribute_name = "horizontalDatumReference"
    horizontal_datum_value_attribute_name = "horizontalDatumValue"
    epoch_attribute_name = "epoch"
    geographic_identifier_attribute_name = "geographicIdentifier"
    vertical_datum_attribute_name = "verticalDatum"
    meta_features_attribute_name = "metaFeatures"
    metadata_attribute_name = "metadata"
    product_specification_attribute_name = "productSpecification"
    issue_time_attribute_name = "issueTime"
    issue_date_attribute_name = "issueDate"

    @property
    def __version__(self) -> int:
        return 1

    def write(self, group_object):
        super().write(group_object)
        # any grids that were had datasets which possible chunking should now be written
        # and we can look through those to get the overall chunking attribute
        # and put that into the GroupF FeatureInformation object
        feat_info = None
        for property_name in self.get_standard_properties():
            if issubclass(self.__getattribute__(property_name + "_type"), GroupFBase):
                feat_info = self.__getattribute__(property_name)
        # we have the GroupF data now, we can look at the names of the FeatureInstances and then search each for its respective chunking
        if feat_info is not None:
            # this will be the names of the feature instances
            # e.g. BathymetryCoverage for S102 or SurfaceCurrent for S111
            for feat_name in feat_info.feature_code:
                # get the associated python name for the feature, e.g. turn SurfaceCurrent into surface_current
                chunking = None
                try:
                    python_name = self.get_standard_properties_mapping()[feat_name]
                except KeyError:
                    python_name = self.get_standard_properties_mapping()[feat_name.decode()]
                # grab the root/SurfaceCurrent data
                feat_container = self.__getattribute__(python_name)
                # now look through all the SurfaceCurrent_01, SurfaceCurrent_02...
                # so find the list object (there really only should be one and it should match the naming but we'll be general here)
                for pattern, list_name in feat_container.get_standard_list_properties().items():
                    try:
                        list_of_features = feat_container.__getattribute__(list_name)
                    except KeyError:  # not initialized
                        list_of_features = []
                    for feat_instance in list_of_features:
                        try:
                            chunking = feat_instance.instance_chunking
                        except:
                            pass
                if chunking is not None:
                    # find the GroupF feature dataset, e.g. /GroupF/SurrfaceCurrent
                    try:
                        groupf_python_name = feat_info.get_standard_properties_mapping()[feat_name]
                    except KeyError:
                        groupf_python_name = feat_info.get_standard_properties_mapping()[feat_name.decode()]
                    # Get the python object
                    feat_info_dataset = feat_info.__getattribute__(groupf_python_name)
                    # set chunking
                    feat_info_dataset.chunking = chunking
                    # and now use HDF5 pathing to write the chunking part back out
                    # in theory whis would work from any level, not just the root
                    relative_hdf5_dataset_path = "/".join(feat_info_dataset._hdf5_path.split("/")[-2:])
                    feat_info_dataset.write_simple_attributes(group_object[relative_hdf5_dataset_path])

    @property
    def feature_information(self) -> GroupFBase:
        return self._attributes[self.feature_information_attribute_name]

    @feature_information.setter
    def feature_information(self, val: GroupFBase):
        self._attributes[self.feature_information_attribute_name] = val

    @property
    @abstractmethod
    def feature_information_type(self):
        raise NotImplementedError()

    def feature_information_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_information = self.feature_information_type()

    @property
    def product_specification(self) -> str:
        return self._attributes[self.product_specification_attribute_name]

    @product_specification.setter
    def product_specification(self, val: str):
        self._attributes[self.product_specification_attribute_name] = val

    @property
    def product_specification_type(self) -> Type[str]:
        return str

    def product_specification_create(self):
        """ Creates a blank, empty or zero value for product_specification"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.product_specification = self.product_specification_type()

    @property
    def issue_time(self) -> datetime.time:
        return self._attributes[self.issue_time_attribute_name]

    @issue_time.setter
    def issue_time(self, val: Union[datetime.time, datetime.datetime, str]):
        self.set_datetime_attribute(val, self.issue_time_attribute_name, self.issue_time_type)

    @property
    def issue_time_type(self) -> Type[datetime.time]:
        return datetime.time

    def issue_time_create(self):
        """ Creates a blank, empty or zero value for issue_time"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.issue_time = self.issue_time_type(0, 0)  # midnight

    @property
    def issue_date(self) -> datetime.date:
        return self._attributes[self.issue_date_attribute_name]

    @issue_date.setter
    def issue_date(self, val: Union[datetime.date, datetime.datetime, str]):
        self.set_datetime_attribute(val, self.issue_date_attribute_name, self.issue_date_type)

    @property
    def issue_date_type(self) -> Type[datetime.date]:
        return datetime.date

    def issue_date_create(self):
        """ Creates a blank, empty or zero value for issue_date"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.issue_date = self.issue_date_type(1970, 1, 1)

    @property
    def horizontal_datum_reference(self) -> str:
        return self._attributes[self.horizontal_datum_reference_attribute_name]

    @horizontal_datum_reference.setter
    def horizontal_datum_reference(self, val: str):
        self._attributes[self.horizontal_datum_reference_attribute_name] = val

    @property
    def horizontal_datum_reference_type(self) -> Type[str]:
        return str

    def horizontal_datum_reference_create(self):
        """ Creates a blank, empty or zero value for horizontal_datum_reference"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_datum_reference = self.horizontal_datum_reference_type()

    @property
    def horizontal_datum_value(self) -> int:
        return self._attributes[self.horizontal_datum_value_attribute_name]

    @horizontal_datum_value.setter
    def horizontal_datum_value(self, val: Union[str, int]):
        val = int(val)
        self._attributes[self.horizontal_datum_value_attribute_name] = val

    @property
    def horizontal_datum_value_type(self) -> Type[int]:
        return int

    def horizontal_datum_value_create(self):
        """ Creates a blank, empty or zero value for horizontal_datum_value"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_datum_value = self.horizontal_datum_value_type()

    @property
    def epoch(self) -> str:
        return self._attributes[self.epoch_attribute_name]

    @epoch.setter
    def epoch(self, val: str):
        self._attributes[self.epoch_attribute_name] = val

    @property
    def epoch_type(self) -> Type[str]:
        return str

    def epoch_create(self):
        """ Creates a blank, empty or zero value for epoch"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.epoch = self.epoch_type()

    @property
    def geographic_identifier(self) -> str:
        return self._attributes[self.geographic_identifier_attribute_name]

    @geographic_identifier.setter
    def geographic_identifier(self, val: str):
        self._attributes[self.geographic_identifier_attribute_name] = val

    @property
    def geographic_identifier_type(self) -> Type[str]:
        return str

    def geographic_identifier_create(self):
        """ Creates a blank, empty or zero value for geographic_identifier"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.geographic_identifier = self.geographic_identifier_type()

    @property
    def metadata(self) -> str:
        return self._attributes[self.metadata_attribute_name]

    @metadata.setter
    def metadata(self, val: str):
        self._attributes[self.metadata_attribute_name] = val

    @property
    def metadata_type(self) -> Type[str]:
        return str

    def metadata_create(self):
        """ Creates a blank, empty or zero value for metadata"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.metadata = self.metadata_type()

    @property
    def vertical_datum(self) -> Enum:
        return self._attributes[self.vertical_datum_attribute_name]

    @vertical_datum.setter
    def vertical_datum(self, val: Union[int, str, VERTICAL_DATUM]):
        self.set_enum_attribute(val, self.vertical_datum_attribute_name, self.vertical_datum_type)
        # if isinstance(val, str):
        #     val = self.vertical_datum_type[val]
        # if isinstance(val , int):
        #     val = self.vertical_datum_type(val)
        # self._attributes[self.vertical_datum_attribute_name] = val

    @property
    def vertical_datum_type(self) -> Type[Enum]:
        return VERTICAL_DATUM

    def vertical_datum_create(self):
        """ Creates a blank, empty or zero value for vertical_datum"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum = self.vertical_datum_type["MLLW"]

    @property
    def meta_features(self) -> str:
        return self._attributes[self.meta_features_attribute_name]

    @meta_features.setter
    def meta_features(self, val: str):
        self._attributes[self.meta_features_attribute_name] = val

    @property
    def meta_features_type(self) -> Type[str]:
        return str

    def meta_features_create(self):
        """ Creates a blank, empty or zero value for meta_features"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.meta_features = self.meta_features_type()


class S100File(S1XXFile):
    PRODUCT_SPECIFICATION = numpy.string_('INT.IHO.S-100.4.0')

    def __init__(self, *args, **kywrds):
        # kywrds['root'] = S100Root
        super().__init__(*args, root=S100Root, **kywrds)
