import os
import subprocess
import datetime
import logging
import functools
from abc import ABC, abstractmethod
from typing import Callable, Iterator, Union, Optional, List, Type

import xml
import xml.etree.ElementTree
import pprint

import numpy
import h5py

try:
    from .. import s1xx
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py.s102"

from ..s1xx import s1xx_sequence, S1xxAttributesBase, S1xxMetadataListBase, S1xxGridsBase, S1XXFile, h5py_string_dtype
from ..s100 import GridCoordinate, DirectPosition, GeographicExtent, GridEnvelope, SequenceRule, VertexPoint, \
    FeatureInformation, FeatureInformationDataset, FeatureContainerDCF2, S100Root, S100Exception, FeatureInstanceDCF2, GroupFBase, \
    CommonPointRule


class S102Exception(S100Exception):
    pass


BATHY_COVERAGE = "BathymetryCoverage"
TRACKING_COVERAGE = "TrackingListCoverage"
DEPTH = "depth"
UNCERTAINTY = "uncertainty"

"""Contains s102 metadata to pass to S102File.

PRODUCT_SPECIFICATION: The product specification used to create this dataset.
HORIZONTAL_DATUM_REFERENCE: Reference to the register from which the horizontal datum value is taken.
DATA_CODING_FORMAT: Reference to the type of S102 product.
INTERPOLATION_TYPE: Interpolation method recommended for evaluation of the S100_GridCoverage.
COMMON_POINT_RULE: The procedure used for evaluating geometric objects that overlap or lie fall on boundaries.
DIMENSION: The dimension of the feature instance.
SEQUENCING_RULE_TYPE: Method to assign values from the sequence of values to the grid coordinates (e.g. "linear").
SEQUENCING_RULE_SCAN_DIRECTION: AxisNames, comma-separated (e.g. "longitude,latitude").
START_SEQUENCE: Starting location of the scan.

"""


# figure 4.4 in section 4.2.1 of v2.0.0 shows an overview of many of the S102 classes used
# Annex B in the S102 spec has an example layout (still determining if it matches the docs)
# S100 doc part 10C has HDF5 layout information as well.
# S100 doc 10C-7 has some representation guidelines


# override the basic S100 spec that says to use an underscore and use a dot instead
class S102MetadataListBase(S1xxMetadataListBase):
    write_format_str = ".%03d"


# # @TODO -- determine if this is old.  The spec seems to describe a one dimensional array or list of points but the values in the grid is a 2 x N x M dataset
# class BathymetryValueRecord(S1xxAttributesBase):
#     """ 4.2.1.1.2.2 and Figure 4.4 of v2.0.0
#     The attribute values has the value type S102_BathymetryValueRecord which is a sequence of value items that
#     shall assign values to the grid points.
#     There are two attributes in the bathymetry value record, depth and uncertainty in the S102_BathymetryValues class.
#     The definition for the depth is defined by the depthCorrectionType attribute in the S102_DataIdentification class.
#     The definition of the type of data in the values record is defined by the verticalUncertaintyType attribute in the
#     S102_DataIdentification class.
#     """
#     depth_attribute_name = "depth"  #: HDF5 naming
#     uncertainty_attribute_name = "uncertainty"  #: HDF5 naming
#
#     @property
#     def __version__(self) -> int:
#         return 1
#
#     @property
#     def __version__(self) -> int:
#         return 1
#
#
#     @property
#     def depth_type(self):
#         return numpy.ndarray
#
#     def depth_create(self):
#         self.depth = self.depth_type([], numpy.float)
#
#     @property
#     def depth(self) -> float:
#         return self._attributes[self.depth_attribute_name]
#
#     @depth.setter
#     def depth(self, val: float):
#         self._attributes[self.depth_attribute_name] = val
#
#
#     @property
#     def uncertainty_type(self):
#         return numpy.ndarray
#
#     def uncertainty_create(self):
#         self.uncertainty = self.uncertainty_type([], numpy.float)
#
#     @property
#     def uncertainty(self) -> float:
#         return self._attributes[self.uncertainty_attribute_name]
#
#     @uncertainty.setter
#     def uncertainty(self, val: float):
#         self._attributes[self.uncertainty_attribute_name] = val
#
#
# class BathymetryValuesList(S102MetadataListBase):
#     """ 4.2.1.1.2 and Figure 4.4 of v2.0.0
#     The class S102_BathymetryValues is related to BathymetryCoverage by a composition relationship in which
#     an ordered sequence of depth values provide data values for each grid cell.
#     The class S102_BathymetryValues inherits from S100_Grid.
#     """
#
#     @property
#     def __version__(self) -> int:
#         return 1
#
#     @property
#     def metadata_name(self) -> str:
#         return "values"
#
#     @property
#     def metadata_type(self) -> type:
#         return BathymetryValueRecord


class BathymetryValues(S1xxGridsBase):
    depth_attribute_name = "depth"  #: HDF5 naming
    uncertainty_attribute_name = "uncertainty"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "values"

    @property
    def depth(self) -> s1xx_sequence:
        """This is the depth array.  For bathymetric gridded data, the dataset includes a two-dimensional array containing both the depth and uncertainty data.
        These dimensions are defined by numPointsLongitudinal and numPointsLatitudinal.
        By knowing the grid origin and the grid spacing, the position of every point in the grid can be computed by simple formulae"""
        return self._attributes[self.depth_attribute_name]

    @depth.setter
    def depth(self, val: s1xx_sequence):
        self._attributes[self.depth_attribute_name] = val

    @property
    def depth_type(self) -> s1xx_sequence:
        return numpy.float32

    def depth_create(self):
        """ Creates a blank, empty or zero value for depth"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.depth = self.depth_type([], numpy.float)

    @property
    def uncertainty(self) -> s1xx_sequence:
        """This is the uncertainty array.  For bathymetric gridded data, the dataset includes a two-dimensional array containing both the depth and uncertainty data.
        These dimensions are defined by numPointsLongitudinal and numPointsLatitudinal.
        By knowing the grid origin and the grid spacing, the position of every point in the grid can be computed by simple formulae"""
        return self._attributes[self.uncertainty_attribute_name]

    @uncertainty.setter
    def uncertainty(self, val: s1xx_sequence):
        self._attributes[self.uncertainty_attribute_name] = val

    @property
    def uncertainty_type(self) -> s1xx_sequence:
        return numpy.float32

    def uncertainty_create(self):
        """ Creates a blank, empty or zero value for uncertainty"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.uncertainty = self.uncertainty_type([], numpy.float)

    def get_write_order(self):
        return [self.depth_attribute_name, self.uncertainty_attribute_name]

    def get_compound_dtype(self):
        return [self.depth_type, self.uncertainty_type]


class BathymetryCoverage(S1xxAttributesBase):
    """ This is the Group.NNN object that contains the grid data in a values dataset and other metadata about the grids.

    4.2.1.1.1 and Figure 4.4 of v2.0.0
    also see section 12.3 and table 12.5

    """

    write_format_str = ".%03d"

    values_attribute_name = "values"  #: HDF5 naming
    minimum_depth_attribute_name = "minimumDepth"  #: HDF5 naming
    maximum_depth_attribute_name = "maximumDepth"  #: HDF5 naming
    maximum_display_scale_attribute_name = "maximumDisplayScale"  #: HDF5 naming
    minimum_display_scale_attribute_name = "minimumDisplayScale"  #: HDF5 naming
    minimum_uncertainty_attribute_name = "minimumUncertainty"  #: HDF5 naming
    maximum_uncertainty_attribute_name = "maximumUncertainty"  #: HDF5 naming
    origin_attribute_name = "origin"  #: HDF5 naming
    offset_vectors_attribute_name = "offsetVectors"  #: HDF5 naming
    dimension_attribute_name = "dimension"  #: HDF5 naming
    axis_names_attribute_name = "axisNames"  #: HDF5 naming
    extent_attribute_name = "extent"  #: HDF5 naming
    sequencing_rule_attribute_name = "sequencingRule"  #: HDF5 naming
    start_sequence_attribute_name = "startSequence"  #: HDF5 naming

    @property
    def values(self) -> BathymetryValues:
        """ The grids for depth and uncertainty.

        4.2.1.1.2.1 S102_BathymetryValues semantics

        The class S102_BathymetryValues is related to BathymetryCoverage by a composition relationship in which an ordered sequence
        of depth values provide data values for each grid cell.
        The class S102_BathymetryValues inherits from S100_Grid

        4.2.1.1.2.2 values

        The attribute values has the value type S102_BathymetryValueRecord which is
         a sequence of value items that shall assign values to the grid points.
        There are two attributes in the bathymetry value record, depth and uncertainty in the S102_BathymetryValues class.
        The definition for the depth is defined by the depthCorrectionType attribute in the S102_DataIdentification class.
        The definition of the type of data in the values record is defined by the verticalUncertaintyType attribute
        in the S102_DataIdentification class
        """
        return self._attributes[self.values_attribute_name]

    @values.setter
    def values(self, val: BathymetryValues):
        self._attributes[self.values_attribute_name] = val

    @property
    def values_type(self) -> Type[BathymetryValues]:
        return BathymetryValues

    def values_create(self):
        """ Creates a blank, empty or zero value for values"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.values = self.values_type()

    @property
    def __version__(self) -> int:
        return 1

    @property
    def minimum_depth_type(self):
        return float

    def minimum_depth_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.minimum_depth = self.minimum_depth_type()

    @property
    def minimum_depth(self) -> float:
        """From 4.2.1.1.1.3,
        The attribute minimumDepth has the value type Real and describes the lower bound of the depth estimate
        for all the depth values in S102_BathymetryValues record.
        This attribute is required. There is no default"""
        return self._attributes[self.minimum_depth_attribute_name]

    @minimum_depth.setter
    def minimum_depth(self, val: float):
        self._attributes[self.minimum_depth_attribute_name] = val

    @property
    def maximum_depth_type(self):
        return float

    def maximum_depth_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.maximum_depth = self.maximum_depth_type()

    @property
    def maximum_depth(self) -> float:
        """From 4.2.1.1.1.4,
        The attribute minimumDepth has the value type Real and describes the lower bound of the depth estimate
        for all the depth values in S102_BathymetryValues record.
        This attribute is required. There is no default
        """
        return self._attributes[self.maximum_depth_attribute_name]

    @maximum_depth.setter
    def maximum_depth(self, val: float):
        self._attributes[self.maximum_depth_attribute_name] = val

    @property
    def minimum_display_scale(self) -> int:
        """From 4.2.1.1.1.2,
        The larger value of the ratio of the linear dimensions of the features of a dataset presented in the display and
        the actual dimensions of the features represented (largest scale) of the scale range of the dataset.
        A list of display scale ranges is available in Figure 11.1, 1st column
        """
        return self._attributes[self.minimum_display_scale_attribute_name]

    @minimum_display_scale.setter
    def minimum_display_scale(self, val: int):
        self._attributes[self.minimum_display_scale_attribute_name] = val

    @property
    def maximum_display_scale_type(self):
        return float

    def maximum_display_scale_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.maximum_display_scale = self.maximum_display_scale_type()

    @property
    def maximum_display_scale(self) -> int:
        """From 4.2.1.1.1.5,
        The smaller value of the ratio of the linear dimensions of the features of a dataset presented in the display and
        the actual dimensions of the features represented (smallest scale) of the scale range of the dataset.
        A list of display scale ranges is available in Table 11.1, 1st column
        """
        return self._attributes[self.maximum_display_scale_attribute_name]

    @maximum_display_scale.setter
    def maximimum_display_scale(self, val: int):
        self._attributes[self.maximum_display_scale_attribute_name] = val

    @property
    def minimum_display_scale_type(self):
        return float

    def minimum_display_scale_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.minimum_display_scale = self.minimum_display_scale_type()

    @property
    def minimum_uncertainty(self) -> float:
        """From 4.2.1.1.1.6,
        The attribute minimumUncertainty has the value type Real and describes the lower bound of the uncertainty of the
        depth estimate for all the depth values in S102_BathymetryValues record.
        This attribute is required. There is no default
        """
        return self._attributes[self.minimum_uncertainty_attribute_name]

    @minimum_uncertainty.setter
    def minimum_uncertainty(self, val: float):
        self._attributes[self.minimum_uncertainty_attribute_name] = val

    @property
    def minimum_uncertainty_type(self):
        return float

    def minimum_uncertainty_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.minimum_uncertainty = self.minimum_uncertainty_type()

    @property
    def maximum_uncertainty(self) -> float:
        """From 4.2.1.1.1.7,
        The attribute minimumUncertainty has the value type Real and describes the lower bound of the uncertainty
        of the depth estimate for all the depth values in S102_BathymetryValues record.
        This attribute is required. There is no default"""
        return self._attributes[self.maximum_uncertainty_attribute_name]

    @maximum_uncertainty.setter
    def maximum_uncertainty(self, val: float):
        self._attributes[self.maximum_uncertainty_attribute_name] = val

    @property
    def maximum_uncertainty_type(self):
        return float

    def maximum_uncertainty_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.maximum_uncertainty = self.maximum_uncertainty_type()

    @property
    def origin(self) -> DirectPosition:
        """From 4.2.1.1.1.8,
        The attribute origin has the value class DirectPosition which is a position that shall locate the origin of the rectified grid
        in the coordinate reference system.
        This attribute is required. There is no default
        """
        return self._attributes[self.origin_attribute_name]

    @origin.setter
    def origin(self, val: DirectPosition):
        self._attributes[self.origin_attribute_name] = val

    @property
    def origin_type(self):
        return DirectPosition

    def origin_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.origin = self.origin_type()

    @property
    def origin_attribute_type(self) -> Type[DirectPosition]:
        return DirectPosition

    @property
    def offset_vectors(self) -> s1xx_sequence:
        """sequence of s102 Vectors  From 4.2.1.1.1.9 in S102 v2.0.0
        The attribute offsetVectors has the value class Sequence<Vector> that shall be a sequence of offset vector elements
        that determine the grid spacing in each direction.
        The data type Vector is specified in ISO/TS 19103. This attribute is required.
        There is no default.
        """
        return self._attributes[self.offset_vectors_attribute_name]

    @offset_vectors.setter
    def offset_vectors(self, val: s1xx_sequence):
        self._attributes[self.offset_vectors_attribute_name] = val

    @property
    def offset_vectors_type(self):
        return numpy.ndarray

    def offset_vectors_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.offset_vectors = self.offset_vectors_type([2, ], numpy.float64)

    @property
    def dimension(self) -> int:
        """From 4.2.1.1.1.10,
        The attribute dimension has the value class Integer that shall identify the dimensionality of the grid.
        The value of the grid dimension in this product specification is 2.
        This value is fixed in this Product Specification and does not need to be encoded
        """
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
        self.dimension = self.dimension_type(2)

    @property
    def axis_names(self) -> s1xx_sequence:
        """sequence of character strings From 4.2.1.1.1.11,
        The attribute axisNames has the value class Sequence<CharacterString> that shall be used to assign names to the grid axis.
        The grid axis names shall be "Latitude" and "Longitude" for unprojected data sets or “Northing” and “Easting” in a projected space
        """
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
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.axis_names = numpy.array(["", ""], dtype=h5py_string_dtype)

    @property
    def extent(self) -> GridEnvelope:
        """From 4.2.1.1.1.12,
        The attribute extent has the value class CV_GridEnvelope that shall contain the extent of the spatial domain of the coverage.
        It uses the value class CV_GridEnvelope which provides the grid coordinate values for the diametrically opposed corners of the grid.
        The default is that this value is derived from the bounding box for the data set or tile in a multi tile data set"""
        return self._attributes[self.extent_attribute_name]

    @extent.setter
    def extent(self, val: GridEnvelope):
        self._attributes[self.extent_attribute_name] = val

    @property
    def extent_type(self) -> Type[GridEnvelope]:
        return GridEnvelope

    def extent_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.extent = self.extent_type()

    @property
    def sequencing_rule(self) -> SequenceRule:
        """From 4.2.1.1.1.13,
        The attribute sequencingRule has the value class CV_SequenceRule (ISO 19123) that shall
        describe how the grid points are ordered for association to the elements of the sequence values.
        The default value is "Linear".
        No other options are allowed.
        (note that for S100: Only the values "linear" (for a simple regular cell size grid) and "Morton" (for a
        Quad Tree Grid) shall be used for data that conforms to this standard.)
        """
        return self._attributes[self.sequencing_rule_attribute_name]

    @sequencing_rule.setter
    def sequencing_rule(self, val: SequenceRule):
        self._attributes[self.sequencing_rule_attribute_name] = val

    @property
    def sequencing_rule_type(self):
        return SequenceRule

    def sequencing_rule_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.sequencing_rule = self.sequencing_rule_type()

    @property
    def start_sequence(self) -> GridCoordinate:
        """ From4.2.1.1.1.14,
        The attribute startSequence has the value class CV_GridCoordinate that shall identify the grid
        point to be associated with the first record in the values sequence.
        The default value is the lower left corner of the grid.
        No other options are allowed.

        Returns
        -------

        """
        return self._attributes[self.start_sequence_attribute_name]

    @start_sequence.setter
    def start_sequence(self, val: GridCoordinate):
        self._attributes[self.start_sequence_attribute_name] = val

    @property
    def start_sequence_type(self):
        return GridCoordinate

    def start_sequence_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.start_sequence = self.start_sequence_type()


class SurfaceCorrectionValues(VertexPoint):
    pass


# this S102_SurfaceCorrectionValues may not be right as the docs refer to S100_VertexPoint,
# but I'm betting they would both change if there ever is an update
class TrackingListValues(SurfaceCorrectionValues):
    """ From 4.2.1.1.10 of S102 v2.0.0
    """
    track_code_attribute_name = "trackCode"  #: HDF5 naming
    list_series_attribute_name = "listSeries"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def track_code(self) -> str:
        """ From 4.2.1.1.10.2,
        The optional attribute trackCode has the value type CharacterString which may contain a text string
        describing the reason for the override of the corresponding depth and uncertainty values in the bathymetry coverage.
        This is a user definable field with values defined in the lineage metadata.

        Returns
        -------

        """
        return self._attributes[self.track_code_attribute_name]

    @track_code.setter
    def track_code(self, val: str):
        self._attributes[self.track_code_attribute_name] = val

    @property
    def track_code_type(self):
        return str

    def track_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.track_code = self.track_code_type()

    @property
    def list_series(self) -> int:
        """From 4.2.1.1.10.3,
        The attribute listSeries has the value type Integer which contains an index number into a list of metadata
        elements describing the reason for the override of the corresponding depth and uncertainty values in the bathymetry coverage.
        """
        return self._attributes[self.list_series_attribute_name]

    @list_series.setter
    def list_series(self, val: int):
        self._attributes[self.list_series_attribute_name] = val

    @property
    def list_series_type(self):
        return int

    def list_series_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.list_series = self.list_series_type()


class TrackingListValuesList(S102MetadataListBase):
    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "point"

    @property
    def metadata_type(self) -> type:
        return TrackingListValues


class TrackingListSetList(S102MetadataListBase):
    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "set"

    @property
    def metadata_type(self) -> type:
        return TrackingListValuesList


class TrackingListCoverage(CommonPointRule, S1xxAttributesBase):
    """ 4.2.1.1.9 and Figure 4.4 of v2.0.0
    commonPointRule is defined to be an S100_PointCoverage with a value of default and it therefore optional.
    a metadata attribute from S100 is allowed but not necessary as well.

    """
    write_format_str = ".%02d"

    domain_extent_attribute_name = "domainExtent"  #: HDF5 naming
    common_point_rule_attribute_name = "commonPointRule"  #: HDF5 naming
    set_attribute_name = "set"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def domain_extent(self) -> S102MetadataListBase:
        return self._attributes[self.domain_extent_attribute_name]

    @domain_extent.setter
    def domain_extent(self, val: S102MetadataListBase):
        self._attributes[self.domain_extent_attribute_name] = val

    @property
    def domain_extent_type(self):
        return GeographicExtent

    def domain_extent_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.domain_extent = self.domain_extent_type()

    def common_point_rule_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.common_point_rule = self.common_point_rule_type["average"]

    @property
    def set_type(self):
        return TrackingListSetList

    def set_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.set = self.set_type()

    @property
    def set(self) -> S102MetadataListBase:
        """
        Returns
        -------
        TrackingListValuesList
            list of TrackingListValues
        """
        return self._attributes[self.set_attribute_name]

    @set.setter
    def set(self, val: S102MetadataListBase):
        self._attributes[self.set_attribute_name] = val

    # @TODO  I don't think this is right, but not sure where I found it
    #
    # geometry_attribute_name = "geometry"  #: HDF5 naming
    #
    # @property
    # def geometry(self) -> S1xxAttributesBase:
    #     return self._attributes[self.geometry_attribute_name]
    #
    # @geometry.setter
    # def geometry(self, val: S1xxAttributesBase):
    #     self._attributes[self.geometry_attribute_name] = val
    #
    # value_attribute_name = "value"  #: HDF5 naming
    #
    # @property
    # def value(self) -> s1xx_sequence:
    #     return self._attributes[self.value_attribute_name]
    #
    # @value.setter
    # def value(self, val: s1xx_sequence):
    #     self._attributes[self.value_attribute_name] = val


class TrackingListGroupList(S102MetadataListBase):
    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "Group"

    @property
    def metadata_type(self) -> type:
        return TrackingListCoverage


class BathymetryGroupList(S102MetadataListBase):
    """ This is the list of Group.NNN that are held as a list.
    Each Group.NNN has a dataset of depth and uncertainty.
    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "Group"

    @property
    def metadata_type(self) -> type:
        return BathymetryCoverage


class TrackingListCoveragesList(S102MetadataListBase):
    """ 4.2.1.1.9 and Figure 4.4 and Table 10.1 of v2.0.0
    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return TRACKING_COVERAGE

    @property
    def metadata_type(self) -> type:
        return TrackingListGroupList


class BathymetryFeatureInstance(FeatureInstanceDCF2):
    """ This will be the BathymetryCoverage.001 element in HDF5.
    It will contain a Group.NNN which will have the "values" dataset of the deptha dn uncertainty.
    """
    bathymetry_group_attribute_name = "Group" + r"[\._]\d+"
    """ Basic template for HDF5 naming of the attribute.  
    Attribute name will be automatically determined based on the list's index of the data. 
    """

    @property
    def bathymetry_group_type(self):
        return BathymetryGroupList

    def bathymetry_group_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.bathymetry_group = self.bathymetry_group_type()

    @property
    def bathymetry_group(self) -> S102MetadataListBase:
        """ The bathymetry data, a list of Bathymetrygroup
        Returns
        -------
        S102MetadataListBase
            Contains a list of BathymetryCoverage objects via the BathymetryCoveragesList class
        """
        return self._attributes[self.bathymetry_group_attribute_name]

    @bathymetry_group.setter
    def bathymetry_group(self, val: S102MetadataListBase):
        self._attributes[self.bathymetry_group_attribute_name] = val


class BathymetryCoveragesList(S102MetadataListBase):
    """ 4.2.1.1.2 and Figure 4.4 and Table 10.1 of v2.0.0
    This is the set of BathymetryCoverage.NN that act like a list here.
    They will contain a list of Groups.NNN as well as other attributes etc.
    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return BATHY_COVERAGE

    @property
    def metadata_type(self) -> Type[BathymetryFeatureInstance]:
        return BathymetryFeatureInstance


class BathymetryContainer(FeatureContainerDCF2):
    """ This is the BathymetryCoverage right off the root of the HDF5 which has possible attributes from S100 spec table 10c-10
    This will hold child groups named BathymetryCoverage.NN
    """
    #: attribute name will be automatically determined based on the containing list's index
    bathymetry_coverage_attribute_name = BATHY_COVERAGE + r"[\._]\d+"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def bathymetry_coverage_type(self):
        return BathymetryCoveragesList

    def bathymetry_coverage_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.bathymetry_coverage = self.bathymetry_coverage_type()

    @property
    def bathymetry_coverage(self) -> S102MetadataListBase:
        """ The bathymetry data, a list of BathymetryCoverage

        Returns
        -------
        S102MetadataListBase
            Contains a list of BathymetryCoverage objects via the BathymetryCoveragesList class
        """
        return self._attributes[self.bathymetry_coverage_attribute_name]

    @bathymetry_coverage.setter
    def bathymetry_coverage(self, val: S102MetadataListBase):
        self._attributes[self.bathymetry_coverage_attribute_name] = val

    def data_coding_format_create(self):
        """ Creates a blank, empty or zero value for data_coding_format"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_coding_format = self.data_coding_format_type(2)  # regular grid

    def dimension_create(self):
        """ Creates a blank, empty or zero value for dimension"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.dimension = self.dimension_type(2)


class TrackingListContainer(FeatureContainerDCF2):
    """
    Table 10.1 of v2.0.0
    """

    tracking_list_coverage_attribute_name = TRACKING_COVERAGE + r"[\._]\d+"

    @property
    def __version__(self) -> int:
        return 1

    def data_coding_format_create(self):
        """ Creates a blank, empty or zero value for data_coding_format"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_coding_format = self.data_coding_format_type(1)  # point set

    @property
    def tracking_list_coverage_type(self):
        return TrackingListCoveragesList

    def tracking_list_coverage_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.tracking_list_coverage = self.tracking_list_coverage_type()

    @property
    def tracking_list_coverage(self) -> S1xxAttributesBase:
        """ The tracking list data, a list of TrackingListCoverage
        Returns
        -------
        S102MetadataListBase
            Contains a list of TrackingListCoverage objects via the TrackingListCoveragesList class
        """
        return self._attributes[self.tracking_list_coverage_attribute_name]

    @tracking_list_coverage.setter
    def tracking_list_coverage(self, val: S1xxAttributesBase):
        self._attributes[self.tracking_list_coverage_attribute_name] = val


class S102FeatureInformation(FeatureInformation):
    """ S102 specifc version of FeatureInformation.
    Sets defaults of uom.name to metres, fillValue to 1000000, upper and lower to 12000, -12000 and closure to closedInterval
    and datatype to H5T_NATIVE_FLOAT.
    The user should set code and name to 'depth' or 'uncertainty' as needed.
    """

    def __init__(self, *args, **kwrds):
        super().__init__(*args, **kwrds)
        self.datatype_create()  # make this first so anyone who tries to set values before the datatype isn't surprised by an exception
        self.initialize_properties()  # initialize everything since there aren't any choices except name and code
        del self.name
        del self.code

    def unit_of_measure_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.unit_of_measure = self.unit_of_measure_type("metres")

    def fill_value_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.fill_value = self.fill_value_type(1000000)

    def datatype_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.datatype = self.datatype_type("H5T_NATIVE_FLOAT")

    def lower_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.lower = self.lower_type(-12000)

    def upper_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.upper = self.upper_type(12000)

    def closure_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.closure = self.closure_type("closedInterval")


class S102FeatureInformationDataset(FeatureInformationDataset, ABC):
    """   In S102, 10.2.1 and table 10.2 and Table 10.1 of v2.0.0

    This is used to describe the BathymetryCoverage and TrackingListCoverage within the GroupF feature listing.
    The features described under GroupF have a matching named entry parallel to GroupF (top level).
    The actual data (depths etc) is stored in the top level element while basic metadata is stored in this element.

    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_type(self) -> Type[S102FeatureInformation]:
        return S102FeatureInformation


class TrackingListCoverageDataset(S102FeatureInformationDataset):
    @property
    def metadata_name(self) -> str:
        return TRACKING_COVERAGE


class BathymetryCoverageDataset(S102FeatureInformationDataset):
    @property
    def metadata_name(self) -> str:
        return BATHY_COVERAGE


class FeatureCodes(GroupFBase):
    """ Table 10.1 and sect 10.2.1 of v2.0.0
    """

    feature_name_attribute_name = "featureName"  #: HDF5 naming
    bathymetry_coverage_dataset_attribute_name = BATHY_COVERAGE
    tracking_list_coverage_attribute_name = TRACKING_COVERAGE

    @property
    def __version__(self) -> int:
        return 1

    @property
    def feature_name_type(self):
        return numpy.array

    def feature_name_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_name = self.feature_name_type([BATHY_COVERAGE, TRACKING_COVERAGE], dtype=h5py_string_dtype)

    @property
    def feature_name(self) -> s1xx_sequence:
        return self._attributes[self.feature_name_attribute_name]

    @feature_name.setter
    def feature_name(self, val: s1xx_sequence):
        self._attributes[self.feature_name_attribute_name] = val

    def feature_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_code = self.feature_code_type([BATHY_COVERAGE, TRACKING_COVERAGE], dtype=h5py_string_dtype)

    @property
    def bathymetry_coverage_dataset_type(self):
        return BathymetryCoverageDataset

    def bathymetry_coverage_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.bathymetry_coverage_dataset = self.bathymetry_coverage_dataset_type()

    @property
    def bathymetry_coverage_dataset(self) -> BathymetryCoverageDataset:
        return self._attributes[self.bathymetry_coverage_dataset_attribute_name]

    @bathymetry_coverage_dataset.setter
    def bathymetry_coverage_dataset(self, val: BathymetryCoverageDataset):
        self._attributes[self.bathymetry_coverage_dataset_attribute_name] = val

    @property
    def tracking_list_coverage_type(self):
        return TrackingListCoverageDataset

    def tracking_list_coverage_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.tracking_list_coverage = self.tracking_list_coverage_type()

    @property
    def tracking_list_coverage(self) -> TrackingListCoverageDataset:
        return self._attributes[self.tracking_list_coverage_attribute_name]

    @tracking_list_coverage.setter
    def tracking_list_coverage(self, val: TrackingListCoverageDataset):
        self._attributes[self.tracking_list_coverage_attribute_name] = val


class S102Root(S100Root):
    """The root group contains a feature information group and N feature containers.
    In S102 there are currently two feature containers which are the 'coverages'  bathymetry and tracking list.
    The coverage names are determined from the matching CoveragesAttributes
    10.2 and Figure 10.1 of v2.0.0
    """
    feature_information_attribute_name = "Group_F"  #: HDF5 naming
    bathymetry_coverage_attribute_name = BATHY_COVERAGE
    tracking_list_coverage_attribute_name = TRACKING_COVERAGE

    @property
    def __version__(self) -> int:
        return 1

    @property
    def feature_information(self) -> FeatureCodes:
        """Feature Information stored in GroupF in the HDF5 using :class:`FeatureCodes`"""
        return self._attributes[self.feature_information_attribute_name]

    @feature_information.setter
    def feature_information(self, val: FeatureCodes):
        self._attributes[self.feature_information_attribute_name] = val

    @property
    def feature_information_type(self):
        return FeatureCodes

    def feature_information_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_information = self.feature_information_type()

    @property
    def bathymetry_coverage(self) -> S1xxAttributesBase:
        """Bathymetry instance stored under the HDF5 root using :class:`BathymetryContainer`
        """
        return self._attributes[self.bathymetry_coverage_attribute_name]

    @property
    def bathymetry_coverage_type(self):
        return BathymetryContainer

    def bathymetry_coverage_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.bathymetry_coverage = self.bathymetry_coverage_type()

    @bathymetry_coverage.setter
    def bathymetry_coverage(self, val: S1xxAttributesBase):
        self._attributes[self.bathymetry_coverage_attribute_name] = val

    @property
    def tracking_list_coverage_type(self):
        return TrackingListContainer

    def tracking_list_coverage_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.tracking_list_coverage = self.tracking_list_coverage_type()

    @property
    def tracking_list_coverage(self) -> S1xxAttributesBase:
        return self._attributes[self.tracking_list_coverage_attribute_name]

    @tracking_list_coverage.setter
    def tracking_list_coverage(self, val: S1xxAttributesBase):
        self._attributes[self.tracking_list_coverage_attribute_name] = val


class S102File(S1XXFile):
    PRODUCT_SPECIFICATION = numpy.string_('INT.IHO.S-102.2.0')
    # these keys allow backward compatibility with NAVO data, the first key is current at time of writing
    top_level_keys = ('BathymetryCoverage', 'S102_Grid', 'S102_BathymetryCoverage')
    tracking_list_top_level = ("TrackingListCoverage",)
    tracking_list_second_level = ("TrackingListCoverage.01",)
    tracking_list_group_level = ("Group.001",)
    second_level_keys = (
        'BathymetryCoverage.01', 'S102_Grid.01', 'S102_BathymetryCoverage.01', 'BathymetryCoverage_01', 'S102_Grid_01', 'S102_BathymetryCoverage_01',)
    group_level_keys = ('Group.001', 'Group_001',)
    value_level_keys = ("values",)
    depth_keys = ("depth", "depths", 'elevation', "elevations", "S102_Elevation")

    def __init__(self, *args, **kywrds):
        super().__init__(*args, root=S102Root, **kywrds)

    def print_overview(self, display_nodes=10):
        depths = self.get_depths()
        print("shape of grid is", depths.shape, "of type", depths.dtype)
        with numpy.printoptions(precision=2, suppress=True, linewidth=200):
            x, y = depths.shape
            r = max(x, y)
            step = int(r / display_nodes)
            print(depths[::step, ::step])

    def print_depth_attributes(self):
        hdf5 = self.get_depth_dataset()
        print(hdf5.attrs)

    def get_depth_dataset(self):
        for k in self.top_level_keys:
            if k in self:
                d = self[k]
                break
        try:
            d
        except NameError:
            raise KeyError(str(self.top_level_keys) + " were not found in " + str(list(self.keys())))

        for k in self.second_level_keys:
            if k in d:
                g = d[k]
                break

        try:
            g
        except NameError:
            raise KeyError(str(self.second_level_keys) + " were not found in " + str(list(d.keys())))

        for k in self.group_level_keys:
            if k in g:
                gp = g[k]
                break

        try:
            gp
        except NameError:
            raise KeyError(str(self.group_level_keys) + " were not found in " + str(list(g.keys())))

        for k in self.value_level_keys:
            if k in gp:
                v = gp[k]
                break
        try:
            v
            return v
        except NameError:
            raise KeyError(str(self.value_level_keys) + " were not found in " + str(list(gp.keys())))

    def get_depths(self):
        v = self.get_depth_dataset()
        # v.dtype
        # dtype([('S102_Elevation', '<f4'), ('S102_Uncertainty', '<f4')])
        for k in self.depth_keys:
            if k in v.dtype.names:
                return v[k]
        raise KeyError(str(self.depth_keys) + " were not found in " + str(list(v.dtype.names)))
