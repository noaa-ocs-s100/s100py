from __future__ import annotations

import os
import pathlib
import pathlib
import shutil
import sys
import subprocess
import datetime
import tempfile
import warnings
import logging
import functools
from abc import ABC, abstractmethod
from typing import Callable, Iterator, Union, Optional, List, Type

import xml
import xml.etree.ElementTree
import pprint
from xml.etree import ElementTree as et

import numpy
import h5py
try:
    from osgeo import gdal, osr
except:
    pass

try:
    from ... import s1xx
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py.s102"

from ...s1xx import s1xx_sequence, S1xxObject, S1xxCollection, S1xxGridsBase, S1XXFile, h5py_string_dtype, h5py_string_comp
from ...s100.v4_0.api import S100File, GridCoordinate, DirectPosition, GeographicExtent, GridEnvelope, SequenceRule, VertexPoint, \
    FeatureInformation, FeatureInformationDataset, FeatureContainerDCF2, S100Root, S100Exception, FeatureInstanceDCF2, GroupFBase, \
    CommonPointRule

EDITION = 2.0
PRODUCT_SPECIFICATION = 'INT.IHO.S-102.2.0'

class S102Exception(S100Exception):
    pass


BATHY_COVERAGE = "BathymetryCoverage"
TRACKING_COVERAGE = "TrackingListCoverage"
DEPTH = "depth"
UNCERTAINTY = "uncertainty"

gco = "{http://www.isotc211.org/2005/gco}"

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
class S102MetadataListBase(S1xxCollection):
    write_format_str = ".%03d"


# # @TODO -- determine if this is old.  The spec seems to describe a one dimensional array or list of points but the values in the grid is a 2 x N x M dataset
# class BathymetryValueRecord(S1xxObject):
#     """ 4.2.1.1.2.2 and Figure 4.4 of v2.0.0
#     The attribute values has the value type S102_BathymetryValueRecord which is a sequence of value items that
#     shall assign values to the grid points.
#     There are two attributes in the bathymetry value record, depth and uncertainty in the S102_BathymetryValues class.
#     The definition for the depth is defined by the depthCorrectionType attribute in the S102_DataIdentification class.
#     The definition of the type of data in the values record is defined by the verticalUncertaintyType attribute in the
#     S102_DataIdentification class.
#     """
#     __depth_hdf_name__ = "depth"  #: HDF5 naming
#     __uncertainty_hdf_name__ = "uncertainty"  #: HDF5 naming
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
#     def __depth_type__(self):
#         return numpy.ndarray
#
#     def depth_create(self):
#         self.depth = self.__depth_type__([], numpy.float64)
#
#     @property
#     def depth(self) -> float:
#         return self._attributes[self.__depth_hdf_name__]
#
#     @depth.setter
#     def depth(self, val: float):
#         self._attributes[self.__depth_hdf_name__] = val
#
#
#     @property
#     def __uncertainty_type__(self):
#         return numpy.ndarray
#
#     def uncertainty_create(self):
#         self.uncertainty = self.__uncertainty_type__([], numpy.float64)
#
#     @property
#     def uncertainty(self) -> float:
#         return self._attributes[self.__uncertainty_hdf_name__]
#
#     @uncertainty.setter
#     def uncertainty(self, val: float):
#         self._attributes[self.__uncertainty_hdf_name__] = val
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
    __depth_hdf_name__ = "depth"  #: HDF5 naming
    __uncertainty_hdf_name__ = "uncertainty"  #: HDF5 naming

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
        return self._attributes[self.__depth_hdf_name__]

    @depth.setter
    def depth(self, val: s1xx_sequence):
        self._attributes[self.__depth_hdf_name__] = val

    @property
    def __depth_type__(self) -> s1xx_sequence:
        return numpy.ndarray

    @property
    def depth_dtype(self) -> Type[float]:
        return numpy.float32

    def depth_create(self):
        """ Creates a blank, empty or zero value for depth"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.depth = self.__depth_type__([], self.depth_dtype)

    @property
    def uncertainty(self) -> s1xx_sequence:
        """This is the uncertainty array.  For bathymetric gridded data, the dataset includes a two-dimensional array containing both the depth and uncertainty data.
        These dimensions are defined by numPointsLongitudinal and numPointsLatitudinal.
        By knowing the grid origin and the grid spacing, the position of every point in the grid can be computed by simple formulae"""
        return self._attributes[self.__uncertainty_hdf_name__]

    @uncertainty.setter
    def uncertainty(self, val: s1xx_sequence):
        self._attributes[self.__uncertainty_hdf_name__] = val

    @property
    def __uncertainty_type__(self) -> s1xx_sequence:
        return numpy.ndarray

    @property
    def uncertainty_dtype(self) -> Type[float]:
        return numpy.float32

    def uncertainty_create(self):
        """ Creates a blank, empty or zero value for uncertainty"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.uncertainty = self.__uncertainty_type__([], self.uncertainty_dtype)

    def get_write_order(self):
        return [self.__depth_hdf_name__, self.__uncertainty_hdf_name__]

    def get_compound_dtype(self):
        return [self.depth_dtype, self.uncertainty_dtype]


class BathymetryCoverageBase(S1xxObject):
    """ This is the Group.NNN object that contains the grid data in a values dataset and other metadata about the grids.

    4.2.1.1.1 and Figure 4.4 of v2.0.0
    also see section 12.3 and table 12.5

    """

    write_format_str = ".%03d"

    __values_hdf_name__ = "values"  #: HDF5 naming
    __minimum_depth_hdf_name__ = "minimumDepth"  #: HDF5 naming
    __maximum_depth_hdf_name__ = "maximumDepth"  #: HDF5 naming
    __minimum_uncertainty_hdf_name__ = "minimumUncertainty"  #: HDF5 naming
    __maximum_uncertainty_hdf_name__ = "maximumUncertainty"  #: HDF5 naming
    __origin_hdf_name__ = "origin"  #: HDF5 naming
    __offset_vectors_hdf_name__ = "offsetVectors"  #: HDF5 naming
    __dimension_hdf_name__ = "dimension"  #: HDF5 naming
    __axis_names_hdf_name__ = "axisNames"  #: HDF5 naming
    __extent_hdf_name__ = "extent"  #: HDF5 naming
    __sequencing_rule_hdf_name__ = "sequencingRule"  #: HDF5 naming
    __start_sequence_hdf_name__ = "startSequence"  #: HDF5 naming

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
        return self._attributes[self.__values_hdf_name__]

    @values.setter
    def values(self, val: BathymetryValues):
        self._attributes[self.__values_hdf_name__] = val

    @property
    def __values_type__(self) -> Type[BathymetryValues]:
        return BathymetryValues

    def values_create(self):
        """ Creates a blank, empty or zero value for values"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.values = self.__values_type__()

    @property
    def __version__(self) -> int:
        return 1

    @property
    def __minimum_depth_type__(self):
        return float

    def minimum_depth_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.minimum_depth = self.__minimum_depth_type__()

    @property
    def minimum_depth(self) -> float:
        """From 4.2.1.1.1.3,
        The attribute minimumDepth has the value type Real and describes the lower bound of the depth estimate
        for all the depth values in S102_BathymetryValues record.
        This attribute is required. There is no default"""
        return self._attributes[self.__minimum_depth_hdf_name__]

    @minimum_depth.setter
    def minimum_depth(self, val: float):
        self._attributes[self.__minimum_depth_hdf_name__] = val

    @property
    def __maximum_depth_type__(self):
        return float

    def maximum_depth_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.maximum_depth = self.__maximum_depth_type__()

    @property
    def maximum_depth(self) -> float:
        """From 4.2.1.1.1.4,
        The attribute minimumDepth has the value type Real and describes the lower bound of the depth estimate
        for all the depth values in S102_BathymetryValues record.
        This attribute is required. There is no default
        """
        return self._attributes[self.__maximum_depth_hdf_name__]

    @maximum_depth.setter
    def maximum_depth(self, val: float):
        self._attributes[self.__maximum_depth_hdf_name__] = val

    @property
    def minimum_uncertainty(self) -> float:
        """From 4.2.1.1.1.6,
        The attribute minimumUncertainty has the value type Real and describes the lower bound of the uncertainty of the
        depth estimate for all the depth values in S102_BathymetryValues record.
        This attribute is required. There is no default
        """
        return self._attributes[self.__minimum_uncertainty_hdf_name__]

    @minimum_uncertainty.setter
    def minimum_uncertainty(self, val: float):
        self._attributes[self.__minimum_uncertainty_hdf_name__] = val

    @property
    def __minimum_uncertainty_type__(self):
        return float

    def minimum_uncertainty_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.minimum_uncertainty = self.__minimum_uncertainty_type__()

    @property
    def maximum_uncertainty(self) -> float:
        """From 4.2.1.1.1.7,
        The attribute minimumUncertainty has the value type Real and describes the lower bound of the uncertainty
        of the depth estimate for all the depth values in S102_BathymetryValues record.
        This attribute is required. There is no default"""
        return self._attributes[self.__maximum_uncertainty_hdf_name__]

    @maximum_uncertainty.setter
    def maximum_uncertainty(self, val: float):
        self._attributes[self.__maximum_uncertainty_hdf_name__] = val

    @property
    def __maximum_uncertainty_type__(self):
        return float

    def maximum_uncertainty_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.maximum_uncertainty = self.__maximum_uncertainty_type__()

    @property
    def origin(self) -> DirectPosition:
        """From 4.2.1.1.1.8,
        The attribute origin has the value class DirectPosition which is a position that shall locate the origin of the rectified grid
        in the coordinate reference system.
        This attribute is required. There is no default
        """
        return self._attributes[self.__origin_hdf_name__]

    @origin.setter
    def origin(self, val: DirectPosition):
        self._attributes[self.__origin_hdf_name__] = val

    @property
    def __origin_type__(self):
        return DirectPosition

    def origin_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.origin = self.__origin_type__()

    @property
    def __origin_attribute_type__(self) -> Type[DirectPosition]:
        return DirectPosition

    @property
    def offset_vectors(self) -> s1xx_sequence:
        """sequence of s102 Vectors  From 4.2.1.1.1.9 in S102 v2.0.0
        The attribute offsetVectors has the value class Sequence<Vector> that shall be a sequence of offset vector elements
        that determine the grid spacing in each direction.
        The data type Vector is specified in ISO/TS 19103. This attribute is required.
        There is no default.
        """
        return self._attributes[self.__offset_vectors_hdf_name__]

    @offset_vectors.setter
    def offset_vectors(self, val: s1xx_sequence):
        self._attributes[self.__offset_vectors_hdf_name__] = val

    @property
    def __offset_vectors_type__(self):
        return numpy.ndarray

    def offset_vectors_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.offset_vectors = self.__offset_vectors_type__([2, ], numpy.float64)

    @property
    def dimension(self) -> int:
        """From 4.2.1.1.1.10,
        The attribute dimension has the value class Integer that shall identify the dimensionality of the grid.
        The value of the grid dimension in this product specification is 2.
        This value is fixed in this Product Specification and does not need to be encoded
        """
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
        self.dimension = self.__dimension_type__(2)

    @property
    def axis_names(self) -> s1xx_sequence:
        """sequence of character strings From 4.2.1.1.1.11,
        The attribute axisNames has the value class Sequence<CharacterString> that shall be used to assign names to the grid axis.
        The grid axis names shall be "Latitude" and "Longitude" for unprojected data sets or “Northing” and “Easting” in a projected space
        """
        return self._attributes[self.__axis_names_hdf_name__]

    @axis_names.setter
    def axis_names(self, val: s1xx_sequence):
        self._attributes[self.__axis_names_hdf_name__] = val

    @property
    def __axis_names_type__(self) -> Type[numpy.ndarray]:
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
        return self._attributes[self.__extent_hdf_name__]

    @extent.setter
    def extent(self, val: GridEnvelope):
        self._attributes[self.__extent_hdf_name__] = val

    @property
    def __extent_type__(self) -> Type[GridEnvelope]:
        return GridEnvelope

    def extent_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.extent = self.__extent_type__()

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
        return self._attributes[self.__sequencing_rule_hdf_name__]

    @sequencing_rule.setter
    def sequencing_rule(self, val: SequenceRule):
        self._attributes[self.__sequencing_rule_hdf_name__] = val

    @property
    def __sequencing_rule_type__(self):
        return SequenceRule

    def sequencing_rule_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.sequencing_rule = self.__sequencing_rule_type__()

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
        return self._attributes[self.__start_sequence_hdf_name__]

    @start_sequence.setter
    def start_sequence(self, val: GridCoordinate):
        self._attributes[self.__start_sequence_hdf_name__] = val

    @property
    def __start_sequence_type__(self):
        return GridCoordinate

    def start_sequence_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.start_sequence = self.__start_sequence_type__()


# mixin uses _attributes from the main class - ignore its warnings
# noinspection PyUnresolvedReferences
class DisplayScaleMixin:
    __maximum_display_scale_hdf_name__ = "maximumDisplayScale"  #: HDF5 naming
    __minimum_display_scale_hdf_name__ = "minimumDisplayScale"  #: HDF5 naming

    @property
    def minimum_display_scale(self) -> int:
        """From 4.2.1.1.1.2,
        The larger value of the ratio of the linear dimensions of the features of a dataset presented in the display and
        the actual dimensions of the features represented (largest scale) of the scale range of the dataset.
        A list of display scale ranges is available in Figure 11.1, 1st column
        """
        return self._attributes[self.__minimum_display_scale_hdf_name__]

    @minimum_display_scale.setter
    def minimum_display_scale(self, val: int):
        self._attributes[self.__minimum_display_scale_hdf_name__] = val

    @property
    def __minimum_display_scale_type__(self):
        return float

    def minimum_display_scale_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.minimum_display_scale = self.__minimum_display_scale_type__()

    @property
    def __maximum_display_scale_type__(self):
        return float

    def maximum_display_scale_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.maximum_display_scale = self.__maximum_display_scale_type__()

    @property
    def maximum_display_scale(self) -> int:
        """From 4.2.1.1.1.5,
        The smaller value of the ratio of the linear dimensions of the features of a dataset presented in the display and
        the actual dimensions of the features represented (smallest scale) of the scale range of the dataset.
        A list of display scale ranges is available in Table 11.1, 1st column
        """
        return self._attributes[self.__maximum_display_scale_hdf_name__]

    @maximum_display_scale.setter
    def maximum_display_scale(self, val: int):
        self._attributes[self.__maximum_display_scale_hdf_name__] = val


class BathymetryCoverage(BathymetryCoverageBase, DisplayScaleMixin):
    pass

class SurfaceCorrectionValues(VertexPoint):
    pass


# this S102_SurfaceCorrectionValues may not be right as the docs refer to S100_VertexPoint,
# but I'm betting they would both change if there ever is an update
class TrackingListValues(SurfaceCorrectionValues):
    """ From 4.2.1.1.10 of S102 v2.0.0
    """
    __track_code_hdf_name__ = "trackCode"  #: HDF5 naming
    __list_series_hdf_name__ = "listSeries"  #: HDF5 naming

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
        return self._attributes[self.__track_code_hdf_name__]

    @track_code.setter
    def track_code(self, val: str):
        self._attributes[self.__track_code_hdf_name__] = val

    @property
    def __track_code_type__(self):
        return str

    def track_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.track_code = self.__track_code_type__()

    @property
    def list_series(self) -> int:
        """From 4.2.1.1.10.3,
        The attribute listSeries has the value type Integer which contains an index number into a list of metadata
        elements describing the reason for the override of the corresponding depth and uncertainty values in the bathymetry coverage.
        """
        return self._attributes[self.__list_series_hdf_name__]

    @list_series.setter
    def list_series(self, val: int):
        self._attributes[self.__list_series_hdf_name__] = val

    @property
    def __list_series_type__(self):
        return int

    def list_series_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.list_series = self.__list_series_type__()


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


class TrackingListCoverage(CommonPointRule, S1xxObject):
    """ 4.2.1.1.9 and Figure 4.4 of v2.0.0
    commonPointRule is defined to be an S100_PointCoverage with a value of default and it therefore optional.
    a metadata attribute from S100 is allowed but not necessary as well.

    """
    write_format_str = ".%02d"

    __domain_extent_hdf_name__ = "domainExtent"  #: HDF5 naming
    __common_point_rule_hdf_name__ = "commonPointRule"  #: HDF5 naming
    __set_hdf_name__ = "set"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def domain_extent(self) -> S102MetadataListBase:
        return self._attributes[self.__domain_extent_hdf_name__]

    @domain_extent.setter
    def domain_extent(self, val: S102MetadataListBase):
        self._attributes[self.__domain_extent_hdf_name__] = val

    @property
    def __domain_extent_type__(self):
        return GeographicExtent

    def domain_extent_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.domain_extent = self.__domain_extent_type__()

    def common_point_rule_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.common_point_rule = self.__common_point_rule_type__["average"]

    @property
    def __set_type__(self):
        return TrackingListSetList

    def set_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.set = self.__set_type__()

    @property
    def set(self) -> S102MetadataListBase:
        """
        Returns
        -------
        TrackingListValuesList
            list of TrackingListValues
        """
        return self._attributes[self.__set_hdf_name__]

    @set.setter
    def set(self, val: S102MetadataListBase):
        self._attributes[self.__set_hdf_name__] = val

    # @TODO  I don't think this is right, but not sure where I found it
    #
    # __geometry_hdf_name__ = "geometry"  #: HDF5 naming
    #
    # @property
    # def geometry(self) -> S1xxObject:
    #     return self._attributes[self.__geometry_hdf_name__]
    #
    # @geometry.setter
    # def geometry(self, val: S1xxObject):
    #     self._attributes[self.__geometry_hdf_name__] = val
    #
    # __value_hdf_name__ = "value"  #: HDF5 naming
    #
    # @property
    # def value(self) -> s1xx_sequence:
    #     return self._attributes[self.__value_hdf_name__]
    #
    # @value.setter
    # def value(self, val: s1xx_sequence):
    #     self._attributes[self.__value_hdf_name__] = val


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
    __bathymetry_group_hdf_name__ = "Group" + r"[\._]\d+"
    """ Basic template for HDF5 naming of the attribute.  
    Attribute name will be automatically determined based on the list's index of the data. 
    """

    @property
    def __bathymetry_group_type__(self):
        return BathymetryGroupList

    def bathymetry_group_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.bathymetry_group = self.__bathymetry_group_type__()

    @property
    def bathymetry_group(self) -> S102MetadataListBase:
        """ The bathymetry data, a list of Bathymetrygroup
        Returns
        -------
        S102MetadataListBase
            Contains a list of BathymetryCoverage objects via the BathymetryCoveragesList class
        """
        return self._attributes[self.__bathymetry_group_hdf_name__]

    @bathymetry_group.setter
    def bathymetry_group(self, val: S102MetadataListBase):
        self._attributes[self.__bathymetry_group_hdf_name__] = val


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
    __bathymetry_coverage_hdf_name__ = BATHY_COVERAGE + r"[\._]\d+"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def __bathymetry_coverage_type__(self):
        return BathymetryCoveragesList

    def bathymetry_coverage_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.bathymetry_coverage = self.__bathymetry_coverage_type__()

    @property
    def bathymetry_coverage(self) -> S102MetadataListBase:
        """ The bathymetry data, a list of BathymetryCoverage

        Returns
        -------
        S102MetadataListBase
            Contains a list of BathymetryCoverage objects via the BathymetryCoveragesList class
        """
        return self._attributes[self.__bathymetry_coverage_hdf_name__]

    @bathymetry_coverage.setter
    def bathymetry_coverage(self, val: S102MetadataListBase):
        self._attributes[self.__bathymetry_coverage_hdf_name__] = val

    def data_coding_format_create(self):
        """ Creates a blank, empty or zero value for data_coding_format"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_coding_format = self.__data_coding_format_type__(2)  # regular grid

    def dimension_create(self):
        """ Creates a blank, empty or zero value for dimension"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.dimension = self.__dimension_type__(2)


class TrackingListContainer(FeatureContainerDCF2):
    """
    Table 10.1 of v2.0.0
    """

    __tracking_list_coverage_hdf_name__ = TRACKING_COVERAGE + r"[\._]\d+"

    @property
    def __version__(self) -> int:
        return 1

    def data_coding_format_create(self):
        """ Creates a blank, empty or zero value for data_coding_format"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_coding_format = self.__data_coding_format_type__(1)  # point set

    @property
    def __tracking_list_coverage_type__(self):
        return TrackingListCoveragesList

    def tracking_list_coverage_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.tracking_list_coverage = self.__tracking_list_coverage_type__()

    @property
    def tracking_list_coverage(self) -> S1xxObject:
        """ The tracking list data, a list of TrackingListCoverage
        Returns
        -------
        S102MetadataListBase
            Contains a list of TrackingListCoverage objects via the TrackingListCoveragesList class
        """
        return self._attributes[self.__tracking_list_coverage_hdf_name__]

    @tracking_list_coverage.setter
    def tracking_list_coverage(self, val: S1xxObject):
        self._attributes[self.__tracking_list_coverage_hdf_name__] = val


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
        self.unit_of_measure = self.__unit_of_measure_type__("metres")

    def fill_value_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.fill_value = self.__fill_value_type__(1000000)

    def datatype_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.datatype = self.__datatype_type__("H5T_NATIVE_FLOAT")

    def lower_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.lower = self.__lower_type__(-12000)

    def upper_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.upper = self.__upper_type__(12000)

    def closure_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.closure = self.__closure_type__("closedInterval")


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


class FeatureCodesBase(GroupFBase):
    """ Table 10.1 and sect 10.2.1 of v2.0.0
    """

    __feature_name_hdf_name__ = "featureName"  #: HDF5 naming
    __bathymetry_coverage_dataset_hdf_name__ = BATHY_COVERAGE
    __tracking_list_coverage_hdf_name__ = TRACKING_COVERAGE

    def feature_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_code = self.__feature_code_type__([BATHY_COVERAGE, TRACKING_COVERAGE], dtype=h5py_string_dtype)

    @property
    def __version__(self) -> int:
        return 1

    @property
    def __bathymetry_coverage_dataset_type__(self):
        return BathymetryCoverageDataset

    def bathymetry_coverage_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.bathymetry_coverage_dataset = self.__bathymetry_coverage_dataset_type__()

    @property
    def bathymetry_coverage_dataset(self) -> BathymetryCoverageDataset:
        return self._attributes[self.__bathymetry_coverage_dataset_hdf_name__]

    @bathymetry_coverage_dataset.setter
    def bathymetry_coverage_dataset(self, val: BathymetryCoverageDataset):
        self._attributes[self.__bathymetry_coverage_dataset_hdf_name__] = val


# mixin uses _attributes from the main class - ignore its errors
# noinspection PyUnresolvedReferences
class FeatureCodesTrackingMixin:
    @property
    def __tracking_list_coverage_type__(self):
        return TrackingListCoverageDataset

    def tracking_list_coverage_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.tracking_list_coverage = self.__tracking_list_coverage_type__()

    @property
    def tracking_list_coverage(self) -> TrackingListCoverageDataset:
        return self._attributes[self.__tracking_list_coverage_hdf_name__]

    @tracking_list_coverage.setter
    def tracking_list_coverage(self, val: TrackingListCoverageDataset):
        self._attributes[self.__tracking_list_coverage_hdf_name__] = val

    @property
    def __feature_name_type__(self):
        return numpy.array

    @property
    def feature_name(self) -> s1xx_sequence:
        return self._attributes[self.__feature_name_hdf_name__]

    @feature_name.setter
    def feature_name(self, val: s1xx_sequence):
        self._attributes[self.__feature_name_hdf_name__] = val

    def feature_name_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_name = self.__feature_name_type__([BATHY_COVERAGE, TRACKING_COVERAGE], dtype=h5py_string_dtype)



class FeatureCodes(FeatureCodesBase, FeatureCodesTrackingMixin):
    pass


class S102Root(S100Root):
    """The root group contains a feature information group and N feature containers.
    In S102 there are currently two feature containers which are the 'coverages'  bathymetry and tracking list.
    The coverage names are determined from the matching CoveragesAttributes
    10.2 and Figure 10.1 of v2.0.0
    """
    __feature_information_hdf_name__ = "Group_F"  #: HDF5 naming
    __bathymetry_coverage_hdf_name__ = BATHY_COVERAGE
    __tracking_list_coverage_hdf_name__ = TRACKING_COVERAGE

    @property
    def __version__(self) -> int:
        return 1

    @property
    def feature_information(self) -> FeatureCodesBase:
        """Feature Information stored in GroupF in the HDF5 using :class:`FeatureCodes`"""
        return self._attributes[self.__feature_information_hdf_name__]

    @feature_information.setter
    def feature_information(self, val: FeatureCodesBase):
        self._attributes[self.__feature_information_hdf_name__] = val

    @property
    def __feature_information_type__(self):
        return FeatureCodes

    def feature_information_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_information = self.__feature_information_type__()

    @property
    def bathymetry_coverage(self) -> S1xxObject:
        """Bathymetry instance stored under the HDF5 root using :class:`BathymetryContainer`
        """
        return self._attributes[self.__bathymetry_coverage_hdf_name__]

    @property
    def __bathymetry_coverage_type__(self):
        return BathymetryContainer

    def bathymetry_coverage_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.bathymetry_coverage = self.__bathymetry_coverage_type__()

    @bathymetry_coverage.setter
    def bathymetry_coverage(self, val: S1xxObject):
        self._attributes[self.__bathymetry_coverage_hdf_name__] = val

    @property
    def __tracking_list_coverage_type__(self):
        return TrackingListContainer

    def tracking_list_coverage_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.tracking_list_coverage = self.__tracking_list_coverage_type__()

    @property
    def tracking_list_coverage(self) -> S1xxObject:
        return self._attributes[self.__tracking_list_coverage_hdf_name__]

    @tracking_list_coverage.setter
    def tracking_list_coverage(self, val: S1xxObject):
        self._attributes[self.__tracking_list_coverage_hdf_name__] = val


class S102File(S100File):
    PRODUCT_SPECIFICATION = PRODUCT_SPECIFICATION
    # these keys allow backward compatibility with NAVO data, the first key is current at time of writing
    top_level_keys = ('BathymetryCoverage', 'S102_Grid', 'S102_BathymetryCoverage')
    tracking_list_top_level = ("TrackingListCoverage",)
    tracking_list_second_level = ("TrackingListCoverage.01",)
    tracking_list_group_level = ("Group.001",)
    second_level_keys = (
        'BathymetryCoverage.001', 'BathymetryCoverage.01', 'S102_Grid.01', 'S102_BathymetryCoverage.01', 'BathymetryCoverage_01', 'S102_Grid_01', 'S102_BathymetryCoverage_01',)
    group_level_keys = ('Group.001', 'Group_001',)
    value_level_keys = ("values",)
    depth_keys = ("depth", "depths", 'elevation', "elevations", "S102_Elevation")

    def __init__(self, name, *args, **kywrds):
        if 'root' not in kywrds:
            kywrds['root'] = S102Root  # inherited classes will specify their own root type
        super().__init__(name, *args, **kywrds)

    @property
    def z_down(self) -> bool:
        return False


    def subdivide(self, path, rows, cols):
        # hp5y does not have the ability to repack the data.
        # This means that when replacing data the file will not shrink or even grow when you'd expect it to shrink.
        # Basically it allocates space on disk and if you delete data and then add it will reuse that storage block but not reduce the space.
        # If the new data is larger than the old data then it adds space at the end and leave empty space in the old location.
        # There is a free "h5repack" utility from hdfgroup with binaries which would remove the empty space.
        # However, we want to be platform independent so we will do an end run.
        #
        # 0) Make sure any changes are on disk
        # 1) Copy the file to a temporary file
        # 2) Open the copy and delete the datasets we need to subdivide
        #    (it deletes it on disk so we can't operate on the original safely).
        # 3) Then make as many copies as needed using the h5py copy method (which should stop the empty space issue).
        # 4) Finally add the bathy+uncertainty data to each sub-file and revise the extents and min/max values.

        # self.write()
        self.flush()
        base_path = pathlib.Path(path).with_suffix("")
        try:
            tmpname = tempfile.mktemp(".h5", "tmp")
            shutil.copy(self.filename, tmpname)
        except FileExistsError:
            tmpname = tempfile.mktemp(".h5", "tmp")
            shutil.copy(self.name, tmpname)
        tmp = self.__class__(tmpname, "r+")
        # delete the bathy and uncertainty
        # @FIXME @TODO add a remove hdf5 method
        del tmp[tmp.root.bathymetry_coverage.bathymetry_coverage[0]._hdf5_path]  # force the data out of hdf5
        del tmp.root.bathymetry_coverage.bathymetry_coverage[0]
        tmp.flush()
        bathy_01 = self.root.bathymetry_coverage.bathymetry_coverage[0]
        bathy_group_object = bathy_01.bathymetry_group[0]
        grid = bathy_group_object.values
        depth_grid = grid.depth
        uncert_grid = grid.uncertainty
        origin = bathy_group_object.origin.coordinate
        res_x = bathy_01.grid_spacing_longitudinal
        res_y = bathy_01.grid_spacing_latitudinal
        res = numpy.array([res_x, res_y])
        # make list of indices to subdivide with, so a size 1000 array divided 3 times gives [0, 333, 666, None]
        row_indices = (numpy.arange(rows) * int(depth_grid.shape[0] / rows)).tolist() + [None]
        col_indices = (numpy.arange(cols) * int(depth_grid.shape[0] / cols)).tolist() + [None]
        no_data = self.root.feature_information.bathymetry_coverage_dataset[0].fill_value
        fnames = []
        for r in range(rows):
            for c in range(cols):
                out_path = base_path.with_suffix(f".{r+1}_{c+1}.h5")
                fnames.append(str(out_path))
                # I haven't figured out how to copy the root to the temp root - it gives errors about "no name"
                raw_out = h5py.File(str(out_path), "w")
                for key in tmp.keys():  # copy all groups+datasets data from the root
                    tmp.copy(tmp[key], raw_out['/'], key)
                for name in tmp.attrs.keys():  # copy the attributes of the root (I'd think there'd be a better way)
                    raw_out.attrs[name] = tmp.attrs[name]  # out.create(n, a, dtype=a.dtype)
                raw_out.close()

                out = self.__class__(str(out_path), "r+")
                out.root.bathymetry_coverage.bathymetry_coverage_create()  # we deleted the bathymetry_coverage above, so make a new container

                start_row = row_indices[r]
                end_row = row_indices[r+1]
                start_col = col_indices[c]
                end_col = col_indices[c+1]
                local_origin = origin + numpy.array([start_col, start_row]) * res
                sub_depth_grid = depth_grid[start_row:end_row, start_col:end_col]
                sub_uncert_grid = uncert_grid[start_row:end_row, start_col:end_col]
                # just set the origin and res for the load_arrays_with_metadata
                # this will reset the coordinates and min/max automatically
                # but use overwrite=False since we want to maintain the original file's spatial reference, dates etc.
                metadata = {"origin": local_origin, "res": res, 'metadataFile': out_path.with_suffix(".xml").name}
                out.load_arrays_with_metadata(sub_depth_grid, sub_uncert_grid, metadata, overwrite=False, nodata_value=no_data)
        tmp.close()
        del tmp
        try:
            os.remove(tmpname)
        except (FileNotFoundError, PermissionError):
            print(f"Failed to remove temp file {tmpname}")
        return fnames

    @staticmethod
    def get_valid_epsg() -> list:
        """
        Create and return the list of valid EPSG codes for S-102 version 2.0.
        """
        valid_epsg = [4326, 5041, 5042]
        valid_epsg += list(numpy.arange(32601, 32660 + 1))
        valid_epsg += list(numpy.arange(32701, 32760 + 1))
        return valid_epsg

    @staticmethod
    def upgrade_in_place(s100_object):
        try:
            spec = s100_object.root.product_specification
        except:
            spec = "unknown"
        if not h5py_string_comp(spec, PRODUCT_SPECIFICATION):
            raise S102Exception(f"Could not upgrade file of type {spec}")

    @classmethod
    def upgrade(cls, src_filename, dest_filename=None, mode='r'):
        if dest_filename is None:
            dest_filename = src_filename
        else:
            shutil.copy(src_filename, dest_filename)
        s100_object = S100File(dest_filename, "r+")
        cls.upgrade_in_place(s100_object)
        s100_object.close()
        del s100_object
        return cls(dest_filename, mode)

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

    # noinspection PyUnboundLocalVariable
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


    @classmethod
    def _get_S102File(cls, output_file):
        """ Small helper function to convert the output_file parameter into a S102File, currently accepting file path as string or S102File instance.
        Could propbably accept h5py.File or other things in the future"""
        if isinstance(output_file, S1XXFile):
            data_file = output_file
        else:  # try everything else -- pathlib, str, tempfile, io.BytesIO
            try:
                data_file = cls(output_file, "w")
            except TypeError as typeerr:
                msg = "Failed to create S102File using {}".format(str(output_file))
                logging.error(msg)
                raise type(typeerr)(msg).with_traceback(sys.exc_info()[2])

        return data_file

    @classmethod
    def create_s102(cls, output_file, overwrite=True) -> S102File:
        data_file = cls._get_S102File(output_file)
        data_file.set_defaults(overwrite=overwrite)
        return data_file

    def _set_tracking_defaults(self, overwrite=True):
        # I'm not sure what to put here, yet
        root = self.root
        tracking_cov = root.feature_information.tracking_list_coverage

        track_info = tracking_cov.append_new_item()  # append(tracking_cov.metadata_type())
        track_info.initialize_properties(True, overwrite=overwrite)
        track_info.code = "X"
        track_info.name = "X"
        track_info.unit_of_measure = "N/A"

        track_info = tracking_cov.append_new_item()
        track_info.initialize_properties(True, overwrite=overwrite)
        track_info.code = "Y"
        track_info.name = "Y"
        track_info.unit_of_measure = "N/A"

        track_info = tracking_cov.append_new_item()
        track_info.initialize_properties(True, overwrite=overwrite)
        track_info.code = "originalValue"
        track_info.name = "Original Value"

        track_info = tracking_cov.append_new_item()
        track_info.initialize_properties(True, overwrite=overwrite)
        track_info.code = "trackCode"
        track_info.name = "Track Code"
        track_info.unit_of_measure = "N/A"

        track_info = tracking_cov.append_new_item()
        track_info.initialize_properties(True, overwrite=overwrite)
        track_info.code = "listSeries"
        track_info.name = "List Series"
        track_info.unit_of_measure = "N/A"

    def set_defaults(self, overwrite=True):
        """ Creates or updates an S102File object.
        Default values are set for any data that don't have options or are mandatory to be filled in the S102 spec.

        Parameters
        ----------
        output_file
            Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
        overwrite
            If updating an existing file then set this option to False in order to retain data (not sure this is needed).

        Returns
        -------
        S102File
            The object created or updated by this function.


        """
        # @fixme @todo -- I think this will overwrite no matter what, need to look into that
        self.create_empty_metadata()  # init the root with a fully filled out empty metadata set
        self._set_bathy_defaults()
        self._set_tracking_defaults()

    def _set_bathy_defaults(self, overwrite=True):
        """ This function initializes the values in more recent versions of the spec to reduce redundant code in later modules
        """
        root = self.root
        bathy_cov_dset = root.feature_information.bathymetry_coverage_dataset
        bathy_depth_info = bathy_cov_dset.append_new_item()  # bathy_cov_dset.append(bathy_cov_dset.metadata_type())
        bathy_depth_info.initialize_properties(True, overwrite=overwrite)
        bathy_depth_info.code = DEPTH
        bathy_depth_info.name = DEPTH
        # these are auto-filled by the api
        # bathy_depth_info.unit_of_measure="metres"
        # bathy_depth_info.fill_value=1000000.0
        # bathy_depth_info.datatype=H5T_NATIVE_FLOAT
        # bathy_depth_info.lower = -12000
        # bathy_depth_info.upper = 12000
        # bathy_depth_info.closure = "closedInterval"

        bathy_uncertainty_info = bathy_cov_dset.append_new_item()
        bathy_uncertainty_info.initialize_properties(True, overwrite=overwrite)
        bathy_uncertainty_info.code = UNCERTAINTY
        bathy_uncertainty_info.name = UNCERTAINTY

        root.bathymetry_coverage.axis_names = numpy.array(["Longitude", "Latitude"])  # row major order means X/longitude first
        root.bathymetry_coverage.sequencing_rule_scan_direction = "Longitude, Latitude"
        root.bathymetry_coverage.common_point_rule = 1  # average
        # root.bathymetry_coverage.data_coding_format = 2  # default
        # root.bathymetry_coverage.dimension = 2  # default value
        root.bathymetry_coverage.interpolation_type = 1  # nearest neighbor
        root.bathymetry_coverage.num_instances = 1  # how many Bathycoverages
        root.bathymetry_coverage.sequencing_rule_type = 1  # linear
        del root.bathymetry_coverage.time_uncertainty

    @classmethod
    def from_arrays(cls, depth_grid: s1xx_sequence, uncert_grid: s1xx_sequence, output_file, nodata_value=None,
                    flip_x: bool = False, flip_y: bool = False, overwrite: bool = True,
                    flip_z: bool = False) -> S102File:  # num_array, or list of lists accepted
        """  Creates or updates an S102File object based on numpy array/h5py datasets.
        Calls :any:`create_s102` then fills in the HDF5 datasets with the supplied depth_grid and uncert_grid.
        Fills the number of points areas and any other appropriate places in the HDF5 file per the S102 spec.

        For most parameters, see S102File.load_arrays

        Parameters
        ----------
        output_file
            Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
        """
        data_file = cls.create_s102(output_file)
        data_file.load_arrays(depth_grid, uncert_grid, nodata_value=nodata_value,
                              flip_x=flip_x, flip_y=flip_y, overwrite=overwrite,
                              flip_z=flip_z)
        return data_file

    def load_arrays(self, depth_grid: s1xx_sequence, uncert_grid: s1xx_sequence, nodata_value=None,
                    flip_x: bool = False, flip_y: bool = False, overwrite: bool = True,
                    flip_z: bool = False):  # num_array, or list of lists accepted

        """  Updates an S102File object based on numpy array/h5py datasets.
        Calls :any:`create_s102` then fills in the HDF5 datasets with the supplied depth_grid and uncert_grid.
        Fills the number of points areas and any other appropriate places in the HDF5 file per the S102 spec.

        Raises an S102Exception if the shapes of the depth and uncertainty (if not None) grids are not equal.

        Parameters
        ----------
        depth_grid
        uncert_grid
            The uncertainty dataset to embed in the object.
            If None then a numpy.zeros array will be created in the appropriate shape to be stored in the file.
        nodata_value
            Value used to denote an empty cell in the grid.  Used in finding the min/max and then converted to the S102 fillValue.
        flip_x
            boolean if the data should be mirrored on x coordinate (i.e. the original grid is right to left)
            Flips are done here so we can implement a chunked read/write to save memory
        flip_y
            boolean if the data should be mirrored on y coordinate (i.e. the original grid is top to bottom)
            Flips are done here so we can implement a chunked read/write to save memory
        overwrite
            If updating an existing file then set this option to False in order to retain data (not sure this is needed).
        flip_z
            boolean if the data should be reversed in z coordinate (i.e. the original grid is upside down)
            Flips are done here so we can implement a chunked read/write to save memory.
            This is after overwrite for backwards compatibility.

        Returns
        -------
        S102File
            The object created or updated by this function.

        """
        # @todo -- Add logic that if the grids are gdal raster bands then read in blocks and use h5py slicing to write in blocks.
        #   Slower but saves resources
        root = self.root
        try:
            bathy_01 = root.bathymetry_coverage.bathymetry_coverage[0]
        except IndexError:
            bathy_01 = root.bathymetry_coverage.bathymetry_coverage.append_new_item()
        bathy_01.initialize_properties(recursively_create_children=True, overwrite=overwrite)

        del bathy_01.grid_spacing_vertical
        del bathy_01.grid_origin_vertical
        del bathy_01.number_of_times
        del bathy_01.time_record_interval
        del bathy_01.date_time_of_last_record
        del bathy_01.date_time_of_first_record
        bathy_01.num_grp = 1

        try:
            bathy_group_object = bathy_01.bathymetry_group[0]
        except IndexError:
            bathy_group_object = bathy_01.bathymetry_group.append_new_item()
        # bathy_group_object.initialize_properties()  # Not creating everything as I'm not sure if the grid attributes should be there

        # @todo @fixme fix here -- row/column order?
        rows, cols = depth_grid.shape
        if uncert_grid is None:
            uncert_grid = numpy.full(depth_grid.shape, nodata_value, dtype=numpy.float32)
        if depth_grid.shape != uncert_grid.shape:
            raise S102Exception("Depth and Uncertainty grids have different shapes")

        bathy_01.num_points_latitudinal = rows
        bathy_01.num_points_longitudinal = cols
        bathy_01.start_sequence = "0,0"
        del bathy_01.num_points_vertical
        del bathy_01.vertical_extent_maximum_z
        del bathy_01.vertical_extent_minimum_z

        bathy_group_object.extent_create()
        bathy_group_object.extent.initialize_properties(True, overwrite=overwrite)
        bathy_group_object.extent.low.coord_values[0:2] = [0, 0]
        bathy_group_object.extent.high.coord_values[0:2] = [rows, cols]

        try:
            uncertainty_max = uncert_grid[uncert_grid != nodata_value].max()
            uncertainty_min = uncert_grid[uncert_grid != nodata_value].min()
        except ValueError:  # an empty uncertainty array (all values == nodata) will cause this
            uncertainty_max = uncertainty_min = nodata_value
        bathy_group_object.minimum_uncertainty = uncertainty_min
        bathy_group_object.maximum_uncertainty = uncertainty_max

        bathy_group_object.dimension = 2

        bathy_group_object.origin_create()
        bathy_group_object.origin.initialize_properties(True, overwrite=overwrite)
        bathy_group_object.origin.dimension = 2

        bathy_group_object.values_create()
        grid = bathy_group_object.values
        # @todo -- need to make sure nodata values are correct,
        #   especially if converting something other than bag which is supposed to have the same nodata value
        # @todo -- Add logic that if the grids are gdal raster bands then read in blocks and use h5py slicing to write in blocks.
        #   Slower but saves resources
        if flip_x:
            depth_grid = numpy.fliplr(depth_grid)
            uncert_grid = numpy.fliplr(uncert_grid)
        if flip_y:
            depth_grid = numpy.flipud(depth_grid)
            uncert_grid = numpy.flipud(uncert_grid)
        if flip_z:
            depth_grid[depth_grid != nodata_value] *= -1

        try:
            depth_max = depth_grid[depth_grid != nodata_value].max()
            depth_min = depth_grid[depth_grid != nodata_value].min()
        except ValueError:  # an empty depth array (all values == nodata) will cause this, subdivide() may cause this or data to be updated later
            depth_min = depth_max = nodata_value
        bathy_group_object.maximum_depth = depth_max
        bathy_group_object.minimum_depth = depth_min

        if nodata_value != root.feature_information.bathymetry_coverage_dataset[0].fill_value:
            depth_grid = numpy.copy(depth_grid)
            depth_grid[depth_grid == nodata_value] = root.feature_information.bathymetry_coverage_dataset[0].fill_value
            uncert_grid = numpy.copy(uncert_grid)
            uncert_grid[uncert_grid == nodata_value] = root.feature_information.bathymetry_coverage_dataset[1].fill_value

        grid.depth = depth_grid
        grid.uncertainty = uncert_grid


    @classmethod
    def from_arrays_with_metadata(cls, depth_grid: s1xx_sequence, uncert_grid: s1xx_sequence, metadata: dict, output_file, nodata_value=None,
                                  overwrite: bool = True, flip_z: bool = False) -> S102File:  # raw arrays and metadata accepted
        """  Creates or updates an S102File object based on numpy array/h5py datasets.
        Calls :any:`create_s102` then fills in the HDF5 datasets with the supplied depth_grid and uncert_grid.
        Fills the number of points areas and any other appropriate places in the HDF5 file per the S102 spec.

        For most parameters, see S102File.load_arrays

        Parameters
        ----------
        output_file
            Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
        """
        data_file = cls.create_s102(output_file)
        data_file.load_arrays_with_metadata(depth_grid, uncert_grid, metadata, nodata_value=nodata_value,
                                  overwrite=overwrite, flip_z=flip_z)
        return data_file

    def load_arrays_with_metadata(self, depth_grid: s1xx_sequence, uncert_grid: s1xx_sequence, metadata: dict, nodata_value=None,
                                  overwrite: bool = True, flip_z: bool = False):  # raw arrays and metadata accepted
        """ Fills or creates an :any:`S102File` from the given arguments.

        Parameters
        ----------
        depth_grid
            a numpy or hdf5 dataset object of the rectangular grid of depths
        uncert_grid
            a numpy or hdf5 dataset object of the rectangular grid of uncertainties, lower left corner is the first point
        metadata
            a dictionary of metadata describing the grids passed in,
            metadata should have the following key/value pairs:
                - "origin": tuple of the position (x,y) or (lon, lat) for the reference corner node.
                    Other corners are calulated from this corner using the resolution and size of the data array.
                - "res": tuple of the resolution (cell size) of each grid cell (x, y).
                    Lower left corner is the first point of both resolutions are positive.
                    If a resolution is negative then the grid will be flipped in that dimension and the origin adjusted accordingly.
                - "horizontalDatumReference": See :any:`S102Root` horizontal_datum_reference, ex: "EPSG".
                    "EPSG" is the default value.
                - "horizontalDatumValue":  The value for the horizontal data such as the EPSG code ex: 32611
                - "epoch":
                - "geographicIdentifier": Location of the data, ex: "Long Beach, CA, USA".
                    An empty string ("") is the default.
                - "issueDate":
                - "metadataFile": File name for the associated discovery metatadata (xml)
        output_file
            Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
        nodata_value
            the "no data" value used in the grids
        overwrite
            if the output_file was an existing S102File then keep any attributes that might have
        Returns
        -------
        S102File

        """
        # @todo - add logic to see if the coordinate system is lower right, if not then need to mirror the arrays or add flags to do that in from_arrays
        res_x, res_y = metadata["res"]
        flip_x = True if res_x < 0 else False
        flip_y = True if res_y < 0 else False

        self.load_arrays(depth_grid, uncert_grid, nodata_value=nodata_value, overwrite=overwrite, flip_x=flip_x, flip_y=flip_y,
                                flip_z=flip_z)

        rows, cols = depth_grid.shape
        corner_x, corner_y = metadata['origin']

        # S-102 is node based, so distance to far corner is res * (n -1)
        opposite_corner_x = corner_x + res_x * (cols - 1)
        opposite_corner_y = corner_y + res_y * (rows - 1)

        minx = min((corner_x, opposite_corner_x))
        maxx = max((corner_x, opposite_corner_x))
        miny = min((corner_y, opposite_corner_y))
        maxy = max((corner_y, opposite_corner_y))

        # now add the additional metadata
        root = self.root
        bathy_01 = root.bathymetry_coverage.bathymetry_coverage[0]
        bathy_group_object = bathy_01.bathymetry_group[0]

        root.east_bound_longitude = maxx
        root.west_bound_longitude = minx
        root.south_bound_latitude = miny
        root.north_bound_latitude = maxy
        bathy_01.east_bound_longitude = maxx
        bathy_01.west_bound_longitude = minx
        bathy_01.south_bound_latitude = miny
        bathy_01.north_bound_latitude = maxy
        bathy_01.grid_origin_latitude = miny

        bathy_01.grid_origin_longitude = minx
        bathy_01.grid_origin_latitude = miny
        bathy_01.grid_spacing_longitudinal = abs(res_x)  # we adjust for negative resolution in the from_arrays
        bathy_01.grid_spacing_latitudinal = abs(res_y)

        bathy_group_object.origin.coordinate = numpy.array([minx, miny])

        # these names are taken from the S100/S102 attribute names
        # but are hard coded here to allow the S102 spec to change but not affect any tools built on these utility functions
        if "horizontalDatumReference" in metadata or overwrite:
            root.horizontal_datum_reference = metadata.get("horizontalDatumReference", "EPSG")
        if "horizontalDatumValue" in metadata or overwrite:
            source_epsg = int(metadata.get("horizontalDatumValue", 0))
            if source_epsg in self.get_valid_epsg():
                root.horizontal_datum_value = source_epsg
            else:
                raise S102Exception(f'The provided EPSG code {source_epsg} is not within the S102 specified values.')
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(root.horizontal_datum_value)
        if srs.IsProjected():
            axes = ["Easting", "Northing"]  # ["Northing", "Easting"]  # row major instead of
        else:
            axes = ["Longitude", "Latitude"]  # ["Latitude", "Longitude"]  # row major instead of

        bathy_group_object.axis_names = numpy.array(axes)  # row major order means X/longitude first
        root.bathymetry_coverage.axis_names = numpy.array(axes)  # row major order means X/longitude first
        root.bathymetry_coverage.sequencing_rule_scan_direction = ", ".join(axes)

        if "epoch" in metadata or overwrite:
            root.epoch = metadata.get("epoch", "")  # e.g. "G1762"  this is the 2013-10-16 WGS84 used by CRS
        if "geographicIdentifier" in metadata or overwrite:
            root.geographic_identifier = metadata.get("geographicIdentifier", "")
        if "issueDate" in metadata or overwrite:
            root.issue_date = metadata.get('issueDate', "")  # datetime.date.today().isoformat()
        if "metadataFile" in metadata or overwrite:
            root.metadata = metadata.get('metadataFile', "")  # datetime.date.today().isoformat()

        self.write()
        self.flush()

    @classmethod
    def from_raster(cls, input_raster, output_file, metadata: dict = None, flip_z=False) -> S102File:  # gdal instance or filename accepted
        """  Fills or creates an :any:`S102File` from the given arguments.

        For most parameters, see :any:`S102File.load_arrays`

        Parameters
        ----------
        output_file
            Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
        """
        data_file = cls.create_s102(output_file)
        data_file.load_gdal(input_raster, metadata=metadata, flip_z=flip_z)
        return data_file

    from_gdal = from_raster  # alias

    def load_gdal(self, input_raster, metadata: dict = None, flip_z=False):  # gdal instance or filename accepted
        """ Fills or creates an :any:`S102File` from the given arguments.

        Parameters
        ----------
        input_raster
            Either a path to a raster file that GDAL can open or a gdal.Dataset object.
        output_file
            Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
        metadata
            A dictionary of metadata describing the grids passed in.
            All the metadata used in :any:`from_from_arrays_with_metadata` can be specified and
            would override the values that would have been populated based on the GDAL data.

            horizontalDatumReference, horizontalDatumValue, origin, res will be determined from GDAL if not otherwise specified.

        Returns
        -------
        S102File

        """
        if metadata is None:
            metadata = {}
        else:
            metadata = metadata.copy()

        if isinstance(input_raster, gdal.Dataset):
            dataset = input_raster
        else:
            dataset = gdal.Open(str(input_raster))

        # @todo @fixme -- transform the coordinate system to a WGS84.  Strictly this may not end up being square, so how do we handle
        #  transform = osr.CoordinateTransformation( src_srs, tgt_srs)
        # Until we have a working datum engine this module should not do datum transformations - GR 20200402
        if "horizontalDatumReference" not in metadata or "horizontalDatumValue" not in metadata:
            metadata["horizontalDatumReference"] = "EPSG"
            sr = osr.SpatialReference(dataset.GetProjection())
            epsg = sr.GetAuthorityCode(None)
            # FIXME: this is likely incorrect. We probably don't want to get the code of the geographic CRS when the CRS is projected
            if epsg is None and sr.IsProjected():
                sr = sr.CloneGeogCS()
            if epsg:
                metadata["horizontalDatumValue"] = int(epsg)
            else:
                if sr.GetAttrValue("GEOGCS") == 'WGS 84':
                    metadata["horizontalDatumValue"] = 4326
                # elif sr.GetAttrValue("GEOGCS") == 'North_American_Datum_1983':
                #    metadata["horizontalDatumValue"] = 4269
                else:
                    raise S102Exception("Projection not understood, was searching for an EPSG code and found " + osr.SpatialReference(
                        dataset.GetProjection()).ExportToWkt())

        if "epoch" not in metadata:
            # @todo We should be able to pull this from the WKT
            pass

        raster_band = dataset.GetRasterBand(1)
        depth_nodata_value = raster_band.GetNoDataValue()
        uncertainty_band = dataset.GetRasterBand(2)

        ulx, dxx, dxy, uly, dyx, dyy = dataset.GetGeoTransform()
        if dxy != 0.0 or dyx != 0.0:
            raise S102Exception("raster is not north up but is rotated, this is not handled at this time")

        if "origin" not in metadata:
            # shift the gdal geotransform corner point to reference the node (pixel is center) rather than cell (pixel is area)
            metadata["origin"] = [ulx + dxx / 2, uly + dyy / 2]
        if "res" not in metadata:
            metadata["res"] = [dxx, dyy]
        self.load_arrays_with_metadata(raster_band.ReadAsArray(), uncertainty_band.ReadAsArray(), metadata,
                                                   nodata_value=depth_nodata_value, flip_z=flip_z)

    @classmethod
    def from_bag(cls, bagfile, output_file, metadata: dict = None) -> S102File:
        """  Fills or creates an :any:`S102File` from the given arguments.

        For most parameters, see :any:`S102File.load_arrays`

        Parameters
        ----------
        output_file
            Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
        """
        data_file = cls.create_s102(output_file)
        data_file.load_bag(bagfile=bagfile, metadata=metadata)
        return data_file

    def load_bag(self, bagfile, metadata: dict = None):
        """
        Parameters
        ----------
        bagfile
            Either a path to a raster file that GDAL can open or a gdal.Dataset object.
        metadata
            Supports the metadata options in :any:`from_from_arrays_with_metadata`.
            In addition, 'resample_resolution' can supplied to use a particular resolution using gdal "MODE=RESAMPLED_GRID"
        Returns
        -------

        """
        # @todo update method docstring for possible metadata fields
        if metadata is None:
            metadata = {}
        else:
            metadata = metadata.copy()

        if isinstance(bagfile, gdal.Dataset):
            bag = bagfile
        else:
            bag = gdal.Open(str(bagfile))

        # check for and resample variable resolution BAG if able
        gdal_metadata = bag.GetMetadata()
        if 'HAS_SUPERGRIDS' in gdal_metadata and gdal_metadata['HAS_SUPERGRIDS'] == 'TRUE':
            bag_filename = bag.GetFileList()[0]
            if "resample_resolution" in metadata:
                res = metadata["resample_resolution"]
                bag = None
                bag = gdal.OpenEx(bag_filename, open_options=['MODE=RESAMPLED_GRID', f'RESX={res}', f'RESY={res}'])
            else:
                warnings.warn(f'No resampling resolution provided for variable resolution bag {bag_filename}.  Using overview resolution.',
                              category=RuntimeWarning)

        # populate the issueDate if possible from a simple string search
        if 'issueDate' not in metadata:
            xml_str = bag.GetMetadata('xml:BAG')[0]
            root = et.fromstring(xml_str)
            elem = root.find(".//" + gco + "Date")
            if elem is not None and elem.text:
                metadata['issueDate'] = elem.text

        self.load_gdal(bag, metadata=metadata, flip_z=self.z_down)

    def to_geotiff(self, output_path: (str, pathlib.Path), creation_options: list=None):
        """ Creates a GeoTIFF file from the S102File object.

        Parameters
        ----------
        output_path
            Full path of the desired geotiff file with extension
        creation_options
            List of GDAL creation options

        Returns
        -------
        None
        """
        # this is a backport from 2.2 (RAT - feature attribute table) so doesn't do much than allow specifying the output path
        #  since there is no feature attribute table in 2.0
        instances = list(self.to_raster_datasets())
        if len(instances) > 1:
            raise NotImplementedError('Only one coverage per S102 file is expected')
        dataset, group_instance, flipx, flipy = instances[0]

        # Output the geotiff with creation options if provided
        if creation_options is None:
            creation_options = []
        tiff_ds = gdal.GetDriverByName('GTiff').CreateCopy(str(output_path), dataset, options=creation_options)

    def to_geotiffs(self, output_path: (str, pathlib.Path), creation_options: list=None):
        # there is only one coverage in an S102 file
        split_path = os.path.split(self.filename)
        filename = os.path.splitext(split_path[1])
        name = os.path.join(str(output_path), f'{filename[0]}.tif')
        self.to_geotiff(name, creation_options=None)
        return [name]

    @property
    def depth(self):
        return self.root.bathymetry_coverage.bathymetry_coverage[0].bathymetry_group[0].values.depth

    @property
    def uncertainty(self):
        return self.root.bathymetry_coverage.bathymetry_coverage[0].bathymetry_group[0].values.uncertainty

    @property
    def grid_origin_latitude(self):
        return self.root.bathymetry_coverage.bathymetry_coverage[0].grid_origin_latitude

    @property
    def grid_origin_longitude(self):
        return self.root.bathymetry_coverage.bathymetry_coverage[0].grid_origin_longitude

    @property
    def grid_spacing_latitudinal(self):
        return self.root.bathymetry_coverage.bathymetry_coverage[0].grid_spacing_latitudinal

    @property
    def grid_spacing_longitudinal(self):
        return self.root.bathymetry_coverage.bathymetry_coverage[0].grid_spacing_longitudinal


# # S102File = S102File_2_0
# def S102File(name, *args, version=2.1, **kwargs):
#     obj = None
#     if version == 2.1:
#         obj = S102File_2_1(name, *args, **kwargs)
#     elif version == 2.0:
#         obj = S102File_2_0(name, *args, **kwargs)
#     elif version < 2.0:
#         raise NotImplementedError("Version 1.x of S102 is not supported")
#     else:
#         raise ValueError(f"Version {version} is not supported in this version of the S100py module")
#     return obj
