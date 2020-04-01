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

try:
    from .. import s1xx
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py.s102"

from ..s1xx import s1xx_sequence, S1xxAttributesBase, S1xxMetadataListBase, S1xxDatasetBase, S1xxGridsBase, S1XXFile
from ..s100 import GridCoordinate, DirectPosition, GeographicBoundingBox, GeographicExtent, GridEnvelope, SequenceRule, VertexPoint, \
    FeatureInformation, FeatureInformationDataset, FeatureContainer, S100Root, S100Exception, FeatureInstanceDCF2


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

def get_valid_epsg() -> list:
    """
    Create and return the list of valid EPSG codes for S-102 version 2.0.
    """
    valid_epsg = [4326, 5041, 5042]
    valid_epsg += list(numpy.arange(32601, 32660 + 1))
    valid_epsg += list(numpy.arange(32701, 32760 + 1))
    return valid_epsg


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
        return self._attributes[self.depth_attribute_name]

    @depth.setter
    def depth(self, val: s1xx_sequence):
        self._attributes[self.depth_attribute_name] = val

    @property
    def depth_type(self) -> s1xx_sequence:
        return numpy.ndarray

    def depth_create(self):
        """ Creates a blank, empty or zero value for depth"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.depth = self.depth_type([], numpy.float)

    @property
    def uncertainty(self) -> s1xx_sequence:
        return self._attributes[self.uncertainty_attribute_name]

    @uncertainty.setter
    def uncertainty(self, val: s1xx_sequence):
        self._attributes[self.uncertainty_attribute_name] = val

    @property
    def uncertainty_type(self) -> s1xx_sequence:
        return numpy.ndarray

    def uncertainty_create(self):
        """ Creates a blank, empty or zero value for uncertainty"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.uncertainty = self.uncertainty_type([], numpy.float)

    def get_write_order(self):
        return [self.depth_attribute_name, self.uncertainty_attribute_name]


# @TODO -- determine where this is supposed to go.
# @TODO -- in some ways it seems to be the 'feature container' and other ways it's like the 'feature instance'
# which are the root/BathymetryCoverage and the root/BathymetryCoverage/BathymetryCoverage.01  groups respectively.
# it looks like NAVO interprets it as the Group.001 datastructure which is also possible.
# I'll go with Group.001 for now
class BathymetryCoverage(S1xxAttributesBase):
    """ 4.2.1.1.1 and Figure 4.4 of v2.0.0
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
        return self._attributes[self.maximum_depth_attribute_name]

    @maximum_depth.setter
    def maximum_depth(self, val: float):
        self._attributes[self.maximum_depth_attribute_name] = val

    @property
    def minimum_display_scale(self) -> int:
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
        """sequence of s102 Vectors  4.2.1.1.1.9 in S102 v2.0.0
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
    def extent(self) -> GridEnvelope:
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
        """ The attribute sequencingRule has the value class CV_SequenceRule (ISO 19123) that shall
        describe how the grid points are ordered for association to the elements of the sequence values.
        The default value is "Linear".
        No other options are allowed.
        (note that for S100: Only the values "linear" (for a simple regular cell size grid) and "Morton" (for a
        Quad Tree Grid) shall be used for data that conforms to this standard.)
        Returns
        -------

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
        """ The attribute startSequence has the value class CV_GridCoordinate that shall identify the grid
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

    # grid_matrix_attribute_name = "gridMatrix"  #: HDF5 naming
    #
    # @property
    # def grid_matrix_type(self):
    #     return BathymetryValuesList
    #
    # def grid_matrix_remove(self):
    #     self._remove_attr(self.grid_matrix_attribute_name)
    #
    # def grid_matrix_create(self):
    #     self.grid_matrix = self.grid_matrix_type()
    #
    # @property
    # def grid_matrix(self) -> S102MetadataListBase:
    #     """
    #     Returns a BathymetryValuesList object
    #     -------
    #
    #     """
    #     return self._attributes[self.grid_matrix_attribute_name]
    #
    # @grid_matrix.setter
    # def grid_matrix(self, val: S102MetadataListBase):
    #     self._attributes[self.grid_matrix_attribute_name] = val


class SurfaceCorrectionValues(VertexPoint):
    pass


# this S102_SurfaceCorrectionValues may not be right as the docs refer to S100_VertexPoint,
# but I'm betting they would both change if there ever is an update
class TrackingListValues(SurfaceCorrectionValues):
    """ 4.2.1.1.10 of v2.0.0
    """
    track_code_attribute_name = "trackCode"  #: HDF5 naming
    list_series_attribute_name = "listSeries"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def track_code(self) -> str:
        """ The optional attribute trackCode has the value type CharacterString which may contain a text string
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
        """ The attribute listSeries has the value type Integer which contains an index number into a list of metadata
        elements describing the reason for the override of the corresponding depth and uncertainty values in the bathymetry coverage.

        Returns
        -------

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


class TrackingListCoverage(S1xxAttributesBase):
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

    @property
    def common_point_rule_type(self):
        return str

    def common_point_rule_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.common_point_rule = self.common_point_rule_type("average")

    @property
    def common_point_rule(self) -> str:
        """ From S100 sec 8-7.1.2 S100_PointCoverage Spatial Model (also see figure 8-21 in edition 4.0.0):
        An S100_Point Coverage is a type of CV_DiscretePointCoverage from ISO 19123. The
        attribute values in the value record for each CV_GeometryValuePair represent values of the
        coverage, such as bathymetric soundings.
        The class S100_Point Coverage (Figure 8-21) represents a set of values, such as
        bathymetric depth values, assigned to a set of arbitrary X,Y points. Each point is identified by
        a horizontal coordinate geometry pair (X,Y) and assigned one or more values as attribute
        values. These values are organized in a record for each point.

        Returns
        -------

        """
        return self._attributes[self.common_point_rule_attribute_name]

    @common_point_rule.setter
    def common_point_rule(self, val: str):
        self._attributes[self.common_point_rule_attribute_name] = val

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


class BathymetryContainer(FeatureContainer):
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


class TrackingListContainer(FeatureContainer):
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
    def metadata_type(self) -> Type[FeatureInformation]:
        return FeatureInformation


class TrackingListCoverageDataset(S102FeatureInformationDataset):
    @property
    def metadata_name(self) -> str:
        return TRACKING_COVERAGE


class BathymetryCoverageDataset(S102FeatureInformationDataset):
    @property
    def metadata_name(self) -> str:
        return BATHY_COVERAGE


class FeatureCodes(S1xxAttributesBase):
    """ Table 10.1 and sect 10.2.1 of v2.0.0
    """

    feature_name_attribute_name = "featureName"  #: HDF5 naming
    feature_code_attribute_name = "featureCode"  #: HDF5 naming
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
        self.feature_name = self.feature_name_type([BATHY_COVERAGE, TRACKING_COVERAGE], dtype='S')

    @property
    def feature_name(self) -> s1xx_sequence:
        return self._attributes[self.feature_name_attribute_name]

    @feature_name.setter
    def feature_name(self, val: s1xx_sequence):
        self._attributes[self.feature_name_attribute_name] = val

    @property
    def feature_code_type(self):
        return numpy.array

    def feature_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_code = self.feature_code_type([BATHY_COVERAGE, TRACKING_COVERAGE], dtype='S')

    @property
    def feature_code(self) -> s1xx_sequence:
        """  This is most likely an error in the NAVO software -- I believe it should be named 'featureName' instead.
        This code should be removed if that is correct.  Otherwise remove the featureName code.
        """
        return self._attributes[self.feature_code_attribute_name]

    @feature_code.setter
    def feature_code(self, val: s1xx_sequence):
        self._attributes[self.feature_code_attribute_name] = val

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


class TilingScheme(S1xxAttributesBase):
    """ 4.2.2, table 4.1 in v2.0.0
    """
    tiling_scheme_type_attribute_name = "tilingSchemeType"  #: HDF5 naming
    domain_extent_attribute_name = "domainExtent"  #: HDF5 naming
    range_type_attribute_name = "rangeType"  #: HDF5 naming
    common_point_rule_attribute_name = "commonPointRule"  #: HDF5 naming
    geometry_attribute_name = "geometry"  #: HDF5 naming
    interpolation_type_attribute_name = "interpolationType"  #: HDF5 naming
    dimension_attribute_name = "dimension"  #: HDF5 naming
    axis_names_attribute_name = "axisNames"  #: HDF5 naming
    origin_attribute_name = "origin"  #: HDF5 naming
    offset_vectors_attribute_name = "offsetVectors"  #: HDF5 naming
    extent_attribute_name = "extent"  #: HDF5 naming
    sequencing_rule_attribute_name = "sequencingRule"  #: HDF5 naming
    start_sequence_attribute_name = "startSequence"  #: HDF5 naming

    @property
    def tiling_scheme_type(self) -> Type[str]:
        return self._attributes[self.tiling_scheme_type_attribute_name]

    @tiling_scheme_type.setter
    def tiling_scheme_type(self, val: str):
        self._attributes[self.tiling_scheme_type_attribute_name] = val


class DiscoveryMetadata(S1xxAttributesBase):
    """ 12.1 of v2.0.0
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        raise NotImplementedError()


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
