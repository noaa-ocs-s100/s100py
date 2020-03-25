import os
import subprocess
import datetime
import logging
import functools
from typing import Callable, Iterator, Union, Optional, List, Type

import xml
import xml.etree.ElementTree
import pprint

import numpy

try:
    from . import s1xx
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py"

from .s1xx import s1xx_sequence, S1XX_Attributes_base, S1XX_MetadataList_base, S1XX_Dataset_base, S1XX_WritesOwnGroup_base, S1XXFile
from .s100 import GridCoordinate, DirectPosition, GeographicBoundingBox, GeographicExtent, GridEnvelope, SequenceRule, VertexPoint
from .s100 import FeatureInformation, S100_FeatureContainer, S100Root, FeatureInstance_Format_2

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

# As an example and to test compatibility a derived attributes class is made for an older NAVO format


# override the basic S100 spec that says to use an underscore and use a dot instead
class S102_MetadataList_base(S1XX_MetadataList_base):
    write_format_str = ".%03d"


# # @TODO -- determine if this is old.  The spec seems to describe a one dimensional array or list of points but the values in the grid is a 2 x N x M dataset
# class BathymetryValueRecord(S1XX_Attributes_base):
#     """ 4.2.1.1.2.2 and Figure 4.4 of v2.0.0
#     The attribute values has the value type S102_BathymetryValueRecord which is a sequence of value items that
#     shall assign values to the grid points.
#     There are two attributes in the bathymetry value record, depth and uncertainty in the S102_BathymetryValues class.
#     The definition for the depth is defined by the depthCorrectionType attribute in the S102_DataIdentification class.
#     The definition of the type of data in the values record is defined by the verticalUncertaintyType attribute in the
#     S102_DataIdentification class.
#     """
#
#     @property
#     def __version__(self) -> int:
#         return 1
#
#     @property
#     def __version__(self) -> int:
#         return 1
#
#     @property
#     def depth_attribute_name(self) -> str:
#         return "depth"
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
#     @property
#     def uncertainty_attribute_name(self) -> str:
#         return "uncertainty"
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
# class BathymetryValuesList(S102_MetadataList_base):
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


class BathymetryValues(S1XX_WritesOwnGroup_base):
    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "values"

    @property
    def depth_attribute_name(self) -> str:
        return "depth"

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
        self.depth = self.depth_type([], numpy.float)

    @property
    def uncertainty_attribute_name(self) -> str:
        return "uncertainty"

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
        self.uncertainty = self.uncertainty_type([], numpy.float)

    def get_write_order(self):
        return [self.depth_attribute_name, self.uncertainty_attribute_name]

    def read(self, group_object, indent=0):
        logging.debug("reading depth/uncertainty matrices")
        self.depth = group_object[self.depth_attribute_name]
        self.uncertainty = group_object[self.uncertainty_attribute_name]

    def write(self, group_object, indent=0):
        # @todo - is there a bug here if some instances are missing attributes leading to a mismatched array?
        """ Write out the dataset using order specified with any extra values as unordered but named at the end.

        Parameters
        ----------
        group_object
            HDF5 object to write into
        indent

        Returns
        -------
        HDF5 dataset created during the write method
        """

        try:
            # First determine the write order of the keys
            logging.debug(indent * "  " + "Writing" + " " + str(self))

            dataset = None

            write_keys = []
            if self.get_write_order():  # @todo I think bathycoverage and trackingcoverage in the feature information may want to be ordered
                write_keys.extend(self.get_write_order())

            # to preserve order of other keys - iterate instead of using set logic
            for key in self._attributes:
                if key not in write_keys:
                    write_keys.append(key)
            # write_keys.extend(set(self._attributes.keys()).difference(write_keys))
            write_array = [self._attributes[key] for key in write_keys]

            # hdf5 needs names to the columns which is done in a record array or structured array.
            # but to create that without specifying type we need to transpose first then call 'fromarrays'

            # numpy.array is coming out with wrong (at least different) shape and fromarrays is working -- not sure why right now.
            # rec_array = numpy.array(write_array, dtype=[(name, 'f4') for name in write_keys])
            rec_array = numpy.core.records.fromarrays(write_array, dtype=[(name, 'f4') for name in write_keys])
            dataset = group_object.create_dataset(self.metadata_name, data=rec_array)
            return dataset
        except Exception as e:
            raise e


# @TODO -- determine where this is supposed to go.
# @TODO -- in some ways it seems to be the 'feature container' and other ways it's like the 'feature instance'
# which are the root/BathymetryCoverage and the root/BathymetryCoverage/BathymetryCoverage.01  groups respectively.
# it looks like NAVO interprets it as the Group.001 datastructure which is also possible.
# I'll go with Group.001 for now
class BathymetryCoverage(S1XX_Attributes_base):
    """ 4.2.1.1.1 and Figure 4.4 of v2.0.0
    also see section 12.3 and table 12.5

    """
    write_format_str = ".%03d"

    @property
    def values_attribute_name(self) -> str:
        return "values"

    @property
    def values(self) -> BathymetryValues:
        return self._attributes[self.values_attribute_name]

    @values.setter
    def values(self, val: BathymetryValues):
        self._attributes[self.values_attribute_name] = val

    @property
    def values_type(self) -> BathymetryValues:
        return BathymetryValues

    def values_create(self):
        """ Creates a blank, empty or zero value for values"""
        self.values = self.values_type()

    @property
    def __version__(self) -> int:
        return 1

    @property
    def minimum_depth_attribute_name(self) -> str:
        return "minimumDepth"

    @property
    def minimum_depth_type(self):
        return float

    def minimum_depth_create(self):
        self.minimum_depth = self.minimum_depth_type()

    @property
    def minimum_depth(self) -> float:
        return self._attributes[self.minimum_depth_attribute_name]

    @minimum_depth.setter
    def minimum_depth(self, val: float):
        self._attributes[self.minimum_depth_attribute_name] = val

    @property
    def maximum_depth_attribute_name(self) -> str:
        return "maximumDepth"

    @property
    def maximum_depth_type(self):
        return float

    def maximum_depth_create(self):
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
    def maximum_display_scale_attribute_name(self) -> str:
        return "maximumDisplayScale"

    @property
    def maximum_display_scale_type(self):
        return float

    def maximum_display_scale_create(self):
        self.maximum_display_scale = self.maximum_display_scale_type()

    @property
    def maximimum_display_scale(self) -> int:
        return self._attributes[self.maximimum_display_scale_attribute_name]

    @maximimum_display_scale.setter
    def maximimum_display_scale(self, val: int):
        self._attributes[self.maximimum_display_scale_attribute_name] = val

    @property
    def minimum_display_scale_attribute_name(self) -> str:
        return "minimumDisplayScale"

    @property
    def minimum_display_scale_type(self):
        return float

    def minimum_display_scale_create(self):
        self.minimum_display_scale = self.minimum_display_scale_type()

    @property
    def minimum_uncertainty(self) -> float:
        return self._attributes[self.minimum_uncertainty_attribute_name]

    @minimum_uncertainty.setter
    def minimum_uncertainty(self, val: float):
        self._attributes[self.minimum_uncertainty_attribute_name] = val

    @property
    def minimum_uncertainty_attribute_name(self) -> str:
        return "minimumUncertainty"

    @property
    def minimum_uncertainty_type(self):
        return float

    def minimum_uncertainty_create(self):
        self.minimum_uncertainty = self.minimum_uncertainty_type()

    @property
    def maximum_uncertainty(self) -> float:
        return self._attributes[self.maximum_uncertainty_attribute_name]

    @maximum_uncertainty.setter
    def maximum_uncertainty(self, val: float):
        self._attributes[self.maximum_uncertainty_attribute_name] = val

    @property
    def maximum_uncertainty_attribute_name(self) -> str:
        return "maximumUncertainty"

    @property
    def maximum_uncertainty_type(self):
        return float

    def maximum_uncertainty_create(self):
        self.maximum_uncertainty = self.maximum_uncertainty_type()

    @property
    def origin(self) -> DirectPosition:
        return self._attributes[self.origin_attribute_name]

    @origin.setter
    def origin(self, val: DirectPosition):
        self._attributes[self.origin_attribute_name] = val

    @property
    def origin_attribute_name(self) -> str:
        return "origin"

    @property
    def origin_type(self):
        return DirectPosition

    def origin_create(self):
        self.origin = self.origin_type()

    @property
    def origin_attribute_type(self) -> Type[str]:
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
    def offset_vectors_attribute_name(self) -> str:
        return "offsetVectors"

    @property
    def offset_vectors_type(self):
        return numpy.ndarray

    def offset_vectors_create(self):
        self.offset_vectors = self.offset_vectors_type([2, ], numpy.float64)

    @property
    def dimension(self) -> int:
        return self._attributes[self.dimension_attribute_name]

    @dimension.setter
    def dimension(self, val: int):
        self._attributes[self.dimension_attribute_name] = val

    @property
    def dimension_attribute_name(self) -> str:
        return "dimension"

    @property
    def dimension_type(self):
        return int

    def dimension_create(self):
        self.dimension = self.dimension_type(2)

    @property
    def axis_names(self) -> s1xx_sequence:
        """sequence of character strings"""
        return self._attributes[self.axis_names_attribute_name]

    @axis_names.setter
    def axis_names(self, val: s1xx_sequence):
        self._attributes[self.axis_names_attribute_name] = val

    @property
    def axis_names_attribute_name(self) -> str:
        return "axisNames"

    @property
    def axis_names_type(self) -> Type[str]:
        return numpy.ndarray

    def axis_names_create(self):
        """ The attribute axisNames has the value class Sequence<CharacterString> that shall be used to assign names to the grid axis.
        The grid axis names shall be "Latitude" and "Longitude" for unprojected data sets or “Northing” and “Easting” in a projected space.

        Returns
        -------

        """
        self.axis_names = self.axis_names_type([2], dtype='S')

    @property
    def extent(self) -> GridEnvelope:
        return self._attributes[self.extent_attribute_name]

    @extent.setter
    def extent(self, val: GridEnvelope):
        self._attributes[self.extent_attribute_name] = val

    @property
    def extent_attribute_name(self) -> str:
        return "extent"

    @property
    def extent_type(self) -> Type[str]:
        return GridEnvelope

    def extent_create(self):
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
    def sequencing_rule_attribute_name(self) -> str:
        return "sequencingRule"

    @property
    def sequencing_rule_type(self):
        return SequenceRule

    def sequencing_rule_create(self):
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
    def start_sequence_attribute_name(self) -> str:
        return "startSequence"

    @property
    def start_sequence_type(self):
        return GridCoordinate

    def start_sequence_create(self):
        self.start_sequence = self.start_sequence_type()

    # @property
    # def grid_matrix_attribute_name(self) -> str:
    #     return "gridMatrix"
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
    # def grid_matrix(self) -> S102_MetadataList_base:
    #     """
    #     Returns a BathymetryValuesList object
    #     -------
    #
    #     """
    #     return self._attributes[self.grid_matrix_attribute_name]
    #
    # @grid_matrix.setter
    # def grid_matrix(self, val: S102_MetadataList_base):
    #     self._attributes[self.grid_matrix_attribute_name] = val


class SurfaceCorrectionValues(VertexPoint):
    pass


# this S102_SurfaceCorrectionValues may not be right as the docs refer to S100_VertexPoint,
# but I'm betting they would both change if there ever is an update
class TrackingListValues(SurfaceCorrectionValues):
    """ 4.2.1.1.10 of v2.0.0
    """

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
    def track_code_attribute_name(self) -> str:
        return "trackCode"

    @property
    def track_code_type(self):
        return str

    def track_code_create(self):
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
    def list_series_attribute_name(self) -> str:
        return "listSeries"

    @property
    def list_series_type(self):
        return int

    def list_series_create(self):
        self.list_series = self.list_series_type()


class TrackingListValues_List(S102_MetadataList_base):
    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "point"

    @property
    def metadata_type(self) -> type:
        return TrackingListValues


class TrackingListSet_List(S102_MetadataList_base):
    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "set"

    @property
    def metadata_type(self) -> type:
        return TrackingListValues_List


class TrackingListCoverage(S1XX_Attributes_base):
    """ 4.2.1.1.9 and Figure 4.4 of v2.0.0
    commonPointRule is defined to be an S100_PointCoverage with a value of default and it therefore optional.
    a metadata attribute from S100 is allowed but not necessary as well.

    """
    write_format_str = ".%02d"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def domain_extent(self) -> S102_MetadataList_base:
        return self._attributes[self.domain_extent_attribute_name]

    @domain_extent.setter
    def domain_extent(self, val: S102_MetadataList_base):
        self._attributes[self.domain_extent_attribute_name] = val

    @property
    def domain_extent_attribute_name(self) -> str:
        return "domainExtent"

    @property
    def domain_extent_type(self):
        return GeographicExtent

    def domain_extent_create(self):
        self.domain_extent = self.domain_extent_type()

    @property
    def common_point_rule_attribute_name(self) -> str:
        return "commonPointRule"

    @property
    def common_point_rule_type(self):
        return str

    def common_point_rule_create(self):
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
    def set_attribute_name(self) -> str:
        return "set"

    @property
    def set_type(self):
        return TrackingListSet_List

    def set_create(self):
        self.set = self.set_type()

    @property
    def set(self) -> S102_MetadataList_base:
        """
        Returns
        -------
        TrackingListValues_List
            list of TrackingListValues
        """
        return self._attributes[self.set_attribute_name]

    @set.setter
    def set(self, val: S102_MetadataList_base):
        self._attributes[self.set_attribute_name] = val

    # @TODO  I don't think this is right, but not sure where I found it
    # @property
    # def geometry_attribute_name(self) -> str:
    #     return "geometry"
    #
    # @property
    # def geometry(self) -> S1XX_Attributes_base:
    #     return self._attributes[self.geometry_attribute_name]
    #
    # @geometry.setter
    # def geometry(self, val: S1XX_Attributes_base):
    #     self._attributes[self.geometry_attribute_name] = val
    #
    # @property
    # def value_attribute_name(self) -> str:
    #     return "value"
    #
    # @property
    # def value(self) -> s1xx_sequence:
    #     return self._attributes[self.value_attribute_name]
    #
    # @value.setter
    # def value(self, val: s1xx_sequence):
    #     self._attributes[self.value_attribute_name] = val


class TrackingListGroup_List(S102_MetadataList_base):
    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "Group"

    @property
    def metadata_type(self) -> type:
        return TrackingListCoverage


class BathymetryGroup_List(S102_MetadataList_base):
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


class TrackingListCoverages_List(S102_MetadataList_base):
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
        return TrackingListGroup_List




class BathymetryFeatureInstance(FeatureInstance_Format_2):
    @property
    def bathymetry_group_attribute_name(self) -> str:
        """ attribute name will be automatically determined based on the array position of the S102_MetadataList
        Returns
        -------
        Basic template for the name of the attribute
        """
        return "Group" + r"\.\d+"

    @property
    def bathymetry_group_type(self):
        return BathymetryGroup_List

    def bathymetry_group_create(self):
        self.bathymetry_group = self.bathymetry_group_type()

    @property
    def bathymetry_group(self) -> S102_MetadataList_base:
        """ The bathymetry data, a list of Bathymetrygroup
        Returns
        -------
        S102_MetadataList_base
            Contains a list of BathymetryCoverage objects via the BathymetryCoverages_List class
        """
        return self._attributes[self.bathymetry_group_attribute_name]

    @bathymetry_group.setter
    def bathymetry_group(self, val: S102_MetadataList_base):
        self._attributes[self.bathymetry_group_attribute_name] = val


class BathymetryCoverages_List(S102_MetadataList_base):
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
    def metadata_type(self) -> Type[type]:
        return BathymetryFeatureInstance




class BathymetryContainer(S100_FeatureContainer):
    """ This is the BathymetryCoverage right off the root of the HDF5 which has possible attributes from S100 spec table 10c-10
    This will hold child groups named BathymetryCoverage.NN
    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def bathymetry_coverage_attribute_name(self) -> str:
        """ attribute name will be automatically determined based on the array position of the S102_MetadataList
        Returns
        -------
        Basic template for the name of the attribute
        """
        return BATHY_COVERAGE + r"\.\d+"

    @property
    def bathymetry_coverage_type(self):
        return BathymetryCoverages_List

    def bathymetry_coverage_create(self):
        self.bathymetry_coverage = self.bathymetry_coverage_type()

    @property
    def bathymetry_coverage(self) -> S102_MetadataList_base:
        """ The bathymetry data, a list of BathymetryCoverage
        Returns
        -------
        S102_MetadataList_base
            Contains a list of BathymetryCoverage objects via the BathymetryCoverages_List class
        """
        return self._attributes[self.bathymetry_coverage_attribute_name]

    @bathymetry_coverage.setter
    def bathymetry_coverage(self, val: S102_MetadataList_base):
        self._attributes[self.bathymetry_coverage_attribute_name] = val

    def data_coding_format_create(self):
        """ Creates a blank, empty or zero value for data_coding_format"""
        self.data_coding_format = self.data_coding_format_type(2)

    def dimension_create(self):
        """ Creates a blank, empty or zero value for dimension"""
        self.dimension = self.dimension_type(2)


class TrackingListContainer(S100_FeatureContainer):
    """
    Table 10.1 of v2.0.0
    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def tracking_list_coverage_attribute_name(self) -> str:
        return TRACKING_COVERAGE + r"\.\d+"

    @property
    def tracking_list_coverage_type(self):
        return TrackingListCoverages_List

    def tracking_list_coverage_create(self):
        self.tracking_list_coverage = self.tracking_list_coverage_type()

    @property
    def tracking_list_coverage(self) -> S1XX_Attributes_base:
        """ The tracking list data, a list of TrackingListCoverage
        Returns
        -------
        S102_MetadataList_base
            Contains a list of TrackingListCoverage objects via the TrackingListCoverages_List class
        """
        return self._attributes[self.tracking_list_coverage_attribute_name]

    @tracking_list_coverage.setter
    def tracking_list_coverage(self, val: S1XX_Attributes_base):
        self._attributes[self.tracking_list_coverage_attribute_name] = val


class FeatureInformation_dataset(S1XX_Dataset_base):
    """   In S102, 10.2.1 and table 10.2 and Table 10.1 of v2.0.0

    This is used to describe the BathymetryCoverage and TrackingListCoverage within the GroupF feature listing.
    The features described under GroupF have a matching named entry parallel to GroupF (top level).
    The actual data (depths etc) is stored in the top level element while basic metadata is stored in this element.

    """
    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_type(self) -> Type[type]:
        return FeatureInformation


class TrackingListCoverage_dataset(FeatureInformation_dataset):
    @property
    def metadata_name(self) -> str:
        return TRACKING_COVERAGE


class BathymetryCoverage_dataset(FeatureInformation_dataset):
    @property
    def metadata_name(self) -> str:
        return BATHY_COVERAGE


class FeatureCodes(S1XX_Attributes_base):
    """ Table 10.1 and sect 10.2.1 of v2.0.0
    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def feature_name_attribute_name(self) -> str:
        return "featureName"

    @property
    def feature_name_type(self):
        return numpy.array

    def feature_name_create(self):
        self.feature_name = self.feature_name_type([BATHY_COVERAGE, TRACKING_COVERAGE], dtype='S')

    @property
    def feature_name(self) -> s1xx_sequence:
        return self._attributes[self.feature_name_attribute_name]

    @feature_name.setter
    def feature_name(self, val: s1xx_sequence):
        self._attributes[self.feature_name_attribute_name] = val

    @property
    def feature_code_attribute_name(self) -> str:
        return "featureCode"

    @property
    def feature_code_type(self):
        return numpy.array

    def feature_code_create(self):
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
    def bathymetry_coverage_dataset_attribute_name(self) -> str:
        return BATHY_COVERAGE

    @property
    def bathymetry_coverage_dataset_type(self):
        return BathymetryCoverage_dataset

    def bathymetry_coverage_dataset_create(self):
        self.bathymetry_coverage_dataset = self.bathymetry_coverage_dataset_type()

    @property
    def bathymetry_coverage_dataset(self) -> BathymetryCoverage_dataset:
        return self._attributes[self.bathymetry_coverage_dataset_attribute_name]

    @bathymetry_coverage_dataset.setter
    def bathymetry_coverage_dataset(self, val: BathymetryCoverage_dataset):
        self._attributes[self.bathymetry_coverage_dataset_attribute_name] = val

    @property
    def tracking_list_coverage_attribute_name(self) -> str:
        return TRACKING_COVERAGE

    @property
    def tracking_list_coverage_type(self):
        return TrackingListCoverage_dataset

    def tracking_list_coverage_create(self):
        self.tracking_list_coverage = self.tracking_list_coverage_type()

    @property
    def tracking_list_coverage(self) -> TrackingListCoverage_dataset:
        return self._attributes[self.tracking_list_coverage_attribute_name]

    @tracking_list_coverage.setter
    def tracking_list_coverage(self, val: TrackingListCoverage_dataset):
        self._attributes[self.tracking_list_coverage_attribute_name] = val


class S102Root(S100Root):
    """The root group contains a feature information group and N feature containers.
    In S102 there are currently two feature containers which are the 'coverages'  bathymetry and tracking list.
    The coverage names are determined from the matching CoveragesAttributes
    10.2 and Figure 10.1 of v2.0.0
    """

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
        self.feature_information = self.feature_information_type()

    @property
    def feature_information_attribute_name(self) -> str:
        return "Group_F"

    @property
    def bathymetry_coverage_attribute_name(self) -> str:
        return BATHY_COVERAGE

    @property
    def bathymetry_coverage(self) -> S1XX_Attributes_base:
        return self._attributes[self.bathymetry_coverage_attribute_name]

    @property
    def bathymetry_coverage_type(self):
        return BathymetryContainer

    def bathymetry_coverage_create(self):
        self.bathymetry_coverage = self.bathymetry_coverage_type()

    @bathymetry_coverage.setter
    def bathymetry_coverage(self, val: S1XX_Attributes_base):
        self._attributes[self.bathymetry_coverage_attribute_name] = val

    @property
    def tracking_list_coverage_attribute_name(self) -> str:
        return TRACKING_COVERAGE

    @property
    def tracking_list_coverage_type(self):
        return TrackingListContainer

    def tracking_list_coverage_create(self):
        self.tracking_list_coverage = self.tracking_list_coverage_type()

    @property
    def tracking_list_coverage(self) -> S1XX_Attributes_base:
        return self._attributes[self.tracking_list_coverage_attribute_name]

    @tracking_list_coverage.setter
    def tracking_list_coverage(self, val: S1XX_Attributes_base):
        self._attributes[self.tracking_list_coverage_attribute_name] = val


class TilingScheme(S1XX_Attributes_base):
    """ 4.2.2, table 4.1 in v2.0.0
    """

    @property
    def tiling_scheme_type(self) -> Type[str]:
        return self._attributes[self.tiling_scheme_type_attribute_name]

    @tiling_scheme_type.setter
    def tiling_scheme_type(self, val: str):
        self._attributes[self.tiling_scheme_type_attribute_name] = val

    @property
    def tiling_scheme_type_attribute_name(self) -> str:
        return "tilingSchemeType"

    @property
    def domain_extent_attribute_name(self) -> str:
        return "domainExtent"

    @property
    def range_type_attribute_name(self) -> str:
        return "rangeType"

    @property
    def common_point_rule_attribute_name(self) -> str:
        return "commonPointRule"

    @property
    def geometry_attribute_name(self) -> str:
        return "geometry"

    @property
    def interpolation_type_attribute_name(self) -> str:
        return "interpolationType"

    @property
    def dimension_attribute_name(self) -> str:
        return "dimension"

    @property
    def axis_names_attribute_name(self) -> str:
        return "axisNames"

    @property
    def origin_attribute_name(self) -> str:
        return "origin"

    @property
    def offset_vectors_attribute_name(self) -> str:
        return "offsetVectors"

    @property
    def extent_attribute_name(self) -> str:
        return "extent"

    @property
    def sequencing_rule_attribute_name(self) -> str:
        return "sequencingRule"

    @property
    def start_sequence_attribute_name(self) -> str:
        return "startSequence"


class DiscoveryMetadata(S1XX_Attributes_base):
    """ 12.1 of v2.0.0
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        raise NotImplementedError()


class S102File(S1XXFile):
    PRODUCT_SPECIFICATION = numpy.string_('INT.IHO.S-102.2.0')
    def __init__(self, *args, **kywrds):
        # kywrds['root'] = S102Root
        super().__init__(*args, root=S102Root, **kywrds)


def make_s102(output_path, elev_raster, uncert_raster, metadata):
    """
    Create an S-102 file from the elevation and uncertainty rasters and a
    metadata dictionary.
    """
    
    sfile = S102File(output_path, "w", driver=None)
    sfile.create_empty_metadata()  # init the root with a fully filled out empty metadata set
    root = sfile.root
    del root.feature_information.feature_code  # Guessing at the right dataset name to keep (keep featureName)
    # root.feature_information.feature_code_remove()
    root.product_specification = "INT.IHO.S-102.2.0.0"
    root.metadata = os.path.splitext(os.path.basename(output_path))[0] + ".xml"
    root.vertical_datum = metadata['from_vert_datum']

    bathy_cov_dset = root.feature_information.bathymetry_coverage_dataset
    bathy_depth_info = bathy_cov_dset.append_new_item()  # bathy_cov_dset.append(bathy_cov_dset.metadata_type())
    bathy_depth_info.initialize_properties(True)
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
    bathy_uncertainty_info.initialize_properties(True)
    bathy_uncertainty_info.code = UNCERTAINTY
    bathy_uncertainty_info.name = UNCERTAINTY

    # I'm not sure what to put here, yet
    tracking_cov = root.feature_information.tracking_list_coverage

    track_info = tracking_cov.append_new_item()  # append(tracking_cov.metadata_type())
    track_info.initialize_properties(True)
    track_info.code = "X"
    track_info.name = "X"
    track_info.unit_of_measure = "N/A"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True)
    track_info.code = "Y"
    track_info.name = "Y"
    track_info.unit_of_measure = "N/A"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True)
    track_info.code = "originalValue"
    track_info.name = "Original Value"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True)
    track_info.code = "trackCode"
    track_info.name = "Track Code"
    track_info.unit_of_measure = "N/A"

    track_info = tracking_cov.append_new_item()
    track_info.initialize_properties(True)
    track_info.code = "listSeries"
    track_info.name = "List Series"
    track_info.unit_of_measure = "N/A"

    root.bathymetry_coverage.axis_names = numpy.array(["longitude", "latitude"])  # row major order means X/longitude first
    root.bathymetry_coverage.common_point_rule = 1  # average
    # root.bathymetry_coverage.data_coding_format = 2  # default
    # root.bathymetry_coverage.dimension = 2  # default value
    root.bathymetry_coverage.sequencing_rule_scan_direction = "Longitude, Latitude"
    root.bathymetry_coverage.interpolation_type = 1  # nearest neighbor
    root.bathymetry_coverage.num_instances = 1  # how many Bathycoverages
    root.bathymetry_coverage.sequencing_rule_type = 1  # linear
    del root.bathymetry_coverage.time_uncertainty

    bathy_01 = root.bathymetry_coverage.bathymetry_coverage.append_new_item()
    bathy_01.initialize_properties(True)

    # @todo @fixme
    print("Need to determine if projected coords or not - assuming UTM right now")
    minx = metadata['lon_min']
    maxx = metadata['lon_max']
    miny = metadata['lat_min']
    maxy = metadata['lat_max']

    minx = min([metadata['bounds'][0][0], metadata['bounds'][1][0]])
    maxx = max([metadata['bounds'][0][0], metadata['bounds'][1][0]])
    miny = min([metadata['bounds'][0][1], metadata['bounds'][1][1]])
    maxy = max([metadata['bounds'][0][1], metadata['bounds'][1][1]])

    root.east_bound_longitude = minx
    root.west_bound_longitude = maxx
    root.south_bound_latitude = miny
    root.north_bound_latitude = maxy
    bathy_01.east_bound_longitude = minx
    bathy_01.west_bound_longitude = maxx
    bathy_01.south_bound_latitude = miny
    bathy_01.north_bound_latitude = maxy
    bathy_01.grid_origin_latitude = miny

    bathy_01.grid_origin_longitude = minx
    bathy_01.grid_origin_latitude = miny
    bathy_01.grid_spacing_longitudinal = metadata["res"][0]
    bathy_01.grid_spacing_latitudinal = metadata["res"][1]
    del bathy_01.grid_spacing_vertical
    del bathy_01.grid_origin_vertical
    bathy_01.num_grp = 1

    bathy_group_object = bathy_01.bathymetry_group.append_new_item()
    # bathy_group_object.initialize_properties()  # Not creating everything as I'm not sure if the grid attributes shoul dbe thereTrue)
    print(bathy_group_object.get_standard_properties())
    # @todo  @FIXME
    print("need to determine if this is degrees/metres and set axisNames accordingly")
    # bathy_group_object.axis_names = numpy.array(["longitude", "latitude"])  # row major order means X/longitude first
    # use default dimension =2

    print("fix here -- row/column order?")
    nx, ny = metadata['shape']

    bathy_01.num_points_latitudinal = ny
    bathy_01.num_points_longitudinal = nx
    del bathy_01.num_points_vertical
    bathy_01.start_sequence = "0,0"
    del bathy_01.vertical_extent_maximum_z
    del bathy_01.vertical_extent_minimum_z

    bathy_group_object.extent_create()
    bathy_group_object.extent.initialize_properties(True)
    bathy_group_object.extent.low.coord_values[0:2] = [0, 0]
    bathy_group_object.extent.high.coord_values[0:2] = [nx, ny]

    nodata_value = metadata['nodata']
    depth_max = elev_raster[elev_raster != nodata_value].max()
    depth_min = elev_raster[elev_raster != nodata_value].min()
    bathy_group_object.maximum_depth = depth_max
    bathy_group_object.minimum_depth = depth_min

    uncertainty_max = uncert_raster[uncert_raster != nodata_value].max()
    uncertainty_min = uncert_raster[uncert_raster != nodata_value].min()
    bathy_group_object.minimum_uncertainty = uncertainty_min
    bathy_group_object.maximum_uncertainty = uncertainty_max
    bathy_group_object.dimension = 2

    bathy_group_object.origin_create()
    bathy_group_object.origin.initialize_properties(True)
    bathy_group_object.origin.dimension = 2
    bathy_group_object.origin.coordinate = numpy.array([minx, miny])

    bathy_group_object.values_create()
    grid = bathy_group_object.values
    # @todo -- need to make sure nodata values are correct, especially if converting something other than bag which is supposed to have the same nodata value
    grid.depth = elev_raster
    grid.uncertainty = uncert_raster

    # @todo @fixme hardcoded stuff....
    root.horizontal_datum_reference = "EPSG"
    root.horizontal_datum_value = metadata['epsg']
    root.epoch = "G1762"  # this is the 2013-10-16 WGS84 used by CRS
    root.geographic_identifier = metadata['geographic_location']
    root.issue_date = metadata['date_stamp']  # datetime.date.today().isoformat()

    sfile.write()