from __future__ import annotations

import os
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
from enum import Enum

import xml
import xml.etree.ElementTree
import pprint
from xml.etree import ElementTree as et

import numpy
import h5py
try:
    from osgeo import gdal, osr
    gdal.UseExceptions()
except:
    pass

try:
    from ... import s1xx, s100
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py.s102"

from ...s1xx import s1xx_sequence, S1xxObject, S1xxCollection, S1xxGridsBase, S1XXFile, h5py_string_dtype, make_enum_dtype, change_attr_type
from ...s100.v5_0.api import S100File, GridCoordinate, DirectPosition, GridEnvelope, SequenceRule, VertexPoint, \
    FeatureInformation, FeatureInformationDataset, FeatureContainerDCF2, S100Root, S100Exception, FeatureInstanceDCF2, GroupFBase, \
    CommonPointRule, FeatureInstanceDCF9, FeatureContainerDCF9, S1xxDatasetBase, InterpolationType, \
    VERTICAL_CS, VERTICAL_DATUM_REFERENCE, VERTICAL_COORDINATE_BASE, VERTICAL_DATUM

from .. import v2_0
from .. import v2_1

EDITION = 2.2
PRODUCT_SPECIFICATION = 'INT.IHO.S-102.2.2'

CHANGELOG = """
v2.1
Removed TrackingList  --  4.2.1.1.8 TrackingListCoverage
Removed min/max display scale -- 4.2.1.1.1.2 and 4.2.1.1.1.5 BathymetryCoverage semantics
Added flip_z parameters in utils since z orientation is going from positive up to positive down
Change FeatureInformation datatype to H5T_FLOAT from H5T_NATIVE_FLOAT - per table 10-3 
featureName and featureCode were both used in 2.0 doc, was corrected to only use featureCode in 2.1

v2.2
Add QualityOfSurvey for RasterAttribute storage.
Revisions to the horizontal and vertical datum attributes at the root level.
Stricter datatypes per S102 (but not S100) spec
"""


class S102Exception(S100Exception):
    pass


BATHY_COVERAGE = "BathymetryCoverage"
QUALITY_OF_SURVEY = "QualityOfSurvey"
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

class BATHYMETRIC_UNCERTAINTY_TYPE(Enum):
    """ Note: while a Vertical Datum can be created with the shorthand aliases, ex: MLWS, the string written and
    returned from the file/S100 object will be the official long name, e.g. "meanLowWaterSprings" etc.
    S100 Part 4a Metadata

    S100 v5 Part 17 Vertical and Sounding Datum
    Added balticSeaChartDatum2000 = 44
    """
    unknown = 0
    rawStandardDeviation = 1
    cUBEStandardDeviation = 2
    productUncertainty = 3
    historicalStandardDeviation = 4

# figure 4.4 in section 4.2.1 of v2.0.0 shows an overview of many of the S102 classes used
# Annex B in the S102 spec has an example layout (still determining if it matches the docs)
# S100 doc part 10C has HDF5 layout information as well.
# S100 doc 10C-7 has some representation guidelines


# override the basic S100 spec that says to use an underscore and use a dot instead
class S102MetadataListBase(S1xxCollection):
    write_format_str = ".%03d"


class BathymetryValues(S1xxGridsBase):
    """ S100 v5.0
        8-6.2.8 Grid cell structure S-100 utilizes the same view of grid cell structure as Section 8.2.2 of ISO 19123.
        The grid data in S-100 grid coverages are nominally situated exactly at the grid points defined by the grid coordinates.
        The grid points are therefore the “sample points.” Data values at a sample point represent measurements over a neighbourhood of the sample point.
        This neighbourhood is assumed to extend a half-cell in each dimension.
        The effect is that the sample space corresponding to each grid point is a cell centred at the grid point.
    """
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


# v2.1 Chagne to .01 from .001
# v2.1 removed min/max display scale mixin
class BathymetryCoverage(S1xxObject):
    """ This is the "Values" Group.NNN object that contains the grid data in a values dataset and other metadata about the grids.
    S100 v4.0 table 10c-18
    4.2.1.1.1 and Figure 4.4 of v2.0.0
    also see section 12.3 and table 12.5

    """
    # Changed to %02d in v2.1
    write_format_str = ".%02d"

    __values_hdf_name__ = "values"  #: HDF5 naming
    __minimum_depth_hdf_name__ = "minimumDepth"  #: HDF5 naming
    __maximum_depth_hdf_name__ = "maximumDepth"  #: HDF5 naming
    __minimum_uncertainty_hdf_name__ = "minimumUncertainty"  #: HDF5 naming
    __maximum_uncertainty_hdf_name__ = "maximumUncertainty"  #: HDF5 naming

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

        S100 v5.0
        8-6.2.8 Grid cell structure S-100 utilizes the same view of grid cell structure as Section 8.2.2 of ISO 19123.
        The grid data in S-100 grid coverages are nominally situated exactly at the grid points defined by the grid coordinates.
        The grid points are therefore the “sample points.” Data values at a sample point represent measurements over a neighbourhood of the sample point.
        This neighbourhood is assumed to extend a half-cell in each dimension.
        The effect is that the sample space corresponding to each grid point is a cell centred at the grid point.
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
        return numpy.float32

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
        return numpy.float32

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
        return numpy.float32

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
        return numpy.float32

    def maximum_uncertainty_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.maximum_uncertainty = self.__maximum_uncertainty_type__()


# v2.1 change to Group_001  "." to "_"
class BathymetryGroupList(S102MetadataListBase):
    """ This is the list of Group.NNN that are held as a list.
    Each Group.NNN has a dataset of depth and uncertainty.
    """
    write_format_str = "_%03d"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "Group"

    @property
    def metadata_type(self) -> type:
        return BathymetryCoverage


class BathymetryFeatureInstance(FeatureInstanceDCF2):
    """ This will be the BathymetryCoverage.001 element in HDF5.
    It will contain a Group.NNN which will have the "values" dataset of the deptha dn uncertainty.
    """
    __bathymetry_group_hdf_name__ = "Group" + r"[\._]\d+"
    """ Basic template for HDF5 naming of the attribute.  
    Attribute name will be automatically determined based on the list's index of the data. 
    """

    @property
    def __num_grp_type__(self) -> Type[int]:
        return numpy.uint8

    def num_grp_create(self):
        self.num_grp = self.__num_grp_type__(1)

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
    write_format_str = ".%02d"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return BATHY_COVERAGE

    @property
    def metadata_type(self) -> Type[BathymetryFeatureInstance]:
        return BathymetryFeatureInstance


# TODO FIXME If this is not a feature attributed grid then it should be DataCodingFormat2
#  (which has an additional interpolation attribute) and not any QualityOfSurvey structures
class BathymetryContainer(FeatureContainerDCF9, InterpolationType):
    """ This is the BathymetryCoverage right off the root of the HDF5 which has possible attributes from S100 spec table 10c-10
    This will hold child groups named BathymetryCoverage.NN
    """
    #: attribute name will be automatically determined based on the containing list's index
    __bathymetry_coverage_hdf_name__ = BATHY_COVERAGE + r"[\._]\d+"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def __num_instances_type__(self) -> Type[int]:
        return numpy.uint8  # S102 reduces this from uint32

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
        self.data_coding_format = self.__data_coding_format_type__(9)  # regular grid

    def dimension_create(self):
        """ Creates a blank, empty or zero value for dimension"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.dimension = self.__dimension_type__(2)


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
        self.datatype = self.__datatype_type__("H5T_FLOAT")

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


class BathymetryCoverageDataset(S102FeatureInformationDataset):
    @property
    def metadata_name(self) -> str:
        return BATHY_COVERAGE


# @TODO can I just derive quality from the bathy?
class QualityOfSurvey_GroupNNN(S1xxObject):
    """ This is the "Values" Group.NNN object that contains the grid data in a values dataset and other metadata about the grids.

    4.2.1.1.1 and Figure 4.4 of v2.0.0
    also see section 12.3 and table 12.5

    """

    write_format_str = ".%02d"

    __values_hdf_name__ = "values"  #: HDF5 naming
    # @TODO are these metadata attributes applicable to the QualityCoverage (feature attribute table) or only the BathymetryCoverage objects?
    #  They come from the S100 spec while S102 says no attributes (did they mean no additional attributes?)
    __offset_vectors_hdf_name__ = "offsetVectors"  #: HDF5 naming

    @property
    def values(self) -> numpy.ndarray:
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
    def values(self, val: numpy.ndarray):
        self._attributes[self.__values_hdf_name__] = val

    # FIXME - is this an ndarray or needs a qualityValues class like BathmetryValues?
    # Values for the BathymetryCoverageBase used BathymetryValues(S1xxGridBase) since it needed a compound array (depth, uncertainty)
    # while the Quality can just use a ndarray to hold the integers.
    @property
    def __values_type__(self) -> Type[numpy.ndarray]:
        return numpy.ndarray

    @property
    def quality_dtype(self) -> Type[int]:
        return numpy.int64

    def values_create(self):
        """ Creates a blank, empty or zero value for values"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.values = self.__values_type__([], self.quality_dtype)

    @property
    def __version__(self) -> int:
        return 1

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


class QualityGroupList(S102MetadataListBase):
    """ This is the list of Group.NNN that are held as a list.
    Each Group.NNN has a dataset of depth and uncertainty.
    """
    write_format_str = "_%03d"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "Group"

    @property
    def metadata_type(self) -> type:
        return QualityOfSurvey_GroupNNN


class QualityFeatureInstance(FeatureInstanceDCF9):
    """ This will be the QualityCoverage.001 element in HDF5.
    It will contain a Group.NNN which will have the "values" dataset of the deptha dn uncertainty.
    """
    __quality_group_hdf_name__ = "Group" + r"[\._]\d+"
    """ Basic template for HDF5 naming of the attribute.  
    Attribute name will be automatically determined based on the list's index of the data. 
    """

    @property
    def __num_grp_type__(self) -> Type[int]:
        return numpy.uint8  # S102 only

    def num_grp_create(self):
        self.num_grp = self.__num_grp_type__(1)

    @property
    def __quality_group_type__(self):
        return QualityGroupList

    def quality_group_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.quality_group = self.__quality_group_type__()

    @property
    def quality_group(self) -> S102MetadataListBase:
        """ The quality data, a list of qualitygroup
        Returns
        -------
        S102MetadataListBase
            Contains a list of QualityOfSurvey_GroupNNN objects via the QualityCoveragesList class
        """
        return self._attributes[self.__quality_group_hdf_name__]

    @quality_group.setter
    def quality_group(self, val: S102MetadataListBase):
        self._attributes[self.__quality_group_hdf_name__] = val

    @property
    def __east_bound_longitude_type__(self):
        return numpy.float32


# S102 adds the interpolationType
class QualityOfSurveyContainer(FeatureContainerDCF9, InterpolationType):
    """ This is the QualityOfSurvey right off the root of the HDF5 which has possible attributes from S100 spec table 10c-10
    This will hold child groups named QualityOfSurvey.NN
    """
    #: attribute name will be automatically determined based on the containing list's index
    __quality_of_survey_hdf_name__ = QUALITY_OF_SURVEY + r"[\._]\d+"
    # featureAttributeTable is inherited via FeatureContainerDCF9 from S100 v5.0
    @property
    def __version__(self) -> int:
        return 1

    @property
    def __num_instances_type__(self) -> Type[int]:
        return numpy.uint8  # S102 reduces this from uint32

    @property
    def __quality_of_survey_type__(self):
        return QualityOfSurveysList

    def quality_of_survey_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.quality_of_survey = self.__quality_of_survey_type__()

    @property
    def quality_of_survey(self) -> S102MetadataListBase:
        """ The quality_of_survey data, a list of QualityOfSurvey

        Returns
        -------
        S102MetadataListBase
            Contains a list of QualityOfSurvey objects via the QualityOfSurveyList class
        """
        return self._attributes[self.__quality_of_survey_hdf_name__]

    @quality_of_survey.setter
    def quality_of_survey(self, val: S102MetadataListBase):
        self._attributes[self.__quality_of_survey_hdf_name__] = val

    def data_coding_format_create(self):
        """ Creates a blank, empty or zero value for data_coding_format"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_coding_format = self.__data_coding_format_type__(9)  # regular grid

    def dimension_create(self):
        """ Creates a blank, empty or zero value for dimension"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.dimension = self.__dimension_type__(2)

    # @TODO fixme reduce/remove these as possible and make compound dataset for feature attribute records
    @property
    def __feature_attribute_table_type__(self) -> Type[FeatureAttributeDataset]:
        return FeatureAttributeDataset

    def feature_attribute_table_create(self):
        """ Creates a blank, empty or zero value for feature_attribute_table"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        # FIXME @TODO -- this has to be an arbitrary number of strings, is this limiting it to two?
        self.feature_attribute_table = self.__feature_attribute_table_type__()


class QualityOfSurveysList(S102MetadataListBase):
    """ 4.2.1.1.2 and Figure 4.4 and Table 10.1 of v2.0.0
    This is the set of QualityOfSurvey.NN that act like a list here.
    They will contain a list of Groups.NNN as well as other attributes etc.
    """
    write_format_str = ".%02d"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return QUALITY_OF_SURVEY

    @property
    def metadata_type(self) -> Type[QualityFeatureInstance]:
        return QualityFeatureInstance


# The FeatureAttributeRecord is a new class in 2.2 and stores the quality values (basically the catzoc flags)
# the FeatureInstance class from the dataCodingFormat will have integer values that refer to the ids in these records.
class FeatureAttributeRecord(S1xxObject):
    """ The records that the integer matrix refer to and the actual attribute values for the quality of coverage.
    From section 10.2.7 and table 10.7
    """
    __id_hdf_name__ = "id"  #: HDF5 naming
    __data_assessment_hdf_name__ = "dataAssessment"  #: HDF5 naming
    __least_depth_of_detected_features_measured_hdf_name__ = "featuresDetected.leastDepthOfDetectedFeaturesMeasured"  #: HDF5 naming
    __significant_features_detected_hdf_name__ = "featuresDetected.significantFeaturesDetected"  #: HDF5 naming
    __size_of_features_detected_hdf_name__ = "featuresDetected.sizeOfFeaturesDetected"  #: HDF5 naming
    __feature_size_var_hdf_name__ = "featureSizeVar"  #: HDF5 naming
    __full_seafloor_coverage_achieved_hdf_name__ = "fullSeafloorCoverageAchieved"  #: HDF5 naming
    __bathy_coverage_hdf_name__ = "bathyCoverage"  #: HDF5 naming
    __uncertainty_fixed_hdf_name__ = "zoneOfConfidence.horizontalPositionUncertainty.uncertaintyFixed"  #: HDF5 naming
    __uncertainty_variable_factor_hdf_name__ = "zoneOfConfidence.horizontalPositionUncertainty.uncertaintyVariableFactor"  #: HDF5 naming
    __date_start_hdf_name__ = "surveyDateRange.dateStart"  #: HDF5 naming
    __date_end_hdf_name__ = "surveyDateRange.dateEnd"  #: HDF5 naming
    __source_survey_id_hdf_name__ = "sourceSurveyID"  #: HDF5 naming
    __survey_authority_hdf_name__ = "surveyAuthority"  #: HDF5 naming
    __bathymetric_uncertainty_type_hdf_name__ = "bathymetricUncertaintyType"  #: HDF5 naming

    @property
    def id(self) -> int:
        return self._attributes[self.__id_hdf_name__]

    @id.setter
    def id(self, val:int):
        self._attributes[self.__id_hdf_name__] = val

    @property
    def __id_type__(self) -> Type[int]:
        return numpy.uint32

    def id_create(self):
        """ Creates a blank, empty or zero value for id"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.id = self.__id_type__()

    @property
    def data_assessment(self) -> int:
        return self._attributes[self.__data_assessment_hdf_name__]

    @data_assessment.setter
    def data_assessment(self, val: int):
        self._attributes[self.__data_assessment_hdf_name__] = val

    @property
    def __data_assessment_type__(self) -> Type[int]:
        return numpy.uint8

    def data_assessment_create(self):
        """ Creates a blank, empty or zero value for data_assessment"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_assessment = self.__data_assessment_type__()

    @property
    def least_depth_of_detected_features_measured(self) -> int:
        return self._attributes[self.__least_depth_of_detected_features_measured_hdf_name__]

    @least_depth_of_detected_features_measured.setter
    def least_depth_of_detected_features_measured(self, val: int):
        self._attributes[self.__least_depth_of_detected_features_measured_hdf_name__] = val

    @property
    def __least_depth_of_detected_features_measured_type__(self) -> Type[int]:
        return numpy.uint8

    def least_depth_of_detected_features_measured_create(self):
        """ Creates a blank, empty or zero value for least_depth_of_detected_features_measured"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.least_depth_of_detected_features_measured = self.__least_depth_of_detected_features_measured_type__()

    @property
    def significant_features_detected(self) -> int:
        return self._attributes[self.__significant_features_detected_hdf_name__]

    @significant_features_detected.setter
    def significant_features_detected(self, val: int):
        self._attributes[self.__significant_features_detected_hdf_name__] = val

    @property
    def __significant_features_detected_type__(self) -> Type[int]:
        return numpy.uint8

    def significant_features_detected_create(self):
        """ Creates a blank, empty or zero value for significant_features_detected"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.significant_features_detected = self.__significant_features_detected_type__()

    @property
    def size_of_features_detected(self) -> float:
        return self._attributes[self.__size_of_features_detected_hdf_name__]

    @size_of_features_detected.setter
    def size_of_features_detected(self, val: float):
        self._attributes[self.__size_of_features_detected_hdf_name__] = val

    @property
    def __size_of_features_detected_type__(self) -> Type[float]:
        return numpy.float32

    def size_of_features_detected_create(self):
        """ Creates a blank, empty or zero value for size_of_features_detected"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.size_of_features_detected = self.__size_of_features_detected_type__()

    @property
    def feature_size_var(self) -> float:
        return self._attributes[self.__feature_size_var_hdf_name__]

    @feature_size_var.setter
    def feature_size_var(self, val: float):
        self._attributes[self.__feature_size_var_hdf_name__] = val

    @property
    def __feature_size_var_type__(self) -> Type[float]:
        return numpy.float32

    def feature_size_var_create(self):
        """ Creates a blank, empty or zero value for feature_size_var"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_size_var = self.__feature_size_var_type__()

    @property
    def full_seafloor_coverage_achieved(self) -> int:
        return self._attributes[self.__full_seafloor_coverage_achieved_hdf_name__]

    @full_seafloor_coverage_achieved.setter
    def full_seafloor_coverage_achieved(self, val: int):
        self._attributes[self.__full_seafloor_coverage_achieved_hdf_name__] = val

    @property
    def __full_seafloor_coverage_achieved_type__(self) -> Type[int]:
        return numpy.uint8

    def full_seafloor_coverage_achieved_create(self):
        """ Creates a blank, empty or zero value for full_seafloor_coverage_achieved"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.full_seafloor_coverage_achieved = self.__full_seafloor_coverage_achieved_type__()

    @property
    def bathy_coverage(self) -> int:
        return self._attributes[self.__bathy_coverage_hdf_name__]

    @bathy_coverage.setter
    def bathy_coverage(self, val: int):
        self._attributes[self.__bathy_coverage_hdf_name__] = val

    @property
    def __bathy_coverage_type__(self) -> Type[int]:
        return numpy.uint8

    def bathy_coverage_create(self):
        """ Creates a blank, empty or zero value for bathy_coverage"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.bathy_coverage = self.__bathy_coverage_type__()

    @property
    def uncertainty_fixed(self) -> float:
        return self._attributes[self.__uncertainty_fixed_hdf_name__]

    @uncertainty_fixed.setter
    def uncertainty_fixed(self, val: float):
        self._attributes[self.__uncertainty_fixed_hdf_name__] = val

    @property
    def __uncertainty_fixed_type__(self) -> Type[float]:
        return numpy.float32

    def uncertainty_fixed_create(self):
        """ Creates a blank, empty or zero value for uncertainty_fixed"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.uncertainty_fixed = self.__uncertainty_fixed_type__()

    @property
    def uncertainty_variable_factor(self) -> float:
        return self._attributes[self.__uncertainty_variable_factor_hdf_name__]

    @uncertainty_variable_factor.setter
    def uncertainty_variable_factor(self, val: float):
        self._attributes[self.__uncertainty_variable_factor_hdf_name__] = val

    @property
    def __uncertainty_variable_factor_type__(self) -> Type[float]:
        return numpy.float32

    def uncertainty_variable_factor_create(self):
        """ Creates a blank, empty or zero value for uncertainty_variable_factor"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.uncertainty_variable_factor = self.__uncertainty_variable_factor_type__()

    @property
    def date_start(self) -> str:
        return self._attributes[self.__date_start_hdf_name__]

    @date_start.setter
    def date_start(self, val: str):
        self._attributes[self.__date_start_hdf_name__] = val

    @property
    def __date_start_type__(self) -> Type[str]:
        return str

    def date_start_create(self):
        """ Creates a blank, empty or zero value for date_start"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.date_start = self.__date_start_type__()

    @property
    def date_end(self) -> str:
        return self._attributes[self.__date_end_hdf_name__]

    @date_end.setter
    def date_end(self, val: str):
        self._attributes[self.__date_end_hdf_name__] = val

    @property
    def __date_end_type__(self) -> Type[str]:
        return str

    def date_end_create(self):
        """ Creates a blank, empty or zero value for date_end"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.date_end = self.__date_end_type__()

    @property
    def source_survey_id(self) -> str:
        return self._attributes[self.__source_survey_id_hdf_name__]

    @source_survey_id.setter
    def source_survey_id(self, val: str):
        self._attributes[self.__source_survey_id_hdf_name__] = val

    @property
    def __source_survey_id_type__(self) -> Type[str]:
        return str

    def source_survey_id_create(self):
        """ Creates a blank, empty or zero value for source_survey_id"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.source_survey_id = self.__source_survey_id_type__()

    @property
    def survey_authority(self) -> str:
        return self._attributes[self.__survey_authority_hdf_name__]

    @survey_authority.setter
    def survey_authority(self, val: str):
        self._attributes[self.__survey_authority_hdf_name__] = val

    @property
    def __survey_authority_type__(self) -> Type[str]:
        return str

    def survey_authority_create(self):
        """ Creates a blank, empty or zero value for survey_authority"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.survey_authority = self.__survey_authority_type__()

    @property
    def bathymetric_uncertainty_type(self) -> BATHYMETRIC_UNCERTAINTY_TYPE:
        return self._attributes[self.__bathymetric_uncertainty_type_hdf_name__]

    @bathymetric_uncertainty_type.setter
    def bathymetric_uncertainty_type(self, val: Union[int, str, BATHYMETRIC_UNCERTAINTY_TYPE]):
        self.set_enum_attribute(val, self.__bathymetric_uncertainty_type_hdf_name__, self.__bathymetric_uncertainty_type_type__)

    @property
    def __bathymetric_uncertainty_type_type__(self) -> Type[BATHYMETRIC_UNCERTAINTY_TYPE]:
        return BATHYMETRIC_UNCERTAINTY_TYPE

    def bathymetric_uncertainty_type_create(self):
        """ Creates a blank, empty or zero value for bathymetric_uncertainty_type
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.bathymetric_uncertainty_type = list(self.bathymetric_uncertainty_type_type)[0]

    @property
    def __version__(self) -> int:
        return 1

    def get_write_order(self):
        return [self.__id_hdf_name__,
                self.__data_assessment_hdf_name__,
                self.__least_depth_of_detected_features_measured_hdf_name__,
                self.__significant_features_detected_hdf_name__,
                self.__size_of_features_detected_hdf_name__,
                self.__feature_size_var_hdf_name__,
                self.__full_seafloor_coverage_achieved_hdf_name__,
                self.__bathy_coverage_hdf_name__,
                self.__uncertainty_fixed_hdf_name__,
                self.__uncertainty_variable_factor_hdf_name__,
                self.__date_start_hdf_name__,
                self.__date_end_hdf_name__,
                self.__source_survey_id_hdf_name__,
                self.__survey_authority_hdf_name__,
                self.__bathymetric_uncertainty_type_hdf_name__,
                ]

    def get_write_dtypes(self):
        expected_items = self.get_standard_properties_mapping()
        dtypes = []
        for key in self.get_write_order():
            use_type = self.__getattribute__("__" + expected_items[key] + "_type__")
            dtypes.append((key, use_type))
        return dtypes

# TODO FIXME - somewhere here needs to be the featureAttributeTable from 10.2.7 (or 10.2.8) and table 10.6 (Table 14) that holds the list of FeatureAttributeRecords
class FeatureAttributeDataset(S1xxDatasetBase):
    """ This class comes from S102 -- 10.2.7 Feature information group.
    This class serves to keep a list of FeatureAttributeRecord objects which will be turned into a compound array
    of strings in the HDF5 file.

    The metadata_name property must be overridden.
    """

    @property
    def metadata_type(self) -> Type[FeatureAttributeRecord]:
        return FeatureAttributeRecord

    # analogous to BathymetryCoverageDataset
    @property
    def metadata_name(self) -> str:
        return "featureAttributeTable"


class QualityCoverageDataset(S102FeatureInformationDataset):
    """ This is for the Group_F information group.
    It is the same data structure as the BathymetryCoverageDataset.
    Adds the QualityOfSurvey from S102 v2.2 section 10.2.2 and table 10.3
    """
    @property
    def metadata_name(self) -> str:
        return QUALITY_OF_SURVEY


class FeatureCodes(GroupFBase):
    """ Table 10.1 and sect 10.2.1 of v2.0.0
    Add the QualityOfSurvey from S102 v2.2 section 10.2.2 and table 10.3
    """

    __feature_name_hdf_name__ = "featureName"  #: HDF5 naming
    __bathymetry_coverage_dataset_hdf_name__ = BATHY_COVERAGE
    __quality_of_survey_dataset_hdf_name__ = QUALITY_OF_SURVEY

    def feature_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_code = self.__feature_code_type__([BATHY_COVERAGE, QUALITY_OF_SURVEY], dtype=h5py_string_dtype)

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

    @property
    def __quality_of_survey_dataset_type__(self):
        return QualityCoverageDataset

    def quality_of_survey_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.quality_of_survey_dataset = self.__quality_of_survey_dataset_type__()

    @property
    def quality_of_survey_dataset(self) -> QualityCoverageDataset:
        return self._attributes[self.__quality_of_survey_dataset_hdf_name__]

    @quality_of_survey_dataset.setter
    def quality_of_survey_dataset(self, val: QualityCoverageDataset):
        self._attributes[self.__quality_of_survey_dataset_hdf_name__] = val


class S102Root(S100Root):
    """The root group contains a feature information group and N feature containers.
    In S102 there are currently two feature containers which are the 'coverages'  bathymetry and tracking list.
    The coverage names are determined from the matching CoveragesAttributes
    10.2 and Figure 10.1 of v2.0.0

    v2.2 Adds the QualityOfSurvey from S102 v2.2 section 10.2.2 and table 10.3
    """
    __feature_information_hdf_name__ = "Group_F"  #: HDF5 naming
    __bathymetry_coverage_hdf_name__ = BATHY_COVERAGE
    __quality_of_survey_hdf_name__ = QUALITY_OF_SURVEY  #: HDF5 naming

    # If the S100 file has a type that doesn't match the product spec (like float32 vs float64) then you can just override the type function.
    # @property
    # def __east_bound_longitude_type__(self):
    #     return numpy.float32
    #
    @property
    def __version__(self) -> int:
        return 1

    @property
    def feature_information(self) -> FeatureCodes:
        """Feature Information stored in GroupF in the HDF5 using :class:`FeatureCodes`"""
        return self._attributes[self.__feature_information_hdf_name__]

    @feature_information.setter
    def feature_information(self, val: FeatureCodes):
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
    def quality_of_survey(self) -> S1xxObject:
        return self._attributes[self.__quality_of_survey_hdf_name__]

    @quality_of_survey.setter
    def quality_of_survey(self, val:S1xxObject):
        self._attributes[self.__quality_of_survey_hdf_name__] = val

    @property
    def __quality_of_survey_type__(self) -> Type[S1xxObject]:
        return QualityOfSurveyContainer

    def quality_of_survey_create(self):
        """ Creates a blank, empty or zero value for quality_of_survey"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.quality_of_survey = self.__quality_of_survey_type__()

    def vertical_cs_create(self):
        """ Sets the vertical_cs to Depth as S102 specifies
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_cs = VERTICAL_CS.Depth

    @property
    def __vertical_datum_type__(self) -> Type[int]:
        return numpy.uint16


class S102File(S100File):
    PRODUCT_SPECIFICATION = PRODUCT_SPECIFICATION
    # these keys allow backward compatibility with NAVO data, the first key is current at time of writing
    top_level_keys = ('BathymetryCoverage', 'S102_Grid', 'S102_BathymetryCoverage')
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
    def z_down(self) -> bool:  # reverse Z direction
        return True


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
        origin = numpy.array([bathy_01.grid_origin_longitude, bathy_01.grid_origin_latitude])
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
        valid_epsg = [4326, 5041, 5042]  # , 6318, 6319 are NAD83
        valid_epsg += list(numpy.arange(32601, 32660 + 1))
        valid_epsg += list(numpy.arange(32701, 32760 + 1))
        return valid_epsg

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
        r = self.root
        del r.horizontal_cs, r.projection_parameter_1, r.projection_parameter_2, r.projection_parameter_3
        del r.projection_parameter_4, r.projection_parameter_5, r.false_easting, r.false_northing
        del r.name_of_horizontal_datum, r.type_of_horizontal_crs, r.horizontal_datum, r.prime_meridian, r.spheriod, r.projection_method
        del r.name_of_horizontal_crs, r.epoch
        self._set_bathy_defaults()
        self._set_quality_defaults()  # added for DataCodingFormat9 only -- need to make this optional

    def _set_quality_defaults(self, overwrite=True):
        root = self.root
        quality_dset = root.feature_information.quality_of_survey_dataset
        quality_info = quality_dset.append_new_item()
        quality_info.initialize_properties(True, overwrite=overwrite)
        quality_info.code = "id"
        quality_info.name = ""
        quality_info.unit_of_measure = ""
        quality_info.fill_value = 0
        quality_info.datatype = "H5T_INTEGER"
        quality_info.lower = 1
        quality_info.upper = ""
        quality_info.closure = "geSemiInterval"

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
        bathy_uncertainty_info.lower = 0
        bathy_uncertainty_info.closure = "gtLeInterval"

        root.bathymetry_coverage.axis_names = numpy.array(["Longitude", "Latitude"])  # row major order means X/longitude first
        root.bathymetry_coverage.sequencing_rule_scan_direction = "Longitude, Latitude"
        root.bathymetry_coverage.common_point_rule = 1  # average
        # root.bathymetry_coverage.data_coding_format = 2  # default
        # root.bathymetry_coverage.dimension = 2  # default value
        root.bathymetry_coverage.interpolation_type = 1  # nearest neighbor
        root.bathymetry_coverage.num_instances = 1  # how many Bathycoverages
        root.bathymetry_coverage.sequencing_rule_type = 1  # linear
        del root.bathymetry_coverage.time_uncertainty
        del root.bathymetry_coverage.feature_attribute_table  # this only goes in the QualityOfSurvey group


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
                    flip_z: bool = False, quality_grid=None):  # num_array, or list of lists accepted

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
        quality_grid
            The quality dataset to embed in the object. Optional for S102, leave as None if not needed.

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

        del root.geographic_identifier
        del root.meta_features
        del root.bathymetry_coverage.data_offset_vector
        del root.bathymetry_coverage.data_offset_code
        del bathy_01.grid_spacing_vertical
        del bathy_01.grid_origin_vertical
        del bathy_01.number_of_times
        del bathy_01.time_record_interval
        del bathy_01.date_time_of_last_record
        del bathy_01.date_time_of_first_record
        bathy_01.num_grp = 1

        if quality_grid is not None:  # I'd like to just do a loop on bathy and quality but they have some different names and might have different values at some point in the future.
            try:
                quality_01 = root.quality_of_survey.quality_of_survey[0]
            except IndexError:
                quality_01 = root.quality_of_survey.quality_of_survey.append_new_item()
            quality_01.initialize_properties(recursively_create_children=True, overwrite=overwrite)

            del root.quality_of_survey.data_offset_vector
            del root.quality_of_survey.data_offset_code
            del quality_01.grid_spacing_vertical
            del quality_01.grid_origin_vertical
            del quality_01.number_of_times
            del quality_01.time_record_interval
            del quality_01.date_time_of_last_record
            del quality_01.date_time_of_first_record
            quality_01.num_grp = 1

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

        if quality_grid is not None:
            try:
                quality_group_object = quality_01.quality_group[0]
            except IndexError:
                quality_group_object = quality_01.quality_group.append_new_item()
            if depth_grid.shape != quality_grid.shape:
                raise S102Exception("Depth and Quality grids have different shapes")

        bathy_01.num_points_latitudinal = rows
        bathy_01.num_points_longitudinal = cols
        bathy_01.start_sequence = "0,0"
        del bathy_01.num_points_vertical
        del bathy_01.vertical_extent_maximum_z
        del bathy_01.vertical_extent_minimum_z

        if quality_grid is not None:
            quality_01.num_points_latitudinal = rows
            quality_01.num_points_longitudinal = cols
            quality_01.start_sequence = "0,0"
            del quality_01.num_points_vertical
            del quality_01.vertical_extent_maximum_z
            del quality_01.vertical_extent_minimum_z

        bathy_01.extent_create()
        bathy_01.extent = numpy.array([[0,0], [rows, cols]], bathy_01.extent_dtype)
        # bathy_01.extent.initialize_properties(True, overwrite=overwrite)
        # bathy_01.extent.low.coord_values[0:2] = [0, 0]
        # bathy_01.extent.high.coord_values[0:2] = [rows, cols]

        bathy_group_object.dimension = 2

        # bathy_group_object.origin_create()
        # bathy_group_object.origin.initialize_properties(True, overwrite=overwrite)
        # bathy_group_object.origin.dimension = 2

        bathy_group_object.values_create()

        grid = bathy_group_object.values
        if quality_grid is not None:
            quality_group_object.values_create()
            quality_grid_values = quality_group_object.values
        # @todo -- need to make sure nodata values are correct,
        #   especially if converting something other than bag which is supposed to have the same nodata value
        # @todo -- Add logic that if the grids are gdal raster bands then read in blocks and use h5py slicing to write in blocks.
        #   Slower but saves resources
        if flip_x:
            depth_grid = numpy.fliplr(depth_grid)
            uncert_grid = numpy.fliplr(uncert_grid)
            if quality_grid is not None:
                quality_grid = numpy.fliplr(quality_grid)
        if flip_y:
            depth_grid = numpy.flipud(depth_grid)
            uncert_grid = numpy.flipud(uncert_grid)
            if quality_grid is not None:
                quality_grid = numpy.flipud(quality_grid)

        # @TODO backport this to 2.1, 2.0
        # change the grids to use the no data values specified in the S102 file
        depth_mask = None
        uncert_mask = None
        fill_value = root.feature_information.bathymetry_coverage_dataset[0].fill_value
        uncert_fill_value = root.feature_information.bathymetry_coverage_dataset[1].fill_value

        # Numpy comparison using NaN doesn't work (must use isnan function)
        # so check for NaN first and convert everything to a non-NaN (S102 uses 1000000)
        if numpy.isnan(nodata_value):  # if the nodata value is nan then use the fill value from the dataset
            if not numpy.isnan(fill_value):  # if fill value isn't also NaN then use it
                depth_mask = numpy.isnan(depth_grid)
                uncert_mask = numpy.isnan(uncert_grid)
                if quality_grid is not None:
                    quality_mask = numpy.isnan(quality_grid)
        # the fill value is nan and nodata wasn't OR
        # Both fills are numbers - make sure they match
        elif numpy.isnan(fill_value) or \
            nodata_value != fill_value:
            depth_mask = depth_grid == nodata_value
            uncert_mask = uncert_grid == nodata_value
            if quality_grid is not None:
                quality_mask = quality_grid == nodata_value
        if depth_mask is not None:
            depth_grid = numpy.copy(depth_grid)
            depth_grid[depth_mask] = fill_value
            uncert_grid = numpy.copy(uncert_grid)
            uncert_grid[uncert_mask] = uncert_fill_value
            if quality_grid is not None:
                quality_grid = numpy.copy(quality_grid)
                quality_grid[quality_mask] = root.feature_information.quality_of_survey_dataset[0].fill_value  # will be zero per spec

        # flip the z values if requested now that we have fixed the possible NaN logic issue
        if flip_z:
            depth_grid[depth_grid != fill_value] *= -1
        # Do this after the fill value is set
        try:
            depth_max = depth_grid[depth_grid != fill_value].max()
            depth_min = depth_grid[depth_grid != fill_value].min()
        except ValueError:  # an empty depth array (all values == nodata) will cause this, subdivide() may cause this or data to be updated later
            depth_min = depth_max = nodata_value
        bathy_group_object.maximum_depth = depth_max
        bathy_group_object.minimum_depth = depth_min

        try:
            uncertainty_max = uncert_grid[uncert_grid != uncert_fill_value].max()
            uncertainty_min = uncert_grid[uncert_grid != uncert_fill_value].min()
        except ValueError:  # an empty uncertainty array (all values == nodata) will cause this
            uncertainty_max = uncertainty_min = nodata_value
        bathy_group_object.minimum_uncertainty = uncertainty_min
        bathy_group_object.maximum_uncertainty = uncertainty_max


        grid.depth = depth_grid
        grid.uncertainty = uncert_grid
        if quality_grid is not None:
            # there are no metadata attributes (min/max) for QualityOfSurvey Group_001 - see 10.2.10
            quality_group_object.values_create()
            try:
                if not numpy.issubdtype(quality_grid.dtype, numpy.integer):
                    quality_grid = quality_grid.astype(numpy.uint32)
            except:
                pass
            quality_group_object.values = quality_grid  # @TODO is this right or do we need to do like depth+uncertainty
            # quality_group_object.values.quality_of_survey = quality_grid


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
                                  overwrite: bool = True, flip_z: bool = False, quality_grid=None):  # raw arrays and metadata accepted
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
                - "horizontalDatumReference": Removed in S102 v2.2
                    See :any:`S102Root` horizontal_datum_reference, ex: "EPSG".
                    "EPSG" is the default value.
                - "horizontalDatumValue":  The value for the horizontal data such as the EPSG code ex: 32611
                - "epoch":
                - "geographicIdentifier": Location of the data, ex: "Long Beach, CA, USA".
                    An empty string ("") is the default.
                - "issueDate":  ISO 8601 date string, ex: "2019-01-01"
                - "issueTime": ISO 8601 time string, ex: "00:00:00"
                - "metadataFile": File name for the associated discovery metatadata (xml)
                = "verticalDatumReference": VERTICAL_DATUM_REFERENCE enumeration value.
                    VERTICAL_DATUM_REFERENCE.s100VerticalDatum is the default
                - "verticalDatum": VERTICAL_DATUM enumeration value or EPSG code based on "verticalDatumReference".
                    VERTICAL_DATUM.MLLW is the default
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
                                flip_z=flip_z, quality_grid=quality_grid)

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

        if "verticalDatumReference" in metadata or "verticalDatum" in metadata or overwrite:
            # root.vertical_cs = VERTICAL_CS.Depth
            root.vertical_coordinate_base = VERTICAL_COORDINATE_BASE.verticalDatum
            root.vertical_datum_reference = metadata.get('verticalDatumReference', VERTICAL_DATUM_REFERENCE.s100VerticalDatum)
            root.vertical_datum = metadata.get("verticalDatum", VERTICAL_DATUM.MLLW)

        # these names are taken from the S100/S102 attribute names
        # but are hard coded here to allow the S102 spec to change but not affect any tools built on these utility functions
        # if "horizontalDatumReference" in metadata or overwrite:
        #     root.horizontal_datum_reference = metadata.get("horizontalDatumReference", "EPSG")
        if "horizontalDatumValue" in metadata or overwrite:
            source_epsg = int(metadata.get("horizontalDatumValue", 0))
            if source_epsg in self.get_valid_epsg():
                root.horizontal_crs = source_epsg
            else:
                raise S102Exception(f'The provided EPSG code {source_epsg} is not within the S102 specified values.')
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(root.horizontal_crs))
        if srs.IsProjected():
            axes = ["Easting", "Northing"]  # ["Northing", "Easting"]  # row major instead of
            wgs = osr.SpatialReference()
            wgs.ImportFromEPSG(4326)  # 4326 is WGS84 geodetic - and S102 specifies WGS84
            transform = osr.CoordinateTransformation(srs, wgs)
            # mytransf = Transformer.from_crs(root.horizontal_crs, CRS.from_epsg(4326), always_xy=True)
            south_lat, west_lon = transform.TransformPoint(minx, miny)[:2]
            north_lat, east_lon = transform.TransformPoint(maxx, maxy)[:2]
        else:
            axes = ["Longitude", "Latitude"]  # ["Latitude", "Longitude"]  # row major instead of
            south_lat, west_lon = miny, minx
            north_lat, east_lon = maxy, maxx

        root.east_bound_longitude = east_lon
        root.west_bound_longitude = west_lon
        root.south_bound_latitude = south_lat
        root.north_bound_latitude = north_lat

        feature_instance_groups = [bathy_01]
        if quality_grid is not None:
            quality_01 = root.quality_of_survey.quality_of_survey[0]
            feature_instance_groups.append(quality_01)
            root.quality_of_survey.axis_names = numpy.array(axes)  # row major order means X/longitude first
            root.quality_of_survey.sequencing_rule_scan_direction = ", ".join(axes)

        # the bathymetry_coverage.01 and quality_of_survey.01 are defined to have the same attributes
        for instance in feature_instance_groups:
            # S102 says this is in the CRS of the data (projected) against S100 which says units of Arc Degrees (lat/lon)
            instance.east_bound_longitude = maxx
            instance.west_bound_longitude = minx
            instance.south_bound_latitude = miny
            instance.north_bound_latitude = maxy

            # S102 says this is in the CRS of the data (projected) while S100 says units of Arc Degrees (lat/lon)
            instance.grid_origin_longitude = minx
            instance.grid_origin_latitude = miny
            instance.grid_spacing_longitudinal = abs(res_x)  # we adjust for negative resolution in the from_arrays
            instance.grid_spacing_latitudinal = abs(res_y)

        root.bathymetry_coverage.axis_names = numpy.array(axes)  # row major order means X/longitude first
        root.bathymetry_coverage.sequencing_rule_scan_direction = ", ".join(axes)

        if "epoch" in metadata or overwrite:
            root.epoch = metadata.get("epoch", "")  # e.g. "G1762"  this is the 2013-10-16 WGS84 used by CRS
            if not root.epoch:  # remove optional field if empty
                del root.epoch
        if "geographicIdentifier" in metadata or overwrite:
            root.geographic_identifier = metadata.get("geographicIdentifier", "")
            if not root.geographic_identifier:  # remove optional field if empty
                del root.geographic_identifier
        if "issueDate" in metadata or overwrite:
            root.issue_date = metadata.get('issueDate', datetime.date.today().isoformat())  # datetime.date.today().isoformat()
        if "issueTime" in metadata or overwrite:
            root.issue_time = metadata.get('issueTime', "")  # datetime.date.today().isoformat()
            if not root.issue_time:  # remove optional field if empty
                del root.issue_time
        if "metadataFile" in metadata or overwrite:
            root.metadata = metadata.get('metadataFile', "")  # datetime.date.today().isoformat()

        self.write()
        self.flush()

    @classmethod
    def from_raster(cls, input_raster, output_file, metadata: dict = None, flip_z=False) -> S102File:  # gdal instance or filename accepted
        """  Fills or creates an :any:`S102File` from the given arguments.
        Assumes that depth is in band 1, uncertainty is in band 2, quality is in band 3.

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
        if "horizontalDatumValue" not in metadata:
            sr = osr.SpatialReference(dataset.GetProjection())
            epsg = sr.GetAuthorityCode(None)
            if epsg is None and sr.IsProjected():
                crs_2d = osr.SpatialReference(dataset.GetProjection())
                # in GDAL 3.2 this was added which may make getting the horizontal CRS more obvious
                crs_2d.DemoteTo2D()
                epsg = crs_2d.GetAuthorityCode(None)
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
        if dataset.RasterCount > 2:
            quality_band = dataset.GetRasterBand(3)
            qual_data = quality_band.ReadAsArray()
        else:
            qual_data = None
        # Fill the QualityOfSurvey table
        if qual_data is not None:
            table = self.root.quality_of_survey.feature_attribute_table
            rat = quality_band.GetDefaultRAT()
            column_map = {rat.GetNameOfCol(ncol):ncol for ncol in range(rat.GetColumnCount())}

            for nrow in range(rat.GetRowCount()):
                rec = table.append_new_item()
                rec.id = rat.GetValueAsInt(nrow, column_map['value'])
                rec.data_assessment = rat.GetValueAsInt(nrow, column_map['data_assessment'])
                rec.least_depth_of_detected_features_measured = rat.GetValueAsInt(nrow, column_map['feature_least_depth'])
                rec.significant_features_detected = rat.GetValueAsInt(nrow, column_map['significant_features'])
                rec.size_of_features_detected = rat.GetValueAsDouble(nrow, column_map['feature_size'])
                rec.feature_size_var = 0  # set to "does not scale with depth" instead of rat.GetValueAsDouble(nrow, column_map[''])
                rec.full_seafloor_coverage_achieved = rat.GetValueAsInt(nrow, column_map['coverage'])
                rec.bathy_coverage = rat.GetValueAsInt(nrow, column_map['bathy_coverage'])
                rec.uncertainty_fixed = rat.GetValueAsDouble(nrow, column_map['horizontal_uncert_fixed'])
                rec.uncertainty_variable_factor = rat.GetValueAsDouble(nrow, column_map['horizontal_uncert_var'])
                rec.date_start = rat.GetValueAsString(nrow, column_map['survey_date_start'])
                rec.date_end = rat.GetValueAsString(nrow, column_map['survey_date_end'])
                rec.source_survey_id = rat.GetValueAsString(nrow, column_map['source_survey_id'])
                rec.survey_authority = rat.GetValueAsString(nrow, column_map['source_institution'])
                rec.bathymetric_uncertainty_type = 0  # use "unknown" instead of rat.GetValueAsInt(nrow, column_map[''])

        self.load_arrays_with_metadata(raster_band.ReadAsArray(), uncertainty_band.ReadAsArray(), metadata,
                                                   nodata_value=depth_nodata_value, flip_z=flip_z, quality_grid=qual_data)

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


    @staticmethod
    def upgrade_in_place(s100_object):
        if s100_object.root.product_specification != PRODUCT_SPECIFICATION:
            v2_1.S102File.upgrade_in_place(s100_object)
        if s100_object.root.product_specification == v2_1.PRODUCT_SPECIFICATION:
            # update product specification
            s100_object.attrs['productSpecification'] = S102File.PRODUCT_SPECIFICATION
            # Update horizontal CRS
            del s100_object.attrs['horizontalDatumReference']
            s100_object.attrs.create('horizontalCRS', s100_object.attrs['horizontalDatumValue'], dtype=numpy.int32)
            del s100_object.attrs['horizontalDatumValue']

            for name in ['westBoundLongitude', 'eastBoundLongitude', 'southBoundLatitude', 'northBoundLatitude']:
                change_attr_type(s100_object, name, numpy.float32)
                change_attr_type(s100_object['BathymetryCoverage']['BathymetryCoverage.01'], name, numpy.float32)
            change_attr_type(s100_object['BathymetryCoverage'], 'dimension', numpy.uint8)
            change_attr_type(s100_object['BathymetryCoverage'], 'horizontalPositionUncertainty', numpy.float32)
            change_attr_type(s100_object['BathymetryCoverage'], 'verticalUncertainty', numpy.float32)
            change_attr_type(s100_object['BathymetryCoverage'], 'numInstances', numpy.uint8)
            try:
                # This was a mistake in s100py created data.
                del s100_object['BathymetryCoverage']['BathymetryCoverage.01'].attrs['extentTypeCode']
            except KeyError:
                pass
            change_attr_type(s100_object['BathymetryCoverage']['BathymetryCoverage.01'], 'numGRP', numpy.uint8)
            # Update the vertical CRS
            s100_object.attrs.create('verticalCS', 6498, dtype=numpy.int32)
            s100_object.attrs.create('verticalCoordinateBase', 2, dtype=make_enum_dtype(VERTICAL_COORDINATE_BASE))
            try:
                # This was a mistake in s100py created data.
                del s100_object['BathymetryCoverage']['BathymetryCoverage.01']['Group_001'].attrs['dimension']
            except KeyError:
                pass
            # changes min/max values to single precision.
            for name in ['maximumDepth', 'maximumUncertainty', 'minimumDepth', 'minimumUncertainty']:
                change_attr_type(s100_object['BathymetryCoverage']['BathymetryCoverage.01']['Group_001'], name, numpy.float32)
            try:
                # This was a mistake in s100py created data.
                del s100_object['BathymetryCoverage']['BathymetryCoverage.01']['Group_001']['extent']
            except KeyError:
                pass

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
        instances = list(self.to_raster_datasets())
        if len(instances) > 1:
            raise NotImplementedError('Only one coverage per S102 file is expected')
        dataset, group_instance, flipx, flipy = instances[0]
        # Add the feature attribute table if it exists, an IndexError will be thrown if it doesn't
        try:
            band_values = self.quality
        except IndexError:
            pass  # no quality band, just produce the two band tif of depth and uncertainty
        else:
            dataset.AddBand(gdal.GDT_Float32)  # add a band for the feature attribute table
            if flipx:  # switch the order of the feature attribute band so it matches the depth and uncertainty bands
                band_values = numpy.fliplr(band_values)
            if flipy:
                band_values = numpy.flipud(band_values)
            # add the feature attribute table band to the memory dataset
            dataset.GetRasterBand(3).WriteArray(band_values)
            dataset.GetRasterBand(3).SetDescription("id")  # @TODO  grab from the GroupF
            dataset.GetRasterBand(3).SetNoDataValue(0)  # @TODO grab from the GroupF
            # Add the metadata descriptions equating to the feature attribute table
            rat = gdal.RasterAttributeTable()
            # read the columns of the feature attribute table and translate to GDAL RAT types
            # each column has a defined type in HDF5 so we can use the first record to determine the type
            for key, val in self.feature_attribute_table[0]._attributes.items():
                try:
                    val = val.item()  # convert numpy types to python types
                except AttributeError:
                    pass
                if isinstance(val, (str, Enum)):
                    col_type = gdal.GFT_String
                elif isinstance(val, int):
                    col_type = gdal.GFT_Integer
                elif isinstance(val, float):
                    col_type = gdal.GFT_Real
                else:
                    raise TypeError(f"Unknown data type ({type(val)} submitted for "
                                    "gdal raster attribute table.")
                usage = gdal.GFU_Generic if key != "id" else gdal.GFU_MinMax
                rat.CreateColumn(key, col_type, usage)
            # allocate space for all the existing rows
            rat.SetRowCount(len(self.feature_attribute_table))
            # iterate through the rows and add the values to the RAT
            # check the type to determine which SetValue method to use
            for row_idx, row in enumerate(self.feature_attribute_table):
                for col_idx, (key, val) in enumerate(row._attributes.items()):
                    if isinstance(val, str):
                        rat.SetValueAsString(row_idx, col_idx, val)
                    elif isinstance(val, Enum):
                        rat.SetValueAsString(row_idx, col_idx, val.name)
                    elif isinstance(val, int):
                        rat.SetValueAsInt(row_idx, col_idx, val)
                    elif isinstance(val, float):
                        rat.SetValueAsDouble(row_idx, col_idx, val)
            dataset.GetRasterBand(3).SetDefaultRAT(rat)
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

    @property
    def quality(self):
        return self.root.quality_of_survey.quality_of_survey[0].quality_group[0].values

    @property
    def feature_attribute_table(self):
        return self.root.quality_of_survey.feature_attribute_table

    def get_feature_attribute_dict(self):
        d = {}
        for rec in self.feature_attribute_table:
            d[rec.id] = rec
        return d

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
