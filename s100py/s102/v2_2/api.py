try:
    from ... import s1xx
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py.s102"

from ..v2_0 import api as v2_0
from ..v2_1 import api as v2_1
# Anything not overridden in this module will use whatever was available in the previous version
from ..v2_1.api import *

EDITION = 2.2

QUALITY_COVERAGE = "QualityOfSurveyCoverage"


# @TODO can I just derive quality from the bathy?


class QualityOfSurveyCoverageBase(S1xxObject):
    """ This is the Group.NNN object that contains the grid data in a values dataset and other metadata about the grids.

    4.2.1.1.1 and Figure 4.4 of v2.0.0
    also see section 12.3 and table 12.5

    """

    write_format_str = ".%03d"

    __values_hdf_name__ = "values"  #: HDF5 naming
    # @TODO are these metadata attributes applicable to the QualityCoverage (feature attribute table) or only the BathymetryCoverage objects?
    #  They come from the S100 spec while S102 says no attributes (did they mean no additional attributes?)
    __origin_hdf_name__ = "origin"  #: HDF5 naming
    __offset_vectors_hdf_name__ = "offsetVectors"  #: HDF5 naming
    __dimension_hdf_name__ = "dimension"  #: HDF5 naming
    __axis_names_hdf_name__ = "axisNames"  #: HDF5 naming
    __extent_hdf_name__ = "extent"  #: HDF5 naming
    __sequencing_rule_hdf_name__ = "sequencingRule"  #: HDF5 naming
    __start_sequence_hdf_name__ = "startSequence"  #: HDF5 naming

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

    # Values for the BathymetryCoverageBase used BathymetryValues(S1xxGridBase) since it needed a compound array (depth, uncertainty)
    # while the Quality can just use a ndarray to hold the integers.
    @property
    def __values_type__(self) -> Type[numpy.ndarray]:
        return numpy.ndarray

    def values_create(self):
        """ Creates a blank, empty or zero value for values"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.values = self.__values_type__()

    @property
    def __version__(self) -> int:
        return 1

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

class QualityOfSurveyCoverage(QualityOfSurveyCoverageBase):  # , DisplayScaleMixin
    pass


class QualityGroupList(S102MetadataListBase):
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
        return QualityOfSurveyCoverage

class QualityFeatureInstance(FeatureInstanceDCF9):
    """ This will be the QualityCoverage.001 element in HDF5.
    It will contain a Group.NNN which will have the "values" dataset of the deptha dn uncertainty.
    """
    __quality_group_hdf_name__ = "Group" + r"[\._]\d+"
    """ Basic template for HDF5 naming of the attribute.  
    Attribute name will be automatically determined based on the list's index of the data. 
    """

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
            Contains a list of QualityOfSurveyCoverage objects via the QualityCoveragesList class
        """
        return self._attributes[self.__quality_group_hdf_name__]

    @quality_group.setter
    def quality_group(self, val: S102MetadataListBase):
        self._attributes[self.__quality_group_hdf_name__] = val


class QualityContainer(FeatureContainerDCF2):
    """ This is the QualityOfSurveyCoverage right off the root of the HDF5 which has possible attributes from S100 spec table 10c-10
    This will hold child groups named QualityCoverage.NN
    """
    #: attribute name will be automatically determined based on the containing list's index
    __quality_coverage_hdf_name__ = QUALITY_COVERAGE + r"[\._]\d+"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def __quality_coverage_type__(self):
        return QualityCoveragesList

    def quality_coverage_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.quality_coverage = self.__quality_coverage_type__()

    @property
    def quality_coverage(self) -> S102MetadataListBase:
        """ The quality data, a list of QualityCoverage

        Returns
        -------
        S102MetadataListBase
            Contains a list of QualityCoverage objects via the QualityCoveragesList class
        """
        return self._attributes[self.__quality_coverage_hdf_name__]

    @quality_coverage.setter
    def quality_coverage(self, val: S102MetadataListBase):
        self._attributes[self.__quality_coverage_hdf_name__] = val

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


class QualityCoveragesList(S102MetadataListBase):
    """ 4.2.1.1.2 and Figure 4.4 and Table 10.1 of v2.0.0
    This is the set of qualityCoverage.NN that act like a list here.
    They will contain a list of Groups.NNN as well as other attributes etc.
    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return BATHY_COVERAGE

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

    @property
    def id(self) -> int:
        return self._attributes[self.__id_hdf_name__]

    @id.setter
    def id(self, val:int):
        self._attributes[self.__id_hdf_name__] = val

    @property
    def __id_type__(self) -> Type[int]:
        return int

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
        return int

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
        return int

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
        return int

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
        return float

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
        return float

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
        return int

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
        return int

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
        return float

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
        return float

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


class FeatureAttributeDataset(S1xxDatasetBase):  # Chunking
    """ This class comes from S102 -- 10.2.7 Feature information group.
    This class serves to keep a list of FeatureAttributeRecord objects which will be turned into a compound array
    of strings in the HDF5 file.

    The metadata_name property must be overridden.
    """

    @property
    def metadata_type(self) -> Type[FeatureAttributeRecord]:
        return FeatureAttributeRecord

# @TODO just put the metadata_name in the FeatureAttributeDataset?
class QualityOfSurveyCoverageDataset(FeatureAttributeDataset):
    # analogous to BathymetryCoverageDataset
    @property
    def metadata_name(self) -> str:
        return QUALITY_COVERAGE


class QualityCoverageDataset(S102FeatureInformationDataset):
    """ This is for the Group_F information group.
    It is the same data structure as the BathymetryCoverageDataset.
    Adds the QualityOfSurveyCoverage from S102 v2.2 section 10.2.2 and table 10.3
    """
    @property
    def metadata_name(self) -> str:
        return QUALITY_COVERAGE

class FeatureCodesQualityMixin:
    """ Table 10.1 and sect 10.2.1 of v2.0.0
    Adds the QualityOfSurveyCoverage from S102 v2.2 section 10.2.2 and table 10.3
    """

    __quality_coverage_dataset_hdf_name__ = QUALITY_COVERAGE

    def feature_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_code = self.__feature_code_type__([BATHY_COVERAGE, QUALITY_COVERAGE], dtype=h5py_string_dtype)

    @property
    def __version__(self) -> int:
        return 1

    @property
    def __quality_coverage_dataset_type__(self):
        return QualityCoverageDataset

    def quality_coverage_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.quality_coverage_dataset = self.__quality_coverage_dataset_type__()

    @property
    def quality_coverage_dataset(self) -> QualityCoverageDataset:
        return self._attributes[self.__quality_coverage_dataset_hdf_name__]

    @quality_coverage_dataset.setter
    def quality_coverage_dataset(self, val: QualityCoverageDataset):
        self._attributes[self.__quality_coverage_dataset_hdf_name__] = val


# FIXME @TODO Is this mixed in with Bathymetry or a separate attribute?  Should have Quality parallel to Bathymetry
class FeatureCodes(FeatureCodesQualityMixin, v2_1.FeatureCodes):
    pass


class S102RootQualityMixin:
    """ Adds the QualityOfSurveyCoverage from S102 v2.2 section 10.2.2 and table 10.3
    """
    __quality_of_survey_coverage_hdf_name__ = QUALITY_COVERAGE  #: HDF5 naming

    @property
    def __feature_information_type__(self):
        return FeatureCodes

    @property
    def quality_of_survey_coverage(self) -> FeatureCodesBase:
        return self._attributes[self.__quality_of_survey_coverage_hdf_name__]

    @quality_of_survey_coverage.setter
    def quality_of_survey_coverage(self, val:FeatureCodesBase):
        self._attributes[self.__quality_of_survey_coverage_hdf_name__] = val

    @property
    def __quality_of_survey_coverage_type__(self) -> Type[FeatureCodesBase]:
        return FeatureCodes

    def quality_of_survey_coverage_create(self):
        """ Creates a blank, empty or zero value for quality_of_survey_coverage"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.quality_of_survey_coverage = self.__quality_of_survey_coverage_type__()

# S100Root was updated in v5.0 so re-inherit here
class S102Root(S102RootQualityMixin, v2_1.S102RootMixin, S100Root):
    """ BathymetryCoverage, QualityOfSurveyCoverage and Group_F should be dataset children of this.
    ProductSpecification along with coordinate reference system info and extents should be an attributes of this.
    """
    pass


class S102File(v2_1.S102File):
    """ Derives from an HDF5 file object and defines the root object of an HDF5 file.
    """
    PRODUCT_SPECIFICATION = 'INT.IHO.S-102.2.2'

    def __init__(self, name, *args, **kywrds):
        if 'root' not in kywrds:
            kywrds['root'] = S102Root  # inherited classes will specify their own root type
        super().__init__(name, *args, **kywrds)

    @staticmethod
    def upgrade_in_place(s100_object):
        if s100_object.root.product_specification != v2_1.S102File.PRODUCT_SPECIFICATION:
            v2_1.S102File.upgrade_in_place(s100_object)
        # FIXME this function is not complete
        raise NotImplementedError('Upgrade not implemented for S102 v2.2')
        if s100_object.root.product_specification == v2_1.S102File.PRODUCT_SPECIFICATION:
            # update product specification
            s100_object.attrs['productSpecification'] = S102File.PRODUCT_SPECIFICATION
            # remove TrackingList
            del s100_object['TrackingListCoverage']
            del s100_object['Group_F']['TrackingListCoverage']
            del s100_object['Group_F']['featureName']
            fc20 = s100_object['Group_F']['featureCode']
            if fc20[0] in (b'BathymetryCoverage', 'BathymetryCoverage'):
                fc21 = fc20[:1]  # keep the bathymetery and delete the trackinglist
            elif fc20[1] in (b'BathymetryCoverage', 'BathymetryCoverage'):
                fc21 = fc20[1:]  # keep the bathymetery and delete the trackinglist
            del s100_object['Group_F']['featureCode']
            s100_object['Group_F'].create_dataset('featureCode', data=fc21)

            # remove display scale and reverse the Z direction
            for top in v2_0.S102File.top_level_keys:
                try:
                    bathy_top = s100_object[top]
                    groupf_bathy = s100_object['Group_F'][top]
                    # get the fill value to use when reversing the Z value
                    fill_val = float(groupf_bathy[0]['fillValue'])
                    # update the datatype definition
                    # groupf_bathy[0]['datatype'] = 'H5T_FLOAT' fails to adjust the file as the groupf_bathy[0] creates a temporary copy
                    # groupf_bathy[0, 'datatype'] = 'H5T_FLOAT' raises a typeError about changing the datatype
                    # copying the data with temp=groupf_bathy[0] then changing values then setting groupf_bathy[0]=temp seems to work
                    # similar to revising the depth values later in this function
                    for nrow in range(len(groupf_bathy)):
                        row = groupf_bathy[nrow]
                        row['datatype'] = 'H5T_FLOAT'
                        if row['name'].lower() in ("uncertainty", b"uncertainty"):
                            row['lower'] = 0
                            row['closure'] = 'gtLeInterval'
                        groupf_bathy[nrow] = row
                    # depth_string = groupf_bathy[0]['code']
                except KeyError:
                    pass
                else:
                    for second in v2_0.S102File.second_level_keys:
                        try:
                            bathy_cov = bathy_top[second]
                        except KeyError:
                            pass
                        else:
                            for group in v2_0.S102File.group_level_keys:
                                try:
                                    bathy_group = bathy_cov[group]
                                except KeyError:
                                    pass
                                else:
                                    try:
                                        del bathy_group[v2_0.DisplayScaleMixin.__maximum_display_scale_hdf_name__]
                                    except KeyError:
                                        pass
                                    try:
                                        del bathy_group[v2_0.DisplayScaleMixin.__minimum_display_scale_hdf_name__]
                                    except KeyError:
                                        pass
                                    try:
                                        depth_uncert = bathy_group['values']
                                        # h5py does not allow editing via fancy slicing se we need to convert to numpy, edit and then put it back
                                        # i.e. depth[depth!=fill_val] *= -1 won't work but doesn't raise an error either
                                        a = numpy.array(depth_uncert['depth'])
                                        a[a != fill_val] *= -1
                                        depth_uncert['depth'] = a
                                    except KeyError:
                                        pass
                                    # standardize with the 2.1 required Group_001
                                    if group != "Group_001":
                                        bathy_cov.move(group, "Group_001")
                            # standardize with the 2.1 required BathymetryCoverage.01
                            if second != "BathymetryCoverage.01":
                                bathy_top.move(second, "BathymetryCoverage.01")


    def _set_bathy_defaults(self, overwrite=True):
        """ This function initializes the values in more recent versions of the spec to reduce redundant code in later modules
        """
        # FIXME this function is not complete
        raise NotImplementedError('Upgrade not implemented for S102 v2.2')
        super()._set_bathy_defaults(overwrite=overwrite)
        root = self.root
        bathy_cov_dset = root.feature_information.bathymetry_coverage_dataset

        bathy_uncertainty_info = bathy_cov_dset[1]
        bathy_uncertainty_info.lower = 0
        bathy_uncertainty_info.closure = "gtLeInterval"


    def set_defaults(self, overwrite=True):  # remove tracking list
        self.create_empty_metadata()  # init the root with a fully filled out empty metadata set
        self._set_bathy_defaults()

    @classmethod
    def from_gdal(cls, input_raster, output_file, metadata: dict = None, flip_z=False) -> S102File:  # gdal instance or filename accepted
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
