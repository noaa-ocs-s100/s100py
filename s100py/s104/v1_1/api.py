import datetime
from typing import Callable, Iterator, Union, Optional, List, Type
from enum import Enum, IntEnum
import numpy
import h5py

from s100py.s1xx import s1xx_sequence, S1xxObject, S1xxCollection, S1xxDatasetBase, S1xxGridsBase, S1XXFile, h5py_string_dtype
from ...v5_0.s100 import S100File, S100Root, S100Exception, FeatureContainerDCF2, FeatureInstanceDCF2, \
    FeatureInformation, FeatureInformationDataset, GroupFBase, VERTICAL_CS, VERTICAL_DATUM_REFERENCE, VERTICAL_DATUM

WATER_LEVEL = "WaterLevel"

FILLVALUE_HEIGHT = -9999.00
FILLVALUE_TREND = 0

"""Contains s104 metadata to pass to S104File.

PRODUCT_SPECIFICATION: The product specification used to create this dataset.
HORIZONTAL_CRS: EPSG code or -1 if user defined
DATA_CODING_FORMAT: Reference to the type of S104 product.
INTERPOLATION_TYPE: Interpolation method recommended for evaluation of the S100_GridCoverage.
COMMON_POINT_RULE: The procedure used for evaluating geometric objects that overlap or lie fall on boundaries.
DIMENSION: The dimension of the feature instance.
SEQUENCING_RULE_TYPE: Method to assign values from the sequence of values to the grid coordinates (e.g. "linear").
SEQUENCING_RULE_SCAN_DIRECTION: AxisNames, comma-separated (e.g. "longitude,latitude").
START_SEQUENCE: Starting location of the scan.

"""


EDITION = 1.1
PRODUCT_SPECIFICATION = 'INT.IHO.S-104.1.1'

CHANGELOG = """
v1.1 
* S100 Edition 5.0 Only
* Enum list changed from TYPE_OF_WATER_LEVEL_DATA to DATA_DYNAMICITY --  Table 12.10 - S104_DataDynamicity
* General metadata "datasetDeliveryInterval" added
* General metadata "trendInterval" added
* Beginning S-100 5.0.0, seaSurface and seaFloor have been added to the S100_VerticalAndSoundingDatum enumeration, 
which makes attribute verticalCoordinateBase redundant. Since it is optional in S-100, S-104 no longer uses it. 
If verticalCoordinateBase is encoded, applications may ignore it.
* The S-100 attribute verticalCoordinateBase is no longer used as of S-111 Edition 1.2 because
its “sea surface” and “sea bottom” values have been added to the vertical datums enumeration (Table 12-8)
* Datatypes sizes have been specified for all attributes
"""


class S104Exception(S100Exception):
    """Raised when input is not S104 compliant"""
    pass


class DATA_DYNAMICITY(Enum):
    """
     S104 v1.1 Table 12.10 - Classification of data according to the
     relationship between the time of its collection, generation, or
     calculation of generation parameters, in relation to the time of
     publication of the dataset.
    """
    observation = 1
    astronomicalPrediction = 2
    analysisOrHybrid = 3
    hydrodynamicHindcast = 4
    hydrodynamicForecast = 5
    observedMinusPredicted = 6
    observedMinusAnalysis = 7
    observedMinusHindcast = 8
    observedMinusForecast = 9


class WaterLevelTrend(IntEnum):
    """Water level trend enumerated constant and returns an int object"""

    Decreasing = 1
    Increasing = 2
    Steady = 3


# noinspection PyAbstractClass
class S104MetadataListBase(S1xxCollection):
    """Define group name format"""
    write_format_str = ".%02d"


class WaterLevelUncertaintyInformation(S1xxObject):
    """S100 code and uncertainty of data values"""
    __name_hdf_name__ = "name"  #: HDF5 naming
    __value_hdf_name__ = "value"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def name(self) -> str:
        """ The plain text name of the data

        Returns:
            str: Name of the data ("waterLevelHeight" or "waterLevelTrend")
        """
        return self._attributes[self.__name_hdf_name__]

    @name.setter
    def name(self, val: str):
        """Incoming value datatype validation"""
        self._attributes[self.__name_hdf_name__] = val

    @property
    def __name_type__(self):
        """Uncertainty name datatype"""
        return str

    def name_create(self):
        """Create empty object"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.name = self.__name_type__()

    @property
    def value(self) -> str:
        """ The uncertainty value"""
        return self._attributes[self.__value_hdf_name__]

    @value.setter
    def value(self, val: int):
        """Incoming value datatype validation"""
        self._attributes[self.__value_hdf_name__] = val

    @property
    def __value_type__(self):
        """Uncertainty value datatype"""
        return str

    def value_create(self):
        """Create empty object"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.value = self.__value_type__()


class WaterLevelUncertaintyDataset(S1xxDatasetBase):
    """Create uncertainty dataset"""

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        """ The plain text name of the dataset"""
        return "uncertainty"

    @property
    def metadata_type(self) -> Type[WaterLevelUncertaintyInformation]:
        """S104 datatype"""
        return WaterLevelUncertaintyInformation


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
        return numpy.float32

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
        return numpy.float32

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

    __geometry_values_hdf_name__ = "geometry_values"

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


class WaterLevelValues(S1xxGridsBase):
    """NNN Group Datasets"""
    __water_level_height_hdf_name__ = "waterLevelHeight"
    __water_level_trend_hdf_name__ = "waterLevelTrend"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        """ The plain text name of the dataset
        Returns:
            str: Name of the dataset
        """
        return "values"

    @property
    def water_level_height(self) -> s1xx_sequence:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.__water_level_height_hdf_name__]

    @water_level_height.setter
    def water_level_height(self, val: s1xx_sequence):
        self._attributes[self.__water_level_height_hdf_name__] = val

    @property
    def __water_level_height_type__(self) -> s1xx_sequence:
        """Define array datatype"""
        return numpy.ndarray

    @property
    def water_level_height_dtype(self) -> Type[float]:
        """Define array datatype"""
        return numpy.float32

    def water_level_height_create(self):
        """ Creates a blank, empty or zero value for water_level_height"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level_height = self.__water_level_height_type__([], self.water_level_height_dtype)

    @property
    def water_level_trend(self) -> WaterLevelTrend:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.__water_level_trend_hdf_name__]

    @water_level_trend.setter
    def water_level_trend(self, val: Union[int, str, WaterLevelTrend]):
        self.set_enum_attribute(val, self.__water_level_trend_hdf_name__, self.__water_level_trend_type__)

    @property
    def __water_level_trend_type__(self) -> s1xx_sequence:
        """Define datatype"""
        return numpy.ndarray

    @property
    def water_level_trend_dtype(self) -> Type[WaterLevelTrend]:
        """Define array datatype"""
        return h5py.enum_dtype(dict([(water_level_trend.name, water_level_trend.value) for water_level_trend in WaterLevelTrend]))

    def water_level_trend_create(self):
        """ Creates a blank, empty or zero value for water_level_trend"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level_trend = self.__water_level_trend_type__([], self.water_level_trend_dtype)

    def get_write_order(self):
        return [self.__water_level_height_hdf_name__, self.__water_level_trend_hdf_name__]

    def get_compound_dtype(self):
        return [self.water_level_height_dtype, self.water_level_trend_dtype]


class WaterLevelGroup(S1xxObject):

    __values_hdf_name__ = "values"
    __time_point_hdf_name__ = "timePoint"

    @property
    def values(self) -> WaterLevelValues:
        """Plain text name of the dataset (e.g values)"""
        return self._attributes[self.__values_hdf_name__]

    @values.setter
    def values(self, val: WaterLevelValues):
        self._attributes[self.__values_hdf_name__] = val

    @property
    def __values_type__(self) -> Type[WaterLevelValues]:
        return WaterLevelValues

    def values_create(self):
        """ Creates a blank, empty or zero value for values"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.values = self.__values_type__()

    @property
    def time_point(self) -> datetime.datetime:
        return self._attributes[self.__time_point_hdf_name__]

    @time_point.setter
    def time_point(self, val: Union[datetime.date, datetime.datetime, str]):
        self.set_datetime_attribute(val, self.__time_point_hdf_name__, self.__time_point_type__)

    @property
    def __time_point_type__(self) -> Type[datetime.datetime]:
        return datetime.datetime

    @property
    def __time_point_repr__(self) -> str:
        return self._attributes[self.__time_point_hdf_name__].strftime("%Y%m%dT%H%M%SZ")

    def time_point_create(self):
        """ Creates a blank, empty or zero value for time_point"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.time_point = self.__time_point_type__(1970, 1, 1)

    @property
    def __version__(self) -> int:
        return 1


class WaterLevelGroupList(S1xxCollection):
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
        return WaterLevelGroup


class WaterLevelFeatureInstance(FeatureInstanceDCF2):
    """ Basic template for the name of the attribute.
    Attribute name will be automatically determined based on the array position
    of the S104_MetadataList
    """
    __water_level_group_hdf_name__ = "Group" + r"[\._]\d+"
    __uncertainty_dataset_hdf_name__ = "uncertainty"
    __data_dynamicity_hdf_name__ = "dataDynamicity"

    @property
    def __water_level_group_type__(self):
        return WaterLevelGroupList

    def water_level_group_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level_group = self.__water_level_group_type__()

    @property
    def water_level_group(self) -> S1xxCollection:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.__water_level_group_hdf_name__]

    @water_level_group.setter
    def water_level_group(self, val: S1xxCollection):
        self._attributes[self.__water_level_group_hdf_name__] = val

    @property
    def uncertainty_dataset(self) -> S1xxDatasetBase:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.__uncertainty_dataset_hdf_name__]

    @uncertainty_dataset.setter
    def uncertainty_dataset(self, val: S1xxDatasetBase):
        self._attributes[self.__uncertainty_dataset_hdf_name__] = val

    @property
    def __uncertainty_dataset_type__(self) -> Type[WaterLevelUncertaintyDataset]:
        return WaterLevelUncertaintyDataset

    def uncertainty_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.uncertainty_dataset = self.__uncertainty_dataset_type__()

    @property
    def __positioning_group_hdf_name__(self) -> str:
        return "Positioning"

    @property
    def positioning_group(self) -> S1xxObject:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.__positioning_group_hdf_name__]

    @positioning_group.setter
    def positioning_group(self, val: S1xxObject):
        self._attributes[self.__positioning_group_hdf_name__] = val

    @property
    def __positioning_group_type__(self):
        """Defines datatype"""
        return PositioningGroup

    def positioning_group_create(self):
        """ Creates a blank, empty or zero value for positioning_group"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.positioning_group = self.__positioning_group_type__()

    @property
    def data_dynamicity(self) -> DATA_DYNAMICITY:
        return self._attributes[self.__data_dynamicity_hdf_name__]

    @data_dynamicity.setter
    def data_dynamicity(self, val: Union[int, str, DATA_DYNAMICITY]):
        self.set_enum_attribute(val, self.__data_dynamicity_hdf_name__, self.__data_dynamicity_type__)

    @property
    def __data_dynamicity_type__(self) -> Type[DATA_DYNAMICITY]:
        return DATA_DYNAMICITY

    def data_dynamicity_create(self):
        """ Creates a value using the first item in the enumeration of data_dynamicity"""
        # Make the enum into a list and take the first value
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_dynamicity = list(self.__data_dynamicity_type__)[0]

    @property
    def __number_of_times_type__(self) -> Type[int]:
        return numpy.uint32

    @property
    def __time_record_interval_type__(self) -> Type[int]:
        return numpy.uint16


class WaterLevelList(S104MetadataListBase):
    """
    This is the set of WaterLevel.NN that act like a list here.
    They will contain a list of Groups.NNN as well as other attributes etc.
    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return WATER_LEVEL

    @property
    def metadata_type(self) -> Type[WaterLevelFeatureInstance]:
        return WaterLevelFeatureInstance


class WaterLevelContainer(FeatureContainerDCF2):
    """ This is the WaterLevel right off the root of the HDF5 which has possible attributes from S100 spec table 10c-10
    This will hold child groups named WaterLevel.NN
    """

    #: Basic template for the name of the attribute
    #: Attribute name will be automatically determined based on the containing list's index
    __water_level_hdf_name__ = WATER_LEVEL + r"[\._]\d+"
    __min_dataset_height_hdf_name__ = "minDatasetHeight"
    __max_dataset_height_hdf_name__ = "maxDatasetHeight"
    __method_water_level_product_hdf_name__ = "methodWaterLevelProduct"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def __water_level_type__(self):
        return WaterLevelList

    def water_level_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level = self.__water_level_type__()

    @property
    def water_level(self) -> S104MetadataListBase:
        """ The water level data, a list of WaterLevel

        Returns:
            S104_MetadataList_base: Contains a list of WaterLevel objects
            via the WaterLevel_List class
        """
        return self._attributes[self.__water_level_hdf_name__]

    @water_level.setter
    def water_level(self, val: S104MetadataListBase):
        self._attributes[self.__water_level_hdf_name__] = val

    @property
    def min_dataset_height(self) -> S1xxObject:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.__min_dataset_height_hdf_name__]

    @min_dataset_height.setter
    def min_dataset_height(self, val: S1xxObject):
        self._attributes[self.__min_dataset_height_hdf_name__] = val

    @property
    def __min_dataset_height_type__(self) -> Type[float]:
        """Defines datatype"""
        return float

    def min_dataset_height_create(self):
        """ Creates a blank, empty or zero value for min_dataset_height"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.min_dataset_height = self.__min_dataset_height_type__()

    @property
    def max_dataset_height(self) -> S1xxObject:
        return self._attributes[self.__max_dataset_height_hdf_name__]

    @max_dataset_height.setter
    def max_dataset_height(self, val: S1xxObject):
        self._attributes[self.__max_dataset_height_hdf_name__] = val

    @property
    def __max_dataset_height_type__(self) -> Type[numpy.float32]:
        """Defines datatype"""
        return numpy.float32

    def max_dataset_height_create(self):
        """ Creates a blank, empty or zero value for max_dataset_height"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.max_dataset_height = self.__max_dataset_height_type__()

    @property
    def method_water_level_product(self) -> S1xxObject:
        return self._attributes[self.__method_water_level_product_hdf_name__]

    @method_water_level_product.setter
    def method_water_level_product(self, val: S1xxObject):
        self._attributes[self.__method_water_level_product_hdf_name__] = val

    @property
    def __method_water_level_product_type__(self) -> Type[str]:
        """Defines datatype"""
        return str

    def method_water_level_product_create(self):
        """ Creates a blank, empty or zero value for method_water_level_product"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.method_water_level_product = self.__method_water_level_product_type__()


class WaterLevelFeatureDataset(FeatureInformationDataset):
    """Create group_f feature dataset"""

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        """S104 feature information dataset name"""
        return WATER_LEVEL

    @property
    def metadata_type(self) -> Type[FeatureInformation]:
        """Feature information base class"""
        return FeatureInformation


class GroupF(GroupFBase):
    """From S100 Table 10c-8 – Components of feature information group"""

    __water_level_feature_dataset_hdf_name__ = WATER_LEVEL

    @property
    def __version__(self) -> int:
        return 1

    def feature_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_code = self.__feature_code_type__([WATER_LEVEL], dtype=h5py_string_dtype)

    @property
    def __water_level_feature_dataset_type__(self):
        return WaterLevelFeatureDataset

    def water_level_feature_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level_feature_dataset = self.__water_level_feature_dataset_type__()

    @property
    def water_level_feature_dataset(self) -> WaterLevelFeatureDataset:
        return self._attributes[self.__water_level_feature_dataset_hdf_name__]

    @water_level_feature_dataset.setter
    def water_level_feature_dataset(self, val: WaterLevelFeatureDataset):
        self._attributes[self.__water_level_feature_dataset_hdf_name__] = val


class S104Root(S100Root):
    """The root group contains a feature information group and N feature containers.
    In S104 there is one feature container 'water level'.
    The coverage names are determined from the matching CoveragesAttributes
    Table 3 of v0.0.7
    """
    __water_level_hdf_name__ = WATER_LEVEL
    __water_level_trend_threshold_hdf_name__ = "waterLevelTrendThreshold"
    __dataset_delivery_interval_hdf_name__ = "datasetDeliveryInterval"
    __trend_interval_hdf_name__ = "trendInterval"
    __issue_time_hdf_name__ = "issueTime"
    __issue_date_hdf_name__ = "issueDate"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def __feature_information_type__(self):
        return GroupF

    @property
    def water_level(self) -> S1xxObject:
        return self._attributes[self.__water_level_hdf_name__]

    @property
    def __water_level_type__(self):
        return WaterLevelContainer

    def water_level_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level = self.__water_level_type__()

    @water_level.setter
    def water_level(self, val: S1xxObject):
        self._attributes[self.__water_level_hdf_name__] = val

    @property
    def water_level_trend_threshold(self) -> S1xxObject:
        return self._attributes[self.__water_level_trend_threshold_hdf_name__]

    @water_level_trend_threshold.setter
    def water_level_trend_threshold(self, val: S1xxObject):
        self._attributes[self.__water_level_trend_threshold_hdf_name__] = val

    @property
    def __water_level_trend_threshold_type__(self) -> Type[numpy.float32]:
        return numpy.float32

    def water_level_trend_threshold_create(self):
        """ Creates a blank, empty or zero value for water level trend threshold"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level_trend_threshold = self.__water_level_trend_threshold_type__()

    @property
    def dataset_delivery_interval(self) -> S1xxObject:
        return self._attributes[self.__dataset_delivery_interval_hdf_name__]

    @dataset_delivery_interval.setter
    def dataset_delivery_interval(self, val: S1xxObject):
        self._attributes[self.__dataset_delivery_interval_hdf_name__] = val

    @property
    def __dataset_delivery_interval_type__(self):
        return str

    def dataset_delivery_interval_create(self):
        """ Creates a blank, empty or zero value for dataset_delivery_interval"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.dataset_delivery_interval = self.__dataset_delivery_interval_type__()

    @property
    def trend_interval(self) -> S1xxObject:
        return self._attributes[self.__trend_interval_hdf_name__]

    @trend_interval.setter
    def trend_interval(self, val: S1xxObject):
        self._attributes[self.__trend_interval_hdf_name__] = val

    @property
    def __trend_interval_type__(self) -> Type[numpy.int32]:
        return numpy.int32

    def trend_interval_create(self):
        """ Creates a blank, empty or zero value for trend interval"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.trend_interval = self.__trend_interval_type__()


    @property
    def __issue_date_repr__(self) -> str:
        return self._attributes[self.__issue_date_hdf_name__].strftime("%Y%m%d")

    @property
    def __issue_time_repr__(self) -> str:
        return self._attributes[self.__issue_time_hdf_name__].strftime('%H%M%SZ')


class DiscoveryMetadata(S1xxObject):
    """ 12.2.6 of v1.0.1"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        raise NotImplementedError()


class S104File(S100File):
    """ HDF5 file object"""
    PRODUCT_SPECIFICATION = 'INT.IHO.S-104.1.1'

    def __init__(self, *args, **kywrds):
        super().__init__(*args, root=S104Root, **kywrds)
