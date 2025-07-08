import datetime
from typing import Callable, Iterator, Union, Optional, List, Type
from enum import Enum, IntEnum
import numpy
import h5py

from ...s1xx import s1xx_sequence, S1xxObject, S1xxCollection, S1xxDatasetBase, S1xxGridsBase, S1XXFile, h5py_string_dtype
from ...s100.v5_2.api import S100File, S100Root, S100Exception, FeatureContainerDCF2, FeatureInstanceDCF2, \
    FeatureInformation, FeatureInformationDataset, GroupFBase

WATER_LEVEL = "WaterLevel"

FILLVALUE_HEIGHT = -9999.0
FILLVALUE_TREND = 0
FILLVALUE_UNCERTAINTY = -1.0

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


EDITION = 2.0
PRODUCT_SPECIFICATION = 'INT.IHO.S-104.2.0'

CHANGELOG = """
v2.0 
* Aligned with S-100 Edition 5.2.0
* Scope reduced to use only the regular grid spatial type
* 104_DataDynamicity classifications hydrodynamicHindcast, observedMinusPredicted, observedMinusAnalysis,
observedMinusHindcast, observedMinusForecast, forecastMinusPredicted no longer allowed (S-104 Ed 2.0 Table 12-10)
* verticalCoordinateBase attribute added (S-104 Ed 2.0 Table 12-1, S-100 Ed 5.2 Table 10c-22)  
the only allowed value is verticalDatum (2), the attribute is now mandatory
* Instance metadata constraints  (Table 12-3) adjusted for Water Level adjustment compatibility
* Restricted maximum length of HDF5 string attributes (clause 12.3)
* Optional uncertainty attribute added to values record, represents the uncertainty at a particular
grid point and may be omitted if the uncertainty is the same at all grid points
* Extended format to include grids with datum jumps (multiple vertical datums)
* Added UTM zones and newer WGS84 realizations
* Removed ISO metadata files
* Group F waterLevelTime has been replaced with uncertainty (S-104 Ed 2.0 Table 10.3)
* Domain extent polygon added to level 3 content (Feature Instance) If and only if the feature
covers an area with a different vertical datum from the root group (Table 10.2)
* Added optional dataOffsetCode enumeration for Offset of data point in cell data (S-104 Ed 2.0 Table 12-2)
Mandatory if data points are at grid cell centres. See S-100 clauses 10c-9.6 and 8-5.2.8.
The allowed values in S-104 are:
1: XMin, YMin (“Lower left”)
5: Barycenter (centroid) of cell 
The default is “Lower left” and this attribute may be omitted if data points are at cell 
lower-left corners. Other cell corners are not allowed.
* Modifications to S100 allowable vertical and sounding datums in (S-104 Ed. 2.0 Table 12-8)
The values seaFloor (47), seaSurface (48), and hydrographicZero (49) are not used in S-104.
NOTE: S-104 Edition 2.0.0 uses only data types which the water level adjustment algorithm 
described in S-98 can process.
"""


class S104Exception(S100Exception):
    """Raised when input is not S104 compliant"""
    pass


class S104UnspecifiedClassException(S100Exception):
    pass


class DATA_DYNAMICITY(Enum):
    """
     S104 v2.0 Table 12.10 - Classification of data according to the
     relationship between the time of its collection, generation, or
     calculation of generation parameters, in relation to the time of
     publication of the dataset.
    """
    observation = 1
    astronomicalPrediction = 2
    analysisOrHybrid = 3
    hydrodynamicForecast = 5


class WaterLevelTrend(IntEnum):
    """
     Water level trend, the tendency of water level to change in a
     particular direction, enumerated constant and returns an int object
    """

    Decreasing = 1
    Increasing = 2
    Steady = 3


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


class WaterLevelDomainExtentPolygonInformation(S1xxObject):
    """S100 code and uncertainty of data values"""
    __longitude_hdf_name__ = "longitude"  #: HDF5 naming
    __latitude_hdf_name__ = "latitude"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def longitude(self) -> str:
        """ The plain text name of the data

        Returns:
            str: `longitude` name
        """
        return self._attributes[self.__longitude_hdf_name__]

    @longitude.setter
    def longitude(self, val: float):
        """Incoming value datatype validation"""
        self._attributes[self.__longitude_hdf_name__] = val

    @property
    def __longitude_type__(self):
        """Longitude datatype"""
        return float

    def longitude_create(self):
        """Create empty object"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.longitude = self.__longitude_type__()

    @property
    def latitude(self) -> str:
        """ The plain text name of the data

        Returns:
            str: `latitude` name
        """
        return self._attributes[self.__latitude_hdf_name__]

    @latitude.setter
    def latitude(self, val: float):
        """Incoming value datatype validation"""
        self._attributes[self.__latitude_hdf_name__] = val

    @property
    def __latitude_type__(self):
        """latitude value datatype"""
        return float

    def latitude_create(self):
        """Create empty object"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.latitude = self.__latitude_type__()


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


class WaterLevelDomainExtentPolygonDataset(S1xxDatasetBase):
    """Create domainExtent.polygon dataset S-100 Ed 5.2
    Table 10c-11 (Optional) – Containing coordinates of
    bounding polygon vertices of the spatial extent
    of the domain of the coverage
    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        """ The plain text name of the dataset"""
        return "domainExtent.polygon"

    @property
    def metadata_type(self) -> Type[WaterLevelDomainExtentPolygonInformation]:
        """S104 datatype"""
        return WaterLevelDomainExtentPolygonInformation


class WaterLevelValues(S1xxGridsBase):
    """Water Level Feature attributes"""
    __water_level_height_hdf_name__ = "waterLevelHeight"
    __water_level_trend_hdf_name__ = "waterLevelTrend"
    __water_level_uncertainty_hdf_name__ = "uncertainty"

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

    @property
    def water_level_uncertainty(self) -> s1xx_sequence:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.__water_level_uncertainty_hdf_name__]

    @water_level_uncertainty.setter
    def water_level_uncertainty(self, val: s1xx_sequence):
        self._attributes[self.__water_level_uncertainty_hdf_name__] = val

    @property
    def __water_level_uncertainty_type__(self) -> s1xx_sequence:
        """Define array datatype"""
        return numpy.ndarray

    @property
    def water_level_uncertainty_dtype(self) -> Type[float]:
        """Define array datatype"""
        return numpy.float32

    def water_level_uncertainty_create(self):
        """ Creates a blank, empty or zero value for water_level_uncertainty"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level_uncertainty = self.__water_level_uncertainty_type__([], self.water_level_uncertainty_dtype)

    def get_write_order(self):
        """Write feature attributes, optionally write uncertainty attribute"""

        feature_attributes = [self.__water_level_height_hdf_name__, self.__water_level_trend_hdf_name__]
        if self.__water_level_uncertainty_hdf_name__ in self._attributes:
            feature_attributes.extend([self.__water_level_uncertainty_hdf_name__])
        return feature_attributes

    def get_compound_dtype(self):
        """Write compound dtype, optionally write uncertainty compound dtype"""
        feature_attributes_dtype = [self.water_level_height_dtype, self.water_level_trend_dtype]
        if self.__water_level_uncertainty_hdf_name__ in self._attributes:
            feature_attributes_dtype.extend([self.water_level_uncertainty_dtype])
        return feature_attributes_dtype


class WaterLevelGroup(S1xxObject):

    __values_hdf_name__ = "values"
    __time_point_hdf_name__ = "timePoint"
    __water_level_trend_threshold_hdf_name__ = "waterLevelTrendThreshold"
    __trend_interval_hdf_name__ = "trendInterval"

    @property
    def values(self) -> WaterLevelValues:
        """Plain text name of the dataset (e.g. values)"""
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
    def trend_interval(self) -> S1xxObject:
        return self._attributes[self.__trend_interval_hdf_name__]

    @trend_interval.setter
    def trend_interval(self, val: S1xxObject):
        self._attributes[self.__trend_interval_hdf_name__] = val

    @property
    def __trend_interval_type__(self) -> Type[int]:
        return numpy.uint32

    def trend_interval_create(self):
        """ Creates a blank, empty or zero value for trend interval"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.trend_interval = self.__trend_interval_type__()

    @property
    def __version__(self) -> int:
        return 1


class WaterLevelGroupList(S1xxCollection):
    """ This is the list of Group.NNN that are held as a list.
    Each Group.NNN has a dataset of water level height, water level trend and optionally uncertainty.
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


class WaterLevelFeatureInstanceBase:
    """ Basic template for the name of the attribute.
    Attribute name will be automatically determined based on the array position
    of the S104_MetadataList
    """
    __water_level_group_hdf_name__ = "Group" + r"[\._]\d+"
    __uncertainty_dataset_hdf_name__ = "uncertainty"
    __data_dynamicity_hdf_name__ = "dataDynamicity"
    __domain_extent_polygon_dataset_hdf_name__ = "domainExtent.polygon"

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
    def domain_extent_polygon_dataset(self) -> S1xxDatasetBase:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.__domain_extent_polygon_dataset_hdf_name__]

    @domain_extent_polygon_dataset.setter
    def domain_extent_polygon_dataset(self, val: S1xxDatasetBase):
        self._attributes[self.__domain_extent_polygon_dataset_hdf_name__] = val

    @property
    def __domain_extent_polygon_dataset_type__(self) -> Type[WaterLevelDomainExtentPolygonDataset]:
        return WaterLevelDomainExtentPolygonDataset

    def domain_extent_polygon_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.domain_extent_polygon_dataset = self.__domain_extent_polygon_dataset_type__()

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


class WaterLevelFeatureInstanceDCF2(FeatureInstanceDCF2, WaterLevelFeatureInstanceBase):
    @property
    def __number_of_times_type__(self) -> Type[int]:
        return numpy.uint32

    @property
    def __time_record_interval_type__(self) -> Type[int]:
        return numpy.uint16


class WaterLevelListBase(S104MetadataListBase):
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


class WaterLevelListDCF2(WaterLevelListBase):
    @property
    def metadata_type(self) -> Type[WaterLevelFeatureInstanceBase]:
        return WaterLevelFeatureInstanceDCF2


class WaterLevelContainerBase:
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
    def __min_dataset_height_type__(self) -> Type[numpy.float32]:
        """Defines datatype"""
        return numpy.float32

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


class WaterLevelContainerDCF2(FeatureContainerDCF2, WaterLevelContainerBase):
    @property
    def __water_level_type__(self):
        return WaterLevelListDCF2


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
    __vertical_datum_epoch_hdf_name__ = "verticalDatumEpoch"
    __vertical_datum_hdf_name__ = "verticalDatum"

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
        return WaterLevelContainerBase

    def water_level_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        print("You should not create a generic surface_current object but manually create"
              "a WaterLevelContainerDCF2")
        raise S104UnspecifiedClassException("You should create WaterLevelContainerDCF2")
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
    def __trend_interval_type__(self) -> Type[int]:
        return numpy.uint32

    def trend_interval_create(self):
        """ Creates a blank, empty or zero value for trend interval"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.trend_interval = self.__trend_interval_type__()

    @property
    def vertical_datum_epoch(self) -> str:
        return self._attributes[self.__vertical_datum_epoch_hdf_name__]

    @vertical_datum_epoch.setter
    def vertical_datum_epoch(self, val: str):
        self._attributes[self.__vertical_datum_epoch_hdf_name__] = val

    @property
    def __vertical_datum_epoch_type__(self) -> Type[str]:
        return str

    def vertical_datum_epoch_create(self):
        """ Creates a blank, empty or zero value for vertical_datum_epoch"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum_epoch = self.__vertical_datum_epoch_type__()

    @property
    def vertical_datum(self) -> Enum:
        val = self._attributes[self.__vertical_datum_hdf_name__]
        return val

    @vertical_datum.setter
    def vertical_datum(self, val: (int, str, VERTICAL_DATUM)):
        self.set_enum_attribute(val, self.__vertical_datum_hdf_name__, VERTICAL_DATUM)

    @property
    def __vertical_datum_type__(self) -> Type[int]:
        return int

    def vertical_datum_create(self):
        """ Creates a blank, empty or zero value for vertical_datum"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum = self.__vertical_datum_type__()

    @property
    def __issue_date_repr__(self) -> str:
        return self._attributes[self.__issue_date_hdf_name__].strftime("%Y%m%d")

    @property
    def __issue_time_repr__(self) -> str:
        return self._attributes[self.__issue_time_hdf_name__].strftime('%H%M%SZ')


class S104File(S100File):
    """ HDF5 file object"""
    PRODUCT_SPECIFICATION = PRODUCT_SPECIFICATION

    @staticmethod
    def make_container_for_dcf(data_coding_format):
        if data_coding_format == 2:
            container = WaterLevelContainerDCF2()
        else:
            raise S104Exception(f"DCF {data_coding_format} not supported in v2.0")
        return container

    @staticmethod
    def set_feature_information_defaults(water_level_feature_dataset):
        """
        Water Level Height is the height of a water surface relative to a
        vertical datum and Water Level Trend is the tendency of water level
        to change in a particular direction. Default values are set for any
        data that don't have options or are mandatory.

        Parameters
        ----------
        water_level_feature_dataset
            water level feature information object

        """

        water_level_height_info = water_level_feature_dataset.append_new_item()
        water_level_height_info.code = "waterLevelHeight"
        water_level_height_info.name = "Water Level Height"
        water_level_height_info.unit_of_measure = "metre"
        water_level_height_info.datatype = "H5T_FLOAT"
        water_level_height_info.fill_value = f"{FILLVALUE_HEIGHT:0.02f}"
        water_level_height_info.lower = "-99.99"
        water_level_height_info.upper = "99.99"
        water_level_height_info.closure = "closedInterval"

        water_level_trend_info = water_level_feature_dataset.append_new_item()
        water_level_trend_info.code = "waterLevelTrend"
        water_level_trend_info.name = "Water Level Trend"
        water_level_trend_info.unit_of_measure = ""
        water_level_trend_info.datatype = "H5T_ENUM"
        water_level_trend_info.fill_value = FILLVALUE_TREND
        water_level_trend_info.lower = ""
        water_level_trend_info.upper = ""
        water_level_trend_info.closure = ""

    @staticmethod
    def set_water_level_uncertainty_defaults(water_level_feature_dataset):
        """
        Estimate characterising the range of values within which the true
        value of a measurement is expected to lie as defined within a particular
        confidence level. It is expressed as a positive value. Default values
        are set for any data that don't have options (Optional Attribute)

        Parameters
        ----------
        water_level_feature_dataset
            water level feature information object

        """

        water_level_uncertainty_info = water_level_feature_dataset.append_new_item()
        water_level_uncertainty_info.code = "uncertainty"
        water_level_uncertainty_info.name = "Uncertainty"
        water_level_uncertainty_info.unit_of_measure = "metre"
        water_level_uncertainty_info.datatype = "H5T_FLOAT"
        water_level_uncertainty_info.fill_value = f"{FILLVALUE_UNCERTAINTY:0.02f}"
        water_level_uncertainty_info.lower = "0.00"
        water_level_uncertainty_info.upper = "99.99"
        water_level_uncertainty_info.closure = "closedInterval"

    def __init__(self, *args, **kywrds):
        super().__init__(*args, root=S104Root, **kywrds)
        # when reading from a file we need to look inside the DataCodingFormat to know what type of object to create
        try:
            container_key = f'/{WATER_LEVEL}'
            dcf = self[container_key].attrs['dataCodingFormat']
        except KeyError:  # must not be reading an existing file or doesn't have data for some reason
            pass
        else:
            self.root.water_level = self.make_container_for_dcf(dcf)
            self.root.water_level.read(self[container_key])

