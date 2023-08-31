import os
from abc import ABC
import datetime
import logging
from typing import Callable, Iterator, Union, Optional, List, Type
from enum import Enum, IntEnum
import numpy
import h5py

from s100py.s1xx import s1xx_sequence, S1xxObject, S1xxCollection, S1xxDatasetBase, S1xxGridsBase, S1XXFile, h5py_string_dtype
from s100py.s100 import S100Root, S100Exception, FeatureContainerDCF2, FeatureInstanceDCF2, FeatureInformation, FeatureInformationDataset, GroupFBase

WATER_LEVEL = "WaterLevel"

FILLVALUE_HEIGHT = -9999.0
FILLVALUE_TREND = 0

TYPE_OF_WATER_LEVEL_DATA = Enum(value="TYPE_OF_WATER_LEVEL_DATA",
                                names=[("Observation", 1),
                                       ("Astronomical prediction", 2),
                                       ("Analysis or hybrid method", 3),
                                       ("Hydrodynamic model hindcast", 4),
                                       ("Hydrodynamic model forecast", 5),
                                       ("Observed minus predicted", 6),
                                       ("Observed minus analysis", 7),
                                       ("Observed minus hindcast", 8),
                                       ("Observed minus forecast", 9),
                                       ("Forecast minus predicted", 10),
                                   ]
                                )


VERTICAL_COORDINATE_BASE = Enum(value="VERTICAL_COORDINATE_BASE",
                                names=[("Sea Surface", 1),
                                       ("Vertical Datum", 2),
                                       ("Sea Bottom", 3)
                                       ]
                                )

VERTICAL_DATUM_REFERENCE = Enum(value="VERTICAL_DATUM_REFERENCE",
                                names=[("S-100 vertical datum", 1),
                                       ("EPSG", 2)
                                       ]
                                )


class S104Exception(S100Exception):
    """Raised when input is not S104 compliant"""
    pass


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
    def time_point(self) -> S1xxObject:
        """Defines the conversion from python naming to HDF5 (S100) naming"""
        return self._attributes[self.__time_point_hdf_name__]

    @time_point.setter
    def time_point(self, val: S1xxObject):
        self._attributes[self.__time_point_hdf_name__] = val

    @property
    def __time_point_type__(self) -> Type[str]:
        """Attribute datatype"""
        return str

    def time_point_create(self):
        """ Creates a blank, empty or zero value for time_point"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.time_point = self.__time_point_type__()

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
    """ Basic template for the name of the attribute
    Attribute name will be automatically determined based on the array position
    of the S104_MetadataList
    """
    __water_level_group_hdf_name__ = "Group" + r"[\._]\d+"
    __uncertainty_dataset_hdf_name__ = "uncertainty"
    __number_of_nodes_hdf_name__ = "numberOfNodes"
    __type_of_water_level_data_hdf_name__ = "typeOfWaterLevelData"

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
    def number_of_nodes(self) -> S1xxObject:
        return self._attributes[self.__number_of_nodes_hdf_name__]

    @number_of_nodes.setter
    def number_of_nodes(self, val: S1xxObject):
        self._attributes[self.__number_of_nodes_hdf_name__] = val

    @property
    def __number_of_nodes_type__(self) -> Type[numpy.int32]:
        return numpy.int32

    def number_of_nodes_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.number_of_nodes = self.__number_of_nodes_type__()

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
    def type_of_water_level_data(self) -> TYPE_OF_WATER_LEVEL_DATA:
        return self._attributes[self.__type_of_water_level_data_hdf_name__]

    @type_of_water_level_data.setter
    def type_of_water_level_data(self, val: Union[int, str, TYPE_OF_WATER_LEVEL_DATA]):
        self.set_enum_attribute(val, self.__type_of_water_level_data_hdf_name__, self.__type_of_water_level_data_type__)

    @property
    def __type_of_water_level_data_type__(self) -> Type[TYPE_OF_WATER_LEVEL_DATA]:
        return TYPE_OF_WATER_LEVEL_DATA

    def type_of_water_level_data_create(self):
        """ Creates a value using the first item in the enumeration of type_of_water_level_data"""
        # Make the enum into a list and take the first value
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.type_of_water_level_data = list(self.__type_of_water_level_data_type__)[0]


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
    __vertical_coordinate_system_hdf_name__ = "verticalCS"
    __vertical_coordinate_base_hdf_name__ = "verticalCoordinateBase"
    __vertical_datum_reference_hdf_name__ = "verticalDatumReference"
    __vertical_datum_epsg_hdf_name__ = "verticalDatum"
    __horizontal_crs_hdf_name__ = "horizontalCRS"

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
    def vertical_coordinate_system(self) -> S1xxObject:
        return self._attributes[self.__vertical_coordinate_system_hdf_name__]

    @vertical_coordinate_system.setter
    def vertical_coordinate_system(self, val: S1xxObject):
        self._attributes[self.__vertical_coordinate_system_hdf_name__] = val

    @property
    def __vertical_coordinate_system_type__(self) -> Type[numpy.int32]:
        """Define S104 datatype"""
        return numpy.int32

    def vertical_coordinate_system_create(self):
        """ Creates a blank, empty or zero value for vertical coordinate system"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_coordinate_system = self.__vertical_coordinate_system_type__()

    @property
    def vertical_coordinate_base(self) -> VERTICAL_COORDINATE_BASE:
        return self._attributes[self.__vertical_coordinate_base_hdf_name__]

    @vertical_coordinate_base.setter
    def vertical_coordinate_base(self, val: Union[int, str, VERTICAL_COORDINATE_BASE]):
        self.set_enum_attribute(val, self.__vertical_coordinate_base_hdf_name__, self.__vertical_coordinate_base_type__)

    @property
    def __vertical_coordinate_base_type__(self) -> Type[VERTICAL_COORDINATE_BASE]:
        """Enumeration data type"""
        return VERTICAL_COORDINATE_BASE

    def vertical_coordinate_base_create(self):
        """ Creates a value using the first item in the enumeration of vertical_coordinate_base"""
        # Make the enum into a list and take the first value
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
        """Defines enumeration datatype"""
        return VERTICAL_DATUM_REFERENCE

    def vertical_datum_reference_create(self):
        """ Creates a value using the first item in the enumeration of vertical_datum_reference"""
        # Make the enum into a list and take the first value
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum_reference = list(self.__vertical_datum_reference_type__)[0]

    @property
    def vertical_datum_epsg(self) -> S1xxObject:
        """EPSG code for vertical datum for verticalDatumReference = 2"""
        return self._attributes[self.__vertical_datum_hdf_name__]

    @vertical_datum_epsg.setter
    def vertical_datum_epsg(self, val: S1xxObject):
        self._attributes[self.__vertical_datum_hdf_name__] = val

    @property
    def __vertical_datum_epsg_type__(self) -> Type[numpy.int32]:
        """Define datatype"""
        return numpy.int32

    def vertical_datum_epsg_create(self):
        """ Creates a blank, empty or zero value for vertical_datum_epsg"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum_epsg = self.__vertical_datum_epsg_type__()

    @property
    def horizontal_crs(self) -> S1xxObject:
        return self._attributes[self.__horizontal_crs_hdf_name__]

    @horizontal_crs.setter
    def horizontal_crs(self, val: S1xxObject):
        self._attributes[self.__horizontal_crs_hdf_name__] = val

    @property
    def __horizontal_crs_type__(self) -> Type[numpy.int32]:
        """Define S104 datatype"""
        return numpy.int32

    def horizontal_crs_create(self):
        """ Creates a blank, empty or zero value for horizontal crs"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_crs = self.__horizontal_crs_type__()


class DiscoveryMetadata(S1xxObject):
    """ 12.2.6 of v1.0.1"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        raise NotImplementedError()


class S104File(S1XXFile):
    """ HDF5 file object"""
    PRODUCT_SPECIFICATION = 'INT.IHO.S-104.0.0'

    def __init__(self, *args, **kywrds):
        super().__init__(*args, root=S104Root, **kywrds)
