import os
from abc import ABC
import datetime
import logging
from typing import Callable, Iterator, Union, Optional, List, Type
from enum import Enum, IntEnum
import numpy
import h5py

from s100py.s1xx import s1xx_sequence, S1xxAttributesBase, S1xxMetadataListBase, S1xxDatasetBase, S1xxGridsBase, S1XXFile, h5py_string_dtype
from s100py.s100 import S100Root, S100Exception, FeatureContainerDCF2, FeatureInstanceDCF2, FeatureInformation, FeatureInformationDataset, GroupFBase

WATER_LEVEL = "WaterLevel"

FILLVALUE_HEIGHT = -9999
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
class S104MetadataListBase(S1xxMetadataListBase):
    """Define group name format"""
    write_format_str = ".%02d"


class WaterLevelUncertaintyInformation(S1xxAttributesBase):
    """S100 code and uncertainty of data values"""
    name_attribute_name = "name"  #: HDF5 naming
    value_attribute_name = "value"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def name(self) -> str:
        """ The plain text name of the data

        Returns:
            str: Name of the data ("waterLevelHeight" or "waterLevelTrend")
        """
        return self._attributes[self.name_attribute_name]

    @name.setter
    def name(self, val: str):
        """Incoming value datatype validation"""
        self._attributes[self.name_attribute_name] = val

    @property
    def name_type(self):
        """Uncertainty name datatype"""
        return str

    def name_create(self):
        """Create empty object"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.name = self.name_type()

    @property
    def value(self) -> str:
        """ The uncertainty value"""
        return self._attributes[self.value_attribute_name]

    @value.setter
    def value(self, val: int):
        """Incoming value datatype validation"""
        self._attributes[self.value_attribute_name] = val

    @property
    def value_type(self):
        """Uncertainty value datatype"""
        return str

    def value_create(self):
        """Create empty object"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.value = self.value_type()


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
    longitude_attribute_name = "longitude"
    latitude_attribute_name = "latitude"

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
        return self._attributes[self.longitude_attribute_name]

    @longitude.setter
    def longitude(self, val: s1xx_sequence):
        """Potential validation or other checks/changes to incoming data"""
        self._attributes[self.longitude_attribute_name] = val

    @property
    def longitude_type(self) -> s1xx_sequence:
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
        self.longitude = self.longitude_type([], self.longitude_dtype)

    @property
    def latitude(self) -> s1xx_sequence:
        """Get the data"""
        return self._attributes[self.latitude_attribute_name]

    @latitude.setter
    def latitude(self, val: s1xx_sequence):
        """Potential validation or other checks/changes to incoming data"""
        self._attributes[self.latitude_attribute_name] = val

    @property
    def latitude_type(self) -> s1xx_sequence:
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
        self.latitude = self.latitude_type([], self.latitude_dtype)

    def get_write_order(self):
        """Specify order of attributes for ordered dict"""
        return [self.longitude_attribute_name, self.latitude_attribute_name]

    def get_compound_dtype(self):
        return [self.longitude_dtype, self.latitude_dtype]


class PositioningGroup(S1xxAttributesBase):

    geometry_values_attribute_name = "geometry_values"

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
        return self._attributes[self.geometry_values_attribute_name]

    @geometry_values.setter
    def geometry_values(self, val: GeometryValuesDataset):
        self._attributes[self.geometry_values_attribute_name] = val

    @property
    def geometry_values_type(self) -> Type[GeometryValuesDataset]:
        """S100 Datatype"""
        return GeometryValuesDataset

    def geometry_values_create(self):
        """ Creates a blank, empty or zero value for geometry_values"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.geometry_values = self.geometry_values_type()


class WaterLevelValues(S1xxGridsBase):
    """NNN Group Datasets"""
    water_level_height_attribute_name = "waterLevelHeight"
    water_level_trend_attribute_name = "waterLevelTrend"

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
        return self._attributes[self.water_level_height_attribute_name]

    @water_level_height.setter
    def water_level_height(self, val: s1xx_sequence):
        self._attributes[self.water_level_height_attribute_name] = val

    @property
    def water_level_height_type(self) -> s1xx_sequence:
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
        self.water_level_height = self.water_level_height_type([], self.water_level_height_dtype)

    @property
    def water_level_trend(self) -> WaterLevelTrend:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.water_level_trend_attribute_name]

    @water_level_trend.setter
    def water_level_trend(self, val: Union[int, str, WaterLevelTrend]):
        self.set_enum_attribute(val, self.water_level_trend_attribute_name, self.water_level_trend_type)

    @property
    def water_level_trend_type(self) -> s1xx_sequence:
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
        self.water_level_trend = self.water_level_trend_type([], self.water_level_trend_dtype)

    def get_write_order(self):
        return [self.water_level_height_attribute_name, self.water_level_trend_attribute_name]

    def get_compound_dtype(self):
        return [self.water_level_height_dtype, self.water_level_trend_dtype]


class WaterLevelGroup(S1xxAttributesBase):

    values_attribute_name = "values"
    time_point_attribute_name = "timePoint"

    @property
    def values(self) -> WaterLevelValues:
        """Plain text name of the dataset (e.g values)"""
        return self._attributes[self.values_attribute_name]

    @values.setter
    def values(self, val: WaterLevelValues):
        self._attributes[self.values_attribute_name] = val

    @property
    def values_type(self) -> Type[WaterLevelValues]:
        return WaterLevelValues

    def values_create(self):
        """ Creates a blank, empty or zero value for values"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.values = self.values_type()

    @property
    def time_point(self) -> S1xxAttributesBase:
        """Defines the conversion from python naming to HDF5 (S100) naming"""
        return self._attributes[self.time_point_attribute_name]

    @time_point.setter
    def time_point(self, val: S1xxAttributesBase):
        self._attributes[self.time_point_attribute_name] = val

    @property
    def time_point_type(self) -> Type[str]:
        """Attribute datatype"""
        return str

    def time_point_create(self):
        """ Creates a blank, empty or zero value for time_point"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.time_point = self.time_point_type()

    @property
    def __version__(self) -> int:
        return 1


class WaterLevelGroupList(S1xxMetadataListBase):
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
    water_level_group_attribute_name = "Group" + r"[\._]\d+"
    uncertainty_dataset_attribute_name = "uncertainty"
    number_of_nodes_attribute_name = "numberOfNodes"

    @property
    def water_level_group_type(self):
        return WaterLevelGroupList

    def water_level_group_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level_group = self.water_level_group_type()

    @property
    def water_level_group(self) -> S1xxMetadataListBase:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.water_level_group_attribute_name]

    @water_level_group.setter
    def water_level_group(self, val: S1xxMetadataListBase):
        self._attributes[self.water_level_group_attribute_name] = val

    @property
    def number_of_nodes(self) -> S1xxAttributesBase:
        return self._attributes[self.number_of_nodes_attribute_name]

    @number_of_nodes.setter
    def number_of_nodes(self, val: S1xxAttributesBase):
        self._attributes[self.number_of_nodes_attribute_name] = val

    @property
    def number_of_nodes_type(self) -> Type[numpy.int32]:
        return numpy.int32

    def number_of_nodes_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.number_of_nodes = self.number_of_nodes_type()

    @property
    def uncertainty_dataset(self) -> S1xxDatasetBase:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.uncertainty_dataset_attribute_name]

    @uncertainty_dataset.setter
    def uncertainty_dataset(self, val: S1xxDatasetBase):
        self._attributes[self.uncertainty_dataset_attribute_name] = val

    @property
    def uncertainty_dataset_type(self) -> Type[WaterLevelUncertaintyDataset]:
        return WaterLevelUncertaintyDataset

    def uncertainty_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.uncertainty_dataset = self.uncertainty_dataset_type()

    @property
    def positioning_group_attribute_name(self) -> str:
        return "Positioning"

    @property
    def positioning_group(self) -> S1xxAttributesBase:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.positioning_group_attribute_name]

    @positioning_group.setter
    def positioning_group(self, val: S1xxAttributesBase):
        self._attributes[self.positioning_group_attribute_name] = val

    @property
    def positioning_group_type(self):
        """Defines datatype"""
        return PositioningGroup

    def positioning_group_create(self):
        """ Creates a blank, empty or zero value for positioning_group"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.positioning_group = self.positioning_group_type()


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
    water_level_attribute_name = WATER_LEVEL + r"[\._]\d+"
    min_dataset_height_attribute_name = "minDatasetHeight"
    max_dataset_height_attribute_name = "maxDatasetHeight"
    method_water_level_product_attribute_name = "methodWaterLevelProduct"
    type_of_water_level_data_attribute_name = "typeOfWaterLevelData"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def water_level_type(self):
        return WaterLevelList

    def water_level_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level = self.water_level_type()

    @property
    def water_level(self) -> S104MetadataListBase:
        """ The water level data, a list of WaterLevel

        Returns:
            S104_MetadataList_base: Contains a list of WaterLevel objects
            via the WaterLevel_List class
        """
        return self._attributes[self.water_level_attribute_name]

    @water_level.setter
    def water_level(self, val: S104MetadataListBase):
        self._attributes[self.water_level_attribute_name] = val

    @property
    def min_dataset_height(self) -> S1xxAttributesBase:
        """Defines the conversion from python naming to HDF5 (S104) naming"""
        return self._attributes[self.min_dataset_height_attribute_name]

    @min_dataset_height.setter
    def min_dataset_height(self, val: S1xxAttributesBase):
        self._attributes[self.min_dataset_height_attribute_name] = val

    @property
    def min_dataset_height_type(self) -> Type[float]:
        """Defines datatype"""
        return float

    def min_dataset_height_create(self):
        """ Creates a blank, empty or zero value for min_dataset_height"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.min_dataset_height = self.min_dataset_height_type()

    @property
    def max_dataset_height(self) -> S1xxAttributesBase:
        return self._attributes[self.max_dataset_height_attribute_name]

    @max_dataset_height.setter
    def max_dataset_height(self, val: S1xxAttributesBase):
        self._attributes[self.max_dataset_height_attribute_name] = val

    @property
    def max_dataset_height_type(self) -> Type[numpy.float32]:
        """Defines datatype"""
        return numpy.float32

    def max_dataset_height_create(self):
        """ Creates a blank, empty or zero value for max_dataset_height"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.max_dataset_height = self.max_dataset_height_type()

    @property
    def method_water_level_product(self) -> S1xxAttributesBase:
        return self._attributes[self.method_water_level_product_attribute_name]

    @method_water_level_product.setter
    def method_water_level_product(self, val: S1xxAttributesBase):
        self._attributes[self.method_water_level_product_attribute_name] = val

    @property
    def method_water_level_product_type(self) -> Type[str]:
        """Defines datatype"""
        return str

    def method_water_level_product_create(self):
        """ Creates a blank, empty or zero value for method_currents_product"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.method_currents_product = self.method_water_level_product_type()

    @property
    def type_of_water_level_data(self) -> TYPE_OF_WATER_LEVEL_DATA:
        return self._attributes[self.type_of_water_level_data_attribute_name]

    @type_of_water_level_data.setter
    def type_of_water_level_data(self, val: Union[int, str, TYPE_OF_WATER_LEVEL_DATA]):
        self.set_enum_attribute(val, self.type_of_water_level_data_attribute_name, self.type_of_water_level_data_type)

    @property
    def type_of_water_level_data_type(self) -> Type[TYPE_OF_WATER_LEVEL_DATA]:
        return TYPE_OF_WATER_LEVEL_DATA

    def type_of_water_level_data_create(self):
        """ Creates a value using the first item in the enumeration of type_of_water_level_data"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.type_of_water_level_data = list(self.type_of_water_level_data_type)[0]


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

    water_level_feature_dataset_attribute_name = WATER_LEVEL

    @property
    def __version__(self) -> int:
        return 1

    def feature_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_code = self.feature_code_type([WATER_LEVEL], dtype=h5py_string_dtype)

    @property
    def water_level_feature_dataset_type(self):
        return WaterLevelFeatureDataset

    def water_level_feature_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level_feature_dataset = self.water_level_feature_dataset_type()

    @property
    def water_level_feature_dataset(self) -> WaterLevelFeatureDataset:
        return self._attributes[self.water_level_feature_dataset_attribute_name]

    @water_level_feature_dataset.setter
    def water_level_feature_dataset(self, val: WaterLevelFeatureDataset):
        self._attributes[self.water_level_feature_dataset_attribute_name] = val


class S104Root(S100Root):
    """The root group contains a feature information group and N feature containers.
    In S104 there is one feature container 'water level'.
    The coverage names are determined from the matching CoveragesAttributes
    Table 3 of v0.0.7
    """
    water_level_attribute_name = WATER_LEVEL
    water_level_trend_threshold_attribute_name = "waterLevelTrendThreshold"
    vertical_coordinate_system_attribute_name = "verticalCS"
    vertical_coordinate_base_attribute_name = "verticalCoordinateBase"
    vertical_datum_reference_attribute_name = "verticalDatumReference"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def feature_information_type(self):
        return GroupF

    @property
    def water_level(self) -> S1xxAttributesBase:
        return self._attributes[self.water_level_attribute_name]

    @property
    def water_level_type(self):
        return WaterLevelContainer

    def water_level_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level = self.water_level_type()

    @water_level.setter
    def water_level(self, val: S1xxAttributesBase):
        self._attributes[self.water_level_attribute_name] = val

    @property
    def water_level_trend_threshold(self) -> S1xxAttributesBase:
        return self._attributes[self.water_level_trend_threshold_attribute_name]

    @water_level_trend_threshold.setter
    def water_level_trend_threshold(self, val: S1xxAttributesBase):
        self._attributes[self.water_level_trend_threshold_attribute_name] = val

    @property
    def water_level_trend_threshold_type(self) -> Type[numpy.float32]:
        return numpy.float32

    def water_level_trend_threshold_create(self):
        """ Creates a blank, empty or zero value for surface_current_depth"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.water_level_trend_threshold = self.water_level_trend_threshold_type()

    @property
    def vertical_coordinate_system(self) -> S1xxAttributesBase:
        return self._attributes[self.vertical_coordinate_system_attribute_name]

    @vertical_coordinate_system.setter
    def vertical_coordinate_system(self, val: S1xxAttributesBase):
        self._attributes[self.vertical_coordinate_system_attribute_name] = val

    @property
    def vertical_coordinate_system_type(self) -> Type[numpy.int32]:
        """Define S104 datatype"""
        return numpy.int32

    def vertical_coordinate_system_create(self):
        """ Creates a blank, empty or zero value for surface_current_depth"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_coordinate_system = self.vertical_coordinate_system_type()

    @property
    def vertical_coordinate_base(self) -> VERTICAL_COORDINATE_BASE:
        return self._attributes[self.vertical_coordinate_base_attribute_name]

    @vertical_coordinate_base.setter
    def vertical_coordinate_base(self, val: Union[int, str, VERTICAL_COORDINATE_BASE]):
        self.set_enum_attribute(val, self.vertical_coordinate_base_attribute_name, self.vertical_coordinate_base_type)

    @property
    def vertical_coordinate_base_type(self) -> Type[VERTICAL_COORDINATE_BASE]:
        """Enumeration data type"""
        return VERTICAL_COORDINATE_BASE

    def vertical_coordinate_base_create(self):
        """ Creates a value using the first item in the enumeration of vertical_coordinate_base"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_coordinate_base = list(self.vertical_coordinate_base_type)[0]

    @property
    def vertical_datum_reference(self) -> VERTICAL_DATUM_REFERENCE:
        return self._attributes[self.vertical_datum_reference_attribute_name]

    @vertical_datum_reference.setter
    def vertical_datum_reference(self, val: Union[int, str, VERTICAL_DATUM_REFERENCE]):
        self.set_enum_attribute(val, self.vertical_datum_reference_attribute_name, self.vertical_datum_reference_type)

    @property
    def vertical_datum_reference_type(self) -> Type[VERTICAL_DATUM_REFERENCE]:
        """Defines enumeration datatype"""
        return VERTICAL_DATUM_REFERENCE

    def vertical_datum_reference_create(self):
        """ Creates a value using the first item in the enumeration of vertical_datum_reference"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum_reference = list(self.vertical_datum_reference_type)[0]


class DiscoveryMetadata(S1xxAttributesBase):
    """ 12.2.6 of v1.0.1"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        raise NotImplementedError()


class S104File(S1XXFile):
    """ HDF5 file object"""
    PRODUCT_SPECIFICATION = 'INT.IHO.S-104.0.0'

    def __init__(self, *args, **kywrds):
        super().__init__(*args, root=S104Root, **kywrds)
