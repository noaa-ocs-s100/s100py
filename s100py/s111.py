import os
from abc import ABC
import datetime
import logging
from typing import Callable, Iterator, Union, Optional, List, Type
from enum import Enum
import numpy

from s100py.s1xx import s1xx_sequence, S1xxAttributesBase, S1xxMetadataListBase, S1xxDatasetBase, S1xxGridsBase, S1XXFile
from s100py.s100 import FeatureContainer, S100Root, FeatureInstanceDCF2, FeatureInformation, FeatureInformationDataset

SURFACE_CURRENT = "SurfaceCurrent"

# Default fill value for NetCDF variables
FILLVALUE = -9999.0

# Default depth in meters
DEFAULT_TARGET_DEPTH = 4.5

TYPE_OF_CURRENT_DATA = Enum(value="TYPE_OF_CURRENT_DATA",
                            names=[("Historical observation (O)", 1),
                                   ("Real-time observation (R)", 2),
                                   ("Astronomical prediction (A)", 3),
                                   ("Analysis or hybrid method (Y)", 4),
                                   ("Hydrodynamic model hindcast (M)", 5),
                                   ("Hydrodynamic model forecast (F)", 6)
                                   ]
                            )

DEPTH_TYPE_INDEX = Enum(value="DEPTH_TYPE_INDEX",
                        names=[("Layer average", 1),
                               ("Sea surface", 2),
                               ("Vertical datum", 3),
                               ("Sea bottom", 4)
                               ]
                        )

"""Contains s111 metadata to pass to S111File.

PRODUCT_SPECIFICATION: The product specification used to create this dataset.
HORIZONTAL_DATUM_REFERENCE: Reference to the register from which the horizontal datum value is taken.
DATA_CODING_FORMAT: Reference to the type of S111 product.
INTERPOLATION_TYPE: Interpolation method recommended for evaluation of the S100_GridCoverage.
COMMON_POINT_RULE: The procedure used for evaluating geometric objects that overlap or lie fall on boundaries.
DIMENSION: The dimension of the feature instance.
SEQUENCING_RULE_TYPE: Method to assign values from the sequence of values to the grid coordinates (e.g. "linear").
SEQUENCING_RULE_SCAN_DIRECTION: AxisNames, comma-separated (e.g. "longitude,latitude").
START_SEQUENCE: Starting location of the scan.

"""


class S111MetadataListBase(S1xxMetadataListBase):
    pass


class SurfaceCurrentUncertaintyInformation(S1xxAttributesBase):
    name_attribute_name = "name"  #: HDF5 naming
    value_attribute_name = "value"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def name(self) -> str:
        """ The plain text name of the data
        Returns
        -------
        str
            Name of the dataset ("surfaceCurrentSpeed" or "surfaceCurrentDirection")
        """
        return self._attributes[self.name_attribute_name]

    @name.setter
    def name(self, val: str):
        self._attributes[self.name_attribute_name] = val

    @property
    def name_type(self):
        return str

    def name_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.name = self.name_type()

    @property
    def value(self) -> str:
        """ The uncertainty value"""
        return self._attributes[self.value_attribute_name]

    @value.setter
    def value(self, val: int):
        self._attributes[self.value_attribute_name] = val

    @property
    def value_type(self):
        return str

    def value_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.value = self.value_type()


class SurfaceCurrentUncertaintyDataset(S1xxDatasetBase):

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "uncertainty"

    @property
    def metadata_type(self) -> Type[SurfaceCurrentUncertaintyInformation]:
        return SurfaceCurrentUncertaintyInformation


class GeometryValuesDataset(S1xxGridsBase):
    longitude_attribute_name = "longitude"
    latitude_attribute_name = "latitude"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "geometryValues"

    @property
    def longitude(self) -> s1xx_sequence:
        return self._attributes[self.longitude_attribute_name]

    @longitude.setter
    def longitude(self, val: s1xx_sequence):
        self._attributes[self.longitude_attribute_name] = val

    @property
    def longitude_type(self) -> s1xx_sequence:
        return numpy.ndarray

    def longitude_create(self):
        """ Creates a blank, empty or zero value for longitude"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.longitude = self.longitude_type([], numpy.float)

    @property
    def latitude(self) -> s1xx_sequence:
        return self._attributes[self.latitude_attribute_name]

    @latitude.setter
    def latitude(self, val: s1xx_sequence):
        self._attributes[self.latitude_attribute_name] = val

    @property
    def latitude_type(self) -> s1xx_sequence:
        return numpy.ndarray

    def latitude_create(self):
        """ Creates a blank, empty or zero value for latitude"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.latitude = self.latitude_type([], numpy.float)

    def get_write_order(self):
        return [self.longitude_attribute_name, self.latitude_attribute_name]


class PositioningGroup(S1xxAttributesBase):

    geometry_values_attribute_name = "geometry_values"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "Positioning"

    @property
    def metadata_type(self) -> type:
        return GeometryValuesDataset

    @property
    def geometry_values(self) -> GeometryValuesDataset:
        return self._attributes[self.geometry_values_attribute_name]

    @geometry_values.setter
    def geometry_values(self, val: GeometryValuesDataset):
        self._attributes[self.geometry_values_attribute_name] = val

    @property
    def geometry_values_type(self) -> Type[GeometryValuesDataset]:
        return GeometryValuesDataset

    def geometry_values_create(self):
        """ Creates a blank, empty or zero value for geometry_values"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.geometry_values = self.geometry_values_type()


class SurfaceCurrentValues(S1xxGridsBase):
    surface_current_speed_attribute_name = "surfaceCurrentSpeed"  #: HDF5 naming
    surface_current_direction_attribute_name = "surfaceCurrentDirection"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "values"

    @property
    def surface_current_speed(self) -> s1xx_sequence:
        return self._attributes[self.surface_current_speed_attribute_name]

    @surface_current_speed.setter
    def surface_current_speed(self, val: s1xx_sequence):
        self._attributes[self.surface_current_speed_attribute_name] = val

    @property
    def surface_current_speed_type(self) -> s1xx_sequence:
        return numpy.ndarray

    def surface_current_speed_create(self):
        """ Creates a blank, empty or zero value for surface_current_speed"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current_speed = self.surface_current_speed_type([], numpy.float)

    @property
    def surface_current_direction(self) -> s1xx_sequence:
        return self._attributes[self.surface_current_direction_attribute_name]

    @surface_current_direction.setter
    def surface_current_direction(self, val: s1xx_sequence):
        self._attributes[self.surface_current_direction_attribute_name] = val

    @property
    def surface_current_direction_type(self) -> s1xx_sequence:
        return numpy.ndarray

    def surface_current_direction_create(self):
        """ Creates a blank, empty or zero value for surface_current_direction"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current_direction = self.surface_current_direction_type([], numpy.float)

    def get_write_order(self):
        return [self.surface_current_speed_attribute_name, self.surface_current_direction_attribute_name]


class SurfaceCurrentGroup(S1xxAttributesBase):
    """ 10.2.5 of v1.0.1
    also see section 12.3 and table 12.5

    """

    values_attribute_name = "values"  #: HDF5 naming
    time_point_attribute_name = "timePoint"  #: HDF5 naming

    @property
    def values(self) -> SurfaceCurrentValues:
        return self._attributes[self.values_attribute_name]

    @values.setter
    def values(self, val: SurfaceCurrentValues):
        self._attributes[self.values_attribute_name] = val

    @property
    def values_type(self) -> Type[SurfaceCurrentValues]:
        return SurfaceCurrentValues

    def values_create(self):
        """ Creates a blank, empty or zero value for values"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.values = self.values_type()

    @property
    def time_point(self) -> S1xxAttributesBase:
        return self._attributes[self.time_point_attribute_name]

    @time_point.setter
    def time_point(self, val: S1xxAttributesBase):
        self._attributes[self.time_point_attribute_name] = val

    @property
    def time_point_type(self) -> Type[str]:
        return str

    def time_point_create(self):
        """ Creates a blank, empty or zero value for time_point"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.time_point = self.time_point_type()

    @property
    def __version__(self) -> int:
        return 1


class SurfaceCurrentGroupList(S1xxMetadataListBase):
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
        return SurfaceCurrentGroup


class SurfaceCurrentFeatureInstance(FeatureInstanceDCF2):
    surface_current_group_attribute_name = "Group" + r"[\._]\d+"
    """ Basic template for the name of the attribute
    Attribute name will be automatically determined based on the array position of the S111_MetadataList
    """

    uncertainty_dataset_attribute_name = "uncertainty"
    number_of_nodes_attribute_name = "numberOfNodes"
    datetime_first_record_attribute_name = "dateTimeOfFirstRecord"
    datetime_last_record_attribute_name = "dateTimeOfLastRecord"

    @property
    def surface_current_group_type(self):
        return SurfaceCurrentGroupList

    def surface_current_group_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current_group = self.surface_current_group_type()

    @property
    def surface_current_group(self) -> S1xxMetadataListBase:
        return self._attributes[self.surface_current_group_attribute_name]

    @surface_current_group.setter
    def surface_current_group(self, val: S1xxMetadataListBase):
        self._attributes[self.surface_current_group_attribute_name] = val

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
        return self._attributes[self.uncertainty_dataset_attribute_name]

    @uncertainty_dataset.setter
    def uncertainty_dataset(self, val: S1xxDatasetBase):
        self._attributes[self.uncertainty_dataset_attribute_name] = val

    @property
    def uncertainty_dataset_type(self) -> Type[SurfaceCurrentUncertaintyDataset]:
        return SurfaceCurrentUncertaintyDataset

    def uncertainty_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.uncertainty_dataset = self.uncertainty_dataset_type()

    @property
    def positioning_group_attribute_name(self) -> str:
        return "Positioning"

    @property
    def positioning_group(self) -> S1xxAttributesBase:
        return self._attributes[self.positioning_group_attribute_name]

    @positioning_group.setter
    def positioning_group(self, val: S1xxAttributesBase):
        self._attributes[self.positioning_group_attribute_name] = val

    @property
    def positioning_group_type(self):
        return PositioningGroup

    def positioning_group_create(self):
        """ Creates a blank, empty or zero value for positioning_group"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.positioning_group = self.positioning_group_type()


class SurfaceCurrentList(S1xxMetadataListBase):
    """ Sect 10.2.4 and Table 12.3 of v1.0.1
    This is the set of SurfaceCurrent.NN that act like a list here.
    They will contain a list of Groups.NNN as well as other attributes etc.
    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return SURFACE_CURRENT

    @property
    def metadata_type(self) -> Type[SurfaceCurrentFeatureInstance]:
        return SurfaceCurrentFeatureInstance


class SurfaceCurrentContainer(FeatureContainer):
    """ This is the SurfaceCurrent right off the root of the HDF5 which has possible attributes from S100 spec table 10c-10
    This will hold child groups named SurfaceCurrent.NN
    """

    #: Basic template for the name of the attribute
    #: Attribute name will be automatically determined based on the containing list's index
    surface_current_attribute_name = SURFACE_CURRENT + r"[\._]\d+"
    min_dataset_current_speed_attribute_name = "minDatasetCurrentSpeed"
    max_dataset_current_speed_attribute_name = "maxDatasetCurrentSpeed"
    method_currents_product_attribute_name = "methodCurrentsProduct"
    type_of_current_data_attribute_name = "typeOfCurrentData"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def surface_current_type(self):
        return SurfaceCurrentList

    def surface_current_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current = self.surface_current_type()

    @property
    def surface_current(self) -> S1xxMetadataListBase:
        """ The surface current data, a list of SurfaceCurrent
        Returns
        -------
        S111_MetadataList_base
            Contains a list of SurfaceCurrent objects via the SurfaceCurrent_List class
        """
        return self._attributes[self.surface_current_attribute_name]

    @surface_current.setter
    def surface_current(self, val: S1xxMetadataListBase):
        self._attributes[self.surface_current_attribute_name] = val

    @property
    def min_dataset_current_speed(self) -> S1xxAttributesBase:
        return self._attributes[self.min_dataset_current_speed_attribute_name]

    @min_dataset_current_speed.setter
    def min_dataset_current_speed(self, val: S1xxAttributesBase):
        self._attributes[self.min_dataset_current_speed_attribute_name] = val

    @property
    def min_dataset_current_speed_type(self) -> Type[float]:
        return float

    def min_dataset_current_speed_create(self):
        """ Creates a blank, empty or zero value for min_dataset_current_speed"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.min_dataset_current_speed = self.min_dataset_current_speed_type()

    @property
    def max_dataset_current_speed(self) -> S1xxAttributesBase:
        return self._attributes[self.max_dataset_current_speed_attribute_name]

    @max_dataset_current_speed.setter
    def max_dataset_current_speed(self, val: S1xxAttributesBase):
        self._attributes[self.max_dataset_current_speed_attribute_name] = val

    @property
    def max_dataset_current_speed_type(self) -> Type[numpy.float32]:
        return numpy.float32

    def max_dataset_current_speed_create(self):
        """ Creates a blank, empty or zero value for max_dataset_current_speed"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.max_dataset_current_speed = self.max_dataset_current_speed_type()

    @property
    def method_currents_product(self) -> S1xxAttributesBase:
        return self._attributes[self.method_currents_product_attribute_name]

    @method_currents_product.setter
    def method_currents_product(self, val: S1xxAttributesBase):
        self._attributes[self.method_currents_product_attribute_name] = val

    @property
    def method_currents_product_type(self) -> Type[str]:
        return str

    def method_currents_product_create(self):
        """ Creates a blank, empty or zero value for method_currents_product"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.method_currents_product = self.method_currents_product_type()

    @property
    def type_of_current_data(self) -> TYPE_OF_CURRENT_DATA:
        return self._attributes[self.type_of_current_data_attribute_name]

    @type_of_current_data.setter
    def type_of_current_data(self, val: Union[int, str, TYPE_OF_CURRENT_DATA]):
        self.set_enum_attribute(val, self.type_of_current_data_attribute_name, self.type_of_current_data_type)

    @property
    def type_of_current_data_type(self) -> Type[TYPE_OF_CURRENT_DATA]:
        return TYPE_OF_CURRENT_DATA

    def type_of_current_data_create(self):
        """ Creates a value using the first item in the enumeration of type_of_current_data"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.type_of_current_data = list(self.type_of_current_data_type)[0]

    def data_coding_format_create(self):
        """ Creates a blank, empty or zero value for data_coding_format"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_coding_format = self.data_coding_format_type(2)  # default


class SurfaceCurrentFeatureDataset(FeatureInformationDataset):

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return SURFACE_CURRENT

    @property
    def metadata_type(self) -> Type[FeatureInformation]:
        return FeatureInformation


class GroupF(S1xxAttributesBase):
    """ Table 10.3 and sect 10.2.2 of v1.0.1
    """
    feature_code_attribute_name = "featureCode"
    surface_current_feature_dataset_attribute_name = SURFACE_CURRENT

    @property
    def __version__(self) -> int:
        return 1

    @property
    def feature_code_type(self):
        return numpy.array

    def feature_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_code = self.feature_code_type([SURFACE_CURRENT], dtype='S')

    @property
    def feature_code(self) -> s1xx_sequence:
        return self._attributes[self.feature_code_attribute_name]

    @feature_code.setter
    def feature_code(self, val: s1xx_sequence):
        self._attributes[self.feature_code_attribute_name] = val

    @property
    def surface_current_feature_dataset_type(self):
        return SurfaceCurrentFeatureDataset

    def surface_current_feature_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current_feature_dataset = self.surface_current_feature_dataset_type()

    @property
    def surface_current_feature_dataset(self) -> SurfaceCurrentFeatureDataset:
        return self._attributes[self.surface_current_feature_dataset_attribute_name]

    @surface_current_feature_dataset.setter
    def surface_current_feature_dataset(self, val: SurfaceCurrentFeatureDataset):
        self._attributes[self.surface_current_feature_dataset_attribute_name] = val


class S111Root(S100Root):
    """The root group contains a feature information group and N feature containers.
    In S111 there is one feature container 'surface current'.
    The coverage names are determined from the matching CoveragesAttributes
    10.2.1, 10.2.2 and Table 12.1 of v1.0.1
    """
    feature_information_attribute_name = "Group_F"
    surface_current_attribute_name = SURFACE_CURRENT
    depth_type_index_attribute_name = "depthTypeIndex"
    surface_current_depth_attribute_name = "surfaceCurrentDepth"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def feature_information(self) -> GroupF:
        return self._attributes[self.feature_information_attribute_name]

    @feature_information.setter
    def feature_information(self, val: GroupF):
        self._attributes[self.feature_information_attribute_name] = val

    @property
    def feature_information_type(self):
        return GroupF

    def feature_information_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_information = self.feature_information_type()

    @property
    def surface_current(self) -> S1xxAttributesBase:
        return self._attributes[self.surface_current_attribute_name]

    @property
    def surface_current_type(self):
        return SurfaceCurrentContainer

    def surface_current_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current = self.surface_current_type()

    @surface_current.setter
    def surface_current(self, val: S1xxAttributesBase):
        self._attributes[self.surface_current_attribute_name] = val

    @property
    def depth_type_index(self) -> DEPTH_TYPE_INDEX:
        return self._attributes[self.depth_type_index_attribute_name]

    @depth_type_index.setter
    def depth_type_index(self, val: Union[int, str, DEPTH_TYPE_INDEX]):
        self.set_enum_attribute(val, self.depth_type_index_attribute_name, self.depth_type_index_type)

    @property
    def depth_type_index_type(self) -> Type[DEPTH_TYPE_INDEX]:
        return DEPTH_TYPE_INDEX

    def depth_type_index_create(self):
        """ Creates a value using the first item in the enumeration of depth_type_index"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.depth_type_index = list(self.depth_type_index_type)[0]

    @property
    def surface_current_depth(self) -> S1xxAttributesBase:
        return self._attributes[self.surface_current_depth_attribute_name]

    @surface_current_depth.setter
    def surface_current_depth(self, val: S1xxAttributesBase):
        self._attributes[self.surface_current_depth_attribute_name] = val

    @property
    def surface_current_depth_type(self) -> Type[numpy.float32]:
        return numpy.float32

    def surface_current_depth_create(self):
        """ Creates a blank, empty or zero value for surface_current_depth"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current_depth = self.surface_current_depth_type()


class DiscoveryMetadata(S1xxAttributesBase):
    """ 12.2.6 of v1.0.1
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        raise NotImplementedError()


class S111File(S1XXFile):
    PRODUCT_SPECIFICATION = numpy.string_('INT.IHO.S-111.1.0')

    def __init__(self, *args, **kywrds):
        super().__init__(*args, root=S111Root, **kywrds)
