from __future__ import annotations

from typing import Callable, Iterator, Union, Optional, List, Type
from enum import Enum
import numpy

try:
    from ... import s1xx
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py.s111"

from ...s1xx import s1xx_sequence, S1xxObject, S1xxCollection, S1xxDatasetBase, S1xxGridsBase, S1XXFile, h5py_string_dtype
from ...s100.v4_0.api import S100File, S100Root, S100Exception, FeatureContainerDCF2, FeatureInstanceDCF2, \
    FeatureContainerDCF3, FeatureInstanceDCF3, FeatureInformation, FeatureInformationDataset, GroupFBase

EDITION = 1.0
PRODUCT_SPECIFICATION = 'INT.IHO.S-111.1.0'

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


class S111Exception(S100Exception):
    pass


class S111UnspecifiedClassException(S100Exception):
    pass


class S111MetadataListBase(S1xxCollection):
    write_format_str = ".%02d"


class SurfaceCurrentUncertaintyInformation(S1xxObject):
    __name_hdf_name__ = "name"  #: HDF5 naming
    __value_hdf_name__ = "value"  #: HDF5 naming

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
        return self._attributes[self.__name_hdf_name__]

    @name.setter
    def name(self, val: str):
        self._attributes[self.__name_hdf_name__] = val

    @property
    def __name_type__(self):
        return str

    def name_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.name = self.__name_type__()

    @property
    def value(self) -> str:
        """ The uncertainty value"""
        return self._attributes[self.__value_hdf_name__]

    @value.setter
    def value(self, val: int):
        self._attributes[self.__value_hdf_name__] = val

    @property
    def __value_type__(self):
        return str

    def value_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.value = self.__value_type__()


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


class SurfaceCurrentValues(S1xxGridsBase):
    __surface_current_speed_hdf_name__ = "surfaceCurrentSpeed"  #: HDF5 naming
    __surface_current_direction_hdf_name__ = "surfaceCurrentDirection"  #: HDF5 naming

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "values"

    @property
    def surface_current_speed(self) -> s1xx_sequence:
        return self._attributes[self.__surface_current_speed_hdf_name__]

    @surface_current_speed.setter
    def surface_current_speed(self, val: s1xx_sequence):
        self._attributes[self.__surface_current_speed_hdf_name__] = val

    @property
    def __surface_current_speed_type__(self) -> s1xx_sequence:
        return numpy.ndarray

    @property
    def surface_current_speed_dtype(self) -> Type[float]:
        return numpy.float32

    def surface_current_speed_create(self):
        """ Creates a blank, empty or zero value for surface_current_speed"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current_speed = self.__surface_current_speed_type__([], self.surface_current_speed_dtype)

    @property
    def surface_current_direction(self) -> s1xx_sequence:
        return self._attributes[self.__surface_current_direction_hdf_name__]

    @surface_current_direction.setter
    def surface_current_direction(self, val: s1xx_sequence):
        self._attributes[self.__surface_current_direction_hdf_name__] = val

    @property
    def __surface_current_direction_type__(self) -> s1xx_sequence:
        return numpy.ndarray

    @property
    def surface_current_direction_dtype(self) -> Type[float]:
        return numpy.float32

    def surface_current_direction_create(self):
        """ Creates a blank, empty or zero value for surface_current_direction"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current_direction = self.__surface_current_direction_type__([], self.surface_current_direction_dtype)

    def get_write_order(self):
        return [self.__surface_current_speed_hdf_name__, self.__surface_current_direction_hdf_name__]

    def get_compound_dtype(self):
        return [self.surface_current_speed_dtype, self.surface_current_direction_dtype]


class SurfaceCurrentGroup(S1xxObject):
    """ 10.2.5 of v1.0.1
    also see section 12.3 and table 12.5

    """
    __values_hdf_name__ = "values"  #: HDF5 naming
    __time_point_hdf_name__ = "timePoint"  #: HDF5 naming

    @property
    def values(self) -> SurfaceCurrentValues:
        return self._attributes[self.__values_hdf_name__]

    @values.setter
    def values(self, val: SurfaceCurrentValues):
        self._attributes[self.__values_hdf_name__] = val

    @property
    def __values_type__(self) -> Type[SurfaceCurrentValues]:
        return SurfaceCurrentValues

    def values_create(self):
        """ Creates a blank, empty or zero value for values"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.values = self.__values_type__()

    @property
    def time_point(self) -> S1xxObject:
        return self._attributes[self.__time_point_hdf_name__]

    @time_point.setter
    def time_point(self, val: S1xxObject):
        self._attributes[self.__time_point_hdf_name__] = val

    @property
    def __time_point_type__(self) -> Type[str]:
        return str

    def time_point_create(self):
        """ Creates a blank, empty or zero value for time_point"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.time_point = self.__time_point_type__()

    @property
    def __version__(self) -> int:
        return 1


class SurfaceCurrentGroupList(S1xxCollection):
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


class SurfaceCurrentFeatureInstanceBase:
    """ Basic template for the name of the attribute will be automatically
    determined based on the array position of the S111_MetadataList
    """

    __surface_current_group_hdf_name__ = "Group" + r"[\._]\d+"
    __uncertainty_dataset_hdf_name__ = "uncertainty"

    @property
    def __surface_current_group_type__(self):
        return SurfaceCurrentGroupList

    def surface_current_group_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current_group = self.__surface_current_group_type__()

    @property
    def surface_current_group(self) -> S1xxCollection:
        return self._attributes[self.__surface_current_group_hdf_name__]

    @surface_current_group.setter
    def surface_current_group(self, val: S1xxCollection):
        self._attributes[self.__surface_current_group_hdf_name__] = val

    @property
    def uncertainty_dataset(self) -> S1xxDatasetBase:
        return self._attributes[self.__uncertainty_dataset_hdf_name__]

    @uncertainty_dataset.setter
    def uncertainty_dataset(self, val: S1xxDatasetBase):
        self._attributes[self.__uncertainty_dataset_hdf_name__] = val

    @property
    def __uncertainty_dataset_type__(self) -> Type[SurfaceCurrentUncertaintyDataset]:
        return SurfaceCurrentUncertaintyDataset

    def uncertainty_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.uncertainty_dataset = self.__uncertainty_dataset_type__()


class SurfaceCurrentFeatureInstanceDCF2(FeatureInstanceDCF2, SurfaceCurrentFeatureInstanceBase):
    pass


class SurfaceCurrentFeatureInstanceDCF3(FeatureInstanceDCF3, SurfaceCurrentFeatureInstanceBase):
    pass


class SurfaceCurrentListBase(S111MetadataListBase):
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


class SurfaceCurrentListDCF2(SurfaceCurrentListBase):
    @property
    def metadata_type(self) -> Type[SurfaceCurrentFeatureInstanceBase]:
        return SurfaceCurrentFeatureInstanceDCF2


class SurfaceCurrentListDCF3(SurfaceCurrentListBase):
    @property
    def metadata_type(self) -> Type[SurfaceCurrentFeatureInstanceBase]:
        return SurfaceCurrentFeatureInstanceDCF3


class SurfaceCurrentContainerBase:
    """ This is the SurfaceCurrent right off the root of the HDF5 which has possible attributes from S100 spec table 10c-10
    This will hold child groups named SurfaceCurrent.NN
    """

    #: Basic template for the name of the attribute
    #: Attribute name will be automatically determined based on the containing list's index
    __surface_current_hdf_name__ = SURFACE_CURRENT + r"[\._]\d+"
    __min_dataset_current_speed_hdf_name__ = "minDatasetCurrentSpeed"
    __max_dataset_current_speed_hdf_name__ = "maxDatasetCurrentSpeed"
    __method_currents_product_hdf_name__ = "methodCurrentsProduct"
    __type_of_current_data_hdf_name__ = "typeOfCurrentData"

    @property
    def __version__(self) -> int:
        return 1

    def surface_current_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current = self.__surface_current_type__()

    @property
    def surface_current(self) -> S111MetadataListBase:
        """ The surface current data, a list of SurfaceCurrent
        Returns
        -------
        S111_MetadataList_base
            Contains a list of SurfaceCurrent objects via the SurfaceCurrent_List class
        """
        return self._attributes[self.__surface_current_hdf_name__]

    @surface_current.setter
    def surface_current(self, val: S111MetadataListBase):
        self._attributes[self.__surface_current_hdf_name__] = val

    @property
    def min_dataset_current_speed(self) -> S1xxObject:
        return self._attributes[self.__min_dataset_current_speed_hdf_name__]

    @min_dataset_current_speed.setter
    def min_dataset_current_speed(self, val: S1xxObject):
        self._attributes[self.__min_dataset_current_speed_hdf_name__] = val

    @property
    def __min_dataset_current_speed_type__(self) -> Type[float]:
        return float

    def min_dataset_current_speed_create(self):
        """ Creates a blank, empty or zero value for min_dataset_current_speed"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.min_dataset_current_speed = self.__min_dataset_current_speed_type__()

    @property
    def max_dataset_current_speed(self) -> S1xxObject:
        return self._attributes[self.__max_dataset_current_speed_hdf_name__]

    @max_dataset_current_speed.setter
    def max_dataset_current_speed(self, val: S1xxObject):
        self._attributes[self.__max_dataset_current_speed_hdf_name__] = val

    @property
    def __max_dataset_current_speed_type__(self) -> Type[numpy.float32]:
        return numpy.float32

    def max_dataset_current_speed_create(self):
        """ Creates a blank, empty or zero value for max_dataset_current_speed"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.max_dataset_current_speed = self.__max_dataset_current_speed_type__()

    @property
    def method_currents_product(self) -> S1xxObject:
        return self._attributes[self.__method_currents_product_hdf_name__]

    @method_currents_product.setter
    def method_currents_product(self, val: S1xxObject):
        self._attributes[self.__method_currents_product_hdf_name__] = val

    @property
    def __method_currents_product_type__(self) -> Type[str]:
        return str

    def method_currents_product_create(self):
        """ Creates a blank, empty or zero value for method_currents_product"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.method_currents_product = self.__method_currents_product_type__()

    @property
    def type_of_current_data(self) -> TYPE_OF_CURRENT_DATA:
        return self._attributes[self.__type_of_current_data_hdf_name__]

    @type_of_current_data.setter
    def type_of_current_data(self, val: Union[int, str, TYPE_OF_CURRENT_DATA]):
        self.set_enum_attribute(val, self.__type_of_current_data_hdf_name__, self.__type_of_current_data_type__)

    @property
    def __type_of_current_data_type__(self) -> Type[TYPE_OF_CURRENT_DATA]:
        return TYPE_OF_CURRENT_DATA

    def type_of_current_data_create(self):
        """ Creates a value using the first item in the enumeration of type_of_current_data"""
        # Make the enum into a list and take the first value
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.type_of_current_data = list(self.__type_of_current_data_type__)[0]


class SurfaceCurrentListDCF2(SurfaceCurrentListBase):
    @property
    def metadata_type(self) -> Type[SurfaceCurrentFeatureInstanceBase]:
        return SurfaceCurrentFeatureInstanceDCF2


class SurfaceCurrentListDCF3(SurfaceCurrentListBase):
    @property
    def metadata_type(self) -> Type[SurfaceCurrentFeatureInstanceBase]:
        return SurfaceCurrentFeatureInstanceDCF3


class SurfaceCurrentContainerDCF2(FeatureContainerDCF2, SurfaceCurrentContainerBase):
    @property
    def __surface_current_type__(self):
        return SurfaceCurrentListDCF2


class SurfaceCurrentContainerDCF3(FeatureContainerDCF3, SurfaceCurrentContainerBase):
    @property
    def __surface_current_type__(self):
        return SurfaceCurrentListDCF3


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


class GroupF(GroupFBase):
    """ Table 10.3 and sect 10.2.2 of v1.0.1
    """
    __surface_current_feature_dataset_hdf_name__ = SURFACE_CURRENT

    @property
    def __version__(self) -> int:
        return 1

    def feature_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_code = self.__feature_code_type__([SURFACE_CURRENT], dtype=h5py_string_dtype)

    @property
    def __surface_current_feature_dataset_type__(self):
        return SurfaceCurrentFeatureDataset

    def surface_current_feature_dataset_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current_feature_dataset = self.__surface_current_feature_dataset_type__()

    @property
    def surface_current_feature_dataset(self) -> SurfaceCurrentFeatureDataset:
        return self._attributes[self.__surface_current_feature_dataset_hdf_name__]

    @surface_current_feature_dataset.setter
    def surface_current_feature_dataset(self, val: SurfaceCurrentFeatureDataset):
        self._attributes[self.__surface_current_feature_dataset_hdf_name__] = val


class S111Root(S100Root):
    """The root group contains a feature information group and N feature containers.
    In S111 there is one feature container 'surface current'.
    The coverage names are determined from the matching CoveragesAttributes
    10.2.1, 10.2.2 and Table 12.1 of v1.0.1
    """
    __surface_current_hdf_name__ = SURFACE_CURRENT
    __depth_type_index_hdf_name__ = "depthTypeIndex"
    __surface_current_depth_hdf_name__ = "surfaceCurrentDepth"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def __feature_information_type__(self):
        return GroupF

    @property
    def surface_current(self) -> S1xxObject:
        return self._attributes[self.__surface_current_hdf_name__]

    @property
    def __surface_current_type__(self):
        return SurfaceCurrentContainerBase

    def surface_current_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        print("You should not create a generic surface_current object but manually create"
              "a SurfaceCurrentContainerDCF2 or SurfaceCurrentContainerDCF3")
        raise S111UnspecifiedClassException("You should create SurfaceCurrentContainerDCFx (x=2,3,8 etc)")
        self.surface_current = self.__surface_current_type__()

    @surface_current.setter
    def surface_current(self, val: S1xxObject):
        self._attributes[self.__surface_current_hdf_name__] = val

    @property
    def depth_type_index(self) -> DEPTH_TYPE_INDEX:
        return self._attributes[self.__depth_type_index_hdf_name__]

    @depth_type_index.setter
    def depth_type_index(self, val: Union[int, str, DEPTH_TYPE_INDEX]):
        self.set_enum_attribute(val, self.__depth_type_index_hdf_name__, self.__depth_type_index_type__)

    @property
    def __depth_type_index_type__(self) -> Type[DEPTH_TYPE_INDEX]:
        return DEPTH_TYPE_INDEX

    def depth_type_index_create(self):
        """ Creates a value using the first item in the enumeration of depth_type_index"""
        # Make the enum into a list and take the first value
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.depth_type_index = list(self.__depth_type_index_type__)[0]

    @property
    def surface_current_depth(self) -> S1xxObject:
        return self._attributes[self.__surface_current_depth_hdf_name__]

    @surface_current_depth.setter
    def surface_current_depth(self, val: S1xxObject):
        self._attributes[self.__surface_current_depth_hdf_name__] = val

    @property
    def __surface_current_depth_type__(self) -> Type[numpy.float32]:
        return numpy.float32

    def surface_current_depth_create(self):
        """ Creates a blank, empty or zero value for surface_current_depth"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.surface_current_depth = self.__surface_current_depth_type__()


class DiscoveryMetadata(S1xxObject):
    """ 12.2.6 of v1.0.1
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        raise NotImplementedError()


class S111File(S100File):

    PRODUCT_SPECIFICATION = 'INT.IHO.S-111.1.0'

    @staticmethod
    def make_container_for_dcf(data_coding_format):
        if data_coding_format == 2:
            container = SurfaceCurrentContainerDCF2()
        elif data_coding_format == 3:
            container = SurfaceCurrentContainerDCF3()
        else:
            raise S111Exception("DCF {} not supported".format(data_coding_format))
        return container

    def __init__(self, *args, **kywrds):
        super().__init__(*args, root=S111Root, **kywrds)
        # when reading from a file we need to look inside the DataCodingFormat to know what type of object to create
        try:
            container_key = f'/{SURFACE_CURRENT}'
            dcf = self[container_key].attrs['dataCodingFormat']
        except KeyError:  # must not be reading an existing file or doesn't have data for some reason
            pass
        else:
            self.root.surface_current = self.make_container_for_dcf(dcf)
            self.root.surface_current.read(self[container_key])
