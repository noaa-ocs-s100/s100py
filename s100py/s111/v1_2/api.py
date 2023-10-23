from __future__ import annotations

from typing import Callable, Iterator, Union, Optional, List, Type
from enum import Enum
import numpy
import datetime

try:
    from ... import s1xx, s100
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py.s111"

from ...s1xx import s1xx_sequence, S1xxObject, S1xxCollection, S1xxDatasetBase, S1xxGridsBase, h5py_string_dtype
from ...s100.v5_0.api import S100File, S100Root, S100Exception, FeatureContainerDCF2, FeatureInstanceDCF2,\
    FeatureContainerDCF3, FeatureInstanceDCF3, FeatureInformation, FeatureInformationDataset, GroupFBase, VERTICAL_CS, \
    VERTICAL_DATUM_REFERENCE, VERTICAL_DATUM

EDITION = 1.2
PRODUCT_SPECIFICATION = 'INT.IHO.S-111.1.2'

CHANGELOG = """
v1.2 
* S100 Edition 5.0
* Enum list changed from TYPE_OF_CURRENT_DATA to DATA_DYNAMICITY --  Table 12.10 - S104_DataDynamicity
* Enum list values changed for DEPTH_TYPE_INDEX --  Table 12.11 -  S111_DepthTypeIndex
* General metadata "horizontalDatumReference" and "horizontalDatumValue removed
* General metadata attributes "horizontalCRS", "typeOfHorizontalCRS" were added -- Table 12.1
* General metadata "datasetDeliveryInterval" added
* The S-100 attribute verticalCoordinateBase is no longer used as of S-111 Edition 1.2 because
its “sea surface” and “sea bottom” values have been added to the vertical datums enumeration (Table 12-8)
* Additional Group F Feature Attribute - surface current time(surfaceCurrentTime) -- A-2 Feature Attributes
* Datatypes sizes have been specified for all attributes
"""


SURFACE_CURRENT = "SurfaceCurrent"

# Default fill value for NetCDF variables
FILLVALUE = -9999.00

# Default depth in meters
DEFAULT_TARGET_DEPTH = 4.5

"""Contains s111 metadata to pass to S111File.

PRODUCT_SPECIFICATION: The product specification used to create this dataset.
HORIZONTAL_CRS: EPSG code or -1 if user defined
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


class DATA_DYNAMICITY(Enum):
    """
     S111 v1.2 Table 12.10 - Classification of data according to the
     relationship between the time of its collection, generation, or
     calculation of generation parameters, in relation to the time of
     publication of the dataset.
    """
    observation = 1
    astronomicalPrediction = 2
    analysisOrHybrid = 3
    hydrodynamicHindcast = 4
    hydrodynamicForecast = 5


class DEPTH_TYPE_INDEX(Enum):
    """
     S111 v1.2 H.6.1 - The vertical location of the current in the water
     column is normally referenced to some vertical datum. In this Product
     Specification, the datum is selectable: it can be the sea surface, the
     sea bottom, or any of 30 standard tidal datums. The coordinate system
     axis is directed upward, so if the level of the current is below the
     datum, the depth will have a negative value. Levels referenced above
     the sea bottom will have a positive value. For a layer average, the
     thickness of the layer is specified as a positive value.
    """
    heightOrDepth = 1
    layerAverage = 2


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
    """ Basic template for the name of the attribute will be
     automatically determined based on the array position of
     the S111_MetadataList
    """

    __surface_current_group_hdf_name__ = "Group" + r"[\._]\d+"
    __uncertainty_dataset_hdf_name__ = "uncertainty"
    __data_dynamicity_hdf_name__ = "dataDynamicity"

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


class SurfaceCurrentFeatureInstanceDCF2(FeatureInstanceDCF2, SurfaceCurrentFeatureInstanceBase):
    @property
    def __number_of_times_type__(self) -> Type[int]:
        return numpy.uint32

    @property
    def __time_record_interval_type__(self) -> Type[int]:
        return numpy.uint16


class SurfaceCurrentFeatureInstanceDCF3(FeatureInstanceDCF3, SurfaceCurrentFeatureInstanceBase):
    @property
    def __number_of_nodes_type__(self) -> Type[int]:
        return numpy.uint32

    @property
    def __number_of_times_type__(self) -> Type[int]:
        return numpy.uint32

    @property
    def __time_record_interval_type__(self) -> Type[int]:
        return numpy.uint16


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
    def __min_dataset_current_speed_type__(self) -> Type[numpy.float64]:
        return numpy.float64

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
    def __max_dataset_current_speed_type__(self) -> Type[numpy.float64]:
        return numpy.float64

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
    def type_of_current_data(self) -> DATA_DYNAMICITY:
        return self._attributes[self.__type_of_current_data_hdf_name__]

    @type_of_current_data.setter
    def type_of_current_data(self, val: Union[int, str, DATA_DYNAMICITY]):
        self.set_enum_attribute(val, self.__type_of_current_data_hdf_name__, self.__type_of_current_data_type__)

    @property
    def __type_of_current_data_type__(self) -> Type[DATA_DYNAMICITY]:
        return DATA_DYNAMICITY

    def type_of_current_data_create(self):
        """ Creates a value using the first item in the enumeration of type_of_current_data"""
        # Make the enum into a list and take the first value
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.type_of_current_data = list(self.__type_of_current_data_type__)[0]


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
    __dataset_delivery_interval_hdf_name__ = "datasetDeliveryInterval"
    __issue_time_hdf_name__ = "issueTime"
    __issue_date_hdf_name__ = "issueDate"

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
    def __issue_date_repr__(self) -> str:
        return self._attributes[self.__issue_date_hdf_name__].strftime("%Y%m%d")

    @property
    def __issue_time_repr__(self) -> str:
        return self._attributes[self.__issue_time_hdf_name__].strftime('%H%M%SZ')


class DiscoveryMetadata(S1xxObject):
    """ 12.2.6 of v1.2.0
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        raise NotImplementedError()


class S111File(S100File):
    PRODUCT_SPECIFICATION = 'INT.IHO.S-111.1.2'

    @staticmethod
    def make_container_for_dcf(data_coding_format):
        if data_coding_format == 2:
            container = SurfaceCurrentContainerDCF2()
        elif data_coding_format == 3:
            container = SurfaceCurrentContainerDCF3()
        else:
            raise S111Exception(f"DCF {data_coding_format} not supported")
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
