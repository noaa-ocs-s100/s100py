"""
S-104 Edition 1.0 utilized features from both S100 Edition 4.0 and 5.0, the Edition 5.0
S100 root object was copied and extended to support Edition 4.0 chunking attributes
"""
import datetime
from typing import Callable, Iterator, Union, Optional, List, Type
from enum import Enum, IntEnum
import numpy
import h5py

from ...s1xx import s1xx_sequence, S1xxObject, S1xxCollection, S1xxDatasetBase, S1xxGridsBase, S1XXFile, \
    h5py_string_dtype, is_sub_class
from ...s100.v4_0.api import S100Exception, FeatureContainerDCF2, FeatureInstanceDCF2, FeatureContainerDCF3, \
    FeatureInstanceDCF3, FeatureInformation, FeatureInformationDataset, GroupFBase, GeographicBoundingBox
from ...s100.v5_0.api import S100File, VERTICAL_DATUM, VERTICAL_DATUM_REFERENCE, VERTICAL_CS, VERTICAL_COORDINATE_BASE,\
    HORIZONTAL_CS, TYPE_OF_HORIZONTAL_CRS, PROJECTION_METHOD

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


EDITION = 1.0
PRODUCT_SPECIFICATION = 'INT.IHO.S-104.1.0'


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


class WaterLevelFeatureInstanceBase:
    """ Basic template for the name of the attribute will be automatically determined
    based on the array position of the S104_MetadataList
    """
    __water_level_group_hdf_name__ = "Group" + r"[\._]\d+"
    __uncertainty_dataset_hdf_name__ = "uncertainty"
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


class WaterLevelFeatureInstanceDCF2(FeatureInstanceDCF2, WaterLevelFeatureInstanceBase):
    pass


class WaterLevelFeatureInstanceDCF3(FeatureInstanceDCF3, WaterLevelFeatureInstanceBase):
    @property
    def __number_of_nodes_type__(self) -> Type[int]:
        return numpy.uint32

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


class WaterLevelListDCF3(WaterLevelListBase):
    @property
    def metadata_type(self) -> Type[WaterLevelFeatureInstanceBase]:
        return WaterLevelFeatureInstanceDCF3


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


class WaterLevelContainerDCF2(FeatureContainerDCF2, WaterLevelContainerBase):
    @property
    def __water_level_type__(self):
        return WaterLevelListDCF2


class WaterLevelContainerDCF3(FeatureContainerDCF3, WaterLevelContainerBase):
    @property
    def __water_level_type__(self):
        return WaterLevelListDCF3


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


class S100Root(GeographicBoundingBox):
    """ From table 10c-6 in S100 v4.0 spec.

    The root of the S100 v5.0 schema.  There are restrictions on many of the CRS attributes based on other attributes.
    For example, the horizontal_cs has different value options based on horizontal_crs.

    """
    __feature_information_hdf_name__ = "Group_F"
    __epoch_hdf_name__ = "epoch"
    __geographic_identifier_hdf_name__ = "geographicIdentifier"
    __vertical_datum_hdf_name__ = "verticalDatum"
    __meta_features_hdf_name__ = "metaFeatures"
    __metadata_hdf_name__ = "metadata"
    __product_specification_hdf_name__ = "productSpecification"
    __issue_time_hdf_name__ = "issueTime"
    __issue_date_hdf_name__ = "issueDate"
    # revised/added in S100 v5.0
    __horizontal_crs_hdf_name__ = "horizontalCRS"  #: HDF5 naming
    __name_of_horizontal_crs_hdf_name__ = "nameOfHorizontalCRS"  #: HDF5 naming
    __type_of_horizontal_crs_hdf_name__ = "typeOfHorizontalCRS"  #: HDF5 naming
    __horizontal_cs_hdf_name__ = "horizontalCS"  #: HDF5 naming
    __horizontal_datum_hdf_name__ = "horizontalDatum"  #: HDF5 naming
    __name_of_horizontal_datum_hdf_name__ = "nameOfHorizontalDatum"  #: HDF5 naming
    __prime_meridian_hdf_name__ = "primeMeridian"  #: HDF5 naming
    __spheriod_hdf_name__ = "spheriod"  #: HDF5 naming
    __projection_method_hdf_name__ = "projectionMethod"  #: HDF5 naming
    __projection_parameter_1_hdf_name__ = "projectionParameter1"  #: HDF5 naming
    __projection_parameter_2_hdf_name__ = "projectionParameter2"  #: HDF5 naming
    __projection_parameter_3_hdf_name__ = "projectionParameter3"  #: HDF5 naming
    __projection_parameter_4_hdf_name__ = "projectionParameter4"  #: HDF5 naming
    __projection_parameter_5_hdf_name__ = "projectionParameter5"  #: HDF5 naming
    __false_easting_hdf_name__ = "falseEasting"  #: HDF5 naming
    __false_northing_hdf_name__ = "falseNorthing"  #: HDF5 naming
    __vertical_cs_hdf_name__ = "verticalCS"  #: HDF5 naming
    __vertical_coordinate_base_hdf_name__ = "verticalCoordinateBase"  #: HDF5 naming
    __vertical_datum_reference_hdf_name__ = "verticalDatumReference"  #: HDF5 naming
    __vertical_datum_hdf_name__ = "verticalDatum"  #: HDF5 naming

    # Removed in v5.0
    # __horizontal_datum_reference_hdf_name__ = "horizontalDatumReference"
    # __horizontal_datum_value_hdf_name__ = "horizontalDatumValue"

    @property
    def __version__(self) -> int:
        return 1

    def write(self, group_object):
        super().write(group_object)
        # any grids that were had datasets which possible chunking should now be written
        # and we can look through those to get the overall chunking attribute
        # and put that into the GroupF FeatureInformation object
        feat_info = None
        for property_name in self.get_standard_properties():
            if is_sub_class(self.__getattribute__("__" + property_name + "_type__"), GroupFBase):
                feat_info = self.__getattribute__(property_name)
        # we have the GroupF data now, we can look at the names of the FeatureInstances and then search each for its respective chunking
        if feat_info is not None:
            # this will be the names of the feature instances
            # e.g. BathymetryCoverage for S102 or SurfaceCurrent for S111
            for feat_name in feat_info.feature_code:
                # get the associated python name for the feature, e.g. turn SurfaceCurrent into surface_current
                chunking = None
                try:
                    python_name = self.get_standard_properties_mapping()[feat_name]
                except KeyError:
                    python_name = self.get_standard_properties_mapping()[feat_name.decode()]
                # grab the root/SurfaceCurrent data
                feat_container = self.__getattribute__(python_name)
                # now look through all the SurfaceCurrent_01, SurfaceCurrent_02...
                # so find the list object (there really only should be one and it should match the naming but we'll be general here)
                for pattern, list_name in feat_container.get_standard_list_properties().items():
                    try:
                        list_of_features = feat_container.__getattribute__(list_name)
                    except KeyError:  # not initialized
                        list_of_features = []
                    for feat_instance in list_of_features:
                        try:
                            chunking = feat_instance.instance_chunking
                        except:
                            pass
                if chunking is not None:
                    # find the GroupF feature dataset, e.g. /GroupF/SurrfaceCurrent
                    try:
                        groupf_python_name = feat_info.get_standard_properties_mapping()[feat_name]
                    except KeyError:
                        groupf_python_name = feat_info.get_standard_properties_mapping()[feat_name.decode()]
                    # Get the python object
                    feat_info_dataset = feat_info.__getattribute__(groupf_python_name)
                    # set chunking
                    feat_info_dataset.chunking = chunking
                    # and now use HDF5 pathing to write the chunking part back out
                    # in theory whis would work from any level, not just the root
                    relative_hdf5_dataset_path = "/".join(feat_info_dataset._hdf5_path.split("/")[-2:])
                    feat_info_dataset.write_simple_attributes(group_object[relative_hdf5_dataset_path])


    @property
    def feature_information(self) -> GroupFBase:
        return self._attributes[self.__feature_information_hdf_name__]

    @feature_information.setter
    def feature_information(self, val: GroupFBase):
        self._attributes[self.__feature_information_hdf_name__] = val

    @property
    def __feature_information_type__(self):
        print("Supposed to override feature_information before using")
        return float
        raise NotImplementedError()

    def feature_information_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_information = self.__feature_information_type__()

    @property
    def product_specification(self) -> str:
        return self._attributes[self.__product_specification_hdf_name__]

    @product_specification.setter
    def product_specification(self, val: str):
        self._attributes[self.__product_specification_hdf_name__] = val

    @property
    def __product_specification_type__(self) -> Type[str]:
        return str

    def product_specification_create(self):
        """ Creates a blank, empty or zero value for product_specification"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.product_specification = self.__product_specification_type__()

    @property
    def issue_time(self) -> datetime.time:
        return self._attributes[self.__issue_time_hdf_name__]

    @issue_time.setter
    def issue_time(self, val: Union[datetime.time, datetime.datetime, str]):
        self.set_datetime_attribute(val, self.__issue_time_hdf_name__, self.__issue_time_type__)

    @property
    def __issue_time_type__(self) -> Type[datetime.time]:
        return datetime.time

    def issue_time_create(self):
        """ Creates a blank, empty or zero value for issue_time"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.issue_time = self.__issue_time_type__(0, 0)  # midnight

    @property
    def issue_date(self) -> datetime.date:
        return self._attributes[self.__issue_date_hdf_name__]

    @issue_date.setter
    def issue_date(self, val: Union[datetime.date, datetime.datetime, str]):
        self.set_datetime_attribute(val, self.__issue_date_hdf_name__, self.__issue_date_type__)

    @property
    def __issue_date_type__(self) -> Type[datetime.date]:
        return datetime.date

    def issue_date_create(self):
        """ Creates a blank, empty or zero value for issue_date"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.issue_date = self.__issue_date_type__(1970, 1, 1)

    @property
    def epoch(self) -> str:
        return self._attributes[self.__epoch_hdf_name__]

    @epoch.setter
    def epoch(self, val: str):
        self._attributes[self.__epoch_hdf_name__] = val

    @property
    def __epoch_type__(self) -> Type[str]:
        return str

    def epoch_create(self):
        """ Creates a blank, empty or zero value for epoch"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.epoch = self.__epoch_type__()

    @property
    def geographic_identifier(self) -> str:
        return self._attributes[self.__geographic_identifier_hdf_name__]

    @geographic_identifier.setter
    def geographic_identifier(self, val: str):
        self._attributes[self.__geographic_identifier_hdf_name__] = val

    @property
    def __geographic_identifier_type__(self) -> Type[str]:
        return str

    def geographic_identifier_create(self):
        """ Creates a blank, empty or zero value for geographic_identifier"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.geographic_identifier = self.__geographic_identifier_type__()

    @property
    def metadata(self) -> str:
        return self._attributes[self.__metadata_hdf_name__]

    @metadata.setter
    def metadata(self, val: str):
        self._attributes[self.__metadata_hdf_name__] = val

    @property
    def __metadata_type__(self) -> Type[str]:
        return str

    def metadata_create(self):
        """ Creates a blank, empty or zero value for metadata"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.metadata = self.__metadata_type__()

    @property
    def vertical_datum(self) -> Enum:
        return self._attributes[self.__vertical_datum_hdf_name__]

    @vertical_datum.setter
    def vertical_datum(self, val: Union[int, str, VERTICAL_DATUM]):
        self.set_enum_attribute(val, self.__vertical_datum_hdf_name__, self.__vertical_datum_type__)
        # if isinstance(val, str):
        #     val = self.__vertical_datum_type__[val]
        # if isinstance(val , int):
        #     val = self.__vertical_datum_type__(val)
        # self._attributes[self.__vertical_datum_hdf_name__] = val

    @property
    def __vertical_datum_type__(self) -> Type[int]:
        return numpy.int32

    def vertical_datum_create(self):
        """ Creates a blank, empty or zero value for vertical_datum"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum = VERTICAL_DATUM["MLLW"]

    @property
    def meta_features(self) -> str:
        return self._attributes[self.__meta_features_hdf_name__]

    @meta_features.setter
    def meta_features(self, val: str):
        self._attributes[self.__meta_features_hdf_name__] = val

    @property
    def __meta_features_type__(self) -> Type[str]:
        return str

    def meta_features_create(self):
        """ Creates a blank, empty or zero value for meta_features"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.meta_features = self.__meta_features_type__()

    @property
    def name_of_horizontal_crs(self) -> str:
        return self._attributes[self.__name_of_horizontal_crs_hdf_name__]

    @name_of_horizontal_crs.setter
    def name_of_horizontal_crs(self, val: str):
        self._attributes[self.__name_of_horizontal_crs_hdf_name__] = val

    @property
    def __name_of_horizontal_crs_type__(self) -> Type[str]:
        return str

    def name_of_horizontal_crs_create(self):
        """ Creates a blank, empty or zero value for name_of_horizontal_crs"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.name_of_horizontal_crs = self.__name_of_horizontal_crs_type__()

    @property
    def type_of_horizontal_crs(self) -> TYPE_OF_HORIZONTAL_CRS:
        return self._attributes[self.__type_of_horizontal_crs_hdf_name__]

    @type_of_horizontal_crs.setter
    def type_of_horizontal_crs(self, val: Union[int, str, TYPE_OF_HORIZONTAL_CRS]):
        self.set_enum_attribute(val, self.__type_of_horizontal_crs_hdf_name__, self.__type_of_horizontal_crs_type__)

    @property
    def __type_of_horizontal_crs_type__(self) -> Type[TYPE_OF_HORIZONTAL_CRS]:
        return TYPE_OF_HORIZONTAL_CRS

    def type_of_horizontal_crs_create(self):
        """ Creates a blank, empty or zero value for type_of_horizontal_crs
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.type_of_horizontal_crs = list(self.__type_of_horizontal_crs_type__)[0]

    @property
    def horizontal_crs(self) -> int:
        return self._attributes[self.__horizontal_crs_hdf_name__]

    @horizontal_crs.setter
    def horizontal_crs(self, val: int):
        self._attributes[self.__horizontal_crs_hdf_name__] = val

    @property
    def __horizontal_crs_type__(self) -> Type[int]:
        return numpy.int32

    def horizontal_crs_create(self):
        """ Creates a blank, empty or zero value for horizontal_crs"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_crs = self.__horizontal_crs_type__()

    # Even though the horizontal_cs is an enum, it is stored as an int per the spec in table 10c-6
    @property
    def horizontal_cs(self) -> int:
        return self._attributes[self.__horizontal_cs_hdf_name__]

    @horizontal_cs.setter
    def horizontal_cs(self, val: Union[int, str, HORIZONTAL_CS]):
        # Use an enumeration to control the values, though not officially an enum it essentially is, throws an exception if not in the enum
        self.set_enum_attribute(val, self.__horizontal_cs_hdf_name__, HORIZONTAL_CS)
        # convert the legal values into an integer
        self._attributes[self.__horizontal_cs_hdf_name__] = self._attributes[self.__horizontal_cs_hdf_name__].value

    @property
    def __horizontal_cs_type__(self) -> Type[int]:
        return numpy.int32

    def horizontal_cs_create(self):
        """ Creates a blank, empty or zero value for horizontal_cs
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_cs = list(HORIZONTAL_CS)[0].value  # this is limited to an enumerated list but they specify it as an int

    @property
    def horizontal_datum(self) -> int:
        return self._attributes[self.__horizontal_datum_hdf_name__]

    @horizontal_datum.setter
    def horizontal_datum(self, val: int):
        self._attributes[self.__horizontal_datum_hdf_name__] = val

    @property
    def __horizontal_datum_type__(self) -> Type[int]:
        return numpy.int32

    def horizontal_datum_create(self):
        """ Creates a blank, empty or zero value for horizontal_datum"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_datum = self.__horizontal_datum_type__()

    @property
    def name_of_horizontal_datum(self) -> str:
        return self._attributes[self.__name_of_horizontal_datum_hdf_name__]

    @name_of_horizontal_datum.setter
    def name_of_horizontal_datum(self, val: str):
        self._attributes[self.__name_of_horizontal_datum_hdf_name__] = val

    @property
    def __name_of_horizontal_datum_type__(self) -> Type[str]:
        return str

    def name_of_horizontal_datum_create(self):
        """ Creates a blank, empty or zero value for name_of_horizontal_datum"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.name_of_horizontal_datum = self.__name_of_horizontal_datum_type__()

    @property
    def prime_meridian(self) -> int:
        return self._attributes[self.__prime_meridian_hdf_name__]

    @prime_meridian.setter
    def prime_meridian(self, val: int):
        self._attributes[self.__prime_meridian_hdf_name__] = val

    @property
    def __prime_meridian_type__(self) -> Type[int]:
        return numpy.int32

    def prime_meridian_create(self):
        """ Creates a blank, empty or zero value for prime_meridian"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.prime_meridian = self.__prime_meridian_type__()

    """NOTE: All latitudes and longitudes of the projection parameters must be given in degrees (south and
west negative). Azimuths are given in degrees. For detailed description of the projection method refer to
the EPSG documentation."""

    @property
    def spheriod(self) -> int:
        return self._attributes[self.__spheriod_hdf_name__]

    @spheriod.setter
    def spheriod(self, val: int):
        self._attributes[self.__spheriod_hdf_name__] = val

    @property
    def __spheriod_type__(self) -> Type[int]:
        return numpy.int32

    def spheriod_create(self):
        """ Creates a blank, empty or zero value for spheriod"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.spheriod = self.__spheriod_type__()

    @property
    def projection_method(self) -> PROJECTION_METHOD:
        return self._attributes[self.__projection_method_hdf_name__]

    @projection_method.setter
    def projection_method(self, val: Union[int, str, PROJECTION_METHOD]):
        # Use an enumeration to control the values, though not officially an enum it essentially is
        self.set_enum_attribute(val, self.__projection_method_hdf_name__, PROJECTION_METHOD)
        # convert the legal values into an integer
        self._attributes[self.__projection_method_hdf_name__] = self._attributes[self.__projection_method_hdf_name__].value

    @property
    def __projection_method_type__(self) -> Type[PROJECTION_METHOD]:
        return numpy.int32  # PROJECTION_METHOD

    def projection_method_create(self):
        """ Creates a blank, empty or zero value for projection_method
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_method = list(PROJECTION_METHOD)[0]

    @property
    def projection_parameter_1(self) -> float:
        return self._attributes[self.__projection_parameter_1_hdf_name__]

    @projection_parameter_1.setter
    def projection_parameter_1(self, val: float):
        self._attributes[self.__projection_parameter_1_hdf_name__] = val

    @property
    def __projection_parameter_1_type__(self) -> Type[float]:
        return numpy.float64

    def projection_parameter_1_create(self):
        """ Creates a blank, empty or zero value for projection_parameter_1"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_parameter_1 = self.__projection_parameter_1_type__()

    @property
    def projection_parameter_2(self) -> float:
        return self._attributes[self.__projection_parameter_2_hdf_name__]

    @projection_parameter_2.setter
    def projection_parameter_2(self, val: float):
        self._attributes[self.__projection_parameter_2_hdf_name__] = val

    @property
    def __projection_parameter_2_type__(self) -> Type[float]:
        return numpy.float64

    def projection_parameter_2_create(self):
        """ Creates a blank, empty or zero value for projection_parameter_2"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_parameter_2 = self.__projection_parameter_2_type__()

    @property
    def projection_parameter_3(self) -> float:
        return self._attributes[self.__projection_parameter_3_hdf_name__]

    @projection_parameter_3.setter
    def projection_parameter_3(self, val: float):
        self._attributes[self.__projection_parameter_3_hdf_name__] = val

    @property
    def __projection_parameter_3_type__(self) -> Type[float]:
        return numpy.float64

    def projection_parameter_3_create(self):
        """ Creates a blank, empty or zero value for projection_parameter_3"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_parameter_3 = self.__projection_parameter_3_type__()

    @property
    def projection_parameter_4(self) -> float:
        return self._attributes[self.__projection_parameter_4_hdf_name__]

    @projection_parameter_4.setter
    def projection_parameter_4(self, val: float):
        self._attributes[self.__projection_parameter_4_hdf_name__] = val

    @property
    def __projection_parameter_4_type__(self) -> Type[float]:
        return numpy.float64

    def projection_parameter_4_create(self):
        """ Creates a blank, empty or zero value for projection_parameter_4"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_parameter_4 = self.__projection_parameter_4_type__()

    @property
    def projection_parameter_5(self) -> float:
        return self._attributes[self.__projection_parameter_5_hdf_name__]

    @projection_parameter_5.setter
    def projection_parameter_5(self, val: float):
        self._attributes[self.__projection_parameter_5_hdf_name__] = val

    @property
    def __projection_parameter_5_type__(self) -> Type[float]:
        return numpy.float64

    def projection_parameter_5_create(self):
        """ Creates a blank, empty or zero value for projection_parameter_5"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_parameter_5 = self.__projection_parameter_5_type__()

    @property
    def false_northing(self) -> float:
        return self._attributes[self.__false_northing_hdf_name__]

    @false_northing.setter
    def false_northing(self, val: float):
        self._attributes[self.__false_northing_hdf_name__] = val

    @property
    def __false_northing_type__(self) -> Type[float]:
        return numpy.float64

    def false_northing_create(self):
        """ Creates a blank, empty or zero value for false_northing"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.false_northing = self.__false_northing_type__()

    @property
    def false_easting(self) -> float:
        return self._attributes[self.__false_easting_hdf_name__]

    @false_easting.setter
    def false_easting(self, val: float):
        self._attributes[self.__false_easting_hdf_name__] = val

    @property
    def __false_easting_type__(self) -> Type[float]:
        return numpy.float64

    def false_easting_create(self):
        """ Creates a blank, empty or zero value for false_easting"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.false_easting = self.__false_easting_type__()

    @property
    def vertical_cs(self) -> VERTICAL_CS:
        return self._attributes[self.__vertical_cs_hdf_name__]

    @vertical_cs.setter
    def vertical_cs(self, val: Union[int, str, VERTICAL_CS]):
        # Use an enumeration to control the values, though not officially an enum it essentially is
        self.set_enum_attribute(val, self.__vertical_cs_hdf_name__, VERTICAL_CS)
        # convert the legal values into an integer
        self._attributes[self.__vertical_cs_hdf_name__] = self._attributes[self.__vertical_cs_hdf_name__].value

    @property
    def __vertical_cs_type__(self) -> Type[VERTICAL_CS]:
        return numpy.int32  # VERTICAL_CS

    def vertical_cs_create(self):
        """ Creates a blank, empty or zero value for vertical_cs
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_cs = list(VERTICAL_CS)[0]

    @property
    def vertical_coordinate_base(self) -> VERTICAL_COORDINATE_BASE:
        return self._attributes[self.__vertical_coordinate_base_hdf_name__]

    @vertical_coordinate_base.setter
    def vertical_coordinate_base(self, val: Union[int, str, VERTICAL_COORDINATE_BASE]):
        self.set_enum_attribute(val, self.__vertical_coordinate_base_hdf_name__, self.__vertical_coordinate_base_type__)

    @property
    def __vertical_coordinate_base_type__(self) -> Type[VERTICAL_COORDINATE_BASE]:
        return VERTICAL_COORDINATE_BASE

    def vertical_coordinate_base_create(self):
        """ Creates a blank, empty or zero value for vertical_coordinate_base
        """
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
        return VERTICAL_DATUM_REFERENCE

    def vertical_datum_reference_create(self):
        """ Creates a blank, empty or zero value for vertical_datum_reference
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum_reference = list(VERTICAL_DATUM_REFERENCE)[0]

    @property
    def vertical_datum(self) -> int:
        return self._attributes[self.__vertical_datum_hdf_name__]

    @vertical_datum.setter
    def vertical_datum(self, val: (int, str, VERTICAL_DATUM)):
        # NOTE: When reading from a file h5py gets attributes alphabetically so we can't rely on vertical_datum_reference being set before vertical_datum
        # verticalDatumReference when 1 is from the enumeration but when 2 is an EPSG code that we don't have a check for.
        try:
            self.vertical_datum_reference  # see if the attribute exists
        except KeyError:
            pass
        else:
            if self.vertical_datum_reference == VERTICAL_DATUM_REFERENCE(1):
                try:
                    self.set_enum_attribute(val, self.__vertical_datum_hdf_name__, VERTICAL_DATUM)
                except S100Exception as e:
                    raise S100Exception(f"When vertical_datum_reference is '1' then vertical_datum must be a value given in the enumeration {VERTICAL_DATUM}, the supplied {val} was not found")
                # convert the enumeration back to an integer
                val = self._attributes[self.__vertical_datum_hdf_name__].value
        self._attributes[self.__vertical_datum_hdf_name__] = val

    @property
    def __vertical_datum_type__(self) -> Type[int]:
        return int

    def vertical_datum_create(self):
        """ Creates a blank, empty or zero value for vertical_datum"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum = self.__vertical_datum_type__()


class S104Root(S100Root):
    """The root group contains a feature information group and N feature containers.
    In S104 there is one feature container 'water level'.
    The coverage names are determined from the matching CoveragesAttributes
    Table 3 of v0.0.7
    """
    __water_level_hdf_name__ = WATER_LEVEL
    __water_level_trend_threshold_hdf_name__ = "waterLevelTrendThreshold"

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


class S104File(S100File):
    """ HDF5 file object"""
    PRODUCT_SPECIFICATION = 'INT.IHO.S-104.1.0'

    @staticmethod
    def make_container_for_dcf(data_coding_format):
        if data_coding_format == 2:
            container = WaterLevelContainerDCF2()
        elif data_coding_format == 3:
            container = WaterLevelContainerDCF3()
        else:
            raise S104Exception(f"DCF {data_coding_format} not supported")
        return container

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

