from itertools import chain
try:
    from ... import s1xx
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py.s100"

from ..v4_0 import s100 as v4_0
# Anything not overridden in this module will use whatever was available in the previous version
from ..v4_0.s100 import *

EDITION = 5.0
# @TODO add support for Positioning group (S100 v5.0 10c-9.10.1, 10c-9.10.2 and table 10c-16) for DataCodingFormats 1,3,4,7,8
#   Basically makes a 1D dataset array of the positions and attributes (like Z).
#   For TINs  DCF=7
#   Optionally has triangles dataset NumTriangles x 3 -- array of the 3 node indices for positions
#   Optionally has adjacency dataset NumTriangles x3 -- array of the indices or the triangles that share an edge with each triangle side
#   adjacency[i][0] = triangle adjacent to the edge specified by triangles[i][0] & triangles[i][1]
#   adjacency[i][1] = triangle adjacent to edge triangles[i][1] & triangles[i][2]
#   adjacency[i][2] = triangle adjacent to edge triangles[i][2] & triangles[i][0]
#   Elements for edges without adjacent triangles are filled with the value -1
#   Applicable only for dataCodingFormat = 7 (TIN), but optional even for TIN.

# Python doesn't allow extending an enum, so we have to create a new one with the additions
# TODO does this work (show all values) in sphinx?
class __VERTICAL_DATUM_additions(Enum):
    """ S100 v5 Part 17 Vertical and Sounding Datum
    Added balticSeaChartDatum2000 = 44
    """
    balticSeaChartDatum2000 = 44
VERTICAL_DATUM = Enum('VERTICAL_DATUM', [(i.name, i.value) for i in chain(v4_0.VERTICAL_DATUM, __VERTICAL_DATUM_additions)])

# Table 10c-4 and Table 10c-23, added values 8 and 9 for Fixed stations - stationwise and Feature oriented regular grid
# @TODO if chaining works above for vertical datum then use for Data Coding Format too
DATA_CODING_FORMAT = Enum(value="DATA_CODING_FORMAT",
                          names=[
                              ('Time series at fixed stations', 1),
                              ('Regularly-gridded arrays', 2),
                              ('Ungeorectified gridded arrays', 3),
                              ('Moving platform', 4),
                              ('Irregular grid', 5),
                              ('Variable cell size', 6),
                              ('TIN', 7),
                              ('Fixed stations stationwise', 8),
                              ('Feature oriented regular grid', 9),
                              # alternate shortcut names that also show up in sphinx, these will be stored with full names including spaces
                              ('TIME', 1),
                              ('REGULAR', 2),
                              ('UNGEORECTIFIED', 3),
                              ('MOVING', 4),
                              ('IRREGULAR', 5),
                              ('VARIABLE', 6),
                              ('STATIONWISE', 8),
                              ('FEATURE_REGULAR', 9)
                          ]
                          )
"""
Sphinx is not interpreting the enum names properly when there are spaces. The correct enum names with spaces are::

  ('Time series at fixed stations', 1),
  ('Regularly-gridded arrays', 2),
  ('Ungeorectified gridded arrays', 3),
  ('Moving platform', 4),
  ('Irregular grid', 5),
  ('Variable cell size', 6),
  ('TIN', 7),
  ('Fixed stations - stationwise', 8),
  ('Feature oriented regular grid', 9),
"""

# Table 10c-24
TYPE_OF_HORIZONTAL_CRS = Enum(value="TYPE_OF_HORIZONTAL_CRS",
                                names=[
                                    ('geodeticCRS2D', 1),
                                    ('projectedCRS', 2),
                              ])
# Table 10c-25
VERTICAL_COORDINATE_BASE = Enum(value="VERTICAL_COORDINATE_BASE",
                                names=[
                                    ('seaSurface', 1),
                                    ('verticalDatum', 2),
                                    ('seaBottom', 3),
                              ])

# Table 10c-26
VERTICAL_DATUM_REFERENCE = Enum(value="VERTICAL_DATUM_REFERENCE",
                                names=[
                                    ('s100VerticalDatum', 1),
                                    ('EPSG', 2),
                                ])

# Table 10c-27
PROJECTION_METHOD = Enum(value="PROJECTION_METHODS",
                            names=[
                                ('Mercator', 9805),
                                ('Transverse Mercator', 9807),
                                ('Oblique Mercator', 9815),
                                ('Hotine Oblique Mercator', 9812),
                                ('Lambert Conic Conformal 1SP', 9801),
                                ('Lambert Conic Conformal 2SP', 9802),
                                ('Oblique Stereographic', 9809),
                                ('Polar Stereographic', 9810),
                                ('Krovak Oblique Conic Conformal', 9819),
                                ('American Polyconic', 9818),
                                ('Albers Equal Area', 9822),
                                ('Lambert Azimuthal Equal Area', 9820),
                          ])

# Table 10c-6
HORIZONTAL_CS = Enum(value="HORIZONTAL_CS",
                        names=[
                            ('Lat Lon', 6422),
                            ('Easting, Northing', 4400),
                            ('Northing, Easting', 4500),
                     ])

#Table 10c-6
VERTICAL_CS = Enum(value="VERTICAL_CS",
                        names=[
                            ('Depth', 6498),
                            ('Height', 6499),
                     ])

DATA_OFFSET_CODE = Enum(value="DATA_OFFSET_CODE",
                        names=[
                            ('XMin, YMin', 1),
                            ('XMax, YMax', 2),
                            ('XMax, YMin', 3),
                            ('XMin, YMax', 4),
                            ('Barycenter', 5),
                            ('Lower left', 1),
                            ('Upper right', 2),
                            ('Lower right', 3),
                            ('Upper left', 4),
                            ('Centroid', 5),
                        ])



class FeatureContainer(v4_0.FeatureContainer):
    """ Feature Container from S100 v5.0 Table 10c-9
    The feature attribute table was added, primarily for S102 to supply extended information about the origin of the data.
    The Data Coding Format must be 9 to use the feature attribute table.
    The Columns will be "id" followed by an arbitrary number of column names.
    Each "id" is an unsigned integer that matches the values in the values record of a data values group.

    Additional columns
    The number of additional columns depends on the information to be transported.  The number must always be greater than zero (>0).
    Table 10c-8 defines that the name of the "feature container group" must always correspond to the
    camelCase name of an object of the type "Feature Type" of the IHO Geospatial Information (GI) Registry, Data Dictionary Register.
    The names of the other columns are based on the attributes of the selected feature type of the Data Dictionary Register.
    Due to the mapping form of the attributes in the HDF5 format, only simple data types (see clause 10c-7) can be used.
    For the mapping of complex data types, a corresponding notation by means of a dot is necessary.
    The respective levels of the complex data type are separated by a dot up to the simple data type;
    for example, "surveyDateRange.dateEnd" & "surveyDateRange.dateStart".
    """
    __feature_attribute_table_hdf_name__ = "featureAttributeTable"  #: HDF5 naming

    @property
    def feature_attribute_table(self) -> s1xx_sequence:
        return self._attributes[self.__feature_attribute_table_hdf_name__]

    @feature_attribute_table.setter
    def feature_attribute_table(self, val: s1xx_sequence):
        self._attributes[self.__feature_attribute_table_hdf_name__] = val

    @property
    def __feature_attribute_table_type__(self) -> Type[numpy.array]:
        return numpy.array

    def feature_attribute_table_create(self):
        """ Creates a blank, empty or zero value for feature_attribute_table"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        # FIXME @TODO -- this has to be an arbitrary number of strings, is this limiting it to two?
        self.feature_attribute_table = self.__feature_attribute_table_type__(["", ""], dtype=h5py_string_dtype)

class DataOffset:
    """ Mixin class for Data Offset from S100 v5.0 Table 10c-10
    Added dataOffsetCode and dataOffsetVector to the featureContainer for DataCodingFormats 2,5,6,9.
    """
    __data_offset_code_hdf_name__ = "dataOffsetCode"  #: HDF5 naming
    __data_offset_vector_hdf_name__ = "dataOffsetVector"  #: HDF5 naming

    @property
    def data_offset_code(self) -> DATA_OFFSET_CODE:
        return self._attributes[self.__data_offset_code_hdf_name__]

    @data_offset_code.setter
    def data_offset_code(self, val: Union[int, str, DATA_OFFSET_CODE]):
        self.set_enum_attribute(val, self.__data_offset_code_hdf_name__, self.__data_offset_code_type__)

    @property
    def __data_offset_code_type__(self) -> Type[DATA_OFFSET_CODE]:
        return DATA_OFFSET_CODE

    def data_offset_code_create(self):
        """ Creates a blank, empty or zero value for data_offset_code
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_offset_code = list(self.data_offset_code_type)[0]

    # FIXME - this is an array not just a float
    @property
    def data_offset_vector(self) -> s1xx_sequence:
        return self._attributes[self.__data_offset_vector_hdf_name__]

    @data_offset_vector.setter
    def data_offset_vector(self, val: s1xx_sequence):
        self._attributes[self.__data_offset_vector_hdf_name__] = val

    @property
    def __data_offset_vector_type__(self) -> Type[numpy.array]:
        return numpy.array

    def data_offset_vector_create(self):
        """ Creates a blank, empty or zero value for data_offset_vector"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.data_offset_vector = self.__data_offset_vector_type__()


class FeatureInstanceDCF8(NumberOfStations, FeatureInstanceBase):
    """ Fixed stations - stationwise from S100 v5.0 Table 10c-12"""

class FeatureInstanceDCF9(FeatureInstanceDCF2):
    pass

class FeatureContainerDCF9(DataOffset, SequencingRule, FeatureContainer):
    """ Feature Oriented Regular Grid from S100 v5.0 Table 10c-12
    Container for Data Coding Format 9 """

    def data_coding_format_create(self):
        self.data_coding_format = self.__data_coding_format_type__(9)

class FeatureContainerDCF2(DataOffset, v4_0.FeatureContainerDCF2):
    pass


class FeatureContainerDCF5(DataOffset, v4_0.FeatureContainerDCF5):
    pass


class FeatureContainerDCF6(DataOffset, v4_0.FeatureContainerDCF6):
    pass


# FIXME @TODO - this FeatureInformationDataset class in v4.0 writes the compound array to the HDF5 file.
#  Do we need a new Compound type that does this instead (we didn't implement any other compounds in v4.0)
#  The FeatureContainer's new feature_attribute_table needs to be stored similarly.  Also be able to store an array of compounds
# class FeatureInformationDataset(S1xxDatasetBase, ABC):  # Chunking
#     """ This class comes from S100 -- 10c-9.5 Feature information group.
#     This class serves to keep a list of FeatureInformation objects which will be turned into a compound array
#     of strings in the HDF5 file.
#
#     The metadata_name property must be overridden.
#     The metadata_type will likely be overridden with a specific subclass for the s100+ spec
#     """
#
#     @property
#     def metadata_type(self) -> Type[FeatureInformation]:
#         return FeatureInformation


class S100Root(v4_0.S100Root):
    """
    The root of the S100 v5.0 schema.  There are restrictions on many of the CRS attributes based on other attributes.
    For example, the horizontal_cs has different value options based on horizontal_crs.
    """
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
        self.type_of_horizontal_crs = list(self.type_of_horizontal_crs_type)[0]

    @property
    def horizontal_crs(self) -> int:
        return self._attributes[self.__horizontal_crs_hdf_name__]

    @horizontal_crs.setter
    def horizontal_crs(self, val: int):
        self._attributes[self.__horizontal_crs_hdf_name__] = val

    @property
    def __horizontal_crs_type__(self) -> Type[int]:
        return int

    def horizontal_crs_create(self):
        """ Creates a blank, empty or zero value for horizontal_crs"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_crs = self.__horizontal_crs_type__()

    @property
    def horizontal_cs(self) -> HORIZONTAL_CS:
        return self._attributes[self.__horizontal_cs_hdf_name__]

    @horizontal_cs.setter
    def horizontal_cs(self, val: Union[int, str, HORIZONTAL_CS]):
        self.set_enum_attribute(val, self.__horizontal_cs_hdf_name__, self.__horizontal_cs_type__)

    @property
    def __horizontal_cs_type__(self) -> Type[HORIZONTAL_CS]:
        return HORIZONTAL_CS

    def horizontal_cs_create(self):
        """ Creates a blank, empty or zero value for horizontal_cs
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.horizontal_cs = list(self.horizontal_cs_type)[0]

    @property
    def horizontal_datum(self) -> int:
        return self._attributes[self.__horizontal_datum_hdf_name__]

    @horizontal_datum.setter
    def horizontal_datum(self, val: int):
        self._attributes[self.__horizontal_datum_hdf_name__] = val

    @property
    def __horizontal_datum_type__(self) -> Type[int]:
        return int

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
        return int

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
        return int

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
        self.set_enum_attribute(val, self.__projection_method_hdf_name__, self.__projection_method_type__)

    @property
    def __projection_method_type__(self) -> Type[PROJECTION_METHOD]:
        return PROJECTION_METHOD

    def projection_method_create(self):
        """ Creates a blank, empty or zero value for projection_method
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.projection_method = list(self.projection_method_type)[0]

    @property
    def projection_parameter_1(self) -> float:
        return self._attributes[self.__projection_parameter_1_hdf_name__]

    @projection_parameter_1.setter
    def projection_parameter_1(self, val: float):
        self._attributes[self.__projection_parameter_1_hdf_name__] = val

    @property
    def __projection_parameter_1_type__(self) -> Type[float]:
        return float

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
        return float

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
        return float

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
        return float

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
        return float

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
        return float

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
        return float

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
        self.set_enum_attribute(val, self.__vertical_cs_hdf_name__, self.__vertical_cs_type__)

    @property
    def __vertical_cs_type__(self) -> Type[VERTICAL_CS]:
        return VERTICAL_CS

    def vertical_cs_create(self):
        """ Creates a blank, empty or zero value for vertical_cs
        """
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_cs = list(self.vertical_cs_type)[0]


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
        self.vertical_coordinate_base = list(self.vertical_coordinate_base_type)[0]

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
        self.vertical_datum_reference = list(self.vertical_datum_reference_type)[0]

    @property
    def vertical_datum(self) -> int:
        return self._attributes[self.__vertical_datum_hdf_name__]

    @vertical_datum.setter
    def vertical_datum(self, val: int):
        self._attributes[self.__vertical_datum_hdf_name__] = val

    @property
    def __vertical_datum_type__(self) -> Type[int]:
        return int

    def vertical_datum_create(self):
        """ Creates a blank, empty or zero value for vertical_datum"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.vertical_datum = self.__vertical_datum_type__()


class FeatureContainerDCF8(FeatureContainer):
    """ Container for Data Coding Format 8 """

    def data_coding_format_create(self):
        self.data_coding_format = self.__data_coding_format_type__(8)



class S100File(v4_0.S100File):
    PRODUCT_SPECIFICATION = 'INT.IHO.S-100.5.0'

    def __init__(self, *args, **kywrds):
        if 'root' not in kywrds:
            kywrds['root'] = S100Root  # inherited classes will specify their own root type
        super().__init__(*args, **kywrds)
