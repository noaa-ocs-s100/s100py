try:
    from ... import s1xx
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py.s100"

from ..v4_0 import s100 as v4_0
# Anything not overridden in this module will use whatever was available in the previous version
from ..v4_0.s100 import *

EDITION = 5.0

# Table 10c-4 and Table 10c-23
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
                              ]
                              )
# Table 10c-25
VERTICAL_COORDINATE_BASE = Enum(value="VERTICAL_COORDINATE_BASE",

# Table 10c-26
VERTICAL_DATUM_REFERENCE = Enum(value="VERTICAL_DATUM_REFERENCE",

# Table 10c-27
PROJECTION_METHODS = Enum(value="PROJECTION_METHODS",

class S100Root(v4_0.S100Root):

    __horizontal_crs_hdf_name__ = "horizontalCRS"
    __name_of_horizontal_crs_hdf_name__ = "nameOfHorizontalCRS"  #: HDF5 naming
    
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
        
        
    __type_of_horizontal_crs_hdf_name__ = "typeOfHorizontalCRS"  #: HDF5 naming
    
    @property
    def type_of_horizontal_crs(self) -> TYPE_OF_HORIZONTAL_CRS:
        """Codes for describing the type of the two-dimensional horizontal CRS
        see Table 10c-24 or :data:`~TYPE_OF_HORIZONTAL_CRS`
        """
        return self._attributes[self.__type_of_horizontal_crs_hdf_name__]
    
    @type_of_horizontal_crs.setter
    def type_of_horizontal_crs(self, val: ):
        self._attributes[self.__type_of_horizontal_crs_hdf_name__] = val
    
    @property
    def __type_of_horizontal_crs_type__(self) -> Type[TYPE_OF_HORIZONTAL_CRS]:
        return TYPE_OF_HORIZONTAL_CRS
    
    def type_of_horizontal_crs_create(self):
        """ Creates a blank, empty or zero value for type_of_horizontal_crs"""
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.type_of_horizontal_crs = self.__type_of_horizontal_crs_type__()
    
    
    horizontalCS
    horizontalDatum
    nameOfHorizontalDatum
    primeMeridian
    """NOTE: All latitudes and longitudes of the projection parameters must be given in degrees (south and
west negative). Azimuths are given in degrees. For detailed description of the projection method refer to
the EPSG documentation."""
    projectionMethod
    projectionParameter1
    projectionParameter2
    projectionParameter3
    projectionParameter4
    projectionParameter5
    falseNorthing
    falseEasting


class S100File(v4_0.S100File):
    PRODUCT_SPECIFICATION = 'INT.IHO.S-100.5.0'

    def __init__(self, *args, **kywrds):
        if 'root' not in kywrds:
            kywrds['root'] = S100Root  # inherited classes will specify their own root type
        super().__init__(*args, **kywrds)
