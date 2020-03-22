"""

An Example API
===============

.. _@property:  https://docs.python.org/3/library/functions.html#property>
.. _Pycharm Live Template: https://www.jetbrains.com/help/pycharm/tutorial-creating-and-applying-live-templates-code-snippets.html

Let's create a new sample api to show the basic pieces of the s100py framework.  We'll lay out an example
data spec and then create the code and classes to make it work.  The basic data types are:

    - basic attribute (string, int, float, enum)
    - named dataset
    - grid dataset
    - named sub group
    - variable number of occurrence subgroup (e.g. Group_01, Group_02...)
    - file

Here is our pretend data spec::

    S999
    |
    + datasetWithNames (type=named dataset)
    | |
    | + attrInt (type=int)
    | |
    | + attrFloat (type=float)
    | |
    | + attrStr (type=string)
    |
    + dataGroup_01, dataGroup_02...  (type=multi-occurrence group)
    | |
    | + nameOfData (type=enumeration)
    | |
    | + dataGrid (type=rectangular grid dataset)
    |
    + myLocationGroup (type=sub group but based on a S100 type, GeographicBoundingBox)
    | |
    | + utmZone (type=int)
    |
    + myFirstObject (type=sub group)
    | |
    | + dataValue (type=string)

Each piece of data from the S100 spec there should be a python friendly (pep8 - "snake case") name to hold the data.
Also there is a read only property that contains the S100 name (which may be camel case and not as human friendly).
There are also two helper functions to help discover the type of the data that should be stored and one to create
an instance of the data with an appropriate default.
"""
# builtings
from enum import Enum
from typing import Callable, Iterator, Union, Optional, List, Type
# 3rd party
import numpy
# custom modules
from s100py import s1xx, s100


class MONTY(Enum):
    spam = 1
    cheese = 2


class MyObject(s1xx.S1XX_Attributes_base):
    @property
    def __version__(self) -> int:
        return 1

    @property
    def data_value_attribute_name(self) -> str:
        return "dataValue"

    @property
    def data_value(self) -> str:
        return self._attributes[self.data_value_attribute_name]

    @data_value.setter
    def data_value(self, val: str):
        self._attributes[self.data_value_attribute_name] = val

    @property
    def data_value_type(self) -> str:
        return str

    def data_value_create(self):
        """ Creates a blank, empty or zero value for data_value"""
        self.data_value = self.data_value_type()


class MyLocation(s100.GeographicBoundingBox):
    empty_zone = 999
    @property
    def __version__(self) -> int:
        return 1

    @property
    def utm_zone_attribute_name(self) -> str:
        return "utmZone"

    @property
    def utm_zone(self) -> int:
        return self._attributes[self.utm_zone_attribute_name]

    @utm_zone.setter
    def utm_zone(self, val: int):
        if isinstance(val, str):
            val = int(val)
        if (val <= 0 or val > 60) and val != self.empty_zone:
            raise Exception("Illegal zone number, must be between 1 and 60")
        self._attributes[self.utm_zone_attribute_name] = val

    @property
    def utm_zone_type(self) -> int:
        return int

    def utm_zone_create(self):
        """ Creates a blank, empty or zero value for utm_zone"""
        self.utm_zone = self.utm_zone_type(self.empty_zone)


class DataGroupObject(s1xx.S1XX_Attributes_base):
    @property
    def __version__(self) -> int:
        return 1

    @property
    def name_of_data_attribute_name(self) -> str:
        return "nameOfData"

    @property
    def name_of_data(self) -> MONTY:
        return self._attributes[self.name_of_data_attribute_name]

    @name_of_data.setter
    def name_of_data(self, val: Union[int, str, MONTY]):
        self.set_enum_attribute(val, self.name_of_data_attribute_name, self.name_of_data_type)

    @property
    def name_of_data_type(self) -> Type[Enum]:
        return MONTY

    def name_of_data_create(self):
        """ Creates an enumerated value of 'spam' (because it's first in the list) """
        self.name_of_data = list(self.name_of_data_type)[0]

    @property
    def data_grid_attribute_name(self) -> str:
        return "dataGrid"

    @property
    def data_grid(self) -> s1xx.s1xx_sequence:
        return self._attributes[self.data_grid_attribute_name]

    @data_grid.setter
    def data_grid(self, val: s1xx.s1xx_sequence):
        self._attributes[self.data_grid_attribute_name] = val

    @property
    def data_grid_type(self) -> s1xx.s1xx_sequence:
        return numpy.ndarray

    def data_grid_create(self):
        """ Creates a blank, empty or zero value for data_grid"""
        self.data_grid = self.data_grid_type([2], numpy.float32)


class DataGroups(s1xx.S1XX_MetadataList_base):
    """ This is the list of dataGroup_NNN that are held as a list.
    Each dataGroup_NNN has a data_grid dataset and name_of_data attribute.
    """

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_name(self) -> str:
        return "dataGroup"

    @property
    def metadata_type(self) -> type:
        return DataGroupObject


class datasetWithNames(s1xx.S1XX_Attributes_base):
    def get_write_order(self):
        return ["attrInt", "attrStr", "attrFloat"]

    @property
    def __version__(self) -> int:
        return 1

    @property
    def attr_int_attribute_name(self) -> str:
        return "attrInt"

    @property
    def attr_int(self) -> int:
        return self._attributes[self.attr_int_attribute_name]

    @attr_int.setter
    def attr_int(self, val: int):
        self._attributes[self.attr_int_attribute_name] = val

    @property
    def attr_int_type(self) -> Type[int]:
        return int

    def attr_int_create(self):
        """ Creates a blank, empty or zero value for attr_int"""
        self.attr_int = self.attr_int_type(55)

    @property
    def attr_float_attribute_name(self) -> str:
        return "attrFloat"

    @property
    def attr_float(self) -> float:
        return self._attributes[self.attr_float_attribute_name]

    @attr_float.setter
    def attr_float(self, val: float):
        self._attributes[self.attr_float_attribute_name] = val

    @property
    def attr_float_type(self) -> Type[float]:
        return float

    def attr_float_create(self):
        """ Creates a blank, empty or zero value for attr_float"""
        self.attr_float = self.attr_float_type(123.4)

    @property
    def attr_str_attribute_name(self) -> str:
        return "attrStr"

    @property
    def attr_str(self) -> str:
        return self._attributes[self.attr_str_attribute_name]

    @attr_str.setter
    def attr_str(self, val: str):
        self._attributes[self.attr_str_attribute_name] = val

    @property
    def attr_str_type(self) -> Type[str]:
        return str

    def attr_str_create(self):
        """ Creates a blank, empty or zero value for attr_str"""
        self.attr_str = self.attr_str_type("used a default string")


class DatasetWithNames_List(s1xx.S1XX_Dataset_base):

    @property
    def __version__(self) -> int:
        return 1

    @property
    def metadata_type(self) -> Type[type]:
        return datasetWithNames

    @property
    def metadata_name(self) -> str:
        return "datasetWithNames"


class S999Root(s1xx.S1XX_Attributes_base):
    @property
    def __version__(self) -> int:
        return 1

    @property
    def dataset_with_names_attribute_name(self) -> str:
        return "datasetWithNames"

    @property
    def dataset_with_names(self) -> DatasetWithNames_List:
        return self._attributes[self.dataset_with_names_attribute_name]

    @dataset_with_names.setter
    def dataset_with_names(self, val: DatasetWithNames_List):
        self._attributes[self.dataset_with_names_attribute_name] = val

    @property
    def dataset_with_names_type(self) -> Type[DatasetWithNames_List]:
        return DatasetWithNames_List

    def dataset_with_names_create(self):
        """ Creates a blank, empty or zero value for dataset_with_names"""
        self.dataset_with_names = self.dataset_with_names_type()

    @property
    def data_group_attribute_name(self) -> str:
        return "dataGroup"

    @property
    def data_group(self) -> DataGroups:
        return self._attributes[self.data_group_attribute_name]

    @data_group.setter
    def data_group(self, val: DataGroups):
        self._attributes[self.data_group_attribute_name] = val

    @property
    def data_group_type(self) -> Type[DataGroups]:
        return DataGroups

    def data_group_create(self):
        """ Creates a blank, empty or zero value for data_group"""
        self.data_group = self.data_group_type()

    @property
    def my_location_group_attribute_name(self) -> str:
        return "myLocationGroup"

    @property
    def my_location_group(self) -> MyLocation:
        return self._attributes[self.my_location_group_attribute_name]

    @my_location_group.setter
    def my_location_group(self, val: MyLocation):
        self._attributes[self.my_location_group_attribute_name] = val

    @property
    def my_location_group_type(self) -> Type[MyLocation]:
        return MyLocation

    def my_location_group_create(self):
        """ Creates a blank, empty or zero value for my_location_group"""
        self.my_location_group = self.my_location_group_type()

    @property
    def my_first_object_attribute_name(self) -> str:
        return "myFirstObject"

    @property
    def my_first_object(self) -> MyObject:
        return self._attributes[self.my_first_object_attribute_name]

    @my_first_object.setter
    def my_first_object(self, val: MyObject):
        self._attributes[self.my_first_object_attribute_name] = val

    @property
    def my_first_object_type(self) -> Type[MyObject]:
        return MyObject

    def my_first_object_create(self):
        """ Creates a blank, empty or zero value for my_first_object"""
        self.my_first_object = self.my_first_object_type()


class S999File(s1xx.S1XXFile):
    PRODUCT_SPECIFICATION = numpy.string_('INT.IHO.S-Fake')

    def __init__(self, *args, **kywrds):
        # kywrds['root'] = S999Root
        super().__init__(*args, root=S999Root, **kywrds)

def test_api():
    import tempfile
    import os

    print("running the sample api")
    h, fstr = tempfile.mkstemp(suffix=".h5", prefix="sample_s999", dir=os.path.split(__file__)[0])
    os.close(h)
    write_to_file = S999File(fstr)
    h, fstr_rewrite = tempfile.mkstemp(suffix=".h5", prefix="sample_revised_s999", dir=os.path.split(__file__)[0])
    os.close(h)

    the_first_object = MyObject()
    the_first_object.data_value = "A sample string"
    write_to_file.root.my_first_object = the_first_object

    a_location = MyLocation()
    a_location.initialize_properties(True)
    a_location.utm_zone = 18
    a_location.east_bound_longitude = 33.5
    del a_location.west_bound_longitude
    write_to_file.root.my_location_group = a_location

    write_to_file.root.data_group_create()  # this makes the DataGroups which is a list container for the DataGroupObject
    data_1 = DataGroupObject()
    data_1.name_of_data = data_1.name_of_data_type["spam"]
    data_1.data_grid = numpy.zeros([2, 5])
    write_to_file.root.data_group.append(data_1)
    del data_1.name_of_data

    data_2 = DataGroupObject()
    data_2.name_of_data = MONTY(2)
    data_2.data_grid = numpy.ones([3, 4])
    write_to_file.root.data_group.append(data_2)

    data_3 = write_to_file.root.data_group.append_new_item()
    data_3.data_grid = numpy.arange(0, 10, .75)
    data_3.name_of_data = "cheese"

    write_to_file.root.data_group.append_new_item().name_of_data = 2

    attrib_datasets = DatasetWithNames_List()  # rather than call create we can also make the list independently and set it later
    attr_1 = datasetWithNames()
    attr_1.initialize_properties(True)

    attr_2 = datasetWithNames()
    attr_2.initialize_properties(True)
    attr_2.attr_str = "A custom string this time"
    attr_2.attr_int = 27
    attr_2.attr_float = 35.0
    write_to_file.root.dataset_with_names = DatasetWithNames_List((attr_1, attr_2))

    write_to_file.write()
    write_to_file.close()

    read_from_file = S999File(fstr)
    assert read_from_file.root.dataset_with_names[1].attr_int == 27
    assert read_from_file.root.dataset_with_names[0].attr_str in b"used a default string"
    assert write_to_file.root.data_group[1].name_of_data == MONTY(2)
    assert write_to_file.root.data_group[2].data_grid[1] == 0.75  # the second element of the range
    assert write_to_file.root.my_location_group.east_bound_longitude == 33.5
    try:
        write_to_file.root.my_location_group.west_bound_longitude
    except:
        pass  # all good, the value shouldn't exist
    else:
        assert write_to_file.root.my_location_group.west_bound_longitude is None  # this should not exist, even as None

    copy_of_file = S999File(fstr_rewrite)
    copy_of_file.root = read_from_file.root
    # this shows how to initialize on creation
    copy_of_file.root.my_location_group = MyLocation(utm_zone=22, east_bound_longitude=11, extra_attr="This shouldn't even be here, but it works")
    copy_of_file.write()

if __name__ == "__main__":
    test_api()