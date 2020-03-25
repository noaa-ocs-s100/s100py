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
    + datasetWithNames (type=compound dataset)
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
import tempfile
import os
from enum import Enum
from typing import Callable, Iterator, Union, Optional, List, Type
# 3rd party
import h5py
import numpy
try:
    import pytest
except:
    pass
# custom modules
from s100py import s1xx, s100


class MONTY(Enum):
    spam = 1
    cheese = 2


class MyObject(s1xx.S1XX_Attributes_base):
    data_value_attribute_name = "dataValue"

    @property
    def __version__(self) -> int:
        return 1


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
        """ This will limit the utm zones to 1 thru 60 but also allow for a special 'empty' zone of 999 """
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
    attr_int_attribute_name = "attrInt"
    attr_float_attribute_name = "attrFloat"
    attr_str_attribute_name = "attrStr"

    def get_write_order(self):
        return [self.attr_int_attribute_name, self.attr_str_attribute_name, self.attr_float_attribute_name]

    @property
    def __version__(self) -> int:
        return 1

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
    dataset_with_names_attribute_name = "datasetWithNames"

    @property
    def __version__(self) -> int:
        return 1

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


try:
    pytest
except NameError:
    pass
else:
    @pytest.fixture(scope="module")
    def filename():
        h, fstr = tempfile.mkstemp(".h5", dir=os.path.split(__file__)[0])
        os.close(h)
        yield fstr
        try:
            os.remove(fstr)
        except (FileNotFoundError, PermissionError):
            print("failed to remove "+fstr)

    @pytest.fixture(scope="module")
    def revised_filename():
        h, fstr = tempfile.mkstemp(".revised.h5", dir=os.path.split(__file__)[0])
        os.close(h)
        yield fstr
        try:
            os.remove(fstr)
        except (FileNotFoundError, PermissionError):
            print("failed to remove "+fstr)


def test_api(filename, revised_filename):

    print("running the sample api")
    write_to_file = S999File(filename)
    if not revised_filename:
        revised_filename = filename + ".revised.h5"

    the_first_object = MyObject()
    the_first_object.data_value = "A sample string"
    my_other_first_object = MyObject(data_value="A sample string")
    assert the_first_object.data_value == my_other_first_object.data_value
    write_to_file.root.my_first_object = the_first_object

    a_location = MyLocation()
    # create values for all child attributes and by passing in True it would recurse all grandchildren and beyond too
    a_location.initialize_properties(True)
    a_location.utm_zone = 18
    a_location.east_bound_longitude = 33.5
    a_location.add_metadata("extraData", 12345)

    print("Show what attributes are held inside our location class, can use either the instance or the class name itself")
    print(a_location.get_standard_properties())
    # ['east_bound_longitude', 'extent_type_code', 'north_bound_latitude', 'south_bound_latitude', 'utm_zone', 'west_bound_longitude']
    print(MyLocation.get_standard_properties())
    # ['east_bound_longitude', 'extent_type_code', 'north_bound_latitude', 'south_bound_latitude', 'utm_zone', 'west_bound_longitude']
    assert a_location.get_standard_properties() == MyLocation.get_standard_properties()

    # the mapping to see the HDF5 naming only works on an instance - not the class unfortunately
    print(a_location.get_standard_properties_mapping())
    # returns --
    # {'eastBoundLongitude': 'east_bound_longitude',
    #  'extentTypeCode': 'extent_type_code',
    #  'northBoundLatitude': 'north_bound_latitude',
    #  'southBoundLatitude': 'south_bound_latitude',
    #  'utmZone': 'utm_zone',
    #  'westBoundLongitude': 'west_bound_longitude'}

    # get_standard_keys will return the names expected from an HDF5 file based on the S100 specs
    print(a_location.get_standard_keys())
    # returns --
    # ['eastBoundLongitude',
    #  'extentTypeCode',
    #  'northBoundLatitude',
    #  'southBoundLatitude',
    #  'utmZone',
    #  'westBoundLongitude']

    # get_all_keys will show all the HDF5 data names held in the object -- this can include non-standard data
    print(a_location.get_all_keys())
    # returns -- (notice the extraData in the list.
    # {'eastBoundLongitude',
    #  'extentTypeCode',
    #  'extraData',
    #  'northBoundLatitude',
    #  'southBoundLatitude',
    #  'utmZone',
    #  'westBoundLongitude'}

    # get rid of west -- it shouldn't show up in the HDF5 file after this
    del a_location.west_bound_longitude
    write_to_file.root.my_location_group = a_location

    write_to_file.root.data_group_create()  # this makes the DataGroups which is a list container for the DataGroupObject

    # Introspect the data group to figure out what it wants without reading the docs :)
    print(write_to_file.root.get_standard_properties())
    # ['data_group', 'dataset_with_names', 'my_first_object', 'my_location_group']
    print(write_to_file.root.data_group.get_standard_properties())
    # []
    print(type(write_to_file.root.data_group.metadata_type()))
    # <class 'sample_api_test.DataGroupObject'>
    print(type(write_to_file.root.data_group.metadata_type()).get_standard_properties())
    # ['data_grid', 'name_of_data']

    data_1 = DataGroupObject()
    data_1.name_of_data = "spam"
    # data_1.name_of_data = data_1.name_of_data_type["spam"]
    data_1.data_grid = numpy.zeros([2, 5])
    write_to_file.root.data_group.append(data_1)
    del data_1.name_of_data

    data_2 = DataGroupObject()
    data_2.name_of_data = 2
    data_2.data_grid = numpy.ones([3, 4])
    write_to_file.root.data_group.append(data_2)

    data_3 = write_to_file.root.data_group.append_new_item()
    data_3.data_grid = numpy.arange(0, 10, .75)
    data_3.name_of_data = data_3.name_of_data_type["cheese"]

    write_to_file.root.data_group.append_new_item().name_of_data = MONTY(2)

    attr_1 = datasetWithNames()
    attr_1.initialize_properties(True)

    attr_2 = datasetWithNames()
    attr_2.initialize_properties(True)
    attr_2.attr_str = "A custom string this time"
    attr_2.attr_int = 27
    attr_2.attr_float = 35.0
    write_to_file.root.dataset_with_names = DatasetWithNames_List((attr_1, attr_2))  # also could have used _create and append/append_new_item

    write_to_file.write()
    write_to_file.close()

    read_from_file = S999File(filename, "r")
    assert read_from_file.root.dataset_with_names[1].attr_int == 27
    assert read_from_file.root.dataset_with_names[0].attr_str in b"used a default string"
    assert read_from_file.root.data_group[1].name_of_data == MONTY(2)
    assert read_from_file.root.data_group[2].data_grid[1] == 0.75  # the second element of the range
    assert read_from_file.root.my_location_group.east_bound_longitude == 33.5
    try:
        read_from_file.root.my_location_group.west_bound_longitude
    except:
        pass  # all good, the value shouldn't exist
    else:
        assert read_from_file.root.my_location_group.west_bound_longitude is None  # this should not exist, even as None

    copy_of_file = S999File(revised_filename)
    copy_of_file.root = read_from_file.root
    # this shows how to initialize on creation
    copy_of_file.root.my_location_group = MyLocation(utm_zone=22, east_bound_longitude=11, extra_attr="This shouldn't even be here, but it works")
    copy_of_file.write()
    copy_of_file.close()


def test_direct_access_just_attributes(revised_filename):
    """ Test/demonstrate accessing an HDF5 file directly on a simple datatype"""
    h5file = h5py.File(revised_filename)
    obj_location = "/myFirstObject"
    new_string = "A revised simple string"
    data = MyObject()
    data.read(h5file[obj_location])
    print(data.data_value)  # prints:  'A sample string'
    data.data_value = new_string
    # overwrite
    data.write(h5file[obj_location])
    # duplicate
    h5file.require_group(obj_location + "_new")
    data.write(h5file[obj_location + "_new"])

    revised_data = MyObject()
    revised_data.read(h5file[obj_location])
    assert revised_data.data_value == data.data_value
    assert revised_data.data_value == new_string

    h5file.close()


def test_direct_access_attr_and_dataset(revised_filename):
    """ Test/demonstrate accessing an HDF5 file directly for enum and dataset (numpy array)"""
    h5file = h5py.File(revised_filename)
    obj_location = "/dataGroup_003"
    new_enum = MONTY['spam']
    new_val = 99
    data = DataGroupObject()
    data.read(h5file[obj_location])
    # print(data.name_of_data, data.data_grid)
    data.name_of_data = new_enum
    data.data_grid[3] = new_val
    # overwrite
    data.write(h5file[obj_location])
    # duplicate -- we will need to create a group for it to write into or put it into an existing group
    h5file.create_group(obj_location + "_new")
    data.write(h5file[obj_location + "_new"])

    revised_data = DataGroupObject()
    revised_data.read(h5file[obj_location])
    assert revised_data.name_of_data == data.name_of_data
    assert revised_data.name_of_data == new_enum
    assert revised_data.name_of_data == new_enum
    assert revised_data.data_grid[3] == new_val

    h5file.close()


def test_direct_access_list(revised_filename):
    """ For a multi-occurrence object, test/demonstrate accessing an HDF5 file directly
     Lists actually name themselves, so we will give it a different parent than the root as a test too"""
    h5file = h5py.File(revised_filename)
    obj_location = "/"  # this is the root of the file
    new_val = 35
    data = DataGroups()
    data_alternate_read = DataGroups()
    # show that the file object is the same as using ["/"]
    data.read(h5file)
    data_alternate_read.read(h5file[obj_location])
    print("number of data groups read", len(data))
    assert len(data) == len(data_alternate_read)
    assert len(data) == 4

    data[2].data_grid[4] = new_val
    # overwrite
    data.write(h5file)
    # duplicate
    h5file.create_group("extra_list")
    data.write(h5file["extra_list"])

    check_data = DataGroups()
    check_data.read(h5file)
    assert check_data[2].data_grid[4] == new_val

    check_other_data = DataGroups()
    check_other_data.read(h5file["extra_list"])
    assert check_other_data[2].data_grid[4] == new_val

    h5file.close()

def test_direct_access_compound_array(revised_filename):
    h5file = h5py.File(revised_filename)
    obj_location = "/datasetWithNames"  # this is the root of the file
    data = DatasetWithNames_List()
    data.read(h5file[obj_location])
    data[0].attr_int = 5
    # overwrite
    # the dataset type needs to create a dataset with it's own name, so we pass in the parent.
    data.write(h5file["/"])

    # duplicate
    data[0].attr_str = "duplicated"
    h5file.require_group("/new_compound_array")
    data.write(h5file["/new_compound_array"])

    h5file.close()

def test_change_names_on_new_data(revised_filename):
    """ This plays some games with the attribute_name.  Because the data is held in a dictionary based on the hdf5 names,
    changing the mapping between python name and HDF5 name can have consequences.  """
    h5file = h5py.File(revised_filename)

    # set up a standard object but store it in a non-standard group
    obj_with_standard_name = MyObject()
    obj_with_standard_name.data_value = "standard"
    assert obj_with_standard_name.data_value_attribute_name == "dataValue"
    h5file.require_group("/test_standard_name")
    obj_with_standard_name.write(h5file["/test_standard_name"])

    # change just the instance's name for HDF5, doing this BEFORE adding data works fine
    obj_with_non_standard_name = MyObject()
    obj_with_non_standard_name.data_value_attribute_name = "Change_instance_name"
    obj_with_non_standard_name.data_value = "Testing just the curreent instance"
    h5file.require_group("/test_instance_names")
    obj_with_non_standard_name.write(h5file["/test_instance_names"])

    assert obj_with_non_standard_name.data_value == "Testing just the curreent instance"

    # Change the class definition, which can be easier if ALL the data you ever want to read uses that different naming
    MyObject.data_value_attribute_name = "Change_all_classes"
    changed_class_obj = MyObject()
    changed_class_obj.data_value = "Change_the_class_itself"
    h5file.require_group("/test_class_names")
    changed_class_obj.write(h5file["/test_class_names"])

    assert changed_class_obj.data_value_attribute_name == "Change_all_classes"
    assert obj_with_non_standard_name.data_value_attribute_name == "Change_instance_name"

    # but watch out, existing data will also get the new name (but the one we changed just the instance of will be unaffected).
    # our standard name data will now have data that is orphaned and adding/changing the data via the api will only use the new names
    obj_with_standard_name.data_value = "still standard?"
    h5file.require_group("/test_standard_whoa")
    obj_with_standard_name.write(h5file["/test_standard_whoa"])
    assert obj_with_standard_name.data_value_attribute_name == "Change_all_classes"
    assert obj_with_non_standard_name.data_value == "Testing just the curreent instance"

    h5file.close()

def test_changing_names_on_existing_data(revised_filename):
    """ Change the attribute names in existing data, this requires re-mapping the old data to the new name or deleting the old."""
    h5file = h5py.File(revised_filename)
    obj_location = "/datasetWithNames"  # this is the root of the file
    data = DatasetWithNames_List()
    data.read(h5file[obj_location])

    # change the names and values of exising data
    for index, compund_arr in enumerate(data):
        del compund_arr.attr_int  #delete the old data before we rename
        compund_arr.attr_int_attribute_name = "changed_individual_int"
        compund_arr.attr_int = (index + 5) * 2

    # change the name in all the classes in existence at once.  This could corrupt other data in memory, in theory!
    old_name = datasetWithNames.attr_float_attribute_name
    datasetWithNames.attr_float_attribute_name = "changed_class_float"
    for index, compund_arr in enumerate(data):
        compund_arr.__delattr__(old_name)
        compund_arr.attr_float = (index + 6) * 3

    h5file.require_group("/compound_array_changed_names")
    data.write(h5file["/compound_array_changed_names"])

    h5file.close()


if __name__ == "__main__":
    filename = "test_sample.h5"
    rev_filename = "test_sample.revised.h5"
    for fname in (filename, rev_filename):
        try:
            os.remove(fname)
        except (FileNotFoundError, PermissionError):
            pass

    test_api(filename, rev_filename)
