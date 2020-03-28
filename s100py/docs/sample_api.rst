An Example API
===============

.. _@property:  https://docs.python.org/3/library/functions.html#property>
.. _Pycharm Live Template: https://www.jetbrains.com/help/pycharm/tutorial-creating-and-applying-live-templates-code-snippets.html

Let's create a new sample api to show the basic pieces of the s100py framework.  We'll lay out an example
data spec and then create the code and classes to make it work.  The basic data types are:

    - basic attribute (string, int, float, enum)
    - named 'compound array' dataset
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
There are also two helpers to help discover the type of the data that should be stored and one to create
an instance of the data with an appropriate default.

Let's start from the bottom of the tree and make a group object that has a string attribute.
This could look like::

    class myFirstObject:
        def __init__(self, val):
            self.data_value=val
            self.data_value_attribute_name = "dataValue"
            self.data_value_type = float
        def data_value_create(self):
            self.data=999

But doing this would require the data types read and write themselves to disk and do other housekeeping.
Also data validation would be harder.

For datatypes that have child data members there is a base class :any:`s1xx.S1XX_Attributes_base`.
This eventually equates to an HDF5 Group object that can have HDF5 attributes, datasets or subgroups.
To allow for qc of the input and output of data by the user or what is being read from disk,
we used properties (`@property`_) instead of plain member data.
They work the same when writing code but allow the API to validate data in addition.

The data is actually held in a 'private' dictionary "_attributes" which is based on the S100 names (not the python style naming).
This allows for any data in an HDF5 file when read to also be stored, so if reading non-standard data or a newer
or older format, there is a better chance that the data will still behave naturally and not be lost.
We also provide type hints to help auto create the docs and when writing code in an IDE or shell that supports typehints.

Here is how the API is actually implemented, notice the name MyObject does not have to match the S100 name
of myFirstObject::

    # Import modules we are going to use below.
    from enum import Enum
    import numpy
    from typing import Callable, Iterator, Union, Optional, List, Type
    from s100py import s1xx, s100

    class MyObject(s1xx.S1XX_Attributes_base):
        """ Create our first data class with properties etc """
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

That would be a lot of typing, but there is a template in :any:`extending_the_api` that makes it much faster
and is even better when used as a `PyCharm Live Template`_.  If using PyCharm just type in the S100 camelcase name
and run the live template and it will automatically make the python style name.  Hit tab and you can specify the datatype
and it will fill it into multiple locations at once for you.

To recap:
    - @property to get the data and do any reformatting needed etc.
    - @property.setter potential validation or other checks/changes to incoming data
    - *_attribute_name which defines the conversion from python naming to HDF5 (S100) naming
    - *_type to help the user of the api know the type to use and for the api to load from disk
    - *_create to make empty objects or supply default values as specified by S100

Now let's try a datatype that has eastBoundLongitude, westBoundLongitude, northBoundLongitude, southBoundLatitude and
utmZone.  The first four attributes are already part of an :any:`s100.GeographicBoundingBox` so let's derive a class
from there.

Use the template for utmZone and notice the attribute will be an int (and in PyCharm you'll be done in an instant).
Let's also add some limits on the zone number in the @property.setter, define an 'empty_zone' and
make 'empty_zone' the default for the utm_zone::

    class MyLocation(s100.GeographicBoundingBox):
        empty_zone = 999  # a way to mark the utm not being set
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
            """ Use 999 by default """
            self.utm_zone = self.utm_zone_type(self.empty_zone)


Next is a multi-occurrence object.  These are groups that S100 says has an integer at the end of it's name, like Group_001.
To store these there is a class that makes them act as python lists, :any:`s1xx.S1XX_MetadataList_base`.
This class needs to know what the acceptable name patterns are for reading/writing the data,
the default is an underscore OR dot followed by one or more integers.
You also have to supply a `@property`_ "metadata_name" and "metadata_type" for the name and type of the data to be held in the list.

But first, our example says that this dataGroup_01 will contain an attribute and a rectangular grid dataset.
We know how to encode an attribute which is a simple string or number but not a dataset.
Actually, a straight rectangular grid is simple, it is just a property that has a numpy array or h5py dataset as it's type.

The other attribute says it's an enumeration.  Let's say the document defines:
    - "spam" = 1
    - "cheese" = 2

Let's encode that as a python enumeration::

    from enum import Enum
    class MONTY(Enum):
        spam = 1
        cheese = 2

Now let's make the class that has the enumeration and the dataset.  The enumeration data doesn't quite follow
the standard template, so there is a second one just for enumerations in :any:`extending_the_api`_ ::

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
            return return numpy.ndarray

        def data_grid_create(self):
            """ Creates a blank, empty or zero value for data_grid"""
            self.data_grid = self.data_grid_type()

Ok, now let's make the list object that will actually have these data groups.  Recall the :any:`s1xx.S1XX_MetadataList_base`
base class::

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

For the last datatype we'll make the compund dataset "datasetWithNames".  This is to encapsulate S100 specs that lay out
data with names, like attributes, but say they belong in a dataset.   The :any:`s1xx.S1XX_Dataset_base`_ takes care of this.
Similar to the List we jsut made above, this class uses a list to keep an arbitrary number of data arrays and read/write
them to HDF%.

For example, the S100 spec Table 10c-8 describes a compound array stored as a dataset which is more naturally used
as a multiple lists of attributes.  Our example will make a datatype to hold three attributes and a datatype that
holds them in a list.  Notice we will implement the get_write_order() to make the HDF5 array be written in the order
we want and not just by name.::

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
            self.attr_int = self.attr_int_type()


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
            self.attr_float = self.attr_float_type()


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
            self.attr_str = self.attr_str_type()

Now we'll wrap this data class inside a :any:`s1xx.S1XX_Dataset_base`_  class so it reads and writes to arrays
and can be accessed as a python list.::

    class DatasetWithNames_List(s1xx.S1XX_Dataset_base):

        @property
        def metadata_type(self) -> Type[type]:
            return datasetWithNames

        @property
        def metadata_name(self) -> str:
            return "datasetWithNames"

The final data class we'll make is make a root object that contains all the datatypes we just made and associate that with a
file object (which is derived from an h5py File).  The root object itself is just another
class derived from :any:`s1xx.S1XX_Attributes_base`.::


    class S999Root(s1xx.S1XX_Attributes_base):
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

The final thing to do is to associate the root data class to a S1XXFile.
The file is derived from a h5py.File object and will accept any of the creation arguments h5py will take.
All we need to do is add a product specification string and add a 'root' keyword. ::

    class S999File(s1xx.S1XXFile):
        PRODUCT_SPECIFICATION = numpy.string_('INT.IHO.S-Fake')

        def __init__(self, *args, **kywrds):
            # kywrds['root'] = S999Root
            super().__init__(*args, root=S999Root **kywrds)

All that is left is :any:`using_sample_api`