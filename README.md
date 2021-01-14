s100py
======
[![Build Status](https://travis-ci.com/noaa-ocs-s100/s100py.svg?branch=master)](https://travis-ci.com/noaa-ocs-s100/s100py)

Python Utilities and API for Working with IHO S-100 Data Formats

Overview
--------

This python package provides utilities for encoding hydrographic
datasets in the International Hydrographic Organization (IHO) S-100
format.

Background
----------

The IHO S-100 standard is a data framework for digital products and
services for hydrographic, maritime, and GIS communities, comprised of
multiple data encoding formats designed for interoperability with
Electronic Navigational Charts (ENCs).

The initial focus of this package is on three of the S-100 encoding
formats:

-   S-102 Bathymetry 
-   S-104 Water Level Information for Surface Navigation
-   S-111 Surface Currents

However, support for additional formats will likely be added in the
future.

For further information about S-100 formats, see the [IHO
website](http://s100.iho.int/S100/).

Using the API
-------------

- Follow [S-102](s100py/s102/README.md) to create an S-102 File
- Follow [S-104](s100py/s104/README.md) to create an S-104 File
- Follow [S-111](s100py/s111/README.md) to create an S-111 File


S100 API Developers
-------------------

S100 is a family of specifications where each data product has it’s own 
extensions to the basic format dictated by S100. S100 files are HDF5 where
the data’s types, names and structure is in the aforementioned specs. 
s100py uses pep8 naming for the api itself.

The S1XXFile class handles the top level data management and is derived from
the h5py.File class. There is a “root” data member for the file object which
is where the S100/HDF5 data is held. Each data specification should create 
it’s own S100 type file object but there is a basic s100.S100File class which
can read the general top level metadata that pertains to all data specifications.

The three primary data types in HDF5 are attributes, datasets and groups. 
s100py encapsulates groups (which can have attributes) with the S1xxAttributesBase 
class. Datasets are derived from either S1xxDatasetBase or S1xxGridsBase. 
S100 also adds groups that have a trailing number and can have an arbitrary 
number of occurrences. These objects are managed with the S1xxMetadataListBase 
class. Classes derived from S1xxDatasetBase and S1xxMetadataListBase will act as 
lists of groups (S1xxAttributesBase)

Each S100 class has data and some methods to determine what data is expected to be
contained therein. The function get_standard_properties() will show what child data
is referenced in the specs. The initialize_properties() method will create default 
values for all expected child data. Using initialize_properties(True) will recurse 
the data spec and create a skeleton of all children and their children. 
The get_standard_properties_mapping() show the HDF5 names and they pythonic pep8 
names that are used with them.

Inside of each S100 data class are data values and functions related to each data 
value. For a “dataname” value there are two properties, dataname_type and 
dataname_attribute_name, and a function dataname_create.

When the File is being written any data that has not been set or initialized will 
be omitted. This way optional data will not appear in the HDF5 file at all.



###Extending the API

- Create a new API from this framework or modify an existing one, use the following 
  classes from s100py.s1xx as well as the ones in s100py.s100<br/><br/>

  * S1xxAttributesBase -> For HDF5 groups with attributes, datasets and sub-groups <br/><br/>  

  * S1xxDatasetBase -> For datasets (numpy arrays)  <br/><br/>
      
  * S1xxMetadataListBase -> For groups that have multiple occurrences using the S100 
    naming scheme (Group_NNN) <br/><br/>

- To make data members for one of the above classes, here is how to make a template
  in PyCharm to speed creation: <br/><br/>
    
  - File-Settings-Editor-LiveTemplates  <br/><br/>
    - Create a new template (hit the plus button on the side) and name it S102 and give
    it a description paste in the code below (@property lines) <br/><br/>
    - Click the Edit Variables and for $attr$ under expression put: snakeCase(SELECTION)
    at the bottom is a line that says “applicable in” and has a hyperlinked word 
    (define or change) – click that and select Python <br/><br/>
    
- To use, highlight the camel case S102 attribute name and press `ctrl-alt-J` 
  (or `ctrl-alt-T`) and select your S102 from the list. It will fill in the 
  names and then put the cursor in the right place to specify the type
    
    - Take the selected text (name from the S102 doc) and:
    
      - Make a read only property returning the S102 HDF5 attribute name
      
      - Make a read/write property to contain the data
      
      - Make a property to tell the type of the data
      
      - make a _create() function to make a default value for the attribute<br/><br/>
    
**Code to put in your live template**
  ```    
      $attr$_attribute_name = return "$SELECTION$"  #: HDF5 naming
      
      @property
      def $attr$(self) -> $type$:
          return self._attributes[self.$attr$_attribute_name]
      
      @$attr$.setter
      def $attr$(self, val: $type$):
          self._attributes[self.$attr$_attribute_name] = val
      
      @property
      def $attr$_type(self) -> Type[$type$]:
          return $type$
      
      def $attr$_create(self):
          """ Creates a blank, empty or zero value for $attr$"""
          # noinspection PyAttributeOutsideInit
          # pylint: disable=attribute-defined-outside-init
          self.$attr$ = self.$attr$_type()
  ```
- For enumeration data types use this template which is very similar. Again, click the 
  Edit Variables and for $attr$ under expression put: snakeCase(SELECTION) at the bottom
  is a line that says “applicable in” and has a hyperlinked word (define or change) <br/><br/>
  - Click and select Python

  ```  
      $attr$_attribute_name = return "$SELECTION$"  #: HDF5 naming
      
      @property
      def $attr$(self) -> $type$:
          return self._attributes[self.$attr$_attribute_name]
      
      @$attr$.setter
      def $attr$(self, val: Union[int, str, $type$]):
          self.set_enum_attribute(val, self.$attr$_attribute_name, self.$attr$_type)
      
      @property
      def $attr$_type(self) -> Type[$type$]:
          return $type$
      
      def $attr$_create(self):
          """ Creates a value using the first item in the enumeration of $attr$"""
          # noinspection PyAttributeOutsideInit
          # pylint: disable=attribute-defined-outside-init
          self.$attr$ = list(self.$attr$_type)[0]
  ```
- This template makes an attribute but specifies the type as well before you run it. Put the
  HDF5 name from the S100+ spec first followed by an arrow (->) then the type. <br/><br/>
    - Click the Edit Variables:
  ```
      for “type” under expression put: regularExpression(SELECTION, “.*->”, “”)
  
      for “HDF5NAME” for expression put: regularExpression(SELECTION, “->.*”, “”)
  
      for “attr” under expression put: snakeCase(regularExpression(SELECTION, “->.*”, “”))
  
      at the bottom is a line that says “applicable in” and has a hyperlinked word 
      (define or change) – click that and select Python template variable dialog

  ```
    
- Types can be basic python types or custom created classes.
  - Ex: eastBoundLongitude->float  
  ```
      $attr$_attribute_name = return "$HDF5NAME$"  #: HDF5 naming
      
      @property
      def $attr$(self) -> $type$:
          return self._attributes[self.$attr$_attribute_name]
      
      @$attr$.setter
      def $attr$(self, val: $type$):
          self._attributes[self.$attr$_attribute_name] = val
      
      @property
      def $attr$_type(self) -> Type[$type$]:
          return $type$
      
      def $attr$_create(self):
          """ Creates a blank, empty or zero value for $attr$
          $SELECTION$
          """
          # noinspection PyAttributeOutsideInit
          # pylint: disable=attribute-defined-outside-init
          self.$attr$ = self.$attr$_type()
  ```

- And finally a similar one for enumerations. Same syntax or HDF5 name followed by the Enumeration name.
    
    - Click the Edit Variables and:
    ```
    for “type” under expression put: regularExpression(SELECTION, “.*->”, “”)

    for “HDF5NAME” for expression put: regularExpression(SELECTION, “->.*”, “”)

    for “attr” under expression put: snakeCase(regularExpression(SELECTION, “->.*”, “”))

    at the bottom is a line that says “applicable in” and has a hyperlinked word 
    (define or change) – click that and select Python
    ```
- If you used the enumeration ‘MONTY’ from the sample api, this would look like dataName->MONTY
  ```
      $attr$_attribute_name = return "$HDF5NAME$"  #: HDF5 naming
      
      @property
      def $attr$(self) -> $type$:
          return self._attributes[self.$attr$_attribute_name]
      
      @$attr$.setter
      def $attr$(self, val: Union[int, str, $type$]):
          self.set_enum_attribute(val, self.$attr$_attribute_name, self.$attr$_type)
      
      @property
      def $attr$_type(self) -> Type[$type$]:
          return $type$
      
      def $attr$_create(self):
          """ Creates a blank, empty or zero value for $attr$
          $SELECTION$
          """
          # noinspection PyAttributeOutsideInit
          # pylint: disable=attribute-defined-outside-init
          self.$attr$ = list(self.$attr$_type)[0]
  ```

###An Example API

- Let’s create a new sample api to show the basic pieces of the s100py framework. We’ll lay out an 
  example data spec and then create the code and classes to make it work. The basic data types are: <br></br>

  - Basic attribute (string, int, float, enum)

  - Named ‘compound array’ dataset
    
  - Grid dataset

  - Named sub group 

  - Variable number of occurrence subgroup (e.g. Group_01, Group_02…)

  - File <br></br>
 
- Here is our pretend data spec:
  ```
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
  ```
- Each piece of data from the S100 spec there should be a python friendly (pep8 - “snake case”) name to hold
  the data. Also there is a read only property that contains the S100 name (which may be camel case and not 
  as human friendly). There are also two helpers to help discover the type of the data that should be stored
  and one to create an instance of the data with an appropriate default. <br></br>


- Let’s start from the bottom of the tree and make a group object that has a string attribute. This could 
  look like:
  ```python
  class myFirstObject:
      def __init__(self, val):
          self.data_value=val
          self.data_value_attribute_name = "dataValue"
          self.data_value_type = float
      def data_value_create(self):
          self.data=999
  ```
- But doing this would require the data types read and write themselves to disk and do other housekeeping. 
  Also data validation would be harder.

- For datatypes that have child data members there is a base class s1xx.S1xxAttributesBase. This eventually
  equates to an HDF5 Group object that can have HDF5 attributes, datasets or subgroups. To allow for qc of 
  the input and output of data by the user or what is being read from disk, we used properties (@property) 
  instead of plain member data. They work the same when writing code but allow the API to validate data in 
  addition. <br></br>

- The data is actually held in a ‘private’ dictionary “_attributes” which is based on the S100 names 
  (not the python style naming). This allows for any data in an HDF5 file when read to also be stored, 
  so if reading non-standard data or a newer or older format, there is a better chance that the data will 
  still behave naturally and not be lost. We also provide type hints to help auto create the docs and when
  writing code in an IDE or shell that supports typehints. <br></br>


- Here is how the API is actually implemented, notice the name MyObject does not have to match the S100 name
  of myFirstObject:

  ```python
  # Import modules we are going to use below.
  from enum import Enum
  import numpy
  from typing import Callable, Iterator, Union, Optional, List, Type
  from s100py import s1xx, s100
  
  class MyObject(s1xx.S1xxAttributesBase):
      """ Create our first data class with properties etc """
      @property
      def __version__(self) -> int:
          return 1
  
      data_value_attribute_name = "dataValue"  #: HDF5 naming
  
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
  ```
- That would be a lot of typing, but there is a template in Extending the API that makes it much faster and 
  is even better when used as a PyCharm Live Template. If using PyCharm just type in the S100 camelcase name 
  and run the live template and it will automatically make the python style name. Hit tab and you can specify 
  the datatype and it will fill it into multiple locations at once for you.<br></br>

- To recap:
  ```
      @property to get the data and do any reformatting needed etc.

      @property.setter potential validation or other checks/changes to incoming data

      *_attribute_name which defines the conversion from python naming to HDF5 (S100) naming

      *_type to help the user of the api know the type to use and for the api to load from disk

      *_create to make empty objects or supply default values as specified by S100
  ```
- Now let’s try a datatype that has eastBoundLongitude, westBoundLongitude, northBoundLongitude, 
  southBoundLatitude and utmZone. The first four attributes are already part of an s100.GeographicBoundingBox
  so let’s derive a class from there. <br></br>

- Use the template for utmZone and notice the attribute will be an int (and in PyCharm you’ll be done in an instant).
  Let’s also add some limits on the zone number in the @property.setter, define an ‘empty_zone’ and make 
  ‘empty_zone’ the default for the utm_zone:
  ```python
  class MyLocation(s100.GeographicBoundingBox):
      empty_zone = 999  # a way to mark the utm not being set
      @property
      def __version__(self) -> int:
          return 1
  
      utm_zone_attribute_name = "utmZone"  #: HDF5 naming
  
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
  ```
- Next is a multi-occurrence object. These are groups that S100 says has an integer at the end of it’s name,
  like Group_001. To store these there is a class that makes them act as python lists, s1xx.S1xxMetadataListBase.
  This class needs to know what the acceptable name patterns are for reading/writing the data, the default is 
  an underscore OR dot followed by one or more integers. You also have to supply a @property “metadata_name” 
  and “metadata_type” for the name and type of the data to be held in the list. <br></br>

- But first, our example says that this dataGroup_01 will contain an attribute and a rectangular grid dataset.
  We know how to encode an attribute which is a simple string or number but not a dataset. Actually, a straight
  rectangular grid is simple, it is just a property that has a numpy array or h5py dataset as it’s type. <br></br>

- The other attribute says it’s an enumeration. Let’s say the document defines:
  ```
    “spam” = 1
    “cheese” = 2
  ```
- Let’s encode that as a python enumeration:
  ```python
  from enum import Enum
  class MONTY(Enum):
      spam = 1
      cheese = 2
  ```
- Now let’s make the class that has the enumeration and the dataset. The enumeration data doesn’t quite follow
  the standard template, so there is a second one just for enumerations in Extending the API
  ```python
  class DataGroupObject(s1xx.S1xxAttributesBase):
      @property
      def __version__(self) -> int:
          return 1
  
      name_of_data_attribute_name = "nameOfData"  #: HDF5 naming
  
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
  
      data_grid_attribute_name = "dataGrid"  #: HDF5 naming
  
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
  ```
- Ok, now let’s make the list object that will actually have these data groups. Recall the s1xx.S1xxMetadataListBase base class:
  ```python
  class DataGroups(s1xx.S1xxMetadataListBase):
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
  ```
-For the last datatype we’ll make the compund dataset “datasetWithNames”. This is to encapsulate S100 specs
that lay out data with names, like attributes, but say they belong in a dataset. The s1xx.S1xxDatasetBase 
takes care of this. Similar to the List we jsut made above, this class uses a list to keep an arbitrary number
of data arrays and read/write them to HDF5. <br></br>

-For example, the S100 spec Table 10c-8 describes a compound array stored as a dataset which is more naturally
used as a multiple lists of attributes. Our example will make a datatype to hold three attributes and a datatype
that holds them in a list. Notice we will implement the get_write_order() to make the HDF5 array be written in 
the order we want and not just by name:
  ```python
   class datasetWithNames(s1xx.S1xxAttributesBase):
        def get_write_order(self):
            return ["attrInt", "attrStr", "attrFloat"]
    
        @property
        def __version__(self) -> int:
            return 1
    
        attr_int_attribute_name = "attrInt"  #: HDF5 naming
    
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
    
    
        attr_float_attribute_name = "attrFloat"  #: HDF5 naming
    
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
    
    
        attr_str_attribute_name = "attrStr"  #: HDF5 naming
    
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
  ```
- Now we’ll wrap this data class inside a s1xx.S1xxDatasetBase class so it reads and writes to arrays 
  and can be accessed as a python list.:
  ```python
  class DatasetWithNames_List(s1xx.S1xxDatasetBase):
  
      @property
      def metadata_type(self) -> Type[type]:
          return datasetWithNames
  
      @property
      def metadata_name(self) -> str:
          return "datasetWithNames"
  ```
- The final data class we’ll make is make a root object that contains all the datatypes we just made and
  associate that with a file object (which is derived from an h5py File). The root object itself is just
  another class derived from s1xx.S1xxAttributesBase.:
  ```python
  class S999Root(s1xx.S1xxAttributesBase):
      dataset_with_names_attribute_name = "datasetWithNames"  #: HDF5 naming
  
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
  
      data_group_attribute_name = "dataGroup"  #: HDF5 naming
  
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
  
      my_location_group_attribute_name = "myLocationGroup"  #: HDF5 naming
  
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
  
      my_first_object_attribute_name = "myFirstObject"  #: HDF5 naming
  
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
  ```
- The final thing to do is to associate the root data class to a S1XXFile. The file is derived from
  a h5py.File object and will accept any of the creation arguments h5py will take. All we need to do
  is add a product specification string and add a ‘root’ keyword.
  ```python
  class S999File(s1xx.S1XXFile):
      PRODUCT_SPECIFICATION = numpy.string_('INT.IHO.S-Fake')
  
      def __init__(self, *args, **kywrds):
          # kywrds['root'] = S999Root
          super().__init__(*args, root=S999Root **kywrds)
    ```
  
###Using Example API

- To use the artificial example API follow the steps below. The code below is also in the “tests/”
  folder in the “sample_api_test.py” file which should be runnable assuming you have a compliant, 
  working python environment. <br></br>

- Open a new file in the local directory

  `write_to_file = S999File("test.sample.h5")` <br></br>

- Create our first, basic object a “MyObject” and give it a string value and add it to the file
  ```
  the_first_object = MyObject()
  the_first_object.data_value = "A sample string"
  write_to_file.root.my_first_object = the_first_object
  ```
- We could have set the attribute via the constructor as well

  `my_other_first_object = MyObject(data_value="A sample string")` <br></br>

- Create a MyLocation object and add it to the file using the root object
  ```
  a_location = MyLocation()
  # create values for all child attributes and by passing in True it would recurse all grandchildren and beyond too
  a_location.initialize_properties(True)
  a_location.utm_zone = 18
  a_location.east_bound_longitude = 33.5
  ```
- Let’s tack on a non-S999 attribute and make up one called “extraData”. This shouldn’t exist per the spec
  but for extensibility it is supported.

  `a_location.add_metadata("extraData", 12345)`

The api supplies some introspection functionality too. You can use a python IDE to see the function names 
but some alternatives are available to see the python data names are and what S100 specifies the names in 
the file and what type of data it should be.

`get_standard_properties() #gives the python names of attributes based on the S100+ specs`
```
print("Show what attributes are held inside our location class, can use either the instance or the class name itself")
print(a_location.get_standard_properties())
print(MyLocation.get_standard_properties())
# both would return --
# ['east_bound_longitude', 'extent_type_code', 'north_bound_latitude', 'south_bound_latitude', 'utm_zone', 'west_bound_longitude']

get_standard_keys() returns a list of the S100+ spec names as will be written to the HDF5 file

# get_standard_keys will return the names expected from an HDF5 file based on the S100 specs
print(a_location.get_standard_keys())
# returns --
# ['eastBoundLongitude',
#  'extentTypeCode',
#  'northBoundLatitude',
#  'southBoundLatitude',
#  'utmZone',
#  'westBoundLongitude']

get_standard_properties_mapping() returns a dictionary of S100+ spec name (seen in HDF5) to python name, basically tying together the two previous functions. The mapping to see the HDF5 naming only works on an instance - not the class unfortunately.

print(a_location.get_standard_properties_mapping())
# returns --
# {'eastBoundLongitude': 'east_bound_longitude',
#  'extentTypeCode': 'extent_type_code',
#  'northBoundLatitude': 'north_bound_latitude',
#  'southBoundLatitude': 'south_bound_latitude',
#  'utmZone': 'utm_zone',
#  'westBoundLongitude': 'west_bound_longitude'}

get_all_keys will show all the HDF5 data names held in the object – this can include non-standard data If there is no additional data then it will be the same as get_standard_keys():

print(a_location.get_all_keys())
# returns -- (notice the extraData in the list.
# {'eastBoundLongitude',
#  'extentTypeCode',
#  'extraData',
#  'northBoundLatitude',
#  'southBoundLatitude',
#  'utmZone',
#  'westBoundLongitude'}
```
- We can use the standard python ‘del’ command to remove a piece of data we don’t want
  ```
  # get rid of west -- it shouldn't show up in the HDF5 file after this
  del a_location.west_bound_longitude
  ```
- Now we will add the location to the file:

  `write_to_file.root.my_location_group = a_location`

- Moving on the the “data_group” which will be a ‘list’ of data with _001 after it’s name in HDF5. <br></br>
- In python we can just treat it as a standard list and not worry about the implementation detail of naming. <br></br>
- We will use the _create() helper to make the list

  `write_to_file.root.data_group_create()`  # this makes the DataGroups which is a list container for the DataGroupObject

- We could also have said `root.data_group = DataGroups()`

- Let’s pretend we are in a prompt and are trying to make this without looking at the docs. 
Use the get_standard_properties() we saw above and we see the ‘data_group’ from the root. 
However when we get the standard properties of the data_group it returns an empty list. This is because
it is a special list class and it holds a list of one type of data which can be seen by checking the
metadata_type(). <br></br>
  ```
  # Introspect the data group to figure out what it wants without reading the docs :)
  print(write_to_file.root.get_standard_properties())
  # ['data_group', 'dataset_with_names', 'my_first_object', 'my_location_group']
  print(write_to_file.root.data_group.get_standard_properties())
  # []
  print(type(write_to_file.root.data_group.metadata_type()))
  # <class 'sample_api_test.DataGroupObject'>
  ```
- Once we realize the data_group is actually a list and what it holds, then we can see what would be 
  inside that data.
  ```
  print(type(write_to_file.root.data_group.metadata_type()).get_standard_properties())
  # ['data_grid', 'name_of_data']
  ```
- So let’s make some DataGroupObjects. The potentially tricky part is that the data_grid will be a HDF5.dataset
while the name_of_data was set to be an enumeration.<br></br>

- Actually both are treated like any other string or numeric attribute. The data_grid (dataset) just needs to 
be a numpy array or hdf5 dataset and it will work. The enumeration can be set using either the strings or 
numbers that the S100+ spec describes. <br></br>

- Being a list we can make an arbitrary amount of them. Let’s plan on three so we can make the enumerations in 
different ways and have different shaped datasets.
  ```
  data_1 = DataGroupObject()
  data_1.name_of_data = "spam"
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
  ```
Ok, let’s make a fourth element for the list. There is a append_new_item() which creates and returns the 
proper dataype. We’ll give it a name by passing in an enumeration value but no grid data. All on one line
no less!!

write_to_file.root.data_group.append_new_item().name_of_data = MONTY(2)

Compound datasets are lists of objects that will eventually be held as HDF5 datasets. They work as lists
of attribute classes in s100py so you never need to know they are datasets in reality.

- Let’s make our first entry and let it take on default values using the initialize_properties().
  ```
  attr_1 = datasetWithNames()
  attr_1.initialize_properties(True)
  ```
- Next we’ll make a second set of attributes
  ```
  attr_2 = datasetWithNames()
  attr_2.initialize_properties(True)
  attr_2.attr_str = "A custom string this time"
  attr_2.attr_int = 27
  attr_2.attr_float = 35.0
  ```
- And we’ll add it to the file by creating a DatasetWithNames_List and passing the attr_1, attr_2 
  to it’s constructor

  `write_to_file.root.dataset_with_names = DatasetWithNames_List((attr_1, attr_2))  # also could have used _create and append/append_new_item`

- Now we’ll save it all to disk
  ```
  write_to_file.write()
  write_to_file.close()
  ```
- Now let’s open that data file and spot check against the values we thought we wrote:
  ```
  read_from_file = S999File("test.sample.h5", "r")
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
  
  Finally, let’s make a copy of the data on disk and modify one of the values. You can use the hdfview app to confirm everything worked.
  
  copy_of_file = S999File("test.rewrite.h5")
  copy_of_file.root = read_from_file.root
  # this shows how to initialize on creation
  copy_of_file.root.my_location_group = MyLocation(utm_zone=22, east_bound_longitude=11, extra_attr="This shouldn't even be here, but it works")
  copy_of_file.write()
  ```

### Advanced Example API Usage

- Using the data and API shown in the previous steps we can look at some advanced usage. <br></br>

- The sample code below is also in the s100py/tests/sample_api_test.py where it should run with pytest. <br></br>

- The S100+ files are built on HDF5 and use h5py for access. You can access data directly and use the api to read
it by passing any api class the appropriate h5py object. <br></br>

- We’ll use a h5py.File object in the samples below, create one like this:

`h5file = h5py.File(filename)`

- So, for our first example let’s read the “MyObject” that would be found in the HDF under the group named
  “myFirstObject”
```
obj_location = "/myFirstObject"
data = MyObject()
data.read(h5file[obj_location])
print(data.data_value)  # prints:  'A sample string' based on the previous examples
```
- We can also edit the data any write it back to that same place, rather than writing an entire file at once.
```
data.data_value = "A revised simple string"
data.write(h5file[obj_location])  # overwrite the data in the HDF5 file
```
- We could even put it in the another place in the file. This wouldn’y be valid S100+ data but you could do
it for compatibility with other software for example.
```
h5file.require_group(obj_location + "_new")  # make a different subgroup called myFirstObject_new
data.write(h5file[obj_location + "_new"])  # and write it to the file
```
- Similarly we can access the enumerations of the datasets from the DataGroupObject we defined previously.
```
data = DataGroupObject()
data.read(h5file["/dataGroup_003"])
print(data.name_of_data, data.data_grid)  # should show original data
```
- Now we can change it and write it back out:
```
data.name_of_data = MONTY['spam']
data.data_grid[3] = 99
data.write(h5file["/dataGroup_003"])  # overwrite the existing data
```
- And make a duplicate record if wherever we want
```python
# duplicate to a new name -- we will need to create a group for it to write into or put it into an existing group
h5file.create_group("/dataGroup_003_new")
data.write(h5file["/dataGroup_003_new"])
```
- It works for entire lists too. Our list was located at the root and unlike the classes above, a list needs to
  start at the parent node. The lists have to name themselves (Group_001, Group_002) so need to start at it’s 
  parent level:
```
data = DataGroups()
data.read(h5file)  # side note: you can use h5file or h5file["/"] interchangeably

data[2].data_grid[4] = 35
data.write(h5file)  # overwrite
```
- We can create a copy of the list, but we can’t make a copy at the root where it already exists, so we have
  to create a subgroup and tell it to write itself there. We will end up with “/extra_list/Group_001”, 
  “/extra_list/Group_002”…
```
# duplicate
h5file.create_group("extra_list")
data.write(h5file["extra_list"])
```
- Compound arrays similarly, when using them you supply it the parent as it need to call create_dataset.
```
data = DatasetWithNames_List()
data.read(h5file["/"])
data[0].attr_int = 5
data.write(h5file["/"])
```
- We can put a copy in a new location but it will create a dataset (named “datasetWithNames” in this case) 
  under the parent location we supply.
```
data[0].attr_str = "duplicated"
h5file.require_group("/new_compound_array")
data.write(h5file["/new_compound_array"])
```
- Now to really abuse the system we can change names of the data but this is dangerous and not recommended. 
  It would both not adhere to the S100+ specs and also potentially be error prone.

- First we’ll change one instance of an object. Remember the MyObject has one string attribute that should be 
  named “dataValue”:
```
# set up a standard object but store it in a non-standard group
obj_with_standard_name = MyObject()
obj_with_standard_name.data_value = "standard"
assert obj_with_standard_name.data_value_attribute_name == "dataValue"
h5file.require_group("/test_standard_name")
obj_with_standard_name.write(h5file["/test_standard_name"])
```
- Now let’s make another copy of MyObject and change data_value_attribute_name which defines the mapping 
  from python name to S100+. Doing this BEFORE adding data works fine.
```
obj_with_non_standard_name = MyObject()
obj_with_non_standard_name.data_value_attribute_name = "Change_instance_name"
obj_with_non_standard_name.data_value = "Testing just the curreent instance"
h5file.require_group("/test_instance_names")
obj_with_non_standard_name.write(h5file["/test_instance_names"])
```
- If you want to get in trouble then you can change the class variable data_value_attribute_name which will 
then affect ALL the future and existing instances of MyObject.
```python
MyObject.data_value_attribute_name = "Change_all_classes"
changed_class_obj = MyObject()
changed_class_obj.data_value = "Change_the_class_itself"
h5file.require_group("/test_class_names")
changed_class_obj.write(h5file["/test_class_names"])
```
- And here is where the weird stuff happens, the obj_with_standard_name we made just above will also now write into
that new location too. Our standard name data will now have data that is orphaned and adding/changing the data via
the api will only use the new names. <br></br>

- This will end up having the old data under the old name and the new data under the new name – definitely not 
what someone probably wants.:
```python
obj_with_standard_name.data_value = "still standard?"
h5file.require_group("/test_standard_whoa")
obj_with_standard_name.write(h5file["/test_standard_whoa"])
```
- But, if you need to change some existing data, you can do it. Changing the attribute names in existing data will 
require re-mapping the old data to the new name or deleting the old data.

- Here we will change some of the items in the compound array. It had attr_int, attr_float and attr_str. First we’ll
  change each instance’s atrr_int naming and delete the old data and set new data.:
```python
data = DatasetWithNames_List()
data.read(h5file["/"])

# change the names and values of exising data
for index, compund_arr in enumerate(data):
    del compund_arr.attr_int  #delete the old data before we rename
    compund_arr.attr_int_attribute_name = "changed_individual_int"
    compund_arr.attr_int = (index + 5) * 2
```
- Then we’ll change the attr_float naming for the whole class (and any other existing data in the processes memory):
```python
old_name = datasetWithNames.attr_float_attribute_name
datasetWithNames.attr_float_attribute_name = "changed_class_float"
for index, compund_arr in enumerate(data):
    compund_arr.__delattr__(old_name)
    compund_arr.attr_float = (index + 6) * 3

h5file.require_group("/compound_array_changed_names")
data.write(h5file["/compound_array_changed_names"])
```


Requirements
------------

This codebase is written for Python 3 and relies on the following python
packages:

-   h5py
-   numpy
-   gdal


Installation
------------

This package requires the GDAL Python bindings be present, so it usually can\'t 
just be installed using `pip install gdal`. We recommend installing GDAL 
either through a package manager (e.g. `conda`, `apt`, `yum`, `pacman`) 
or by compiling from scratch. [Miniconda](https://docs.conda.io/en/latest/miniconda.html) 
is probably the easiest method.

Once `gdal` has been installed, s100py can be installed using `pip`:

```bash
pip install s100py
```

Authors
-------

-   Barry Gallagher, <barry.gallagher@noaa.gov>
-   Erin Nagel (UCAR), <erin.nagel@noaa.gov>
-   Jason Greenlaw (ERT), <jason.greenlaw@noaa.gov>
-   Glen Rice, <glen.rice@noaa.gov>

License
-------

This work, as a whole, falls under Creative Commons Zero (see
[LICENSE](LICENSE)).

Disclaimer
----------

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an 'as is' basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.
